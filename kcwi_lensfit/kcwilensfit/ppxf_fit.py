"""Stellar-continuum fitting with pPXF (Cappellari & Emsellem 2004; Cappellari 2017, 2023).

What this module adds on top of calling ``ppxf`` directly -- each item a real way a lens-galaxy
redshift or sigma goes wrong:

* **One log-rebinning of the joint blue+red spectrum.**  Rebinning each arm separately and gluing
  them violates pPXF's requirement that ``velscale = c*dln(lam)`` is *exactly* constant.
* **Frame-exact wavelengths.**  The galaxy is put in the template library's native frame (air for
  E-MILES, vacuum for FSPS/GALAXEV/XSL) by converting the wavelength *array*, not by a constant
  factor, so no spurious ~85 km/s velocity offset (and no residual dispersion of it).
* **Resolution matching in the right direction.**  ``sps_lib`` can only *broaden* templates.  When
  the galaxy's (rest-frame) resolution is finer than a template's native one -- the usual case for
  z > 0.5 KCWI data against 2.5 A E-MILES -- the *galaxy* is degraded to the template floor (with
  the noise reduced accordingly) so the fitted sigma is not inflated.
* **A free relative red/blue flux scale** solved in closed form between pPXF iterations.
* **Multiplicative polynomials** by default (they cannot erase or fake line depth).
* Bootstrap uncertainties and degree/library sweeps for the systematic error.
"""
from __future__ import annotations

import io
import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from .extraction import ExtractedSpectrum
from .masks import (BALMER, FORBIDDEN_NEBULAR, in_windows, line_windows, merge_windows)
from .templates import Library, UV_GAS_LINES, resolve_library
from .wavelength import C_KMS, air_to_vac, vac_to_air, velocity_to_z

try:
    from ppxf.ppxf import ppxf as _ppxf
    import ppxf.ppxf_util as util
    import ppxf.sps_util as sps_util
    HAVE_PPXF = True
except ImportError:                                           # pragma: no cover
    HAVE_PPXF = False

ProgressCB = Optional[Callable[[float, str], None]]
Window = Tuple[float, float]


# ---------------------------------------------------------------------------
# template cache (cropping the library to the needed rest range makes sps_lib ~8x faster)
# ---------------------------------------------------------------------------
class _TemplateStore:
    def __init__(self):
        self._raw: Dict[str, dict] = {}
        self._crop: Dict[tuple, bytes] = {}

    def raw(self, path: str) -> dict:
        if path not in self._raw:
            d = np.load(path)
            self._raw[path] = {k: d[k] for k in d.files}
        return self._raw[path]

    def native_fwhm(self, path: str, lam: np.ndarray) -> np.ndarray:
        r = self.raw(path)
        return np.interp(lam, r["lam"], r["fwhm"])

    def cropped_file(self, path: str, lam_lo: float, lam_hi: float, age_range=None, metal_range=None,
                     stride: int = 1) -> io.BytesIO:
        """In-memory .npz limited to the needed rest range and (optionally) an age/metallicity subset.

        ``sps_lib`` smooths *every* template before applying its own age/metal cuts, so cutting
        first (and cropping in wavelength) is several times faster.
        """
        key = (path, round(lam_lo / 100.0), round(lam_hi / 100.0), age_range, metal_range, stride)
        if key not in self._crop:
            r = self.raw(path)
            lo, hi = 0.93 * lam_lo, 1.07 * lam_hi
            m = (r["lam"] >= lo) & (r["lam"] <= hi)
            if m.sum() < 50:
                m = np.ones_like(r["lam"], bool)
            T, ages, metals, masses = r["templates"][m], r["ages"], r["metals"], r["masses"]
            ia = np.arange(ages.size)
            im = np.arange(metals.size)
            if age_range is not None:
                sel = np.flatnonzero((ages >= age_range[0]) & (ages <= age_range[1]))
                ia = sel if sel.size else ia
            if metal_range is not None:
                sel = np.flatnonzero((metals >= metal_range[0]) & (metals <= metal_range[1]))
                im = sel if sel.size else im
            if stride > 1:
                ia, im = ia[::stride], im[::stride]
            T = T[:, ia][:, :, im]
            buf = io.BytesIO()
            np.savez(buf, templates=T, masses=masses[np.ix_(ia, im)], ages=ages[ia], metals=metals[im],
                     lam=r["lam"][m], fwhm=r["fwhm"][m])
            self._crop[key] = buf.getvalue()
        return io.BytesIO(self._crop[key])


_STORE = _TemplateStore()


# ---------------------------------------------------------------------------
@dataclass
class FitSpectrum:
    """A joint (or single-arm) spectrum on one logarithmic wavelength grid, ready for pPXF."""
    lam: np.ndarray               # (n,) wavelength in the library frame, log-spaced
    flux: np.ndarray
    noise: np.ndarray
    fwhm: np.ndarray              # (n,) instrumental FWHM [A], observed frame
    good_base: np.ndarray         # (n,) bool: finite, outside gaps and static masks
    velscale: float
    is_red: np.ndarray            # (n,) bool
    cutover: float                # A, air
    frame: str                    # 'air' | 'vacuum'
    arms: Tuple[str, ...]
    windows_air: List[Window]


def build_fit_spectrum(spec: ExtractedSpectrum, arms: Sequence[str], fwhm: Dict[str, float], frame: str,
                       windows_air: Sequence[Window], scale_r: float = 1.0,
                       velscale: Optional[float] = None) -> FitSpectrum:
    """Concatenate arms at one cut-over wavelength and log-rebin ONCE."""
    parts = []
    for arm in arms:
        a = spec.arms[arm]
        g = a.good
        if g.sum() < 20:
            raise ValueError(f"'{spec.id}': arm '{arm}' has fewer than 20 usable pixels")
        parts.append((arm, a.wave[g], a.flux[g], a.var[g]))
    parts.sort(key=lambda p: p[1].mean())
    if len(parts) == 2:
        (_, wb, fb, vb), (_, wr, fr, vr) = parts
        lo, hi = wr.min(), wb.max()
        mid = 0.5 * (lo + hi)               # centre of the overlap (or of the gap, if there is one)
        kb, kr = wb < mid, wr >= mid
        wave = np.concatenate([wb[kb], wr[kr]])
        flux = np.concatenate([fb[kb], scale_r * fr[kr]])
        var = np.concatenate([vb[kb], scale_r ** 2 * vr[kr]])
        is_red = np.concatenate([np.zeros(kb.sum(), bool), np.ones(kr.sum(), bool)])
        names = [p[0] for p in parts]
        fw = np.where(is_red, fwhm[names[1]], fwhm[names[0]])
    else:
        arm, wave, flux, var = parts[0]
        mid = wave.max() + 1.0 if arm == "blue" else wave.min() - 1.0
        is_red = np.full(wave.shape, arm == "red")
        fw = np.full(wave.shape, fwhm[arm])
    wave_fit = air_to_vac(wave) if frame == "vacuum" else wave
    stacked = np.column_stack([flux, var, fw, is_red.astype(float)])
    out, ln_lam, vs = util.log_rebin(wave_fit, stacked, velscale=velscale)
    lam = np.exp(ln_lam)
    noise = np.sqrt(np.clip(out[:, 1], 1e-30, None))
    # pixels that fall inside a data gap (> 3 native pixels) carry interpolated nonsense
    dl = np.diff(wave_fit)
    gaps = np.flatnonzero(dl > 3.0 * np.median(dl))
    good = np.isfinite(out[:, 0]) & np.isfinite(noise)
    for gi in gaps:
        good &= ~((lam > wave_fit[gi]) & (lam < wave_fit[gi + 1]))
    win_f = [(float(air_to_vac(lo)), float(air_to_vac(hi))) if frame == "vacuum" else (lo, hi)
             for lo, hi in windows_air]
    good &= ~in_windows(lam, win_f)
    return FitSpectrum(lam=lam, flux=out[:, 0], noise=noise, fwhm=out[:, 2], good_base=good, velscale=float(vs),
                       is_red=out[:, 3] > 0.5, cutover=float(mid), frame=frame, arms=tuple(arms),
                       windows_air=list(windows_air))


def snr_per_angstrom(fs: FitSpectrum, good: Optional[np.ndarray] = None) -> float:
    g = fs.good_base if good is None else good
    g = g & (fs.noise > 0)
    if g.sum() < 5:
        return float("nan")
    dl = np.median(np.diff(fs.lam)[g[:-1]]) if g[:-1].sum() else np.median(np.diff(fs.lam))
    return float(np.median(fs.flux[g] / fs.noise[g]) / np.sqrt(dl))


# ---------------------------------------------------------------------------
def _nebular_windows_obs(z: float, kinds: str, half_kms: float, frame: str, wmin: float, wmax: float
                         ) -> List[Window]:
    names = {"none": [], "forbidden": FORBIDDEN_NEBULAR, "all": FORBIDDEN_NEBULAR + BALMER}[kinds]
    win = line_windows(names, z, half_kms, wmin, wmax)          # air
    if frame == "vacuum":
        win = [(float(air_to_vac(a)), float(air_to_vac(b))) for a, b in win]
    return win


@dataclass
class FitSettings:
    """Everything one pPXF call needs besides the data."""
    library: Library
    moments: int = 2
    degree: int = -1
    mdegree: int = 8
    age_range: Optional[Tuple[float, float]] = (0.5, 14.0)
    metal_range: Optional[Tuple[float, float]] = (-1.0, 0.3)
    sigma_init: float = 200.0
    vel_bound: float = 1500.0
    mask_nebular: str = "forbidden"
    nebular_half_kms: float = 500.0
    has_gas: bool = False
    sigma_fixed: Optional[float] = None
    rest_range: Optional[Tuple[float, float]] = None   # restrict fit to this rest window [A]


def _degrade_galaxy_to_template(lam_rest, flux, noise, fwhm_rest, tem_fwhm):
    """Broaden the galaxy where it is sharper than the template's native resolution."""
    under = fwhm_rest < tem_fwhm
    frac = float(np.mean(under))
    if frac <= 0.01:
        return flux, noise, fwhm_rest, 0.0
    extra_sig = np.sqrt(np.clip(tem_fwhm ** 2 - fwhm_rest ** 2, 0, None)) / np.sqrt(4 * np.log(4))
    flux2 = util.varsmooth(lam_rest, flux, extra_sig)
    # NOTE: smoothing lowers the per-pixel scatter but makes neighbouring pixels correlated, so the
    # information content is NOT increased.  The noise array is therefore deliberately left
    # unchanged (conservative): reducing it would let chi2 count correlated pixels as independent
    # and make every error bar too small.
    return flux2, noise, np.maximum(fwhm_rest, tem_fwhm), frac


def fit_at_z(fs: FitSpectrum, z_guess: float, st: FitSettings, extra_bad: Sequence[Window] = (),
             quiet: bool = True) -> dict:
    """One pPXF call with the redshift guess ``z_guess`` fixed in the template mapping.

    The returned ``z`` is ``(1+z_guess) exp(V/c) - 1`` where V is pPXF's fitted velocity about the
    guess (so the answer is independent of the guess as long as |V| stays small).
    """
    if not HAVE_PPXF:
        raise ImportError("ppxf is not installed:  pip install ppxf")
    lam_rest = fs.lam / (1.0 + z_guess)
    fwhm_rest = fs.fwhm / (1.0 + z_guess)

    good = fs.good_base.copy()
    neb = _nebular_windows_obs(z_guess, st.mask_nebular, st.nebular_half_kms, fs.frame, fs.lam[0], fs.lam[-1])
    good &= ~in_windows(fs.lam, neb)
    if extra_bad:
        ew = [(float(air_to_vac(a)), float(air_to_vac(b))) if fs.frame == "vacuum" else (a, b) for a, b in extra_bad]
        good &= ~in_windows(fs.lam, ew)
    rest_ignored = False
    if st.rest_range is not None:
        win = (lam_rest >= st.rest_range[0]) & (lam_rest <= st.rest_range[1])
        if (good & win).sum() >= 200:
            good &= win
        else:
            rest_ignored = True        # the data do not reach the requested rest window: use everything

    lib = st.library
    # crop the library to the rest-frame range actually needed (big speed-up for sps_lib)
    cf = _STORE.cropped_file(lib.path, float(lam_rest.min()), float(lam_rest.max()),
                             None if lib.is_empirical else st.age_range, None if lib.is_empirical else st.metal_range)
    tem_fwhm = _STORE.native_fwhm(lib.path, lam_rest)
    flux, noise = fs.flux, fs.noise
    flux, noise, fwhm_rest_eff, deg_frac = _degrade_galaxy_to_template(lam_rest, flux, noise, fwhm_rest, tem_fwhm)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sps = sps_util.sps_lib(cf, fs.velscale, {"lam": lam_rest, "fwhm": fwhm_rest_eff})
    stars = sps.templates.reshape(sps.templates.shape[0], -1)
    if stars.shape[1] == 0:
        raise ValueError("No templates left after the age/metallicity cuts -- widen the ranges.")
    # pPXF needs the template array to span (at least) the galaxy's rest-frame range.  A library with
    # limited coverage (e.g. the empirical kc96 spectra, 1235-9945 A) cannot at high z, so pad the template
    # grid with zero-flux pixels (constant step in ln lambda) and mask the galaxy pixels the library
    # does not actually cover.
    ln_t = np.asarray(sps.ln_lam_temp, float)
    d_ln = (ln_t[-1] - ln_t[0]) / (ln_t.size - 1)
    lo_cov, hi_cov = np.exp(ln_t[0]), np.exp(ln_t[-1])
    good &= (lam_rest >= lo_cov * 1.001) & (lam_rest <= hi_cov / 1.001)
    n_lo = max(0, int(np.ceil((ln_t[0] - np.log(lam_rest.min() * 0.995)) / d_ln)))
    n_hi = max(0, int(np.ceil((np.log(lam_rest.max() * 1.005) - ln_t[-1]) / d_ln)))
    lam_temp = sps.lam_temp
    if n_lo or n_hi:
        stars = np.vstack([np.zeros((n_lo, stars.shape[1])), stars, np.zeros((n_hi, stars.shape[1]))])
        ln_t = ln_t[0] + d_ln * (np.arange(stars.shape[0]) - n_lo)
        lam_temp = np.exp(ln_t)
    templates = stars
    component = np.zeros(stars.shape[1], int)
    moments = st.moments
    start = [0.0, st.sigma_init if st.sigma_fixed is None else st.sigma_fixed]
    gas_names = None
    if st.has_gas:
        gt, gas_names, _ = util.emission_lines(ln_t, [lam_rest.min(), lam_rest.max()],
                                               float(np.mean(fwhm_rest_eff)), vacuum=(fs.frame == "vacuum"))
        if len(gas_names) == 0:
            names = [n for n, w in UV_GAS_LINES.items() if lam_rest.min() < w < lam_rest.max()]
            if not names:
                raise ValueError("has_gas requested but no gas line falls in the fitted rest range")
            gas_names = names
            gt = util.gaussian(ln_t, np.array([UV_GAS_LINES[n] for n in names]),
                               float(np.mean(fwhm_rest_eff)), pixel=True)
        templates = np.column_stack([stars, gt])
        component = np.concatenate([component, np.ones(gt.shape[1], int)])
        moments = [st.moments, 2]
        start = [[0.0, st.sigma_init], [0.0, st.sigma_init]]

    fixed = bounds = None
    if st.sigma_fixed is not None and not st.has_gas:
        fixed = [0, 1] + [0] * (st.moments - 2)
        bounds = [[-st.vel_bound, st.vel_bound], [st.sigma_fixed - 1e-3, st.sigma_fixed + 1e-3]] + \
                 [[-0.3, 0.3]] * (st.moments - 2)
        start = [0.0, st.sigma_fixed]
    elif not st.has_gas:
        bounds = [[-st.vel_bound, st.vel_bound], [20.0, 600.0]] + [[-0.3, 0.3]] * (st.moments - 2)
    gp = np.flatnonzero(good)
    if gp.size < 100:
        raise ValueError(f"only {gp.size} unmasked pixels at z={z_guess:.4f}")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pp = _ppxf(templates, flux, noise, fs.velscale, start=start, moments=moments, degree=st.degree,
                   mdegree=st.mdegree, lam=lam_rest, lam_temp=lam_temp,
                   component=component if st.has_gas else 0, goodpixels=gp, quiet=quiet, plot=False,
                   fixed=fixed, bounds=bounds)
    sol = np.atleast_2d(pp.sol)
    err = np.atleast_2d(pp.error)
    chi2 = float(pp.chi2)
    errs = err[0] * np.sqrt(max(chi2, 1.0))
    v, sig = float(sol[0][0]), float(sol[0][1])
    z = float((1.0 + z_guess) * np.exp(v / C_KMS) - 1.0)
    out = dict(sigma_at_bound=bool(st.sigma_fixed is None and (sig >= 595.0 or sig <= 25.0)), z=z, v=v, sigma=sig, v_err=float(errs[0]), sigma_err=float(errs[1]) if len(errs) > 1 else np.nan,
               z_err=float((1.0 + z) * errs[0] / C_KMS), chi2=chi2, pp=pp, z_guess=z_guess,
               at_bound=bool(abs(abs(v) - st.vel_bound) < 5.0), n_good=int(gp.size), goodpixels=gp,
               lam_rest=lam_rest, flux=flux, noise=noise, templates=templates, component=component,
               sps=sps, lam_temp=lam_temp, settings=st, degraded_frac=deg_frac, fs=fs, gas_names=gas_names, rest_range_ignored=rest_ignored,
               moments_arg=moments)
    if st.moments >= 4 and sol.shape[1] >= 4:
        out["h3"], out["h4"] = float(sol[0][2]), float(sol[0][3])
    if st.has_gas:
        out["z_gas"] = float((1.0 + z_guess) * np.exp(sol[1][0] / C_KMS) - 1.0)
    return out


def fit_joint(spec: ExtractedSpectrum, arms: Sequence[str], fwhm: Dict[str, float], frame: str,
              windows_air: Sequence[Window], z_guess: float, st: FitSettings, fit_scale: bool = True,
              extra_bad: Sequence[Window] = (), n_iter: int = 6, tol: float = 2e-3, scale_init: float = 1.0,
              velscale: Optional[float] = None) -> dict:
    """Iterate (rescale red -> pPXF -> update scale) until the red/blue factor settles."""
    scale = scale_init
    res = None
    for it in range(n_iter if (fit_scale and len(arms) == 2) else 1):
        fs = build_fit_spectrum(spec, arms, fwhm, frame, windows_air, scale_r=scale, velscale=velscale)
        res = fit_at_z(fs, z_guess, st, extra_bad=extra_bad)
        if not (fit_scale and len(arms) == 2):
            break
        pp = res["pp"]
        sel = fs.is_red & np.isin(np.arange(fs.lam.size), res["goodpixels"])
        if sel.sum() < 20 or np.sum(pp.weights > 1e-8) == 0:
            break
        model_r = pp.bestfit[sel]
        data_r = res["flux"][sel]
        w = 1.0 / res["noise"][sel] ** 2
        new = scale * np.sum(w * data_r * model_r) / np.sum(w * data_r ** 2) if np.sum(w * data_r ** 2) > 0 else scale
        # `data_r` already contains `scale`; the least-squares update for the multiplier is
        # sum(w d m)/sum(w d^2) applied on top of the current scale.
        if not (np.isfinite(new) and new > 1e-4):
            break
        d = abs(new - scale) / scale
        scale = float(new)
        if d < tol:
            break
        velscale = fs.velscale
    res["scale_r"] = float(scale)
    return res


# ---------------------------------------------------------------------------
# redshift scan
# ---------------------------------------------------------------------------
def _highpass(x: np.ndarray, size: int) -> np.ndarray:
    from scipy.ndimage import uniform_filter1d
    return x - uniform_filter1d(x, size=size, axis=-1, mode="nearest")


def scan_continuum_z(spec: ExtractedSpectrum, arms, fwhm, frame, windows_air, z_grid: np.ndarray,
                     st: FitSettings, extra_bad: Sequence[Window] = (), progress: ProgressCB = None,
                     highpass_kms: float = 15000.0, stride: int = 2) -> dict:
    """Fast chi2(z) scan: high-pass-filtered non-negative linear template fit at each trial z.

    For every trial redshift the (broadened) templates are shifted onto the galaxy's log-lambda grid
    by interpolation, both galaxy and templates are high-pass filtered (removing continuum shape,
    reddening and any blue/red flux mismatch on scales > ``highpass_kms``), and the non-negative
    combination minimising chi2 is found by NNLS.  No polynomials, no non-linear optimisation:
    ~1 ms per trial instead of ~1 s, so the whole z range is scanned in seconds.  The lowest minima
    are then refined with full pPXF fits by the caller.  sigma is held fixed (``st.sigma_fixed``).
    """
    from scipy.optimize import nnls
    fs = build_fit_spectrum(spec, arms, fwhm, frame, windows_air)
    sigma = st.sigma_fixed if st.sigma_fixed is not None else 220.0
    lnl = np.log(fs.lam)
    n = lnl.size
    size = max(5, int(highpass_kms / fs.velscale) | 1)
    base_good = fs.good_base.copy()
    if extra_bad:
        ew = [(float(air_to_vac(a)), float(air_to_vac(b))) if frame == "vacuum" else (a, b) for a, b in extra_bad]
        base_good &= ~in_windows(fs.lam, ew)
    # fill masked/NaN pixels for the filter, then high-pass the galaxy once
    g = np.where(np.isfinite(fs.flux), fs.flux, 0.0)
    if (~base_good).any() and base_good.sum() > 10:
        g = np.interp(lnl, lnl[base_good], g[base_good])
    g_hp = _highpass(g, size)
    noise = np.where(fs.noise > 0, fs.noise, np.nanmedian(fs.noise))
    chi2 = np.full(z_grid.shape, np.nan)
    ngood = np.zeros(z_grid.shape, int)
    lib = st.library
    chunk_w = 0.05
    centers = np.arange(z_grid.min(), z_grid.max() + chunk_w, chunk_w)
    done = 0
    for zc in centers:
        in_chunk = np.flatnonzero((z_grid >= zc - chunk_w / 2) & (z_grid < zc + chunk_w / 2))
        if in_chunk.size == 0:
            continue
        zc_eff = float(0.5 * (z_grid[in_chunk[0]] + z_grid[in_chunk[-1]]))
        lam_rest = fs.lam / (1.0 + zc_eff)
        fwhm_rest = fs.fwhm / (1.0 + zc_eff)
        tem_fwhm = _STORE.native_fwhm(lib.path, lam_rest)
        cf = _STORE.cropped_file(lib.path, float(lam_rest.min() * 0.97), float(lam_rest.max() * 1.03),
                                 None if lib.is_empirical else st.age_range,
                                 None if lib.is_empirical else st.metal_range, 1 if lib.is_empirical else stride)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sps = sps_util.sps_lib(cf, fs.velscale, {"lam": lam_rest, "fwhm": np.maximum(fwhm_rest, tem_fwhm)})
        T = sps.templates.reshape(sps.templates.shape[0], -1)            # (m, nt)
        # fixed LOSVD: Gaussian broadening by sigma in log-lambda pixels
        from scipy.ndimage import gaussian_filter1d
        T = gaussian_filter1d(T, max(sigma / fs.velscale, 0.3), axis=0, mode="nearest")
        ln_t = sps.ln_lam_temp
        for i in in_chunk:
            z = float(z_grid[i])
            q = lnl - np.log(1.0 + z)
            if q[0] < ln_t[0] or q[-1] > ln_t[-1]:
                inside = (q >= ln_t[0]) & (q <= ln_t[-1])
            else:
                inside = np.ones(n, bool)
            Ts = np.array([np.interp(q, ln_t, T[:, k]) for k in range(T.shape[1])])   # (nt, n)
            Ts_hp = _highpass(Ts, size)
            good = base_good & inside
            neb = _nebular_windows_obs(z, st.mask_nebular, st.nebular_half_kms, frame, fs.lam[0], fs.lam[-1])
            good &= ~in_windows(fs.lam, neb)
            if st.rest_range is not None:
                good &= (fs.lam / (1.0 + z) >= st.rest_range[0]) & (fs.lam / (1.0 + z) <= st.rest_range[1])
            if good.sum() < 200:
                continue
            A = (Ts_hp[:, good] / noise[good]).T
            w, rn = nnls(A, g_hp[good] / noise[good])
            chi2[i] = rn ** 2 / good.sum()
            ngood[i] = int(good.sum())
        done += in_chunk.size
        if progress:
            progress(done / len(z_grid), f"continuum scan z={zc_eff:.3f}")
    return dict(z=z_grid, chi2=chi2, z_fit=z_grid.copy(), n_good=ngood, sigma_fixed=sigma)


def find_chi2_candidates(z, chi2, z_fit, min_sep_kms: float = 1500.0, n_max: int = 8) -> List[dict]:
    """Local minima of chi2(z), ranked.

    ``contrast`` = median(chi2 over the scan) / chi2 at the minimum: the factor by which this
    redshift beats a typical wrong one (1 = no better than chance).  Unlike a MAD-based sigma it
    stays meaningful when chi2(z) is skewed or bimodal, e.g. after restricting the rest window.
    """
    ok = np.isfinite(chi2)
    if ok.sum() < 3:
        return []
    med = float(np.nanmedian(chi2))
    order = np.argsort(np.where(ok, chi2, np.inf))
    out = []
    for k in order:
        if not ok[k]:
            continue
        zk = z[k]
        if any(abs(C_KMS * (zk - o["z"]) / (1 + zk)) < min_sep_kms for o in out):
            continue
        out.append(dict(z=float(z_fit[k] if np.isfinite(z_fit[k]) else zk), z_grid=float(zk), chi2=float(chi2[k]),
                        contrast=float(med / max(chi2[k], 1e-12))))
        if len(out) >= n_max:
            break
    return out


# ---------------------------------------------------------------------------
# uncertainties
# ---------------------------------------------------------------------------
def _local_scatter(resid: np.ndarray, good: np.ndarray, floor: float, window: int = 75) -> np.ndarray:
    """Running robust (MAD) scatter of the residuals, computed over GOOD pixels only (masked
    stretches would otherwise drive the local MAD to zero) and floored at ``floor``."""
    from scipy.ndimage import median_filter
    idx = np.flatnonzero(good)
    r = resid[idx]
    med = median_filter(r, size=min(window, max(5, r.size // 3) | 1), mode="nearest")
    mad = median_filter(np.abs(r - med), size=min(window, max(5, r.size // 3) | 1), mode="nearest")
    sc = np.interp(np.arange(resid.size), idx, 1.4826 * mad)
    return np.clip(sc, floor, None)


def bootstrap_fit(res: dict, n: int = 100, method: str = "wild", seed: Optional[int] = None,
                  progress: ProgressCB = None) -> dict:
    """Empirical errors by re-fitting perturbed copies of the best-fit model.

    ``wild``: model + eps_i * residual_i with eps = +/-1 (Rademacher); preserves heteroscedastic,
    non-Gaussian residual structure and template-mismatch, and the local residual scatter is used
    as the noise so each re-fit is properly weighted.  ``mc``: model + Gaussian(noise).
    """
    rng = np.random.default_rng(seed)
    pp = res["pp"]
    st: FitSettings = res["settings"]
    gp = res["goodpixels"]
    best = pp.bestfit.copy()
    gal = res["flux"]
    resid = gal - best
    noise = res["noise"].copy()
    if method == "wild":
        gmask = np.zeros(resid.size, bool)
        gmask[gp] = True
        noise_b = _local_scatter(resid, gmask, 0.3 * float(np.median(noise[gp])))
    else:
        noise_b = noise
    v0 = [res["v"], res["sigma"]] + [0.0] * (st.moments - 2)
    start = [[res["v"], res["sigma"]]] * 2 if st.has_gas else v0
    vs, ss = [], []
    for i in range(n):
        if method == "wild":
            eps = rng.choice([-1.0, 1.0], size=resid.size)
            g2 = best + eps * resid
        else:
            g2 = best + rng.normal(0.0, noise_b)
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                p2 = _ppxf(res["templates"], g2, noise_b, res["fs"].velscale, start=start,
                           moments=res["moments_arg"], degree=st.degree, mdegree=st.mdegree,
                           lam=res["lam_rest"], lam_temp=res.get("lam_temp", res["sps"].lam_temp),
                           component=res["component"] if st.has_gas else 0, goodpixels=gp, quiet=True,
                           plot=False, bounds=None if st.has_gas else
                           [[-st.vel_bound, st.vel_bound], [20.0, 600.0]] + [[-0.3, 0.3]] * (st.moments - 2))
            s = np.atleast_2d(p2.sol)[0]
            vs.append(float(s[0]))
            ss.append(float(s[1]))
        except Exception:
            pass
        if progress and (i % 5 == 0 or i == n - 1):
            progress((i + 1) / n, f"bootstrap {i + 1}/{n}")
    vs, ss = np.array(vs), np.array(ss)
    zs = (1 + res["z_guess"]) * np.exp(vs / C_KMS) - 1 if vs.size else vs
    ok = vs.size >= max(5, n // 2)
    return dict(v_samples=vs, sigma_samples=ss, z_samples=zs, n_success=int(vs.size), n_requested=n,
                z_err=float(np.std(zs, ddof=1)) if ok else float("nan"),
                v_err=float(np.std(vs, ddof=1)) if ok else float("nan"),
                sigma_err=float(np.std(ss, ddof=1)) if ok else float("nan"), ok=bool(ok), method=method)


def systematic_sweep(spec, arms, fwhm, windows_air, z_guess, st: FitSettings, frame_of: Callable[[Library], str],
                     degrees: Sequence[int] = (), libraries: Sequence[str] = (), fit_scale=True,
                     extra_bad: Sequence[Window] = (), scale_init: float = 1.0, custom=("", "air"),
                     progress: ProgressCB = None) -> dict:
    """Re-fit with different multiplicative-polynomial degrees and/or template libraries.

    Each library is fitted in ITS OWN native wavelength frame.  Returns every (label, V, sigma, z)
    triplet and the scatter, which the caller adds to the systematic error budget.
    """
    rows = []
    jobs = []
    for d in degrees:
        jobs.append((f"mdegree={d}", st.library, {"mdegree": int(d)}))
    for name in libraries:
        try:
            lib = resolve_library(name, *custom)
        except Exception as e:
            rows.append(dict(label=f"lib={name}", error=str(e)))
            continue
        extra = {"library": lib}
        if lib.is_empirical:
            extra.update(age_range=None, metal_range=None)
        jobs.append((f"lib={name}", lib, extra))
    for k, (label, lib, over) in enumerate(jobs):
        s2 = FitSettings(**{**st.__dict__, **over})
        try:
            r = fit_joint(spec, arms, fwhm, frame_of(s2.library), windows_air, z_guess, s2, fit_scale=fit_scale,
                          extra_bad=extra_bad, scale_init=scale_init)
            rows.append(dict(label=label, v=r["v"], sigma=r["sigma"], z=r["z"], chi2=r["chi2"]))
        except Exception as e:
            rows.append(dict(label=label, error=str(e)))
        if progress:
            progress((k + 1) / max(len(jobs), 1), f"sweep: {label}")
    good = [r for r in rows if "error" not in r]
    out = dict(rows=rows)
    if len(good) >= 2:
        # library fits carry their own V_sys zero-point (template/frame differences): remove the
        # per-library mean velocity before judging the scatter of SIGMA; V scatter is reported raw.
        sig = np.array([r["sigma"] for r in good]); z = np.array([r["z"] for r in good])
        out.update(sigma_sys=float(np.std(sig)), z_sys=float(np.std(z)), sigma_range=(float(sig.min()), float(sig.max())),
                   z_range=(float(z.min()), float(z.max())))
    return out
