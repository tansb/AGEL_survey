"""Blind redshift search and line-based redshift measurement.

Method (matched filtering, i.e. maximum likelihood with a Gaussian line model)
------------------------------------------------------------------------------
The spectrum is continuum-subtracted with a running median and its noise re-scaled to the *local*
empirical scatter (so a mis-estimated variance cannot manufacture detections).  At each trial
redshift z and for every catalogue line (a *multiplet* with nominal component ratios,
:mod:`kcwilensfit.lines`) the best-fit amplitude A and its error follow from linear least squares
against the template T(lam):

    A = sum(w r T)/sum(w T^2),     sigma_A = (sum w T^2)^(-1/2),     s = A / sigma_A

Emission lines need ``s > 0`` and absorption lines ``s < 0``.  Under the null hypothesis each
covered line contributes ``max(s_signed, 0)^2`` whose mean is 1/2, so the redshift score is

    S(z) = sum over covered lines [ max(s_signed, 0)^2 - 1/2 ]

-- the null expectation is subtracted line by line so redshifts at which *more* lines happen to
fall inside the wavelength coverage are not favoured over redshifts with fewer (a bias of the naive
"sum of S/N" score, and the reason a single-line detection at high z can lose to a spurious
low-z multi-line alias).  S(z) is therefore roughly a (positive-part) chi-square improvement over
"no lines at all".

After a candidate is chosen the lines are fitted *jointly* with a common redshift and line width,
and each line is also fitted alone, so inconsistent identifications show up as scatter.  Lines of
different physical origin (nebular / stellar vs interstellar-outflow absorption vs resonant
emission) are allowed their own velocity offsets relative to the systemic redshift.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d
from scipy.optimize import least_squares, nnls
from scipy.signal import find_peaks

from .extraction import ExtractedSpectrum
from .lines import (ABSORPTION_KINDS, EMISSION_KINDS, KIND_OFFSET_KMS, LINES, LINES_BY_NAME, Line,
                    default_line_names, observed_wavelengths_air)
from .masks import in_windows
from .wavelength import C_KMS, vac_to_air

ProgressCB = Optional[Callable[[float, str], None]]
Window = Tuple[float, float]
FWHM2SIG = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


@dataclass
class LineSpectrum:
    wave: np.ndarray        # air, ascending
    flux: np.ndarray
    cont: np.ndarray
    resid: np.ndarray       # flux - cont
    noise: np.ndarray       # empirically rescaled 1-sigma
    good: np.ndarray        # bool
    fwhm: np.ndarray        # instrumental FWHM per pixel [A]
    cutover: float
    noise_scale: float      # median(empirical/nominal noise) -- diagnostic

    @property
    def n(self):
        return self.wave.size


def prepare_line_spectrum(spec: ExtractedSpectrum, arms: Sequence[str], fwhm: Dict[str, float],
                          windows_air: Sequence[Window], cont_window_A: float = 150.0) -> LineSpectrum:
    """Splice the arms at one wavelength (centre of their overlap), subtract the continuum, rescale noise."""
    parts = []
    for arm in arms:
        a = spec.arms[arm]
        parts.append((arm, a.wave, a.flux, a.var))
    parts.sort(key=lambda p: p[1].mean())
    if len(parts) == 2:
        (_, wb, fb, vb), (_, wr, fr, vr) = parts
        mid = 0.5 * (wr.min() + wb.max())
        kb, kr = wb < mid, wr >= mid
        wave = np.concatenate([wb[kb], wr[kr]])
        flux = np.concatenate([fb[kb], fr[kr]])
        var = np.concatenate([vb[kb], vr[kr]])
        fw = np.concatenate([np.full(kb.sum(), fwhm[parts[0][0]]), np.full(kr.sum(), fwhm[parts[1][0]])])
    else:
        arm, wave, flux, var = parts[0]
        mid = wave.max() + 1 if arm == "blue" else wave.min() - 1
        fw = np.full(wave.shape, fwhm[arm])
    good = np.isfinite(flux) & np.isfinite(var) & (var > 0) & ~in_windows(wave, windows_air)
    dl = float(np.median(np.diff(wave)))
    # continuum: median filter over good pixels (masked pixels filled by interpolation)
    f_fill = np.interp(wave, wave[good], flux[good]) if good.sum() > 10 else np.nan_to_num(flux)
    w_pix = max(5, int(cont_window_A / dl) | 1)
    cont = median_filter(f_fill, size=w_pix, mode="nearest")
    cont = uniform_filter1d(cont, size=max(3, w_pix // 4), mode="nearest")
    resid = np.where(good, flux - cont, 0.0)
    sig = np.sqrt(np.where(good, var, np.nanmedian(var[good]) if good.any() else 1.0))
    # local empirical scatter of the residuals relative to the nominal noise (robust -> immune to lines)
    ratio = np.where(good, np.abs(resid) / sig, np.nan)
    wl = max(11, int(300.0 / dl) | 1)
    r_fill = np.where(np.isfinite(ratio), ratio, np.nanmedian(ratio))
    s_loc = median_filter(r_fill, size=wl, mode="nearest") / 0.6745
    s_loc = np.clip(s_loc, 0.7, 6.0)
    noise = sig * s_loc
    return LineSpectrum(wave=wave, flux=flux, cont=cont, resid=resid, noise=noise, good=good, fwhm=fw,
                        cutover=float(mid), noise_scale=float(np.nanmedian(s_loc[good])) if good.any() else 1.0)


# ---------------------------------------------------------------------------
def _comp_arrays(line: Line, z: np.ndarray):
    """Observed air wavelengths (nz, ncomp) and relative strengths (ncomp,)."""
    lam_vac = np.array([c[0] for c in line.comps])
    rel = np.array([c[1] for c in line.comps], float)
    mu = vac_to_air(lam_vac[None, :] * (1.0 + np.asarray(z)[:, None]))
    return mu, rel


def _sig_obs(ls: LineSpectrum, mu: np.ndarray, sigma_v: float) -> np.ndarray:
    fw = np.interp(mu, ls.wave, ls.fwhm)
    return np.sqrt((fw * FWHM2SIG) ** 2 + (mu * sigma_v / C_KMS) ** 2)


def line_snr_grid(ls: LineSpectrum, line: Line, z: np.ndarray, sigma_v: float = 100.0
                  ) -> Tuple[np.ndarray, np.ndarray]:
    """Signed S/N ``s(z)`` of one multiplet's best-fit amplitude, and a coverage flag, vectorised in z."""
    z = np.atleast_1d(np.asarray(z, float))
    mu, rel = _comp_arrays(line, z)                         # (nz, nc)
    sg = _sig_obs(ls, mu, sigma_v)
    dl = float(np.median(np.diff(ls.wave)))
    smax = float(sg.max())
    m = int(np.ceil(5.0 * smax / dl)) + 1
    lo_px = np.interp(mu.min(axis=1), ls.wave, np.arange(ls.n))
    hi_px = np.interp(mu.max(axis=1), ls.wave, np.arange(ls.n))
    span = int(np.ceil(np.max(hi_px - lo_px))) if len(z) else 0
    P = span + 2 * m + 1
    start = np.floor(lo_px).astype(int) - m
    idx = start[:, None] + np.arange(P)[None, :]
    inb = (idx >= 0) & (idx < ls.n)
    idc = np.clip(idx, 0, ls.n - 1)
    lam_p = ls.wave[idc]
    wgt = np.where(inb & ls.good[idc], 1.0 / ls.noise[idc] ** 2, 0.0)
    r = ls.resid[idc]
    T = np.zeros_like(lam_p)
    for j in range(mu.shape[1]):
        T += rel[j] * np.exp(-0.5 * ((lam_p - mu[:, j:j + 1]) / sg[:, j:j + 1]) ** 2)
    # ignore components outside the data range entirely
    in_range = (mu.min(axis=1) > ls.wave[0]) & (mu.max(axis=1) < ls.wave[-1])
    N = np.sum(wgt * r * T, axis=1)
    D = np.sum(wgt * T * T, axis=1)
    s = np.where(D > 1e-30, N / np.sqrt(np.where(D > 1e-30, D, 1.0)), 0.0)
    # coverage: enough unmasked weight under the template core
    w_med = float(np.median(1.0 / ls.noise[ls.good] ** 2)) if ls.good.any() else 1.0
    ideal = np.sum((T > 0.1) * w_med, axis=1)                   # weight if nothing were masked
    cover = in_range & (D > 0) & (D > 0.15 * ideal)
    return s, cover


def scan_lines(ls: LineSpectrum, line_names: Sequence[str], z_grid: np.ndarray, sigma_v: float = 100.0,
               detect_snr: float = 3.0, progress: ProgressCB = None) -> dict:
    """Redshift score S(z), number of detected lines, and per-line signed S/N (nz x nlines)."""
    lines = [LINES_BY_NAME[n] for n in line_names if n in LINES_BY_NAME]
    nz = len(z_grid)
    S = np.zeros(nz)
    ndet = np.zeros(nz, int)
    snr = np.zeros((nz, len(lines)))
    cov = np.zeros((nz, len(lines)), bool)
    chunk = 2000
    for k, line in enumerate(lines):
        for a in range(0, nz, chunk):
            s, c = line_snr_grid(ls, line, z_grid[a:a + chunk], sigma_v)
            snr[a:a + chunk, k] = s
            cov[a:a + chunk, k] = c
        if progress:
            progress((k + 1) / len(lines), f"line scan: {line.name}")
    for k, line in enumerate(lines):
        sgn = snr[:, k] if line.is_emission else -snr[:, k]
        contrib = np.where(cov[:, k], np.clip(sgn, 0, None) ** 2 - 0.5, 0.0)
        S += contrib
        ndet += (cov[:, k] & (sgn >= detect_snr)).astype(int)
    return dict(z=z_grid, score=S, n_det=ndet, snr=snr, covered=cov, lines=[l.name for l in lines])


def null_score_threshold(ls: LineSpectrum, line_names: Sequence[str], z_grid: np.ndarray, sigma_v: float,
                         n_sim: int = 20, quantile: float = 95.0, block: int = 7, seed: int = 1,
                         progress: ProgressCB = None) -> float:
    """Look-elsewhere-corrected significance threshold for the maximum of S(z).

    Scanning ~10^4 trial redshifts x ~45 lines guarantees chance 3-sigma 'detections' in pure noise.
    The null distribution of max_z S(z) is estimated by re-running the *whole* scan on copies of the
    residual spectrum whose signs are randomly flipped in blocks of ``block`` pixels (destroying any
    line, preserving the pixel noise and its local correlation).  A candidate is only called
    significant if its score exceeds the ``quantile`` percentile of these maxima (false-alarm rate
    ~ 1 - quantile/100 over the whole search).
    """
    rng = np.random.default_rng(seed)
    n = ls.n
    maxima = []
    for k in range(n_sim):
        signs = np.repeat(rng.choice([-1.0, 1.0], size=n // block + 1), block)[:n]
        sim = LineSpectrum(wave=ls.wave, flux=ls.flux, cont=ls.cont, resid=ls.resid * signs, noise=ls.noise,
                           good=ls.good, fwhm=ls.fwhm, cutover=ls.cutover, noise_scale=ls.noise_scale)
        sc = scan_lines(sim, line_names, z_grid, sigma_v)
        maxima.append(float(sc["score"].max()))
        if progress:
            progress((k + 1) / n_sim, f"null simulations {k + 1}/{n_sim}")
    return float(np.percentile(maxima, quantile))


def find_line_candidates(scan: dict, min_sep_kms: float = 1500.0, n_max: int = 10, ls_signed=None) -> List[dict]:
    z, S = scan["z"], scan["score"]
    if len(z) < 3:
        return []
    step = float(np.median(np.diff(z)))
    dist = max(1, int(min_sep_kms / C_KMS * (1.0 + np.median(z)) / step))
    pk, _ = find_peaks(S, distance=dist)
    order = pk[np.argsort(S[pk])[::-1]][:n_max]
    out = []
    names = scan["lines"]
    for k in order:
        det = []
        for j, n in enumerate(names):
            ln = LINES_BY_NAME[n]
            s = scan["snr"][k, j] if ln.is_emission else -scan["snr"][k, j]
            if scan["covered"][k, j] and s >= 2.0:
                det.append((n, float(s)))
        det.sort(key=lambda t: -t[1])
        out.append(dict(z=float(z[k]), score=float(S[k]), n_det=int(scan["n_det"][k]), lines=det))
    return out


# ---------------------------------------------------------------------------
# joint refinement
# ---------------------------------------------------------------------------
def _gauss_template(wave, line: Line, z: float, sigma_v: float, ls: LineSpectrum):
    mu = vac_to_air(np.array([c[0] for c in line.comps]) * (1.0 + z))
    rel = np.array([c[1] for c in line.comps], float)
    sg = _sig_obs(ls, mu[None, :], sigma_v)[0]
    T = np.zeros_like(wave)
    for j in range(len(mu)):
        T += rel[j] * np.exp(-0.5 * ((wave - mu[j]) / sg[j]) ** 2)
    return T, mu, sg


def refine_lines(ls: LineSpectrum, z0: float, line_names: Sequence[str], sigma_v0: float = 100.0,
                 select_snr: float = 2.0, sigma_range=(30.0, 500.0), z_halfwidth_kms: float = 1200.0,
                 lines_force: Optional[Sequence[str]] = None) -> Optional[dict]:
    """Joint Gaussian fit of the lines significant at ``z0`` with common z and line width."""
    cand = [LINES_BY_NAME[n] for n in line_names if n in LINES_BY_NAME]
    sel_all = []
    zz = np.array([z0])
    for ln in cand:
        s, c = line_snr_grid(ls, ln, zz, sigma_v0)
        sgn = float(s[0]) if ln.is_emission else -float(s[0])
        if c[0] and sgn >= select_snr:
            sel_all.append((ln, sgn))
    if lines_force:
        for n in lines_force:
            if n in LINES_BY_NAME and all(n != l.name for l, _ in sel_all):
                sel_all.append((LINES_BY_NAME[n], 0.0))
    if not sel_all:
        return None
    # drop weaker duplicates of the same observed feature (e.g. an emission/absorption pair): keep sign-consistent best
    sel_all.sort(key=lambda t: -t[1])
    lines = [l for l, _ in sel_all]
    kinds = [l.kind for l in lines]
    ref_kinds = ("neb", "stel_abs")
    have_ref = any(k in ref_kinds for k in kinds)
    if not have_ref:
        # no systemic tracer: the strongest line's class defines z, other classes get offsets
        ref_kinds = (kinds[0],)
    classes = sorted(set(k for k in kinds if k not in ref_kinds))
    # per-line pixel windows (+- 8 sigma around each component) and local baselines
    wave = ls.wave
    pix_sets = []
    for ln in lines:
        T, mu, sg = _gauss_template(wave, ln, z0, sigma_v0, ls)
        win = np.zeros(wave.size, bool)
        for m_, s_ in zip(mu, sg):
            win |= np.abs(wave - m_) < 8.0 * max(s_, 0.8)
        pix_sets.append(win & ls.good)
    allpix = np.zeros(wave.size, bool)
    for p in pix_sets:
        allpix |= p
    if allpix.sum() < 8:
        return None
    xs = wave[allpix]
    r = ls.resid[allpix]
    w = 1.0 / ls.noise[allpix]
    # nonlinear parameters: dz-velocity about z0, sigma_v, class offsets
    names = ["v_sys", "sigma_v"] + [f"dv_{c}" for c in classes]
    p0 = [0.0, float(np.clip(sigma_v0, *sigma_range))] + [KIND_OFFSET_KMS[c][0] for c in classes]
    lo = [-z_halfwidth_kms, sigma_range[0]] + [-900.0] * len(classes)
    hi = [z_halfwidth_kms, sigma_range[1]] + [900.0] * len(classes)

    pix_index = np.flatnonzero(allpix)
    loc = [np.isin(pix_index, np.flatnonzero(p)) for p in pix_sets]    # which fitted pixels each line owns

    def build(theta):
        v_sys, sig = theta[0], theta[1]
        off = {c: theta[2 + i] for i, c in enumerate(classes)}
        cols, sgn = [], []
        for ln, mask in zip(lines, loc):
            v = v_sys + (off.get(ln.kind, 0.0) if ln.kind not in ref_kinds else 0.0)
            zl = (1.0 + z0) * np.exp(v / C_KMS) - 1.0
            T, _, _ = _gauss_template(xs, ln, zl, sig, ls)
            col = np.where(mask, T, 0.0)
            cols.append(col)
            sgn.append(1.0 if ln.is_emission else -1.0)
        # local baselines (constant + slope) per line window, free sign
        base = []
        for mask in loc:
            x0 = xs[mask].mean() if mask.any() else 0.0
            for deg in (0, 1):
                b = np.where(mask, ((xs - x0) / 20.0) ** deg, 0.0)
                base += [b, -b]
        return np.array(cols).T, np.array(sgn), np.array(base).T

    def solve(theta):
        A, sg_, B = build(theta)
        M = np.hstack([A * sg_[None, :], B]) * w[:, None]
        coef, _ = nnls(M, r * w)
        return M, coef

    def resid_fn(theta):
        M, coef = solve(theta)
        return M @ coef - r * w

    res = least_squares(resid_fn, p0, bounds=(lo, hi), x_scale=[100.0, 50.0] + [100.0] * len(classes))
    theta = res.x
    M, coef = solve(theta)
    nlin = len(lines)
    amps = coef[:nlin]
    # uncertainties from the projected-residual Jacobian, scaled by reduced chi2 if > 1
    J = res.jac
    dof = max(int(allpix.sum()) - len(theta) - nlin - 4 * len(lines), 1)
    chi2 = float(np.sum(res.fun ** 2))
    chi2r = chi2 / dof
    try:
        cov = np.linalg.inv(J.T @ J) * max(chi2r, 1.0)
        perr = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        perr = np.full(len(theta), np.nan)
    z_fit = float((1.0 + z0) * np.exp(theta[0] / C_KMS) - 1.0)
    z_err = float((1.0 + z_fit) * perr[0] / C_KMS)
    # per-line S/N (amplitude / error) from the linear solve at the optimum, and single-line redshifts
    ATA = M.T @ M
    try:
        Cl = np.linalg.pinv(ATA) * max(chi2r, 1.0)
        a_err = np.sqrt(np.clip(np.diag(Cl)[:nlin], 0, None))
    except Exception:
        a_err = np.full(nlin, np.nan)
    rows = []
    for k, ln in enumerate(lines):
        snr_k = amps[k] / a_err[k] if a_err[k] > 0 else np.nan
        zj, zj_err = _single_line_z(ls, ln, z_fit, float(theta[1]), xs, pix_index, loc[k])
        rows.append(dict(name=ln.name, kind=ln.kind, emission=ln.is_emission, amp=float(amps[k]),
                         amp_err=float(a_err[k]), snr=float(snr_k), z_line=zj, z_line_err=zj_err,
                         lam_obs=float(observed_wavelengths_air(ln, z_fit)[int(np.argmax([c[1] for c in ln.comps]))])))
    detected = [r_ for r_ in rows if np.isfinite(r_["snr"]) and r_["snr"] >= 3.0]
    # model for plotting
    A, sg_, B = build(theta)
    model = (A * sg_[None, :]) @ amps
    model_full = np.zeros(wave.size)
    model_full[allpix] = model
    # consistency is judged within the reference class only (classes may legitimately differ by outflow offsets)
    ok_rows = [r_ for r_ in rows if np.isfinite(r_["z_line"]) and np.isfinite(r_["z_line_err"])
               and r_["z_line_err"] > 0 and r_["snr"] >= 2.0 and r_["kind"] in ref_kinds]
    zl = np.array([r_["z_line"] for r_ in ok_rows])
    ze = np.array([r_["z_line_err"] for r_ in ok_rows])
    chi2_lines = float(np.sum(((zl - np.sum(zl / ze ** 2) / np.sum(1 / ze ** 2)) / ze) ** 2)) if len(zl) >= 2 else np.nan
    ref_lines = [r_ for r_ in rows if r_["kind"] in ("neb", "stel_abs") and r_["snr"] >= 3.0]
    out = dict(z=z_fit, z_err=z_err, sigma_v=float(theta[1]), sigma_v_err=float(perr[1]),
               class_offsets={c: (float(theta[2 + i]), float(perr[2 + i])) for i, c in enumerate(classes)},
               lines=rows, n_det=len(detected), n_ref_det=len(ref_lines), chi2r=chi2r,
               z_consistency_chi2=chi2_lines, n_consistency=int(len(zl)),
               model=model_full, model_pix=allpix, have_ref=have_ref,
               kinds_used=sorted(set(r_["kind"] for r_ in detected)))
    # systematic: lines that do not trace the systemic velocity directly
    if detected:
        kinds_det = set(r_["kind"] for r_ in detected)
        sys_ref = kinds_det & {"neb", "stel_abs"}
        if sys_ref:
            sys_kms = 30.0
        else:
            rk = ref_kinds[0]
            sys_kms = KIND_OFFSET_KMS[rk][1]
            off = KIND_OFFSET_KMS[rk][0]
            out["typical_offset_kms"] = off
            out["z_sys_estimate"] = float((1 + z_fit) * np.exp(-off / C_KMS) - 1.0)
            out["reference_kind"] = rk
        out["sys_err_kms"] = float(sys_kms)
        out["z_sys_err"] = float((1 + z_fit) * sys_kms / C_KMS)
    return out


def _single_line_z(ls, ln, z_ref, sigma_v, xs, pix_index, own_mask) -> Tuple[float, float]:
    """Redshift from this multiplet alone (amplitude + local baseline linear, z non-linear)."""
    m = own_mask
    if m.sum() < 6:
        return float("nan"), float("nan")
    x, r, w = xs[m], ls.resid[pix_index[m]], 1.0 / ls.noise[pix_index[m]]
    x0 = x.mean()
    sign = 1.0 if ln.is_emission else -1.0

    def res(p):
        zl = (1.0 + z_ref) * np.exp(p[0] / C_KMS) - 1.0
        T, _, _ = _gauss_template(x, ln, zl, sigma_v, ls)
        M = np.column_stack([sign * T, np.ones_like(x), (x - x0) / 20.0, -np.ones_like(x), -(x - x0) / 20.0]) * w[:, None]
        coef, _ = nnls(M, r * w)
        return M @ coef - r * w
    try:
        sol = least_squares(res, [0.0], bounds=([-600.0], [600.0]), x_scale=[50.0])
        J = sol.jac
        dof = max(m.sum() - 4, 1)
        chi2r = float(np.sum(sol.fun ** 2)) / dof
        err = np.sqrt(max(chi2r, 1.0) / max(float(J[:, 0] @ J[:, 0]), 1e-30))
        zl = (1.0 + z_ref) * np.exp(sol.x[0] / C_KMS) - 1.0
        return float(zl), float((1 + zl) * err / C_KMS)
    except Exception:
        return float("nan"), float("nan")


def refine_peak(ls: LineSpectrum, line_names: Sequence[str], z_center: float, half_kms: float = 600.0,
                n: int = 41, sigma_v: float = 100.0) -> float:
    """Maximise S(z) on a fine grid around a coarse peak."""
    dz = (1 + z_center) * half_kms / C_KMS
    zg = np.linspace(z_center - dz, z_center + dz, n)
    sc = scan_lines(ls, line_names, zg, sigma_v)
    return float(zg[int(np.argmax(sc["score"]))])


def identification_table(ls: LineSpectrum, lam_feature: float, sign: str = "emission", sigma_v: float = 100.0,
                         min_cover: int = 1) -> List[dict]:
    """For a SINGLE isolated feature, list every catalogue line that could produce it.

    For each identification the redshift implies where *other* lines must appear; the S/N measured
    there (and how many are covered but missing) separates the plausible from the ruled-out options.
    """
    out = []
    cands = [l for l in LINES if (l.is_emission == (sign == "emission"))]
    strong_other = [l for l in LINES if l.is_emission == (sign == "emission")]
    for ln in cands:
        # redshift putting the strongest component at lam_feature
        k = int(np.argmax([c[1] for c in ln.comps]))
        from .wavelength import air_to_vac
        z = float(air_to_vac(lam_feature) / ln.comps[k][0] - 1.0)
        if z < 0:
            continue
        sn, cv = [], []
        for o in strong_other:
            if o.name == ln.name:
                continue
            s, c = line_snr_grid(ls, o, np.array([z]), sigma_v)
            if c[0]:
                sn.append((o.name, float(s[0] if o.is_emission else -s[0])))
        missing = [n for n, s in sn if s < 1.0]
        present = [(n, s) for n, s in sn if s >= 3.0]
        out.append(dict(line=ln.name, z=z, n_other_covered=len(sn), n_other_detected=len(present),
                        others_detected=present, others_missing=missing[:6]))
    out.sort(key=lambda d: (-d["n_other_detected"], d["n_other_covered"]))
    return out
