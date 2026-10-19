"""Spectral extraction: weighted aperture sums and linear spectral deblending.

Both routines return, per arm, a flux and a **variance** spectrum in the same units as the cube
(after the optional surface-brightness -> flux conversion).  Variances include

* the per-spaxel noise (variance cube or sky-box scatter),
* the noise-correlation inflation measured from the sky box (:mod:`kcwilensfit.noise`),
* the uncertainty of the sky-background estimate itself (fully correlated across spaxels),
* and, for the deblend, the blending covariance -- two overlapping objects cannot both be
  measured as precisely as either alone, and the error bars say so.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter1d, median_filter

from .cubeio import Arm
from .noise import SkyModel, spaxel_variance


@dataclass
class ArmSpectrum:
    wave: np.ndarray
    flux: np.ndarray
    var: np.ndarray
    npix: float = 0.0              # aperture area [spaxels] (or effective area for weighted extractions)

    def copy(self) -> "ArmSpectrum":
        return ArmSpectrum(self.wave.copy(), self.flux.copy(), self.var.copy(), self.npix)

    @property
    def good(self) -> np.ndarray:
        return np.isfinite(self.flux) & np.isfinite(self.var) & (self.var > 0)

    def median_snr(self) -> float:
        g = self.good
        if g.sum() == 0:
            return float("nan")
        return float(np.nanmedian(self.flux[g] / np.sqrt(self.var[g])))

    def snr_per_angstrom(self, wmin: float = None, wmax: float = None) -> float:
        """Median S/N per Angstrom (pixel S/N scaled by sqrt(1/dispersion))."""
        g = self.good
        if wmin is not None:
            g &= self.wave >= wmin
        if wmax is not None:
            g &= self.wave <= wmax
        if g.sum() < 5:
            return float("nan")
        disp = float(np.median(np.diff(self.wave)))
        return float(np.nanmedian(self.flux[g] / np.sqrt(self.var[g])) / np.sqrt(disp))


@dataclass
class ExtractedSpectrum:
    id: str
    role: str
    group: str
    arms: Dict[str, ArmSpectrum]
    meta: dict = field(default_factory=dict)

    def copy(self) -> "ExtractedSpectrum":
        return ExtractedSpectrum(self.id, self.role, self.group, {k: v.copy() for k, v in self.arms.items()},
                                 dict(self.meta))


# ---------------------------------------------------------------------------
def _pixels(mask: np.ndarray, valid: np.ndarray, thresh: float = 0.01):
    sel = (mask > thresh) & valid
    ys, xs = np.nonzero(sel)
    return (ys, xs), mask[ys, xs].astype(np.float64)


def extract_sum(arm: Arm, sky: SkyModel, mask: np.ndarray, weight_image: Optional[np.ndarray] = None,
                subtract_background: bool = True, correct_correlation: bool = True) -> ArmSpectrum:
    """Weighted aperture sum.

    ``S(lam) = N * sum_i w_i (f_i - bg) / sum_i w_i``  with  ``N = sum_i m_i`` (the aperture area,
    ``m`` = fractional coverage) and ``w_i = m_i`` (uniform) or ``m_i * I_i`` (intensity weighted).
    Spaxels that are NaN at a given wavelength are dropped from that channel's weighted mean, and
    the mean is rescaled to the full aperture area ``N`` so the flux scale is wavelength-independent.
    """
    idx, m = _pixels(mask, arm.valid)
    if m.size == 0:
        raise ValueError("Aperture covers no valid spaxels -- check its position/size.")
    N = float(m.sum())
    w = m.copy()
    if weight_image is not None:
        wi = np.clip(weight_image[idx], 0, None)
        wi = np.where(np.isfinite(wi), wi, 0.0)
        if wi.sum() > 0:
            w = m * wi
    f = arm.flux[:, idx[0], idx[1]].astype(np.float64)
    if subtract_background:
        f = f - sky.bg[:, None]
    v = spaxel_variance(arm, sky, idx)
    ok = np.isfinite(f) & np.isfinite(v)
    W = np.where(ok, w[None, :], 0.0)
    wsum = W.sum(axis=1)
    safe = np.where(wsum > 0, wsum, np.nan)
    flux = N * np.nansum(np.where(ok, f, 0.0) * W, axis=1) / safe
    var = N ** 2 * np.nansum(np.where(ok, v, 0.0) * W ** 2, axis=1) / safe ** 2
    neff = float(w.sum() ** 2 / np.sum(w ** 2))
    if correct_correlation:
        var = var * sky.correlation_factor(neff) ** 2
    if subtract_background:
        var = var + N ** 2 * (sky.bg_var if sky.bg_var is not None else 0.0)
    flux = np.where(wsum > 0, flux, np.nan)
    var = np.where(wsum > 0, var, np.nan)
    return ArmSpectrum(arm.wave.copy(), flux, var, npix=N)


def deblend(arm: Arm, sky: SkyModel, templates: Sequence[np.ndarray], ap_masks: Sequence[np.ndarray],
            region: np.ndarray, subtract_background: bool = True, fit_background: bool = False,
            correct_correlation: bool = True) -> Tuple[List[ArmSpectrum], dict]:
    """Solve ``data(x, lam) = sum_k c_k(lam) P_k(x) [+ b(lam)]`` for all wavelengths.

    Parameters
    ----------
    templates : K spatial templates (full-grid, float, unit sum) -- the objects' light profiles
    ap_masks  : K fractional aperture masks; object k's reported spectrum is the part of its
                light that falls inside *its own* aperture, ``c_k * sum(P_k * mask_k)``
    region    : bool mask of the spaxels entering the fit (apertures + padding)

    Returns the K spectra and a diagnostics dict (residual image, template condition number, ...).
    """
    K = len(templates)
    idx = np.nonzero(region & arm.valid)
    npix = idx[0].size
    if npix < K + 3:
        raise ValueError("Deblend region contains too few valid spaxels for the number of components.")
    A = np.array([t[idx] for t in templates]).T.astype(np.float64)          # (npix, K)
    ncol = K + (1 if fit_background else 0)
    if fit_background:
        A = np.hstack([A, np.full((npix, 1), 1.0 / npix)])
    D = arm.flux[:, idx[0], idx[1]].astype(np.float64).T                     # (npix, nw)
    if subtract_background:
        D = D - sky.bg[None, :]
    V = spaxel_variance(arm, sky, idx).T                                     # (npix, nw)
    nw = D.shape[1]
    ok = np.isfinite(D) & np.isfinite(V) & (V > 0)
    uniform = sky.source != "variance_cube"
    coef = np.full((ncol, nw), np.nan)
    var_c = np.full((ncol, nw), np.nan)       # variance of coefficients from per-spaxel noise
    gsum = np.full((ncol, nw), np.nan)        # sensitivity of each coefficient to a flat bg error
    cond_max = 0.0

    if uniform:
        # weights equal across spaxels at each wavelength -> solve once per NaN pattern
        pats = np.packbits(ok, axis=0)
        _, inv = np.unique(pats.T, axis=0, return_inverse=True)
        inv = np.asarray(inv).ravel()
        for g in np.unique(inv):
            lam_idx = np.flatnonzero(inv == g)
            rows = ok[:, lam_idx[0]]
            if rows.sum() < ncol + 2:
                continue
            Ap = A[rows]
            ATA = Ap.T @ Ap
            try:
                C = np.linalg.inv(ATA)
            except np.linalg.LinAlgError:
                C = np.linalg.pinv(ATA)
            cond_max = max(cond_max, float(np.linalg.cond(ATA)))
            Pinv = C @ Ap.T                                                  # (ncol, nrows)
            coef[:, lam_idx] = Pinv @ D[rows][:, lam_idx]
            sig2 = sky.noise[lam_idx] ** 2
            var_c[:, lam_idx] = np.diag(C)[:, None] * sig2[None, :]
            gsum[:, lam_idx] = Pinv.sum(axis=1)[:, None]
    else:
        for j in range(nw):
            rows = ok[:, j]
            if rows.sum() < ncol + 2:
                continue
            Ap = A[rows]
            Wj = 1.0 / V[rows, j]
            ATA = (Ap * Wj[:, None]).T @ Ap
            try:
                C = np.linalg.inv(ATA)
            except np.linalg.LinAlgError:
                C = np.linalg.pinv(ATA)
            Pinv = C @ (Ap * Wj[:, None]).T
            coef[:, j] = Pinv @ D[rows, j]
            var_c[:, j] = np.diag(C)
            gsum[:, j] = Pinv.sum(axis=1)
            if j == 0:
                cond_max = float(np.linalg.cond(ATA))

    out: List[ArmSpectrum] = []
    for k in range(K):
        nk = float(np.sum(templates[k] * ap_masks[k]))
        neff = float(1.0 / np.sum(templates[k][idx] ** 2))
        fac = sky.correlation_factor(neff) if correct_correlation else 1.0
        flux = coef[k] * nk
        var = var_c[k] * fac ** 2
        if subtract_background and sky.bg_var is not None:
            var = var + (gsum[k] ** 2) * sky.bg_var
        var = var * nk ** 2
        out.append(ArmSpectrum(arm.wave.copy(), flux, var, npix=float(ap_masks[k].sum())))

    # residual image: median over wavelength of the unexplained flux, for the QA plot
    model = A @ np.where(np.isfinite(coef), coef, 0.0)
    resid = np.where(ok, D - model, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        rimg = np.nanmedian(resid, axis=1)
    resid_img = np.full(region.shape, np.nan)
    resid_img[idx] = rimg
    diag = dict(cond=cond_max, resid_img=resid_img, n_pix=int(npix))
    return out, diag


# ---------------------------------------------------------------------------
def inverse_variance_combine(specs: Sequence[ArmSpectrum]) -> ArmSpectrum:
    """Pixel-by-pixel inverse-variance combination of spectra on a common wavelength grid.

    Valid for images of one lensed source: lensing is achromatic, so every image has the same
    redshift.  (Source variability + time delays change the *shape* between images but not z.)
    """
    fl = np.array([s.flux for s in specs])
    va = np.array([s.var for s in specs])
    w = np.where(np.isfinite(va) & (va > 0) & np.isfinite(fl), 1.0 / va, 0.0)
    ws = w.sum(axis=0)
    safe = np.where(ws > 0, ws, 1.0)
    flux = np.where(ws > 0, (w * np.where(np.isfinite(fl), fl, 0)).sum(axis=0) / safe, np.nan)
    var = np.where(ws > 0, 1.0 / safe, np.nan)
    return ArmSpectrum(specs[0].wave.copy(), flux, var, npix=float(sum(s.npix for s in specs)))


def reject_spikes(spec: ArmSpectrum, nsigma: float = 7.0, half: int = 4) -> int:
    """Mask isolated single-pixel positive spikes (cosmic rays / hot pixels) in place.

    A spike must exceed ``nsigma`` above a running median of its neighbours *and* both neighbours
    must be unremarkable (< 3 sigma): a genuine emission line is resolved by the instrument line
    spread (FWHM >~ 2 pixels), so it cannot be a single bright pixel.  Returns the number masked.
    """
    f, v = spec.flux, spec.var
    g = spec.good
    if g.sum() < 4 * half:
        return 0
    fi = np.where(g, f, np.nanmedian(f[g]))
    med = median_filter(fi, size=2 * half + 1, mode="nearest")
    sig = np.sqrt(np.where(g, v, np.nanmedian(v[g])))
    r = (fi - med) / sig
    spike = (r > nsigma) & g
    spike[1:-1] &= (r[:-2] < 3.0) & (r[2:] < 3.0)
    spike[0] = spike[-1] = False
    n = int(spike.sum())
    if n:
        spec.flux[spike] = np.nan
        spec.var[spike] = np.nan
    return n


def smooth_for_display(spec: ArmSpectrum, sigma_A: float) -> ArmSpectrum:
    """Gaussian-smoothed COPY for plotting only (variance propagated with the squared kernel)."""
    if sigma_A <= 0:
        return spec
    disp = float(np.median(np.diff(spec.wave)))
    s = max(sigma_A / disp, 0.3)
    f = np.where(np.isfinite(spec.flux), spec.flux, 0.0)
    wt = np.isfinite(spec.flux).astype(float)
    num = gaussian_filter1d(f, s, mode="nearest")
    den = gaussian_filter1d(wt, s, mode="nearest")
    flux = np.where(den > 0.2, num / np.where(den > 0, den, 1), np.nan)
    # variance of a normalised Gaussian smooth: sum k^2 var
    k = np.exp(-0.5 * (np.arange(-int(4 * s) - 1, int(4 * s) + 2) / s) ** 2)
    k /= k.sum()
    vv = np.where(np.isfinite(spec.var), spec.var, 0.0)
    var = np.convolve(vv, k ** 2, mode="same")
    return ArmSpectrum(spec.wave, flux, np.where(np.isfinite(flux), var, np.nan), spec.npix)


def parse_subtraction_rules(text: str):
    """'target: contaminator [scale|auto]' lines -> [(target, contaminator, scale or None)].

    ``None`` means 'estimate the scale automatically'.  Blank lines and '#' comments are ignored.
    """
    rules = []
    for n, raw in enumerate(str(text or "").splitlines(), 1):
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if ":" not in line:
            raise ValueError(f"manual subtraction line {n} ('{raw.strip()}'): expected 'target: contaminator [scale]'")
        tgt, rest = [t.strip() for t in line.split(":", 1)]
        parts = rest.replace(",", " ").split()
        if not tgt or not parts:
            raise ValueError(f"manual subtraction line {n}: need a target and a contaminator")
        con = parts[0]
        scale = None
        if len(parts) > 1 and parts[1].lower() != "auto":
            try:
                scale = float(parts[1])
            except ValueError:
                raise ValueError(f"manual subtraction line {n}: scale must be a number or 'auto', got '{parts[1]}'")
        if tgt == con:
            raise ValueError(f"manual subtraction line {n}: a region cannot be subtracted from itself")
        rules.append((tgt, con, scale))
    return rules


def subtract_spectrum(target: "ArmSpectrum", contaminator: "ArmSpectrum", scale: float) -> None:
    """target -= scale * contaminator, in place; the contaminator's variance is propagated."""
    target.flux = target.flux - scale * contaminator.flux
    target.var = target.var + (scale ** 2) * contaminator.var
