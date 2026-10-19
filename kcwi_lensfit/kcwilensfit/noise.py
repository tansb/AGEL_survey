"""Sky background and noise model.

The DRP cubes carry (at best) a variance cube that ignores (i) the covariance created when the
slicer data are resampled onto a square grid and (ii) any residual sky-subtraction error.  This
module therefore measures both **empirically** from a source-free box:

* per-wavelength *background level* = median over the box (robust to a stray faint star);
* per-wavelength *noise* = 1.4826 * MAD over the box (a variance cube, if supplied, is used for
  the spatial/spectral shape and re-normalised by the measured scatter);
* the **noise-correlation exponent** alpha from block sums of 1x1, 2x2, 3x3 ... spaxels:
  ``sigma_N = sigma_1 * N**alpha`` (alpha = 0.5 for independent pixels, -> 1 for fully
  correlated).  Aperture errors are inflated by ``N**(alpha-0.5)`` -- without this, errors on a
  summed spectrum are underestimated, sometimes by a factor of 2-3.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from .cubeio import Arm


@dataclass
class SkyModel:
    box: Tuple[int, int, int, int]        # (y0, y1, x0, x1), end-exclusive
    bg: np.ndarray                        # (nw,) per-spaxel background level to subtract
    noise: np.ndarray                     # (nw,) per-spaxel 1-sigma noise
    n_pix: int                            # spaxels in the box
    alpha: float = 0.5                    # noise-correlation exponent
    ratio1: float = 1.0                   # measured/assumed noise at N=1 (variance-cube calibration)
    source: str = "sky_box"               # 'sky_box' | 'variance_cube'
    bg_var: Optional[np.ndarray] = None   # (nw,) variance of the background estimate itself
    scaling: Optional[dict] = None        # diagnostic data for plotting

    def correlation_factor(self, npix: float) -> float:
        """Multiplicative inflation of a SUM's sigma over N independent-noise spaxels."""
        n = max(float(npix), 1.0)
        return float(self.ratio1 * n ** (self.alpha - 0.5))


# ---------------------------------------------------------------------------
def suggest_sky_box(img: np.ndarray, valid: np.ndarray, size: int = 12,
                    avoid: Optional[np.ndarray] = None) -> Tuple[int, int, int, int]:
    """Faintest, flattest fully-valid ``size x size`` box of a white-light image.

    Score = mean + 3*std of the box (so a faint but structured region -- e.g. the wing of a lens
    galaxy -- loses to a flat dark one).  Boxes touching the footprint edge (1 pixel erosion) or an
    ``avoid`` mask are excluded because the outermost spaxels have the lowest throughput.
    """
    from scipy.ndimage import binary_erosion
    ny, nx = img.shape
    size = int(max(3, min(size, ny - 2, nx - 2)))
    ok = binary_erosion(valid & np.isfinite(img), iterations=1)
    if avoid is not None:
        ok &= ~avoid
    f = np.where(ok, img, 0.0).astype(float)
    c1 = np.cumsum(np.cumsum(np.pad(f, ((1, 0), (1, 0))), axis=0), axis=1)
    c2 = np.cumsum(np.cumsum(np.pad(f * f, ((1, 0), (1, 0))), axis=0), axis=1)
    cn = np.cumsum(np.cumsum(np.pad(ok.astype(float), ((1, 0), (1, 0))), axis=0), axis=1)

    def boxsum(c):
        return c[size:, size:] - c[:-size, size:] - c[size:, :-size] + c[:-size, :-size]

    n = boxsum(cn)
    full = n >= size * size - 0.5
    if not full.any():
        # relax: allow boxes with >= 90% valid pixels
        full = n >= 0.9 * size * size
        if not full.any():
            raise ValueError("No source-free sky box of the requested size fits inside the valid field.")
    s1, s2 = boxsum(c1), boxsum(c2)
    nn = np.where(n > 0, n, 1.0)
    mean = s1 / nn
    var = np.clip(s2 / nn - mean ** 2, 0, None)
    score = np.where(full, mean + 3.0 * np.sqrt(var), np.inf)
    y0, x0 = np.unravel_index(np.argmin(score), score.shape)
    return int(y0), int(y0 + size), int(x0), int(x0 + size)


def parse_box(text: str) -> Tuple[int, int, int, int]:
    parts = [int(float(t)) for t in str(text).replace(";", ",").split(",") if t.strip()]
    if len(parts) != 4:
        raise ValueError(f"sky box must be 'y0,y1,x0,x1', got {text!r}")
    y0, y1, x0, x1 = parts
    if y1 <= y0 or x1 <= x0:
        raise ValueError(f"sky box limits are empty: {text!r}")
    return y0, y1, x0, x1


def _block_variance_ratio(sub: np.ndarray, expected_var: np.ndarray, s: int) -> float:
    """Mean over wavelength of  Var[block sums]/E[var]  for s x s blocks (tiled, all offsets)."""
    nw, h, w = sub.shape
    ratios = []
    for oy in range(min(s, 2)):
        for ox in range(min(s, 2)):
            hh, ww = (h - oy) // s, (w - ox) // s
            if hh * ww < 4:
                continue
            blk = sub[:, oy:oy + hh * s, ox:ox + ww * s].reshape(nw, hh, s, ww, s)
            sums = np.nansum(blk, axis=(2, 4)).reshape(nw, -1)            # (nw, ntiles)
            ev = np.nansum(expected_var[:, oy:oy + hh * s, ox:ox + ww * s].reshape(nw, hh, s, ww, s),
                           axis=(2, 4)).reshape(nw, -1).mean(axis=1)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                v = np.nanvar(sums, axis=1, ddof=1)
            r = v / np.where(ev > 0, ev, np.nan)
            ratios.append(r)
    if not ratios:
        return np.nan
    r = np.nanmean(np.vstack(ratios), axis=0)
    r = r[np.isfinite(r)]
    if r.size == 0:
        return np.nan
    lo, hi = np.percentile(r, [10, 90])           # reject sky-line / edge channels
    r = r[(r >= lo) & (r <= hi)]
    return float(np.mean(r)) if r.size else np.nan


def estimate_noise_scaling(arm: Arm, box, bg: np.ndarray, noise: np.ndarray,
                           use_var_cube: bool) -> Tuple[float, float, dict]:
    """Return (alpha, ratio1, diagnostics) from block sums inside the sky box."""
    y0, y1, x0, x1 = box
    sub = arm.flux[:, y0:y1, x0:x1].astype(np.float64) - bg[:, None, None]
    if use_var_cube and arm.var is not None:
        ev = arm.var[:, y0:y1, x0:x1].astype(np.float64)
        ev = np.where(np.isfinite(ev), ev, np.nanmedian(ev))
    else:
        ev = np.broadcast_to((noise ** 2)[:, None, None], sub.shape).copy()
    h, w = sub.shape[1:]
    smax = max(1, min(h, w) // 3)
    sizes, ratios = [], []
    for s in range(1, smax + 1):
        r = _block_variance_ratio(sub, ev, s)
        if np.isfinite(r) and r > 0:
            sizes.append(s * s)
            ratios.append(r)
    diag = dict(npix=np.array(sizes), var_ratio=np.array(ratios))
    if len(sizes) < 2:
        r1 = float(np.sqrt(ratios[0])) if ratios else 1.0
        return 0.5, r1, diag
    logn, logr = np.log(sizes), np.log(ratios)
    slope, icpt = np.polyfit(logn, logr, 1)       # ratio_N^2 = r1^2 * N^(2 alpha - 1)
    alpha = float(np.clip((slope + 1.0) / 2.0, 0.5, 1.0))
    ratio1 = float(np.sqrt(np.exp(icpt)))
    if not use_var_cube:
        ratio1 = 1.0                              # sigma already measured from this very scatter
    diag.update(alpha=alpha, ratio1=ratio1, slope=float(slope))
    return alpha, ratio1, diag


def estimate_sky(arm: Arm, box, noise_source: str = "auto", stat: str = "median",
                 correct_correlation: bool = True) -> SkyModel:
    """Measure background, noise and noise correlation inside ``box`` for one arm."""
    y0, y1, x0, x1 = box
    ny, nx = arm.flux.shape[1:]
    if not (0 <= y0 < y1 <= ny and 0 <= x0 < x1 <= nx):
        raise ValueError(f"sky box {box} outside the {ny}x{nx} cube")
    sub = arm.flux[:, y0:y1, x0:x1].reshape(arm.flux.shape[0], -1)
    n_valid = np.isfinite(sub).sum(axis=1)
    if np.median(n_valid) < 6:
        raise ValueError("Sky box contains fewer than 6 valid spaxels -- it probably lies outside the "
                         "IFU footprint.")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        bg = np.nanmedian(sub, axis=1) if stat == "median" else np.nanmean(sub, axis=1)
        mad = np.nanmedian(np.abs(sub - np.nanmedian(sub, axis=1)[:, None]), axis=1)
    noise = 1.4826 * mad
    bg = np.where(np.isfinite(bg), bg, 0.0)
    good = np.isfinite(noise) & (noise > 0)
    if good.sum() < 10:
        raise ValueError("Cannot estimate the noise from this sky box (all-constant data?).")
    noise = np.interp(np.arange(noise.size), np.flatnonzero(good), noise[good])

    use_var = noise_source == "variance_cube" or (noise_source == "auto" and arm.var is not None)
    if use_var and arm.var is None:
        raise ValueError(f"noise_source='variance_cube' but the {arm.name} arm has no variance cube")
    alpha, ratio1, diag = (0.5, 1.0, {})
    if correct_correlation:
        try:
            alpha, ratio1, diag = estimate_noise_scaling(arm, box, bg, noise, use_var)
        except Exception:                                  # tiny box etc.: fall back to independent
            alpha, ratio1, diag = 0.5, 1.0, {}
    n_pix = int(np.median(n_valid))
    model = SkyModel(box=tuple(box), bg=bg, noise=noise, n_pix=n_pix, alpha=alpha, ratio1=ratio1,
                     source="variance_cube" if use_var else "sky_box", scaling=diag)
    # variance of the background *estimate* (fully correlated across every spaxel that uses it)
    fac = model.correlation_factor(n_pix) if correct_correlation else 1.0
    model.bg_var = (fac * noise) ** 2 / max(n_pix, 1) ** 2 * (np.pi / 2 if stat == "median" else 1.0)
    return model


def spaxel_variance(arm: Arm, sky: SkyModel, idx: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    """Per-spaxel variance (nw, npix) for the spaxels ``idx = (ys, xs)``."""
    if sky.source == "variance_cube" and arm.var is not None:
        v = arm.var[:, idx[0], idx[1]].astype(np.float64)
        # fill isolated NaNs with the channel's typical variance so one bad pixel does not poison a sum
        med = np.nanmedian(v, axis=1)
        bad = ~np.isfinite(v)
        if bad.any():
            v[bad] = np.broadcast_to(med[:, None], v.shape)[bad]
        return v      # sky.ratio1 (variance-cube calibration) is applied at aperture level
    return np.broadcast_to((sky.noise ** 2)[:, None], (sky.noise.size, idx[0].size)).copy()
