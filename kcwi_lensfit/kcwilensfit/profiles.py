"""Spatial light-profile models used to *deblend* overlapping objects.

In a lens system the deflector's halo, the lensed arcs and (for compound lenses) a second lens or
an intermediate-redshift source all overlap on the IFU at KCWI's seeing-limited resolution.
Rather than summing spaxels inside an aperture and hoping the contamination is small, we model
every object's *spatial shape* once (from a white-light image) and then solve for each object's
*spectrum* at every wavelength by linear least squares -- exactly how crowded-field IFU
extraction (e.g. PSF-fitting codes) is done.  Shapes:

* ``sersic`` -- elliptical Sersic profile convolved with the (Gaussian) seeing: lens galaxies;
* ``psf``    -- the seeing disc itself: point-like lensed images or stars;
* ``mask``   -- the user's aperture as a fractional mask: arcs, whose true shape is unknown.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import least_squares, nnls

from .apertures import Aperture, aperture_mask

OVERSAMPLE = 3


def _bn(n: float) -> float:
    """Sersic b_n (Ciotti & Bertin 1999 asymptotic expansion; accurate to <1e-4 for n>0.36)."""
    return 2.0 * n - 1.0 / 3.0 + 4.0 / (405.0 * n) + 46.0 / (25515.0 * n ** 2)


def render_profile(shape: Tuple[int, int], kind: str, x0: float, y0: float, psf_sigma: float,
                   re: float = 3.0, n: float = 2.0, q: float = 1.0, pa_deg: float = 0.0,
                   offset: Tuple[int, int] = (0, 0), oversample: int = OVERSAMPLE) -> np.ndarray:
    """Render a PSF-convolved, pixel-integrated profile on an ``(ny, nx)`` window.

    ``x0, y0`` are in the *window's* pixel frame; ``offset`` is unused except for bookkeeping.
    The result is normalised to unit sum over the window (conditioning of the linear solve).
    """
    ny, nx = shape
    s = oversample
    xs = (np.arange(nx * s) + 0.5) / s - 0.5
    ys = (np.arange(ny * s) + 0.5) / s - 0.5
    X, Y = np.meshgrid(xs, ys)
    dx, dy = X - x0, Y - y0
    if kind == "psf":
        sig = max(psf_sigma, 0.3)
        sup = np.exp(-0.5 * (dx * dx + dy * dy) / sig ** 2)
        img = sup.reshape(ny, s, nx, s).mean(axis=(1, 3))
    elif kind == "sersic":
        th = np.deg2rad(pa_deg)
        u = dx * np.cos(th) + dy * np.sin(th)
        v = -dx * np.sin(th) + dy * np.cos(th)
        r = np.sqrt(u * u + (v / max(q, 0.05)) ** 2)
        sup = np.exp(-_bn(n) * ((r / re) ** (1.0 / n) - 1.0))
        if psf_sigma > 0.05:
            sup = gaussian_filter(sup, sigma=psf_sigma * s, mode="constant")
        img = sup.reshape(ny, s, nx, s).mean(axis=(1, 3))
    else:
        raise ValueError(kind)
    tot = img.sum()
    return (img / tot).astype(np.float64) if tot > 0 else img.astype(np.float64)


@dataclass
class Component:
    """One spatial template, always expressed in the pixel frame of ``frame_arm``."""
    ap: Aperture
    kind: str                                  # 'sersic' | 'psf' | 'mask'
    params: Dict[str, float] = field(default_factory=dict)
    frame_arm: str = "blue"
    psf_sigma: float = 1.0

    def render(self, shape: Tuple[int, int], dx: float = 0.0, dy: float = 0.0, scale: float = 1.0
               ) -> np.ndarray:
        """Template on the FULL grid ``shape`` of some arm, shifted by (dx, dy) from the fit frame.

        The template is normalised to unit sum *over the whole grid* (negative-free, float64).
        """
        ny, nx = shape
        if self.kind == "mask":
            m = aperture_mask(shape, self.ap.moved(dx, dy, scale)).astype(np.float64)
            tot = m.sum()
            return m / tot if tot > 0 else m
        p = self.params
        # render on a window containing the object to keep the cost bounded
        pad = int(np.ceil(6 * max(self.psf_sigma, 1.0) + 3 * p.get("re", 1.0) * scale))
        pad = max(10, min(pad, 70))
        cx, cy = p["x0"] + dx, p["y0"] + dy
        x_lo, x_hi = max(0, int(cx) - pad), min(nx, int(cx) + pad + 1)
        y_lo, y_hi = max(0, int(cy) - pad), min(ny, int(cy) + pad + 1)
        out = np.zeros(shape, np.float64)
        if x_hi <= x_lo or y_hi <= y_lo:
            return out
        sub = render_profile((y_hi - y_lo, x_hi - x_lo), self.kind, cx - x_lo, cy - y_lo,
                             self.psf_sigma * scale, re=p.get("re", 3.0) * scale, n=p.get("n", 2.0),
                             q=p.get("q", 1.0), pa_deg=p.get("pa", 0.0))
        out[y_lo:y_hi, x_lo:x_hi] = sub
        tot = out.sum()
        return out / tot if tot > 0 else out


@dataclass
class ProfileFit:
    components: List[Component]
    model: np.ndarray                 # best-fit white-light model on the full grid
    amplitudes: np.ndarray
    chi2red: float
    window: Tuple[int, int, int, int]  # y0, y1, x0, x1
    flags: List[str] = field(default_factory=list)
    background: float = 0.0
    success: bool = True


def fit_profiles(image: np.ndarray, valid: np.ndarray, apertures: Sequence[Aperture], psf_sigma: float,
                 sigma_img: float, frame_arm: str = "blue", pad: int = 8, fit_background: bool = False,
                 fit_psf: bool = False, model_floor: float = 0.03) -> ProfileFit:
    """Fit all spatial components to ONE (white-light) image simultaneously.

    Variable projection: the non-linear shape parameters (centre, r_e, n, q, PA) are optimised by
    ``least_squares`` while the non-negative amplitudes (and optional flat background) are
    re-solved by NNLS inside every residual evaluation -- far better conditioned than fitting
    amplitudes and shapes together.
    """
    ny, nx = image.shape
    ext = [(ap.x - np.hypot(*ap.semi_axes) - pad, ap.x + np.hypot(*ap.semi_axes) + pad,
            ap.y - np.hypot(*ap.semi_axes) - pad, ap.y + np.hypot(*ap.semi_axes) + pad) for ap in apertures]
    x_lo = max(0, int(np.floor(min(e[0] for e in ext))))
    x_hi = min(nx, int(np.ceil(max(e[1] for e in ext))) + 1)
    y_lo = max(0, int(np.floor(min(e[2] for e in ext))))
    y_hi = min(ny, int(np.ceil(max(e[3] for e in ext))) + 1)
    win = (y_lo, y_hi, x_lo, x_hi)
    wshape = (y_hi - y_lo, x_hi - x_lo)
    data = image[y_lo:y_hi, x_lo:x_hi].astype(float)
    ok = np.isfinite(data) & valid[y_lo:y_hi, x_lo:x_hi]
    d = np.where(ok, data, 0.0)
    sig = np.sqrt(sigma_img ** 2 + (model_floor * np.clip(d, 0, None)) ** 2)
    w = np.where(ok, 1.0 / sig, 0.0).ravel()

    comps: List[Component] = []
    p0, lo, hi, slices = [], [], [], []
    for ap in apertures:
        kind = ap.effective_profile
        c = Component(ap=ap, kind=kind, frame_arm=frame_arm, psf_sigma=psf_sigma)
        k0 = len(p0)
        if kind == "sersic":
            ra, rb = ap.semi_axes
            re0 = float(np.clip(0.6 * max(ra, rb), 1.0, 20.0))
            q0 = float(np.clip(min(ra, rb) / max(ra, rb), 0.4, 1.0))
            p0 += [ap.x - x_lo, ap.y - y_lo, np.log(re0), 3.0, q0, ap.angle]
            lo += [ap.x - x_lo - 2.5, ap.y - y_lo - 2.5, np.log(0.6), 0.8, 0.25, -360.0]
            hi += [ap.x - x_lo + 2.5, ap.y - y_lo + 2.5, np.log(30.0), 6.0, 1.0, 360.0]
            slices.append(slice(k0, k0 + 6))
        elif kind == "psf":
            p0 += [ap.x - x_lo, ap.y - y_lo]
            lo += [ap.x - x_lo - 2.0, ap.y - y_lo - 2.0]
            hi += [ap.x - x_lo + 2.0, ap.y - y_lo + 2.0]
            slices.append(slice(k0, k0 + 2))
        elif kind == "mask":
            slices.append(slice(k0, k0))
        else:
            raise ValueError(f"unknown profile kind {kind!r} for aperture {ap.id}")
        comps.append(c)
    if fit_psf:
        p0.append(psf_sigma); lo.append(0.4 * psf_sigma); hi.append(2.5 * psf_sigma)

    def templates(theta):
        ps = theta[-1] if fit_psf else psf_sigma
        cols = []
        for c, sl in zip(comps, slices):
            if c.kind == "mask":
                m = aperture_mask(wshape, c.ap.moved(-x_lo, -y_lo)).astype(float)
                t = m / m.sum() if m.sum() > 0 else m
            elif c.kind == "psf":
                x0, y0 = theta[sl]
                t = render_profile(wshape, "psf", x0, y0, ps)
            else:
                x0, y0, lre, n, q, pa = theta[sl]
                t = render_profile(wshape, "sersic", x0, y0, ps, re=np.exp(lre), n=n, q=q, pa_deg=pa)
            cols.append(t.ravel())
        return np.array(cols).T                                  # (npix, K)

    def solve(theta):
        A = templates(theta)
        cols = [A]
        if fit_background:
            one = np.ones((A.shape[0], 1)) / A.shape[0]
            cols += [one, -one]
        M = np.hstack(cols) * w[:, None]
        coef, _ = nnls(M, d.ravel() * w)
        return A, coef, M

    def resid(theta):
        _, coef, M = solve(theta)
        return M @ coef - d.ravel() * w

    flags: List[str] = []
    success = True
    if len(p0) > 0:
        try:
            res = least_squares(resid, p0, bounds=(lo, hi), x_scale="jac", max_nfev=120, ftol=1e-6, xtol=1e-6)
            theta = res.x
            success = bool(res.success)
        except Exception as e:                                     # pragma: no cover
            flags.append(f"profile optimisation failed ({e}); using initial guesses")
            theta, success = np.array(p0, float), False
    else:
        theta = np.array([], float)
    A, coef, M = solve(theta)
    amp = coef[:A.shape[1]]
    bkg = float(coef[A.shape[1]] - coef[A.shape[1] + 1]) / A.shape[0] if fit_background else 0.0

    # finalise components (convert back to the full-image pixel frame)
    ps = theta[-1] if fit_psf else psf_sigma
    for c, sl in zip(comps, slices):
        c.psf_sigma = float(ps)
        if c.kind == "sersic":
            x0, y0, lre, n, q, pa = theta[sl]
            c.params = dict(x0=float(x0 + x_lo), y0=float(y0 + y_lo), re=float(np.exp(lre)),
                            n=float(n), q=float(q), pa=float(pa % 180.0))
            for nm, v, l, h in zip(("x0", "y0", "ln_re", "n", "q", "pa"), theta[sl], lo[sl.start:sl.stop],
                                   hi[sl.start:sl.stop]):
                if nm in ("x0", "y0", "ln_re", "n", "q") and (abs(v - l) < 1e-3 or abs(v - h) < 1e-3):
                    flags.append(f"{c.ap.id}: {nm.replace('ln_re','r_e')} pinned at its bound -- profile poorly constrained")
        elif c.kind == "psf":
            x0, y0 = theta[sl]
            c.params = dict(x0=float(x0 + x_lo), y0=float(y0 + y_lo))
    # best-fit model on the full grid
    model = np.zeros_like(image, dtype=float)
    model[y_lo:y_hi, x_lo:x_hi] = (A @ amp).reshape(wshape) + bkg
    n_free = max(int(ok.sum()) - len(theta) - len(amp), 1)
    chi2 = float(np.sum(resid(theta) ** 2)) / n_free
    # template degeneracy check
    if A.shape[1] > 1:
        An = A / np.linalg.norm(A, axis=0)
        cs = An.T @ An
        np.fill_diagonal(cs, 0)
        i, j = np.unravel_index(np.argmax(cs), cs.shape)
        if cs[i, j] > 0.97:
            flags.append(f"templates of '{comps[i].ap.id}' and '{comps[j].ap.id}' are nearly identical "
                         f"(cosine similarity {cs[i, j]:.3f}) -- their spectra cannot be separated")
    return ProfileFit(components=comps, model=model, amplitudes=amp, chi2red=chi2, window=win,
                      flags=flags, background=bkg, success=success)
