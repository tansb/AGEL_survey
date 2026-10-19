"""Blue <-> red spatial registration.

KCWI and KCRM sit behind the same image slicer after the dichroic, so a given spaxel index sees
the same patch of sky in both arms *by instrument design*: the physically expected offset is
zero.  Header WCS zero-points are set by separate astrometric solutions and can be arcseconds
off, and image cross-correlation conflates "where the flux is" with "where the geometry is" (a
blue star-forming arc and a red quiescent core do not trace the same structure).  The default is
therefore ``identity``; the other estimates are always *computed* and reported so you can check
them, and the GUI overlays red contours on the blue image for a by-eye test.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np
from scipy.signal import fftconvolve

from .cubeio import Arm


@dataclass
class Registration:
    dx: float = 0.0          # red x = blue x + dx
    dy: float = 0.0
    scale: float = 1.0       # red pixels per blue pixel along both axes (1/pixscale ratio)
    mode: str = "identity"
    note: str = ""
    diagnostics: dict = field(default_factory=dict)

    def to_red(self, ap):
        """Aperture in the red arm's pixel frame."""
        return ap.moved(self.dx, self.dy, self.scale)


def measure_xcorr_offset(img_a: np.ndarray, img_b: np.ndarray, max_shift: int = 12
                         ) -> Tuple[float, float, float]:
    """Sub-pixel (dx, dy) such that features of img_b sit at (x+dx, y+dy) relative to img_a.

    Returns (dx, dy, sharpness) where sharpness = peak/zero-lag correlation.  Values near 1 mean a
    broad, structureless correlation surface -- the offset is then NOT a measurement.
    """
    a = np.where(np.isfinite(img_a), img_a, 0.0).astype(float)
    b = np.where(np.isfinite(img_b), img_b, 0.0).astype(float)
    a = a - np.median(a[np.isfinite(img_a)]) if np.isfinite(img_a).any() else a
    b = b - np.median(b[np.isfinite(img_b)]) if np.isfinite(img_b).any() else b
    a[~np.isfinite(img_a)] = 0
    b[~np.isfinite(img_b)] = 0
    a /= (a.std() + 1e-30)
    b /= (b.std() + 1e-30)
    corr = fftconvolve(a, b[::-1, ::-1], mode="same")
    ny, nx = corr.shape
    cy, cx = ny // 2, nx // 2
    ms = int(min(max_shift, cy - 2, cx - 2))
    sub = corr[cy - ms:cy + ms + 1, cx - ms:cx + ms + 1]
    py, px = np.unravel_index(np.argmax(sub), sub.shape)
    # parabolic sub-pixel refinement
    def refine(c_m, c_0, c_p):
        d = c_m - 2 * c_0 + c_p
        return 0.0 if d == 0 else 0.5 * (c_m - c_p) / d
    fy = refine(sub[py - 1, px], sub[py, px], sub[py + 1, px]) if 0 < py < sub.shape[0] - 1 else 0.0
    fx = refine(sub[py, px - 1], sub[py, px], sub[py, px + 1]) if 0 < px < sub.shape[1] - 1 else 0.0
    dy, dx = py - ms + fy, px - ms + fx
    zero = sub[ms, ms]
    sharp = float(sub[py, px] / zero) if zero > 0 else float("nan")
    return float(dx), float(dy), sharp


def wcs_offset(blue: Arm, red: Arm) -> Optional[Tuple[float, float]]:
    """(dx, dy) in pixels implied by the two header WCSs, evaluated at the blue footprint centre."""
    if blue.wcs is None or red.wcs is None:
        return None
    ys, xs = np.nonzero(blue.valid)
    cx, cy = float(np.mean(xs)), float(np.mean(ys))
    sky = blue.wcs.pixel_to_world(cx, cy)
    rx, ry = red.wcs.world_to_pixel(sky)
    return float(rx - cx), float(ry - cy)


def build_registration(blue: Optional[Arm], red: Optional[Arm], wl_blue, wl_red, mode: str = "identity",
                       manual: Tuple[float, float] = (0.0, 0.0)) -> Registration:
    if blue is None or red is None:
        return Registration(mode="identity", note="single-arm data: no registration needed")
    scale = blue.pixscale[1] / red.pixscale[1] if red.pixscale[1] > 0 else 1.0
    diag = {}
    wo = wcs_offset(blue, red)
    if wo is not None:
        diag["wcs_dx"], diag["wcs_dy"] = wo
    try:
        dx, dy, sharp = measure_xcorr_offset(wl_blue, wl_red)
        diag.update(xcorr_dx=dx, xcorr_dy=dy, xcorr_sharpness=sharp)
    except Exception:
        pass
    if mode == "identity":
        return Registration(0.0, 0.0, scale, "identity",
                            "same pixel index in both arms (instrument-design expectation)", diag)
    if mode == "wcs":
        if wo is None:
            raise ValueError("arm_registration='wcs' needs a celestial WCS in both cubes")
        return Registration(wo[0], wo[1], scale, "wcs", "from header WCS", diag)
    if mode == "xcorr":
        if "xcorr_dx" not in diag:
            raise ValueError("cross-correlation registration failed")
        note = "from white-light cross-correlation"
        if diag.get("xcorr_sharpness", 0) < 1.15:
            note += " (WARNING: correlation surface is flat -- not a trustworthy measurement)"
        return Registration(diag["xcorr_dx"], diag["xcorr_dy"], scale, "xcorr", note, diag)
    if mode == "manual":
        return Registration(float(manual[0]), float(manual[1]), scale, "manual", "user-specified", diag)
    raise ValueError(f"unknown registration mode {mode!r}")
