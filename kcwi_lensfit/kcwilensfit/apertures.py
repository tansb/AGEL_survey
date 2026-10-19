"""Aperture geometry: mask rendering, DS9 region import, automatic peak detection.

Conventions
-----------
* Pixel coordinates are **0-indexed**, ``x`` along the last (column) axis, ``y`` along rows; pixel
  centres are at integer coordinates.  DS9 'image'/'physical' coordinates are 1-indexed and are
  converted on import.
* Apertures are defined once, in the *reference* (blue) cube's pixels, and mapped to the red arm
  through the arm registration (see :mod:`kcwilensfit.registration`).
* Masks are **fractional** (sub-pixel coverage in [0, 1]) rather than binary, so a half-covered
  edge spaxel contributes half its flux -- smooth, unbiased aperture photometry.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field, replace
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .schema import APERTURE_DEFAULTS, normalise_aperture


@dataclass
class Aperture:
    id: str
    role: str = "source"
    group: str = ""
    shape: str = "circle"
    x: float = 0.0
    y: float = 0.0
    a: float = 3.0
    b: float = 3.0
    angle: float = 0.0
    profile: str = "auto"
    mode: str = "auto"
    z_known: Optional[float] = None
    z_min: Optional[float] = None
    z_max: Optional[float] = None
    # optional sky definition (region files in fk5/icrs): mapped to each cube through ITS OWN WCS
    ra: Optional[float] = None          # deg
    dec: Optional[float] = None         # deg
    a_as: Optional[float] = None        # arcsec (radius / semi-major / width)
    b_as: Optional[float] = None        # arcsec (semi-minor / height)

    @property
    def has_sky(self) -> bool:
        return self.ra is not None and self.dec is not None

    @classmethod
    def from_dict(cls, d: dict) -> "Aperture":
        d = normalise_aperture(d)
        return cls(id=d["id"], role=d["role"], group=d["group"], shape=d["shape"], x=d["x"], y=d["y"],
                   a=d["a"], b=d["b"], angle=d["angle"], profile=d["profile"], mode=d["mode"],
                   z_known=None if d["z_known"] == "" else float(d["z_known"]),
                   z_min=None if d["z_min"] == "" else float(d["z_min"]),
                   z_max=None if d["z_max"] == "" else float(d["z_max"]),
                   ra=None if d.get("ra", "") == "" else float(d["ra"]),
                   dec=None if d.get("dec", "") == "" else float(d["dec"]),
                   a_as=None if d.get("a_as", "") == "" else float(d["a_as"]),
                   b_as=None if d.get("b_as", "") == "" else float(d["b_as"]))

    def to_dict(self) -> dict:
        d = dict(id=self.id, role=self.role, group=self.group, shape=self.shape, x=self.x, y=self.y,
                 a=self.a, b=self.b, angle=self.angle, profile=self.profile, mode=self.mode,
                 z_known="" if self.z_known is None else self.z_known,
                 z_min="" if self.z_min is None else self.z_min,
                 z_max="" if self.z_max is None else self.z_max,
                 ra="" if self.ra is None else self.ra, dec="" if self.dec is None else self.dec,
                 a_as="" if self.a_as is None else self.a_as, b_as="" if self.b_as is None else self.b_as)
        return d

    @property
    def semi_axes(self) -> Tuple[float, float]:
        if self.shape == "circle":
            return self.a, self.a
        if self.shape == "box":
            return self.a / 2.0, self.b / 2.0
        return self.a, self.b

    @property
    def effective_profile(self) -> str:
        if self.profile != "auto":
            return self.profile
        return "sersic" if self.role == "deflector" else "mask"

    def moved(self, dx: float = 0.0, dy: float = 0.0, scale: float = 1.0) -> "Aperture":
        """Copy translated by (dx, dy) pixels with sizes multiplied by ``scale``."""
        return replace(self, x=self.x + dx, y=self.y + dy, a=self.a * scale, b=self.b * scale)


def _axes_sky_orientation(wcs):
    """(theta_x, parity): sky angle (deg, in the East-North basis) of the pixel +x axis, and
    +1/-1 according to whether rotating x -> y is counter-clockwise/clockwise on the sky."""
    cx, cy = (wcs.wcs.crpix[0] - 1, wcs.wcs.crpix[1] - 1)
    a = wcs.pixel_to_world(cx, cy)
    bx = wcs.pixel_to_world(cx + 1.0, cy)
    by = wcs.pixel_to_world(cx, cy + 1.0)
    cosd = np.cos(np.deg2rad(a.dec.deg))
    vx = np.array([(bx.ra.deg - a.ra.deg) * cosd, bx.dec.deg - a.dec.deg])
    vy = np.array([(by.ra.deg - a.ra.deg) * cosd, by.dec.deg - a.dec.deg])
    thx = float(np.degrees(np.arctan2(vx[1], vx[0])))
    parity = 1.0 if (vx[0] * vy[1] - vx[1] * vy[0]) > 0 else -1.0
    return thx, parity


def aperture_in_arm(ap: Aperture, arm_wcs, arm_pixscale: float, ref_wcs=None) -> Aperture:
    """Map a sky-defined aperture onto one cube's pixel grid with THAT cube's WCS.

    Centre -> world_to_pixel; sizes -> arcsec / that cube's pixel scale; angle -> rotated by the
    difference between the reference and this cube's pixel-axis orientation on the sky.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    px, py = arm_wcs.world_to_pixel(SkyCoord(ap.ra * u.deg, ap.dec * u.deg, frame="icrs"))
    a = abs(ap.a_as) / arm_pixscale
    b = abs(ap.b_as if ap.b_as is not None else ap.a_as) / arm_pixscale
    ang = ap.angle
    if ref_wcs is not None and ap.shape != "circle":
        thr, pr = _axes_sky_orientation(ref_wcs)
        tha, pa = _axes_sky_orientation(arm_wcs)
        # the same SKY direction: theta_x + parity*phi is invariant
        sky = thr + pr * ap.angle
        ang = pa * (sky - tha)
    return replace(ap, x=float(px), y=float(py), a=float(a), b=float(b), angle=float(ang))


def aperture_mask(shape: Tuple[int, int], ap: Aperture, oversample: int = 5) -> np.ndarray:
    """Fractional-coverage mask (float32, values in [0, 1]) of ``ap`` on an (ny, nx) grid."""
    ny, nx = shape
    mask = np.zeros(shape, np.float32)
    ra, rb = (abs(r) for r in ap.semi_axes)      # DS9 tolerates negative sizes; they are just magnitudes
    ext = np.hypot(ra, rb) + 1.5
    x0, x1 = max(0, int(np.floor(ap.x - ext))), min(nx, int(np.ceil(ap.x + ext)) + 1)
    y0, y1 = max(0, int(np.floor(ap.y - ext))), min(ny, int(np.ceil(ap.y + ext)) + 1)
    if x1 <= x0 or y1 <= y0:
        return mask
    s = oversample
    # sub-pixel sample positions (pixel centres at integers; each pixel spans +-0.5)
    sub = (np.arange(s) + 0.5) / s - 0.5
    xs = (np.arange(x0, x1)[:, None] + sub[None, :]).ravel()
    ys = (np.arange(y0, y1)[:, None] + sub[None, :]).ravel()
    X, Y = np.meshgrid(xs, ys)
    th = np.deg2rad(ap.angle)
    dx, dy = X - ap.x, Y - ap.y
    u = dx * np.cos(th) + dy * np.sin(th)          # along the (rotated) first axis
    v = -dx * np.sin(th) + dy * np.cos(th)
    if ap.shape == "box":
        inside = (np.abs(u) <= ra) & (np.abs(v) <= rb)
    else:
        inside = (u / ra) ** 2 + (v / rb) ** 2 <= 1.0
    cov = inside.reshape(y1 - y0, s, x1 - x0, s).mean(axis=(1, 3))
    mask[y0:y1, x0:x1] = cov
    return mask


def pixel_count(mask: np.ndarray) -> float:
    return float(np.nansum(mask))


# ---------------------------------------------------------------------------
# DS9 region files
# ---------------------------------------------------------------------------
_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"


def _parse_angle_or_len(tok: str, pixscale_arcsec: float, is_sky: bool) -> float:
    """Length token with optional unit -> pixels."""
    tok = tok.strip()
    if tok.endswith('"'):
        return float(tok[:-1]) / pixscale_arcsec
    if tok.endswith("'"):
        return float(tok[:-1]) * 60.0 / pixscale_arcsec
    if tok.endswith("d"):
        return float(tok[:-1]) * 3600.0 / pixscale_arcsec
    if tok.endswith("r"):
        return float(tok[:-1]) * 206264.806 / pixscale_arcsec
    if tok.endswith("i") or tok.endswith("p"):
        return float(tok[:-1])
    v = float(tok)
    return v / pixscale_arcsec if is_sky else v     # DS9 default sky length unit is arcsec


def parse_ds9_regions(path: str, wcs=None, pixscale_arcsec: float = 0.3) -> Tuple[List[Aperture], List[str]]:
    """Parse circle / ellipse / box regions into :class:`Aperture` objects.

    Supports the 'image'/'physical' pixel systems (1-indexed in DS9 -> converted to 0-indexed) and
    'fk5/icrs/j2000' sky systems (needs ``wcs``).  Returns (apertures, warnings).  Region text labels
    (``# text={name}``) become aperture ids; labels containing 'defl'/'lens' mark deflectors.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    warnings: List[str] = []
    system = "physical"
    out: List[Aperture] = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("global"):
                continue
            label = None
            if "#" in line:
                line, _, comment = line.partition("#")
                line = line.strip()
                m = re.search(r"text=\{([^}]*)\}", comment)
                if m:
                    label = m.group(1).strip()
                if not line:
                    continue
            low = line.lower()
            if low in ("image", "physical", "fk5", "icrs", "j2000", "wcs", "galactic", "fk4", "ecliptic"):
                system = low
                continue
            m = re.match(r"^([+-]?)\s*(circle|ellipse|box|annulus|polygon|point)\s*\((.*)\)", line, re.I)
            if not m:
                warnings.append(f"skipped unsupported region line: {line[:60]!r}")
                continue
            neg, kind, body = m.group(1), m.group(2).lower(), m.group(3)
            if neg == "-":
                warnings.append("exclusion regions are ignored")
                continue
            if kind not in ("circle", "ellipse", "box"):
                warnings.append(f"region type '{kind}' is not supported (use circle/ellipse/box)")
                continue
            toks = [t.strip() for t in body.split(",")]
            is_sky = system in ("fk5", "icrs", "j2000", "wcs")
            if is_sky:
                if wcs is None:
                    warnings.append("sky-coordinate region but the cube has no WCS; skipped")
                    continue
                unit = (u.hourangle, u.deg) if (":" in toks[0] or "h" in toks[0]) else (u.deg, u.deg)
                sc = SkyCoord(toks[0], toks[1], unit=unit, frame="icrs")
                px, py = wcs.world_to_pixel(sc)
                x, y = float(px), float(py)
            else:
                x, y = float(toks[0]) - 1.0, float(toks[1]) - 1.0
            ang = 0.0
            sky_info = {}
            if is_sky:
                sky_info = dict(ra=float(sc.icrs.ra.deg), dec=float(sc.icrs.dec.deg))
            if kind == "circle":
                r = _parse_angle_or_len(toks[2], pixscale_arcsec, is_sky)
                ap = Aperture(id="", shape="circle", x=x, y=y, a=abs(r), b=abs(r))
            elif kind == "ellipse":
                a = _parse_angle_or_len(toks[2], pixscale_arcsec, is_sky)
                b = _parse_angle_or_len(toks[3], pixscale_arcsec, is_sky)
                ang = float(toks[4]) if len(toks) > 4 else 0.0
                ap = Aperture(id="", shape="ellipse", x=x, y=y, a=abs(a), b=abs(b), angle=ang)
            else:
                wdt = _parse_angle_or_len(toks[2], pixscale_arcsec, is_sky)
                hgt = _parse_angle_or_len(toks[3], pixscale_arcsec, is_sky)
                ang = float(toks[4]) if len(toks) > 4 else 0.0
                ap = Aperture(id="", shape="box", x=x, y=y, a=abs(wdt), b=abs(hgt), angle=ang)
            if sky_info:
                ps = pixscale_arcsec
                ap.ra, ap.dec = sky_info["ra"], sky_info["dec"]
                ap.a_as = ap.a * ps
                ap.b_as = ap.b * ps
            if abs(ang % 360.0) < 1e-3 or abs(ang % 360.0 - 360.0) < 1e-3:
                ap.angle = 0.0
            n = len(out) + 1
            ap.id = label or f"region_{n}"
            if label and re.search(r"defl|lens|lrg", label, re.I):
                ap.role = "deflector"
            out.append(ap)
    # unique ids
    seen: Dict[str, int] = {}
    for ap in out:
        if ap.id in seen:
            seen[ap.id] += 1
            ap.id = f"{ap.id}_{seen[ap.id]}"
        else:
            seen[ap.id] = 1
    if system in ("physical", "image"):
        warnings.append("region file uses image/physical coordinates: assumed 1-indexed DS9 pixels "
                        "of the BLUE cube (converted to 0-indexed)")
    return out, warnings


def detect_peaks(img: np.ndarray, valid: np.ndarray, noise: Optional[float] = None, nsigma: float = 5.0,
                 smooth_sigma: float = 1.0, max_peaks: int = 8, min_sep: float = 3.0,
                 default_radius: float = 2.5) -> List[Aperture]:
    """Starting-point apertures from local maxima of the (lightly smoothed) white-light image.

    The brightest peak is labelled 'deflector' and the rest 'source'; this is only a convenient
    seed for the aperture table -- the roles are *your* decision (a bright companion or a star is
    easy to mislabel automatically).
    """
    from scipy.ndimage import gaussian_filter, maximum_filter, binary_erosion
    f = np.where(valid & np.isfinite(img), img, np.nanmedian(img[valid]))
    sm = gaussian_filter(f, smooth_sigma)
    base = np.median(sm[valid])
    if noise is None:
        noise = 1.4826 * np.median(np.abs(sm[valid] - base))
    thr = base + nsigma * max(noise, 1e-30)
    size = max(3, int(round(2 * min_sep)) | 1)
    peaks = (sm == maximum_filter(sm, size=size)) & (sm > thr) & binary_erosion(valid, iterations=2)
    ys, xs = np.nonzero(peaks)
    order = np.argsort(sm[ys, xs])[::-1][:max_peaks]
    out = []
    for k, i in enumerate(order):
        out.append(Aperture(id="deflector" if k == 0 else f"src{k}", role="deflector" if k == 0 else "source",
                            shape="circle", x=float(xs[i]), y=float(ys[i]),
                            a=default_radius, b=default_radius))
    return out


def overlap_fraction(m1: np.ndarray, m2: np.ndarray) -> float:
    """Fraction of the smaller aperture's area that lies inside the other."""
    inter = float(np.sum(np.minimum(m1, m2)))
    smaller = float(min(np.sum(m1), np.sum(m2)))
    return inter / smaller if smaller > 0 else 0.0
