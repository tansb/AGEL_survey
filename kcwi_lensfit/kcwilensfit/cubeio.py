"""Loading KCWI / KCRM cubes into a uniform in-memory representation.

Handled here (each is a classic way for a redshift to go quietly wrong):

* the wavelength axis is rebuilt from the header and its **frame (air/vacuum)** is tracked;
  everything downstream works in air, matching the E-MILES/pPXF convention;
* spaxels outside the IFU footprint are stored as exact zeros by some DRP stacks (or NaN by
  others) -- both are turned into NaN so they can never be mistaken for "very quiet sky";
* the usable bandpass is trimmed using the header ``WAVGOOD0/1`` keywords (or user limits);
* surface-brightness cubes (``BUNIT`` per arcsec^2) can be converted to flux per spaxel;
* the instrumental resolution is identified from grating + slicer and converted to an FWHM.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from . import constants as K
from .wavelength import air_to_vac, apply_velocity_correction, vac_to_air


@dataclass
class Arm:
    """One spectrograph arm (blue = KCWI, red = KCRM) loaded as a cube."""
    name: str
    path: str
    flux: np.ndarray                      # (nw, ny, nx) float32; NaN where no data
    var: Optional[np.ndarray]             # same shape or None (user-supplied variance cube)
    header: fits.Header
    wcs: Optional[WCS]                    # celestial 2-D WCS (None if the header has none)
    wave: np.ndarray                      # AIR wavelengths [A] after any velocity correction
    native_frame: str                     # frame of the wavelength axis in the FILE
    fwhm_A: float                         # instrumental FWHM [A] at the central wavelength
    pixscale: Tuple[float, float]         # arcsec per pixel (y, x)
    valid: np.ndarray                     # (ny, nx) bool spaxel footprint
    unit_scale: float = 1.0               # factor multiplied into flux (area conversion)
    info: List[str] = field(default_factory=list)  # human-readable provenance notes
    grating: str = ""
    slicer: str = ""
    resolving_power: float = 0.0
    vcorr_header: Optional[float] = None
    vcorr_applied_kms: float = 0.0

    @property
    def shape(self):
        return self.flux.shape

    @property
    def pix_area(self) -> float:
        return float(self.pixscale[0] * self.pixscale[1])

    @property
    def central_wave(self) -> float:
        return 0.5 * (self.wave.min() + self.wave.max())


# ---------------------------------------------------------------------------
def _wave_axis(header, path: str = "") -> np.ndarray:
    missing = [k for k in ("CRVAL3",) if k not in header]
    if "CD3_3" not in header and "CDELT3" not in header:
        missing.append("CD3_3/CDELT3")
    if missing:
        raise ValueError(
            f"'{os.path.basename(path)}' is not a spectral cube: its header has no {', '.join(missing)} (the wavelength "
            f"axis keywords). It is a 3-D array with {header.get('NAXIS3')} planes -- e.g. a multi-band image such as a "
            "DECaLS cutout or a '_shifted' alignment file. Choose the reduced KCWI/KCRM cube (e.g. '*_icubes*.fits' or "
            "'*_medium_combined.fits').")
    n = int(header["NAXIS3"])
    crval = float(header["CRVAL3"])
    cdelt = header.get("CD3_3", header.get("CDELT3"))
    crpix = float(header.get("CRPIX3", 1.0))
    wave = crval + float(cdelt) * (np.arange(n) + 1 - crpix)
    unit = str(header.get("CUNIT3", "Angstrom")).strip().lower()
    if unit in ("nm", "nanometer", "nanometers"):
        wave = wave * 10.0
    elif unit in ("m", "meter", "meters"):
        wave = wave * 1e10
    elif unit in ("um", "micron", "microns", "micrometer"):
        wave = wave * 1e4
    return wave


def _frame_from_header(header, arm_name: str, requested: str) -> Tuple[str, str]:
    """Return (frame, how) where frame is 'air' or 'vacuum'."""
    if requested in ("air", "vacuum"):
        return requested, "set by user"
    ctype = str(header.get("CTYPE3", "")).strip().upper()
    if ctype.startswith("AWAV"):
        return "air", "CTYPE3 = AWAV"
    if ctype.startswith("WAVE"):
        return "vacuum", "CTYPE3 = WAVE"
    default = "air" if arm_name == "blue" else "vacuum"
    return default, f"CTYPE3 absent/unrecognised ('{ctype}'); assumed {default} for the {arm_name} arm"


def _tokens(path: str) -> set:
    return set(re.split(r"[^A-Za-z0-9]+", os.path.basename(path).upper()))


def detect_grating(header, path, arm_name: str) -> Tuple[Optional[str], str]:
    allowed = K.BLUE_GRATINGS if arm_name == "blue" else K.RED_GRATINGS
    for key in K.GRATING_KEYS[arm_name]:
        if key in header:
            val = str(header[key]).strip().upper()
            for g in allowed:
                if val == g or g in val.split():
                    return g, f"header {key}={header[key]!r}"
    toks = _tokens(path)
    for g in allowed:
        if g in toks:
            return g, "file-name token"
    return None, "not found"


def detect_slicer(header, path) -> Tuple[Optional[str], str]:
    for key in K.SLICER_KEYS:
        if key in header:
            val = str(header[key]).strip().upper()
            for s in K.SLICER_NAMES:
                if s.upper() in val:
                    return s, f"header {key}={header[key]!r}"
    toks = _tokens(path)
    for s in K.SLICER_NAMES:
        if s.upper() in toks:
            return s, "file-name token"
    return None, "not found"


# ---------------------------------------------------------------------------
def load_arm(name: str, path: str, cfg: dict, var_path: str = "",
             sibling_header: Optional[fits.Header] = None) -> Arm:
    """Load one arm.

    Parameters
    ----------
    cfg : the pipeline configuration (frames, trims, FWHM overrides, velocity correction ...)
    sibling_header : header of the *other* arm, used as a fallback for set-up keywords
        (the Keck DRP writes both BGRATNAM and RGRATNAM into the red cube header, for instance).
    """
    path = os.path.expanduser(path)
    info: List[str] = []
    with fits.open(path, memmap=False) as hdul:
        cands = [h for h in hdul if h.data is not None and h.data.ndim == 3]
        if not cands:
            raise ValueError(f"{path}: no 3-D data cube found in any HDU")
        # prefer the HDU that actually has a wavelength axis, then the one with the most planes
        cands.sort(key=lambda h: ("CRVAL3" in h.header, h.data.shape[0]), reverse=True)
        hdu = cands[0]
        header = hdu.header.copy()
        if "RADESYS" not in header and len(hdul) > 1:
            header.extend(hdul[0].header, update=False, unique=True)
        flux = np.asarray(hdu.data, dtype=np.float32)
    nw, ny, nx = flux.shape

    # ----- wavelength axis & frame ------------------------------------------------
    wave_native = _wave_axis(header, path)
    frame, how = _frame_from_header(header, name, cfg.get(f"{name}_wave_frame", "auto"))
    wave_air = vac_to_air(wave_native) if frame == "vacuum" else wave_native.copy()
    info.append(f"wavelength axis {wave_native[0]:.1f}-{wave_native[-1]:.1f} A, native frame = {frame} ({how})"
                + ("; converted vacuum->air" if frame == "vacuum" else ""))

    # ----- velocity-frame correction ---------------------------------------------
    vcorr_header = header.get("VCORR")
    vcorr_type = str(header.get("VCORRTYP", "")).strip()
    vcorr_applied = 0.0
    mode = cfg.get("velocity_correction", "none")
    if vcorr_header is not None:
        info.append(f"header VCORR = {float(vcorr_header):+.2f} km/s ({vcorr_type or 'type unspecified'}); "
                    "the Keck DRP normally applies this already")
    if mode == "header":
        if vcorr_header is None:
            raise ValueError(f"velocity_correction='header' but the {name} cube has no VCORR keyword")
        vcorr_applied = float(vcorr_header)
    elif mode == "manual":
        vcorr_applied = float(cfg.get("manual_vcorr_kms", 0.0))
    if vcorr_applied != 0.0:
        wave_air = apply_velocity_correction(wave_air, vcorr_applied)
        info.append(f"applied an additional velocity correction of {vcorr_applied:+.2f} km/s to the wavelengths")

    # ----- WCS / pixel scale -------------------------------------------------------
    wcs2d = None
    try:
        w = WCS(header).celestial
        if w.naxis == 2 and w.has_celestial:
            wcs2d = w
            sc = proj_plane_pixel_scales(w) * 3600.0       # (x, y) arcsec/pix
            pixscale = (float(sc[1]), float(sc[0]))
        else:
            raise ValueError
    except Exception:
        pxscl = header.get("PXSCL")
        pixscale = (float(pxscl) * 3600.0,) * 2 if pxscl else (0.3, 0.3)
        info.append(f"no usable celestial WCS; pixel scale taken as {pixscale[0]:.3f} arcsec")

    # ----- valid footprint; zeros/NaN outside it -> NaN --------------------------------
    good = np.isfinite(flux) & (flux != 0)
    valid = good.mean(axis=0) > 0.5
    if valid.sum() == 0:
        raise ValueError(f"{path}: no spaxel has valid (finite, non-zero) data")
    flux[:, ~valid] = np.nan
    flux[~np.isfinite(flux)] = np.nan
    info.append(f"valid footprint: {int(valid.sum())} of {ny*nx} spaxels")

    # ----- optional variance cube ------------------------------------------------------
    var = None
    if var_path:
        with fits.open(os.path.expanduser(var_path), memmap=False) as vh:
            vhdu = next((h for h in vh if h.data is not None and h.data.ndim == 3), None)
            var = np.asarray(vhdu.data, dtype=np.float32)
        if var.shape != flux.shape:
            raise ValueError(f"{name} variance cube shape {var.shape} != flux cube shape {flux.shape}")
        var[:, ~valid] = np.nan
        var[~np.isfinite(var) | (var <= 0)] = np.nan
        info.append("variance cube loaded")

    # ----- bandpass trim -----------------------------------------------------------------
    wmin = float(cfg.get(f"{name}_wave_min", 0.0) or 0.0)
    wmax = float(cfg.get(f"{name}_wave_max", 0.0) or 0.0)
    src = "user"
    if wmin <= 0 and "WAVGOOD0" in header:
        wmin = float(header["WAVGOOD0"]); src = "header WAVGOOD0"
    if wmax <= 0 and "WAVGOOD1" in header:
        wmax = float(header["WAVGOOD1"])
    # WAVGOOD* are in the file's native frame; compare in native frame
    wn = wave_native
    sel = np.ones(nw, bool)
    if wmin > 0:
        sel &= wn >= wmin
    if wmax > 0:
        sel &= wn <= wmax
    # drop channels with no data in the footprint
    chan_ok = np.isfinite(flux[:, valid]).mean(axis=1) > 0.5
    sel &= chan_ok
    if sel.sum() < 50:
        raise ValueError(f"{name}: fewer than 50 usable wavelength channels after trimming")
    if not sel.all():
        flux = flux[sel]
        if var is not None:
            var = var[sel]
        wave_air = wave_air[sel]
        wave_native = wave_native[sel]
        info.append(f"trimmed to {wave_native[0]:.1f}-{wave_native[-1]:.1f} A ({int(sel.sum())} channels; limits from {src})")

    # ----- units: surface brightness -> flux per spaxel -------------------------------------
    unit_scale = 1.0
    bunit = str(header.get("BUNIT", "")).strip()
    if cfg.get("convert_sb_to_flux", True) and "arcsec" in bunit.lower():
        area = pixscale[0] * pixscale[1]
        flux *= np.float32(area)
        if var is not None:
            var *= np.float32(area ** 2)
        unit_scale = area
        info.append(f"BUNIT '{bunit}' is per arcsec^2 -> multiplied by spaxel area {area:.4f} arcsec^2")

    # ----- instrumental resolution -------------------------------------------------------------
    grating, g_how = detect_grating(header, path, name)
    if grating is None and sibling_header is not None:
        grating, g_how = detect_grating(sibling_header, "", name)
        g_how = "sibling cube " + g_how
    slicer, s_how = detect_slicer(header, path)
    if slicer is None and sibling_header is not None:
        slicer, s_how = detect_slicer(sibling_header, "")
    g_cfg = cfg.get(f"{name}_grating", "auto")
    if g_cfg != "auto":
        grating, g_how = g_cfg, "set by user"
    if cfg.get("slicer", "auto") != "auto":
        slicer, s_how = cfg["slicer"], "set by user"
    assumed = []
    if grating is None:
        grating = K.DEFAULT_GRATING[name]; assumed.append(f"grating->{grating}")
    if slicer is None:
        slicer = K.DEFAULT_SLICER; assumed.append(f"slicer->{slicer}")
    R = float(K.GRATING_R_TABLE[grating][slicer])
    lam_c = 0.5 * (wave_air.min() + wave_air.max())
    fwhm = lam_c / R
    fwhm_user = float(cfg.get(f"{name}_fwhm_A", 0.0) or 0.0)
    if fwhm_user > 0:
        fwhm = fwhm_user
        info.append(f"instrumental FWHM = {fwhm:.3f} A (user value)")
    else:
        info.append(f"instrumental FWHM = {lam_c:.0f}/{R:.0f} = {fwhm:.3f} A  "
                    f"(grating {grating} [{g_how}], slicer {slicer} [{s_how}]"
                    + (f"; ASSUMED defaults: {', '.join(assumed)}" if assumed else "") + ")")

    return Arm(name=name, path=path, flux=flux, var=var, header=header, wcs=wcs2d, wave=wave_air,
               native_frame=frame, fwhm_A=float(fwhm), pixscale=pixscale, valid=valid,
               unit_scale=unit_scale, info=info, grating=grating, slicer=slicer,
               resolving_power=R, vcorr_header=(float(vcorr_header) if vcorr_header is not None else None),
               vcorr_applied_kms=vcorr_applied)


def white_light(arm: Arm, wmin: Optional[float] = None, wmax: Optional[float] = None,
                trim_frac: float = 0.02) -> np.ndarray:
    """Median-over-wavelength image (NaN outside the footprint)."""
    w = arm.wave
    lo = w.min() + trim_frac * (w.max() - w.min()) if wmin is None else wmin
    hi = w.max() - trim_frac * (w.max() - w.min()) if wmax is None else wmax
    sel = (w >= lo) & (w <= hi)
    if sel.sum() == 0:
        sel = np.ones_like(w, bool)
    import warnings
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        img = np.nanmedian(arm.flux[sel], axis=0)
    img[~arm.valid] = np.nan
    return img
