"""Synthetic compound-lens KCWI/KCRM cubes with known redshifts, for tests and the ``demo`` command.

The mock mimics the real data products closely enough to exercise every code path:

* blue cube on an **air** axis (``CTYPE3 = AWAV``), red cube on a **vacuum** axis (``WAVE``),
  with a blue/red overlap, 1 A sampling, 0.3" spaxels, a ``BUNIT`` per arcsec^2 and an
  all-zero border outside the IFU footprint (as in real DRP stacks);
* a **double-source-plane (compound) lens**: a massive early-type deflector (z_d), a
  star-forming source at z_1 with two images (emission lines: [O II], Mg II, C III]) and a
  higher-redshift source at z_2 seen through the same field (Lyman-alpha emission plus
  interstellar absorption lines, blueshifted by 150 km/s as in real outflowing galaxies);
* seeing-convolved Sersic / arc-like spatial profiles that overlap, correlated noise with
  sky-line structure, and a known velocity dispersion for the deflector.

Everything is deterministic given ``seed``.
"""
from __future__ import annotations

import json
import os
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import numpy as np
from astropy.io import fits
from scipy.ndimage import gaussian_filter

from .lines import LINES_BY_NAME
from .wavelength import C_KMS, air_to_vac, vac_to_air

FWHM2SIG = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


def _gauss(wave, mu, fwhm_A, amp):
    s = fwhm_A * FWHM2SIG
    return amp * np.exp(-0.5 * ((wave - mu) / s) ** 2)


def _deflector_rest_spectrum(sigma_kms: float, inst_fwhm_A: float = 0.0):
    """High-resolution rest-frame (air) old-population spectrum, broadened by the LOSVD.

    Optical part (> 3540 A) from the X-shooter library (0.8 A FWHM -- much sharper than KCWI, so
    the mock is close to an *intrinsic* spectrum and sigma is recoverable without template-LSF
    double counting); UV part from E-MILES (coarse, only used when the fit window includes it).
    """
    from . import templates as T
    source = "xsl"
    try:
        def pick(name, ages=(8.0, 12.0)):
            d = np.load(T.sps_file(name))
            im = int(np.argmin(abs(d["metals"] - 0.0)))
            ia = [int(np.argmin(abs(d["ages"] - a))) for a in ages]
            return d["lam"], 0.5 * (d["templates"][:, ia[0], im] + d["templates"][:, ia[1], im])
        lx, fx = pick("xsl")
        le, fe = pick("emiles")
        mx = (lx > 3540) & (lx < 9500)
        me = (le > 1700) & (le <= 3540)
        j = (le > 3540) & (le < 3700)
        jx = (lx > 3540) & (lx < 3700)
        scale = np.median(fx[jx]) / np.median(fe[j])
        lam = np.concatenate([le[me], lx[mx]])
        flux = np.concatenate([fe[me] * scale, fx[mx]])
    except Exception:
        lam, flux = T._read_kc96("kc96_elliptical")
        source = "kc96_elliptical"
    flux = flux / np.median(flux[(lam > 4000) & (lam < 5500)])
    lnl = np.linspace(np.log(lam.min()), np.log(lam.max()), lam.size * 2)
    f = np.interp(lnl, np.log(lam), flux)
    dln = lnl[1] - lnl[0]
    f = gaussian_filter(f, (sigma_kms / C_KMS) / dln, mode="nearest")
    return np.exp(lnl), f, source


def _source_continuum(kind: str):
    from . import templates as T
    lam, f = T._read_kc96(kind)
    return lam, f / np.median(f[(lam > 3000) & (lam < 6000)])


def build_truth(z_defl=0.62, z_src1=1.30, z_src2=2.09, sigma_defl=250.0):
    return dict(z_defl=z_defl, z_src1=z_src1, z_src2=z_src2, sigma_defl=sigma_defl)


def make_demo_cubes(outdir: str, seed: int = 7, ny: int = 80, nx: int = 80, truth: Optional[dict] = None,
                    deflector_snr: float = 15.0, arc_scale: float = 0.02, noise_level: float = 1.0,
                    write: bool = True) -> Dict[str, object]:
    """Create ``demo_blue.fits`` / ``demo_red.fits`` (+ truth.json, config) in ``outdir``.

    Returns a dict with the file paths, the truth, and a ready-to-use configuration.
    """
    rng = np.random.default_rng(seed)
    t = build_truth() if truth is None else truth
    pix = 0.3                                         # arcsec / spaxel
    # --- wavelength axes ---------------------------------------------------------------
    wave_b = 3300.0 + np.arange(2400)                 # air, 3300-5699
    wave_r_vac = 5600.0 + np.arange(3300)             # vacuum, 5600-8899
    wave_r_air = vac_to_air(wave_r_vac)
    fwhm_b, fwhm_r = 2.5, 6.5                         # instrumental FWHM [A]

    # --- spatial layout (pixel coordinates, 0-indexed, x = column) --------------------------
    cx, cy = nx / 2 - 0.5, ny / 2 - 0.5
    seeing_sig = 1.0 / pix * FWHM2SIG                 # 1.0" FWHM
    yy, xx = np.mgrid[0:ny, 0:nx]

    def sersic(x0, y0, re, n, q, pa):
        th = np.deg2rad(pa)
        u = (xx - x0) * np.cos(th) + (yy - y0) * np.sin(th)
        v = -(xx - x0) * np.sin(th) + (yy - y0) * np.cos(th)
        r = np.sqrt(u ** 2 + (v / q) ** 2)
        bn = 2 * n - 1 / 3 + 4 / (405 * n)
        return np.exp(-bn * ((r / re) ** (1 / n) - 1))

    def arc(x0, y0, length, width, pa):
        th = np.deg2rad(pa)
        u = (xx - x0) * np.cos(th) + (yy - y0) * np.sin(th)
        v = -(xx - x0) * np.sin(th) + (yy - y0) * np.cos(th)
        return np.exp(-0.5 * (u / length) ** 2 - 0.5 * (v / width) ** 2)

    P_defl = gaussian_filter(sersic(cx, cy, 5.0, 4.0, 0.8, 30.0), seeing_sig)
    P_s1a = gaussian_filter(arc(cx + 8.5, cy + 3.0, 3.0, 1.0, 70.0), seeing_sig)     # image A of source 1
    P_s1b = gaussian_filter(arc(cx - 8.0, cy - 4.5, 3.0, 1.0, 70.0), seeing_sig)     # image B of source 1
    P_s2 = gaussian_filter(arc(cx + 1.0, cy - 8.5, 2.0, 1.2, 0.0), seeing_sig)       # source 2 image
    for P in (P_defl, P_s1a, P_s1b, P_s2):
        P /= P.sum()

    # --- spectra ----------------------------------------------------------------------------
    def inst_smooth(wave, f, fwhm):
        d = np.median(np.diff(wave))
        return gaussian_filter(f, fwhm * FWHM2SIG / d, mode="nearest")

    lam_d, f_d, tsrc = _deflector_rest_spectrum(t["sigma_defl"])

    def deflector_on(wave_air):
        # observed air wavelength -> rest air wavelength (approx. air-air scaling used by pPXF/E-MILES)
        return np.interp(wave_air / (1 + t["z_defl"]), lam_d, f_d, left=0, right=0)

    lam_s1, f_s1 = _source_continuum("kc96_starb4")
    lam_s2, f_s2 = _source_continuum("kc96_starb1")

    def source1_on(wave_air):
        # z_src1: emission lines referenced to the *vacuum* rest wavelengths -> observed air
        cont = np.interp(air_to_vac(wave_air) / (1 + t["z_src1"]), lam_s1, f_s1, left=0, right=0)
        sp = 0.35 * cont
        for name, amp in (("[O II] 3727", 7.0), ("Mg II 2799 (em)", 1.5), ("C III] 1908", 1.8),
                          ("[Ne III] 3869", 1.0)):
            ln = LINES_BY_NAME[name]
            for lam_v, rel in ln.comps:
                mu = float(vac_to_air(lam_v * (1 + t["z_src1"])))
                sp += _gauss(wave_air, mu, np.hypot(6.0, 2.355 * 80.0 / C_KMS * mu), amp * rel * 1.0)
        return sp

    def source2_on(wave_air):
        cont = np.interp(air_to_vac(wave_air) / (1 + t["z_src2"]), lam_s2, f_s2, left=0, right=0)
        sp = 0.5 * cont
        # Lyman-alpha, redshifted by +200 km/s
        z_lya = (1 + t["z_src2"]) * (1 + 200.0 / C_KMS) - 1
        mu = float(vac_to_air(1215.67 * (1 + z_lya)))
        sp += _gauss(wave_air, mu, 6.0, 1.2)
        # interstellar absorption, blueshifted by 150 km/s
        z_abs = (1 + t["z_src2"]) * (1 - 150.0 / C_KMS) - 1
        for name in ("Si II 1260", "O I 1302", "C II 1334", "Si IV 1394", "Si II 1527", "Fe II 1608", "Al II 1670",
                     "Fe II 2344", "Fe II 2382", "Fe II 2600", "Mg II 2799 (abs)"):
            ln = LINES_BY_NAME[name]
            for lam_v, rel in ln.comps:
                mu = float(vac_to_air(lam_v * (1 + z_abs)))
                sp -= _gauss(wave_air, mu, np.hypot(5.0, 2.355 * 100.0 / C_KMS * mu), 0.45 * rel * cont.max() * 0.5)
        return sp

    out = {}
    cubes = {}
    for arm, wave_air, fwhm in (("blue", wave_b, fwhm_b), ("red", wave_r_air, fwhm_r)):
        S = {k: inst_smooth(wave_air, fn(wave_air), fwhm) for k, fn in
             (("d", deflector_on), ("s1", source1_on), ("s2", source2_on))}
        nw = wave_air.size
        # amplitude: deflector normalised so its summed S/N matches the request
        amp_d = 1.0
        cube = (S["d"][:, None, None] * P_defl[None] * amp_d
                + arc_scale * S["s1"][:, None, None] * (P_s1a + 1.4 * P_s1b)[None]
                + arc_scale * 0.5 * S["s2"][:, None, None] * P_s2[None])
        cubes[arm] = (cube, S)
    # noise: per-arm sigma chosen so the deflector aperture-sum S/N per pixel hits the request
    ap = (np.hypot(xx - cx, yy - cy) < 5)
    npix = ap.sum()
    corr = 1.6                                    # noise-correlation inflation of a ~80-spaxel sum
    sigma0 = {}
    for arm, snr_t in (("blue", deflector_snr), ("red", 2.0 * deflector_snr)):
        cube = cubes[arm][0]
        nwa = cube.shape[0]
        sel = slice(int(0.25 * nwa), int(0.75 * nwa))
        ap_flux = np.median(cube[sel][:, ap].sum(axis=1))
        sigma0[arm] = ap_flux / snr_t / (np.sqrt(npix) * corr)
    for arm, wave_air in (("blue", wave_b), ("red", wave_r_air)):
        cube, S = cubes[arm]
        nw = wave_air.size
        shape_noise = np.ones(nw)
        for mu, w_, a_ in ((5577.3, 3.0, 4.0), (6300.3, 3.0, 3.0), (6363.8, 3.0, 1.2), (5890.0, 3.0, 2.5)):
            shape_noise += _gauss(wave_air, mu, w_ * 2.355, a_)
        shape_noise += np.where(wave_air > 7600, (wave_air - 7600) / 1500.0, 0.0)   # redder = noisier
        noise = rng.normal(size=cube.shape)
        noise = gaussian_filter(noise, (0, 0.55, 0.55)) / 0.62                         # covariance between spaxels
        noise *= (noise_level * sigma0[arm] * shape_noise)[:, None, None]
        # footprint: valid disc, zeros outside (real DRP stacks)
        data = (cube + noise).astype(np.float32)
        yv, xv = np.mgrid[0:ny, 0:nx]
        foot = np.hypot(xv - cx, yv - cy) < 0.5 * min(ny, nx) - 4
        data[:, ~foot] = 0.0
        # surface-brightness units: flux = data / area
        bunit = "FLAM16/arcsec^2" if arm == "blue" else "1e-16 erg/s/cm2/arcsec2/Angstrom"
        data = data / np.float32(pix * pix)
        h = fits.Header()
        h["NAXIS"] = 3
        h["CTYPE1"], h["CTYPE2"], h["CTYPE3"] = "RA---TAN", "DEC--TAN", ("AWAV" if arm == "blue" else "WAVE")
        h["CUNIT1"] = h["CUNIT2"] = "deg"; h["CUNIT3"] = "Angstrom"
        h["CRVAL1"], h["CRVAL2"] = 150.0, 2.0
        h["CRPIX1"], h["CRPIX2"] = cx + 1, cy + 1
        h["CD1_1"], h["CD2_2"] = -pix / 3600.0, pix / 3600.0
        h["CD1_2"] = h["CD2_1"] = 0.0
        h["CRVAL3"] = float(wave_b[0] if arm == "blue" else wave_r_vac[0])
        h["CRPIX3"] = 1.0
        h["CD3_3"] = 1.0
        h["BUNIT"] = bunit
        h["RADESYS"] = "ICRS"
        h["INSTRUME"] = "KCWI"
        h["IFUNAM"] = "Medium"
        if arm == "blue":
            h["BGRATNAM"] = "BL"
        else:
            h["RGRATNAM"] = "RL"
            h["BGRATNAM"] = "BL"
            h["WAVGOOD0"] = 5620.0
            h["WAVGOOD1"] = 8890.0
        h["VCORR"] = -12.3
        h["VCORRTYP"] = "heliocentric"
        cubes[arm] = (data, h)
    paths = {}
    if write:
        os.makedirs(outdir, exist_ok=True)
        for arm in ("blue", "red"):
            p = os.path.join(outdir, f"demo_{arm}.fits")
            fits.PrimaryHDU(cubes[arm][0], cubes[arm][1]).writeto(p, overwrite=True)
            paths[arm] = p
        with open(os.path.join(outdir, "truth.json"), "w") as fh:
            json.dump(dict(t, deflector_template=tsrc), fh, indent=2)
    # aperture layout matching the mock geometry (blue pixel coordinates)
    apertures = [
        dict(id="deflector", role="deflector", shape="circle", x=cx, y=cy, a=5.0, b=5.0, profile="sersic"),
        dict(id="A", role="source", group="S1", shape="ellipse", x=cx + 8.5, y=cy + 3.0, a=4.5, b=2.0, angle=70.0),
        dict(id="B", role="source", group="S1", shape="ellipse", x=cx - 8.0, y=cy - 4.5, a=4.5, b=2.0, angle=70.0),
        dict(id="C", role="source", group="S2", shape="ellipse", x=cx + 1.0, y=cy - 8.5, a=3.5, b=2.0, angle=0.0),
    ]
    return dict(paths=paths, truth=t, apertures=apertures, cubes=cubes, center=(cx, cy), template_source=tsrc)
