"""Single source of truth for every user-facing setting.

The GUI builds its input forms from :data:`FIELDS`, the command line loads/saves the same
keys from JSON, and the docs are generated from it -- so a default or description can
never drift between the three.  Each :class:`FieldSpec` carries a physical explanation of
*why* the default is what it is.
"""
from __future__ import annotations

import copy
import json
import os
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Sequence

from .constants import (BLUE_GRATINGS, LIBRARY_CHOICES, RED_GRATINGS, SLICER_NAMES)

STAGES = [
    ("inputs", "1. Data & instrument"),
    ("noise", "2. Noise & sky"),
    ("apertures", "3. Apertures"),
    ("extraction", "4. Extraction"),
    ("fitting", "5. Redshift & kinematics"),
    ("output", "6. Output"),
]


@dataclass
class FieldSpec:
    key: str
    label: str
    kind: str                 # 'str','path','dir','int','float','bool','choice','text'
    default: Any
    help: str
    stage: str
    choices: Optional[Sequence[str]] = None
    group: str = ""           # sub-heading inside a stage
    advanced: bool = False    # rendered in a collapsible "advanced" block
    minimum: Optional[float] = None
    maximum: Optional[float] = None
    unit: str = ""
    blank_default: bool = False   # numeric field shown EMPTY (not 0) when it holds its 'use the default' value
    filetypes: Optional[Sequence[str]] = None  # for 'path' fields


F = FieldSpec
FIELDS: List[FieldSpec] = [
    # ------------------------------------------------------------------ 1. DATA
    F("target_name", "Target name", "str", "lens_system",
      "Label used in output file names and plot titles (e.g. 'AGEL012441').", "inputs"),
    F("blue_cube", "Blue (KCWI) intensity cube", "path", "",
      "Reduced KCWI data cube (FITS, 3-D: wavelength x y x x). Leave empty to use the red arm only. "
      "A cube from the Keck DRP '_icubes' (intensity) product is expected; NaNs or all-zero "
      "spaxels are treated as 'no data'.", "inputs", filetypes=["*.fits", "*.fits.gz"]),
    F("red_cube", "Red (KCRM) intensity cube", "path", "",
      "Reduced KCRM data cube. Leave empty to use the blue arm only. Both arms are fitted "
      "jointly (one redshift, one velocity dispersion) when both are given.", "inputs",
      filetypes=["*.fits", "*.fits.gz"]),
    F("blue_var_cube", "Blue variance cube (optional)", "path", "",
      "Variance cube ('_vcubes') matching the blue flux cube. If given it is used as the noise "
      "model; otherwise the noise is measured from a source-free sky box (stage 2).", "inputs",
      advanced=True, filetypes=["*.fits", "*.fits.gz"]),
    F("red_var_cube", "Red variance cube (optional)", "path", "",
      "Variance cube matching the red flux cube (see above).", "inputs", advanced=True,
      filetypes=["*.fits", "*.fits.gz"]),
    F("output_dir", "Output directory", "dir", "kcwilensfit_output",
      "Everything the pipeline writes (spectra, plots, tables, the settings actually used) goes here.",
      "inputs"),

    # ---------------------------------------------- 1b. INSTRUMENT (same tab as Data)
    F("blue_wave_frame", "Blue wavelength frame", "choice", "auto",
      "Whether the blue cube's wavelength axis is in AIR or VACUUM. 'auto' reads CTYPE3 "
      "(AWAV = air, WAVE = vacuum) and assumes air if absent. Getting this wrong shifts every "
      "line by ~85 km/s (dz ~ 3e-4) and goes straight into the redshift.", "inputs",
      choices=["auto", "air", "vacuum"], group="Wavelength frame"),
    F("red_wave_frame", "Red wavelength frame", "choice", "auto",
      "Same for the red cube. 'auto' reads CTYPE3 and assumes VACUUM if absent (the KCRM DRP "
      "convention).", "inputs", choices=["auto", "air", "vacuum"], group="Wavelength frame"),
    F("velocity_correction", "Velocity-frame correction", "choice", "none",
      "'none': use wavelengths as they are in the file (the Keck DRP normally already writes "
      "heliocentric-corrected cubes; check the VCORR/VCORRTYP keywords printed in the log). "
      "'header': ADD the header VCORR on top (only if the file is NOT already corrected). "
      "'manual': add the value set under Advanced. Double-correcting shifts z by ~1e-4.", "inputs",
      choices=["none", "header", "manual"], group="Wavelength frame"),
    F("manual_vcorr_kms", "Manual velocity correction [km/s]", "float", 0.0,
      "Used only if the correction mode is 'manual'. Positive values make redshifts larger.",
      "inputs", unit="km/s", group="Wavelength frame"),
    F("blue_wave_min", "Blue min wavelength [A]", "float", 0.0,
      "Trim the blue cube to this observed wavelength. 0 = automatic (header WAVGOOD0 if present, "
      "else the first valid channel). Edges of the bandpass are low-throughput and noisy.",
      "inputs", unit="A", group="Usable wavelength range", advanced=True),
    F("blue_wave_max", "Blue max wavelength [A]", "float", 0.0,
      "As above for the red end of the blue cube. 0 = automatic. The dichroic cross-over "
      "region near 5600 A is usually the least reliable part of the blue arm.",
      "inputs", unit="A", group="Usable wavelength range", advanced=True),
    F("red_wave_min", "Red min wavelength [A]", "float", 0.0,
      "Trim the red cube (0 = header WAVGOOD0 if present, else first valid channel).",
      "inputs", unit="A", group="Usable wavelength range", advanced=True),
    F("red_wave_max", "Red max wavelength [A]", "float", 0.0,
      "Trim the red cube (0 = header WAVGOOD1 if present, else last valid channel).",
      "inputs", unit="A", group="Usable wavelength range", advanced=True),
    F("blue_grating", "Blue grating", "choice", "auto", "KCWI grating. 'auto' reads the header / file "
      "name; if neither identifies it, BL is assumed.", "inputs", choices=["auto", *BLUE_GRATINGS],
      group="Spectral resolution"),
    F("red_grating", "Red grating", "choice", "auto", "KCRM grating ('auto' as above; RL assumed if "
      "undetected).", "inputs", choices=["auto", *RED_GRATINGS], group="Spectral resolution"),
    F("slicer", "Image slicer", "choice", "auto", "IFU slicer (Large/Medium/Small). Sets the "
      "resolving power R together with the grating. 'auto' reads IFUNAM; Medium assumed otherwise.",
      "inputs", choices=["auto", *SLICER_NAMES], group="Spectral resolution"),
    F("blue_fwhm_A", "Blue instrumental FWHM [A]", "float", 0.0,
      "Instrumental line-spread FWHM of the blue arm at its central wavelength. LEAVE EMPTY to use the "
      "nominal value lambda_c/R from the grating+slicer table. Enter a measured value (arc lamp / "
      "sky lines) if you have one -- it directly sets the velocity-dispersion zero point "
      "(sigma_inst = FWHM/2.355).", "inputs", unit="A", group="Spectral resolution", blank_default=True),
    F("red_fwhm_A", "Red instrumental FWHM [A]", "float", 0.0,
      "As above for the red arm (empty = nominal table value).", "inputs", unit="A", blank_default=True,
      group="Spectral resolution"),
    F("convert_sb_to_flux", "Convert surface brightness to flux per spaxel", "bool", True,
      "If the cube BUNIT is per arcsec^2 (the DRP 'icubes'), multiply by the spaxel area so "
      "extracted spectra are true fluxes (1e-16 erg/s/cm2/A) and summing spaxels is physically "
      "meaningful. Has no effect on redshifts.", "inputs", group="Units"),
    F("arm_registration", "Blue/red spatial registration", "choice", "wcs",
      "How a spaxel in the blue cube maps to the red cube. 'identity': same (x,y) index -- correct "
      "when both arms were resampled to one grid from the same IFU exposure (the KCWI+KCRM design: "
      "one telescope focus, one image slicer, split by a dichroic). 'wcs': follow the header WCS. "
      "'xcorr': measure from white-light images. 'manual': use dx, dy (under Advanced). Header WCS "
      "zero-points can be off by arcsec, so verify in the registration check (Data tab) before trusting 'wcs'.",
      "inputs", choices=["identity", "wcs", "xcorr", "manual"], group="Registration"),
    F("arm_dx", "Manual dx [pix]", "float", 0.0, "Red x = blue x + dx (only for 'manual').",
      "inputs", unit="pix", group="Registration", advanced=True),
    F("arm_dy", "Manual dy [pix]", "float", 0.0, "Red y = blue y + dy (only for 'manual').",
      "inputs", unit="pix", group="Registration", advanced=True),

    # ------------------------------------------------------------------ 3. NOISE
    F("sky_box_mode", "Sky-region selection", "choice", "auto",
      "'auto' picks the faintest, flattest source-free box inside the valid field of view. "
      "'manual' uses the box given below. The box supplies (i) the residual background level to "
      "subtract and (ii) the empirical noise -- it must contain NO source or lens light.",
      "noise", choices=["auto", "manual"], group="Sky box"),
    F("sky_box_size", "Auto sky-box size [pix]", "int", 12,
      "Side of the square box searched for in 'auto' mode. Bigger = better statistics but harder "
      "to keep source-free in a small field.", "noise", unit="pix", minimum=4, group="Sky box"),
    F("sky_box_blue", "Manual sky box, blue (y0,y1,x0,x1)", "str", "",
      "Pixel limits, 0-indexed, end-exclusive, as 'y0,y1,x0,x1'. Ignored in 'auto' mode.", "noise",
      group="Sky box"),
    F("sky_box_red", "Manual sky box, red (y0,y1,x0,x1)", "str", "",
      "Leave empty to reuse the blue box for the red arm (valid with identity registration).",
      "noise", group="Sky box"),
    F("subtract_background", "Subtract residual background", "bool", True,
      "Subtract the sky-box median spectrum from every spaxel. Removes small sky-subtraction "
      "residuals common to the whole field. Turn OFF if the box might contain diffuse lens or "
      "cluster light (it would then over-subtract real signal).", "noise", group="Background"),
    F("noise_source", "Noise model", "choice", "auto",
      "'sky_box': per-wavelength scatter of the sky box (assumes spatially uniform noise). "
      "'variance_cube': use the supplied variance cubes. 'auto': variance cubes if provided, "
      "else sky box. Neither contains source photon noise, so bright apertures are slightly "
      "under-estimated -- pPXF errors are rescaled by sqrt(chi2/dof) to compensate.", "noise",
      choices=["auto", "sky_box", "variance_cube"], group="Noise"),
    F("correct_noise_correlation", "Correct for correlated spaxel noise", "bool", True,
      "Resampling the KCWI slices onto a square grid correlates neighbouring spaxels, so the noise "
      "of a SUM over N spaxels grows faster than sqrt(N). The pipeline measures this from block "
      "sums in the sky box (sigma_N = sigma_1 N^alpha, alpha in [0.5,1]) and inflates aperture "
      "errors accordingly. Skipping it makes every error too optimistic.", "noise",
      group="Noise"),

    # ------------------------------------------------------------- 5. EXTRACTION
    F("sky_regions_use_wcs", "Map sky-coordinate regions with each cube's WCS", "bool", True,
      "Regions loaded from a DS9 file in fk5/icrs (WCS) coordinates are converted to pixels "
      "separately for the blue and the red cube using EACH cube's own header WCS (centre, size and "
      "rotation), and the registered shapes actually used are drawn on every plot. Turn OFF to "
      "treat the regions as defined once on the blue grid and shift them with the blue/red "
      "registration (Data tab) instead (safer if the header WCS zero-points disagree).",
      "apertures"),
    F("extraction_method", "Extraction method", "choice", "deblend",
      "'deblend': simultaneous spatial-profile (PSF-convolved Sersic / point / arc-mask) fit of ALL "
      "apertures at every wavelength -- spectral deblending, the right tool when a deflector and "
      "its arcs (or two lens galaxies) overlap on the sky. Errors include the blending "
      "covariance. 'sum': plain (optionally flux-weighted) aperture sums with no contamination "
      "correction -- adequate only for well-separated, isolated objects.", "extraction",
      choices=["deblend", "sum"]),
    F("weighting", "Weighting for 'sum' extraction", "choice", "uniform",
      "'uniform': plain sum. 'intensity': weight spaxels by the white-light flux (the "
      "S/N-optimal, co-add-then-fit aperture of Cappellari+2006; suppresses faint contaminating "
      "arc light in a deflector aperture).", "extraction", choices=["uniform", "intensity"]),
    F("seeing_fwhm_arcsec", "Seeing / PSF FWHM [arcsec]", "float", 1.0,
      "Gaussian PSF FWHM used to build the spatial profiles for deblending. KCWI is seeing-limited "
      "(typically 0.7-1.3\"); take it from the header or a field star. The Sersic fit "
      "is degenerate with the PSF, so a wrong value mainly biases the profile r_e.", "extraction",
      unit="arcsec", minimum=0.1, maximum=5.0),
    F("profile_band_A", "Band used to fit profiles [A]", "str", "",
      "Observed wavelength range 'lo-hi' of the image used to fit the spatial profiles. Empty = "
      "the redder half of the red arm (or blue arm if only blue is given), where the deflector "
      "dominates and the blue star-forming arcs are faintest.", "extraction"),
    F("deblend_pad_pix", "Deblend window padding [pix]", "int", 8,
      "Margin added around the aperture bounding box for the spatial fit region. The "
      "profile wings outside the apertures are what pin down the underlying lens light.",
      "extraction", unit="pix", minimum=2, advanced=True),
    F("deblend_fit_background", "Fit a flat background in deblend", "bool", False,
      "Add a free constant to the deblend model at every wavelength. Usually unnecessary after "
      "the sky-box subtraction and it can be degenerate with extended lens light.",
      "extraction", advanced=True),
    F("combine_groups", "Combine images of the same source", "bool", True,
      "Apertures that share a 'group' label (multiple images of ONE lensed source) are "
      "inverse-variance combined into an extra spectrum '<group>' for a higher-S/N redshift. "
      "Valid because lensing is achromatic -- all images share one redshift -- but only label "
      "apertures with the same group when you are confident they are images of the same source.",
      "extraction"),
    F("manual_subtractions", "Manual subtraction of one region from another", "text", "",
      "Remove light you KNOW belongs to another region, one rule per line: 'target: contaminator "
      "[scale]'. Example: 'imgA: deflector 0.3' subtracts 0.3 x the deflector's spectrum from imgA's; "
      "'imgA: deflector auto' (or no scale) estimates the fraction from the contaminator's fitted "
      "spatial profile at the target aperture. Works with both extraction methods (with 'deblend' the "
      "overlap is already modelled, so use it only for residuals). The subtracted spectrum's noise is "
      "added to the target's variance. Uses the contaminator's spectrum as extracted, i.e. before "
      "any other manual subtraction.", "extraction", group="Contamination"),
    F("reject_spikes", "Reject single-pixel spikes", "bool", True,
      "Mask isolated one-pixel positive outliers (cosmic-ray / hot-pixel residuals) in the "
      "extracted spectra. Real emission lines are resolved by the instrument (>= ~2 pixels wide), "
      "so a lone bright pixel cannot be a line.", "extraction"),
    F("display_smooth_A", "Display smoothing [A]", "float", 0.0,
      "Gaussian smoothing (sigma, Angstrom) applied to plots ONLY. Fits always use unsmoothed "
      "data: smoothing correlates neighbouring pixels and fakes S/N.", "extraction", unit="A",
      minimum=0.0),

    # ---------------------------------------------------------------- 6. FITTING
    F("fit_arms", "Arms to fit", "choice", "both",
      "'both': blue+red fitted jointly with a single redshift/sigma and one free relative flux "
      "scale for red (only the ratio between arms is identifiable). Otherwise a single arm.",
      "fitting", choices=["both", "blue", "red"], group="Blue / red arms"),
    F("fit_red_scale", "Fit relative red/blue flux scale", "bool", True,
      "(Only when both arms are fitted.) Fit one free multiplicative factor for the red arm relative to the blue (iterated "
      "in closed form). Corrects flux-calibration differences between the arms.", "fitting",
      group="Blue / red arms"),
    F("crosscheck_arms", "Blue-only vs red-only consistency", "bool", True,
      "(Only when both arms are fitted.) Re-fit each arm alone: independent determinations from the two cameras (different "
      "detectors, calibrations and systematics) must agree on z.", "fitting", group="Blue / red arms"),
    F("mask_airglow", "Mask night-sky emission lines", "bool", True,
      "[O I] 5577/6300/6364 and Na D airglow: residuals here are atmospheric, not astrophysical "
      "at any redshift.", "fitting", group="Spectral masks"),
    F("mask_oh_forest", "Mask the OH forest (> 8750 A)", "bool", True,
      "The crowded red OH Meinel bands, at the very edge of the usable KCRM range.", "fitting",
      group="Spectral masks"),
    F("mask_light_pollution", "Mask Hg/Na light-pollution lines", "bool", True,
      "Street-lamp lines; weak on Maunakea but cheap insurance.", "fitting", group="Spectral masks"),
    F("mask_telluric", "Mask telluric absorption bands", "bool", True,
      "O2 B (6860-6890), H2O (7150-7350) and O2 A (7580-7700) A. No telluric correction is "
      "applied, so these bands are excluded from fits.", "fitting", group="Spectral masks"),
    F("extra_mask_windows", "Extra masked windows [A]", "text", "",
      "Additional observed-frame windows to exclude, one 'lo-hi' per line (e.g. a cosmic-ray "
      "residual or a bad column).", "fitting", group="Spectral masks"),
    F("mask_other_planes_lines", "Mask other planes' lines in each fit", "bool", True,
      "Compound-lens essential: once a redshift is secured for ANY aperture, its emission lines "
      "are masked in the continuum fits of the other apertures (e.g. a z=1.3 source's [O II]/Mg II "
      "emission lands inside the deflector's absorption window and biases sigma if left in). Never "
      "applied to a redshift SEARCH (a region may itself lie at another plane's redshift), and never "
      "to a region that turns out to be at the same redshift (see the same-plane tolerance).", "fitting",
      group="Spectral masks"),
    F("other_line_mask_kms", "Line-mask half-width [km/s]", "float", 600.0,
      "Half-width of the mask around each masked emission line (covers line width + velocity "
      "uncertainty of the redshift).", "fitting", unit="km/s", minimum=100.0, group="Spectral masks",
      advanced=True),
    F("same_plane_tol_kms", "Same-plane tolerance [km/s]", "float", 1500.0,
      "When fitting a region at (or near) redshift z, other regions whose redshift lies within this "
      "velocity of z are treated as the SAME source plane and their emission lines are NOT masked "
      "(they would otherwise erase this region's own lines). Also marks ungrouped regions with "
      "agreeing redshifts as 'possibly one source' in the report.", "fitting", unit="km/s", minimum=100.0,
      group="Spectral masks"),
    F("sps_library", "Stellar template library", "choice", "emiles",
      "Template library for continuum fits. 'emiles' (R~2.5 A FWHM, air) is the best match to KCWI "
      "resolution for old stellar populations; 'fsps'/'galaxev'/'xsl' are alternatives; 'kc96_*' "
      "are empirical galaxy spectra (starburst/spiral/elliptical, 5 A) useful for blue "
      "star-forming sources; 'custom' uses your own pPXF-format .npz. Downloaded on first use.",
      "fitting", choices=LIBRARY_CHOICES, group="Continuum fit (pPXF)"),
    F("source_sps_library", "Template library for sources", "choice", "kc96_starburst",
      "Library used for continuum fits of 'source'/'other' apertures. Lensed sources are usually "
      "blue star-forming galaxies whose spectra are better matched by empirical starburst templates "
      "than by old-population SPS models.", "fitting", choices=LIBRARY_CHOICES,
      group="Continuum fit (pPXF)"),
    F("custom_template_path", "Custom template file (.npz)", "path", "",
      "pPXF sps_lib-format .npz (templates, lam, fwhm, ages, metals, masses). Used if the library is "
      "'custom'.", "fitting", group="Continuum fit (pPXF)", advanced=True, filetypes=["*.npz"]),
    F("custom_template_frame", "Custom template frame", "choice", "air",
      "Wavelength frame of the custom file.", "fitting", choices=["air", "vacuum"],
      group="Continuum fit (pPXF)", advanced=True),
    F("age_range", "Template age range [Gyr]", "str", "0.5-14.0",
      "Restrict SPS templates to ages 'lo-hi' Gyr (ignored for the empirical kc96 libraries). "
      "Rarely needs changing: the default drops very young stars, which are degenerate with the "
      "continuum polynomial for a passive deflector; widen it only for a star-forming deflector.", "fitting", group="Continuum fit (pPXF)"),
    F("metal_range", "Template [M/H] range", "str", "-1.0-0.3",
      "Restrict SPS templates to [M/H] 'lo-hi' (ignored for the empirical kc96 libraries). Rarely needs "
      "changing: massive lens galaxies are metal-rich and this removes only very sub-solar templates "
      "that add degeneracy.", "fitting", group="Continuum fit (pPXF)"),
    F("moments", "Kinematic moments", "choice", "2",
      "2 = (V, sigma). 4 adds Gauss-Hermite h3, h4 (needs S/N >~ 30 per pixel to be useful).",
      "fitting", choices=["2", "4"], group="Continuum fit (pPXF)"),
    F("mdegree", "Multiplicative polynomial degree", "int", 8,
      "Order of the multiplicative polynomial. It absorbs the blue/red flux-calibration mismatch "
      "and slow reddening/continuum mismatch WITHOUT distorting line depths (unlike additive "
      "polynomials), so it is the safer choice when sigma matters. Raise it for a very wide "
      "wavelength range.", "fitting", minimum=0, group="Continuum fit (pPXF)"),
    F("degree", "Additive polynomial degree", "int", -1,
      "Order of the additive polynomial (-1 = none). Additive terms can fake or hide line depth; "
      "use sparingly and test the effect with the robustness sweep.", "fitting", minimum=-1,
      group="Continuum fit (pPXF)"),
    F("sigma_init", "Initial sigma [km/s]", "float", 200.0,
      "Starting velocity dispersion for the non-linear fit.", "fitting", unit="km/s", minimum=20.0,
      group="Continuum fit (pPXF)", advanced=True),
    F("vel_bound_kms", "Velocity bound [km/s]", "float", 1500.0,
      "|V| is restricted to +/- this about the chosen redshift guess. A result sitting on the bound "
      "is flagged as unconverged.", "fitting", unit="km/s", group="Continuum fit (pPXF)",
      advanced=True),
    F("mask_nebular", "Mask nebular emission in continuum fits", "choice", "forbidden",
      "Lines at the aperture's OWN redshift that stellar templates cannot reproduce. 'forbidden': "
      "[O II], [Ne III], [O III], [N II], [S II] (Balmer lines are kept because for a passive "
      "galaxy they are stellar ABSORPTION that constrains sigma). 'all': also Balmer. 'none'.",
      "fitting", choices=["forbidden", "all", "none"], group="Continuum fit (pPXF)"),
    F("fit_rest_range", "Rest-frame fitting window for deflectors [A]", "str", "3540-8950",
      "Rest-frame wavelength range used in stellar-continuum fits of 'deflector' apertures (empty = "
      "everything). The default is the range where E-MILES has its native ~2.5 A resolution; in the "
      "rest-UV the libraries are much coarser than KCWI, so including it forces the galaxy to be "
      "degraded and makes sigma depend on the template. Redshifts are insensitive to this choice; "
      "velocity dispersions are not.", "fitting", group="Continuum fit (pPXF)"),
    F("deflector_z_range", "Deflector z search range", "str", "0.05-1.2",
      "Redshift range 'lo-hi' scanned for aperture role 'deflector'. Strong-lens deflectors are "
      "massive early-type galaxies, mostly 0.1 < z < 1; Ca H&K must lie in the data (z < 1.27 "
      "with KCRM to 8940 A).", "fitting", group="Redshift search"),
    F("source_z_range", "Source z search range", "str", "0.2-3.5",
      "Redshift range scanned for aperture role 'source'. Physically a lensed source must be "
      "BEHIND the lens, so results at or below a measured deflector redshift are flagged as "
      "impossible.", "fitting", group="Redshift search"),
    F("z_step_continuum", "z step, continuum scan", "float", 0.002,
      "Spacing of the chi2(z) scan (velocity dispersion held fixed, continuum removed by high-pass "
      "filtering). Must be smaller than the width of the narrowest strong features (sigma ~ "
      "250 km/s -> dz ~ 0.0015 at z~0.7); the best minima are then refined with full pPXF fits.",
      "fitting", minimum=1e-4, group="Redshift search"),
    F("z_step_lines", "z step, line search", "float", 0.0005,
      "Spacing of the matched-filter line scan. Lines are ~150-300 km/s wide, so steps of "
      "dz ~5e-4 (150 km/s) sample them adequately.", "fitting", minimum=1e-5, group="Redshift search"),
    F("scan_sigma_kms", "sigma(s) held fixed in the scan [km/s]", "str", "150, 350",
      "Velocity dispersion(s) assumed while scanning for the redshift. Give ONE value for a single "
      "scan, or several (default '150, 350') to run a scan for each: a candidate present at every "
      "value is stable against the assumption, one that appears only for some is flagged. Holding sigma "
      "fixed stops a noisy fit from smearing the template over the noise; a wrong value mainly "
      "reduces contrast and rarely shifts the minimum. Massive lens galaxies span roughly 150-450 km/s.",
      "fitting", group="Redshift search"),
    F("line_sets", "Line search uses", "choice", "both",
      "Which kinds of features the blind line search scores: 'emission' (nebular + resonant), "
      "'absorption' (stellar + ISM) or 'both'. Source arcs of star-forming galaxies often show "
      "only one or the other.", "fitting", choices=["both", "emission", "absorption"],
      group="Line search"),
    F("line_sigma_kms", "Intrinsic line width [km/s]", "float", 100.0,
      "Assumed intrinsic Gaussian sigma of searched lines, added in quadrature to the "
      "instrumental width. 100 km/s suits compact star-forming galaxies; use 250+ for AGN or "
      "massive galaxies.", "fitting", unit="km/s", minimum=10.0, group="Line search"),
    F("line_continuum_window_A", "Continuum filter width [A]", "float", 150.0,
      "Width of the running-median filter that removes the continuum before the matched-filter "
      "search. Must exceed ~3x the widest feature of interest.", "fitting", unit="A",
      minimum=30.0, group="Line search", advanced=True),
    F("detect_snr", "Per-line detection threshold [sigma]", "float", 3.0,
      "A line counts as 'detected' when its fitted amplitude is at least this many sigma.",
      "fitting", minimum=1.0, group="Line search"),
    F("n_null_sims", "Noise simulations for significance", "int", 20,
      "Number of sign-flipped noise realisations of the whole redshift scan used to measure how large "
      "the best score gets in pure noise (look-elsewhere correction: a scan of ~10^4 redshifts x ~45 "
      "lines produces chance 3-sigma peaks). A candidate must beat the 95th percentile of these maxima "
      "to be called significant. 0 disables the test.", "fitting", minimum=0, group="Line search"),
    F("min_lines_secure", "Lines needed for a secure z", "int", 2,
      "Independent line features that must be detected before a redshift is called 'secure'. "
      "A single line is never proof of a redshift (an isolated emission feature can be [O II], "
      "Mg II, H-beta, [O III], H-alpha...); the report lists the competing identifications.",
      "fitting", minimum=1, group="Line search"),
    F("n_bootstrap", "Bootstrap realisations", "int", 100,
      "Number of re-fits used to measure the redshift and sigma uncertainties (0 = use pPXF's "
      "formal error only, which is usually optimistic). 100 gives ~7% precision on the error; "
      "use 500+ for publication numbers.", "fitting", minimum=0, group="Uncertainties"),
    F("bootstrap_method", "Bootstrap method", "choice", "mc",
      "'mc' (default): re-draw Gaussian noise from the noise array around the best-fit model and "
      "re-fit -- the error then follows the measured noise level. 'wild': multiply the fit residuals "
      "by random +/-1 (Rademacher), scaled by their local scatter -- also captures template "
      "mismatch, but shrinks if the noise array is over-estimated.", "fitting", choices=["wild", "mc"], group="Uncertainties"),
    F("sweep_degrees", "Polynomial-degree sweep", "str", "6,8,10,12",
      "Comma list of multiplicative degrees re-fit after the main fit; the scatter in V and sigma "
      "is reported as a SYSTEMATIC error (empty = off). Continuum modelling is usually the "
      "largest systematic in a sigma measurement.", "fitting", group="Robustness"),
    F("sweep_libraries", "Template-library sweep", "str", "",
      "Comma list of libraries (e.g. 'emiles,fsps') re-fit with the same settings; the scatter "
      "in V and sigma is added to the systematic budget (empty = off). Each library is "
      "compared in its own native wavelength frame.", "fitting", group="Robustness"),
    F("snr_min_kinematics", "Min S/N per A for sigma", "float", 10.0,
      "Below this median S/N per Angstrom the velocity dispersion is reported but flagged "
      "unreliable (the redshift can still be good; sigma needs far more S/N).", "fitting",
      unit="1/A", minimum=1.0, group="Robustness"),
    F("auto_accept", "Auto-select best candidate (headless)", "bool", True,
      "In non-interactive runs, adopt the best-scoring candidate redshift that is not "
      "geometrically impossible for its role. In the GUI you choose from the candidate list.",
      "fitting", group="General"),

    # ------------------------------------------------------------------ 7. OUTPUT
    F("save_plots", "Save diagnostic plots", "bool", True, "Write PNG figures for every stage.",
      "output"),
    F("save_fits_bundle", "Save FITS spectra bundle", "bool", True,
      "Write all extracted spectra into one multi-extension FITS file.", "output"),
    F("save_ascii_spectra", "Save ASCII spectra", "bool", True,
      "Write one text/ECSV file per aperture and arm.", "output"),
]

# ---------------------------------------------------------------------------
# What a user actually needs to look at.  Everything NOT listed here is an "advanced" setting: it is
# still honoured from the config file, still documented, and shown behind "Show advanced settings"
# in the GUI, but its default is right for the large majority of KCWI/KCRM lens data, so the main
# forms stay short.  (E.g. the template age/metallicity cuts only matter for SPS libraries and only
# at the few-km/s level for a passive deflector, so they are advanced.)
# ---------------------------------------------------------------------------
ESSENTIAL_KEYS = {
    "target_name", "blue_cube", "red_cube", "output_dir",                                    # data
    "velocity_correction", "blue_fwhm_A", "red_fwhm_A", "arm_registration",                  # instrument
    "sky_box_mode", "sky_box_size", "sky_box_blue", "sky_box_red", "subtract_background",   # noise
    "sky_regions_use_wcs",                                                                   # apertures
    "extraction_method", "seeing_fwhm_arcsec", "combine_groups",                             # extraction
    "fit_arms", "fit_red_scale", "crosscheck_arms", "sps_library", "source_sps_library", "fit_rest_range", "deflector_z_range", "source_z_range", "scan_sigma_kms",
    "min_lines_secure", "n_bootstrap", "mask_other_planes_lines", "extra_mask_windows",     # fitting
    "save_plots",                                                                            # output
}
for _f in FIELDS:
    _f.advanced = _f.key not in ESSENTIAL_KEYS

FIELD_BY_KEY: Dict[str, FieldSpec] = {f.key: f for f in FIELDS}

# ----------------------------------------------------------------------------
# Aperture table
# ----------------------------------------------------------------------------
APERTURE_ROLES = ["deflector", "source", "other"]
APERTURE_SHAPES = ["circle", "ellipse", "box"]
APERTURE_PROFILES = ["auto", "sersic", "psf", "mask"]
APERTURE_MODES = ["auto", "continuum", "lines", "both", "skip"]

APERTURE_COLUMNS = [
    # key, label, default, help
    ("id", "id", "", "Unique name (e.g. 'defl', 'A', 'B', 'src2')."),
    ("role", "role", "source", "deflector = foreground lens galaxy (stellar-continuum fit, kinematics); "
     "source = lensed background image (line/continuum search; must be BEHIND every deflector); "
     "other = companion/interloper (fitted but exempt from the lensing-geometry check)."),
    ("group", "group", "", "Label shared by all images of ONE lensed source (e.g. 'S1'). Different "
     "source planes of a compound lens get different labels."),
    ("shape", "shape", "circle", "circle (radius a), ellipse (semi-axes a, b, angle), box (width a, height b, angle)."),
    ("x", "x [pix]", 0.0, "Centre, 0-indexed pixel in the blue cube (column)."),
    ("y", "y [pix]", 0.0, "Centre, 0-indexed pixel in the blue cube (row)."),
    ("a", "a [pix]", 3.0, "Radius / semi-major axis / box width."),
    ("b", "b [pix]", 3.0, "Semi-minor axis / box height (ignored for circles)."),
    ("angle", "angle [deg]", 0.0, "Position angle, counter-clockwise from +x."),
    ("profile", "profile", "auto", "Spatial model for deblending: sersic (galaxy), psf (point-like image), "
     "mask (arc: the aperture itself). auto: deflector->sersic, others->mask."),
    ("mode", "fit mode", "auto", "auto: deflector->continuum, source->lines; or force continuum / lines / both / skip."),
    ("z_known", "z (if known)", "", "Skip the search and refine around this redshift."),
    ("z_min", "z min", "", "Override the role's search range (lower bound)."),
    ("z_max", "z max", "", "Override the role's search range (upper bound)."),
    ("ra", "RA [deg]", "", "Sky centre (filled automatically from fk5/icrs region files)."),
    ("dec", "Dec [deg]", "", "Sky centre (filled automatically from fk5/icrs region files)."),
    ("a_as", "a [arcsec]", "", "Size on the sky (filled automatically from sky region files)."),
    ("b_as", "b [arcsec]", "", "Second size on the sky (sky region files)."),
]
APERTURE_DEFAULTS = {c[0]: c[2] for c in APERTURE_COLUMNS}


def defaults() -> Dict[str, Any]:
    """Fresh default configuration dictionary (including an empty aperture list)."""
    cfg = {f.key: copy.deepcopy(f.default) for f in FIELDS}
    cfg["apertures"] = []
    cfg["region_file"] = ""
    return cfg


def coerce(spec: FieldSpec, value: Any) -> Any:
    """Convert a GUI/JSON value to the field's Python type (raises ValueError if impossible)."""
    if spec.kind == "bool":
        if isinstance(value, str):
            return value.strip().lower() in ("1", "true", "yes", "y", "on")
        return bool(value)
    if spec.kind == "int":
        return int(float(value))
    if spec.kind == "float":
        return float(value)
    if spec.kind == "choice":
        v = str(value)
        if spec.choices and v not in spec.choices:
            raise ValueError(f"{spec.key}: '{v}' not one of {list(spec.choices)}")
        return v
    return "" if value is None else str(value)


def validate(cfg: Dict[str, Any]) -> List[str]:
    """Return a list of human-readable problems (empty list = OK)."""
    problems = []
    for f in FIELDS:
        if f.key not in cfg:
            continue
        try:
            v = coerce(f, cfg[f.key])
        except (ValueError, TypeError) as e:
            problems.append(f"{f.label}: {e}")
            continue
        if f.kind in ("int", "float"):
            if f.minimum is not None and v < f.minimum:
                problems.append(f"{f.label}: must be >= {f.minimum}")
            if f.maximum is not None and v > f.maximum:
                problems.append(f"{f.label}: must be <= {f.maximum}")
    if not cfg.get("blue_cube") and not cfg.get("red_cube"):
        problems.append("Provide at least one of the blue or red cubes.")
    for key in ("blue_cube", "red_cube", "blue_var_cube", "red_var_cube"):
        p = cfg.get(key)
        if p and not os.path.isfile(os.path.expanduser(p)):
            problems.append(f"{FIELD_BY_KEY[key].label}: file not found: {p}")
    for key in ("deflector_z_range", "source_z_range", "age_range", "metal_range"):
        try:
            parse_range(cfg.get(key, FIELD_BY_KEY[key].default))
        except Exception:
            problems.append(f"{FIELD_BY_KEY[key].label}: expected 'lo-hi' numbers")
    ids = [a.get("id") for a in cfg.get("apertures", [])]
    if len(ids) != len(set(ids)):
        problems.append("Aperture ids must be unique.")
    for a in cfg.get("apertures", []):
        if not str(a.get("id", "")).strip():
            problems.append("Every aperture needs an id.")
        if a.get("role") not in APERTURE_ROLES:
            problems.append(f"Aperture {a.get('id')}: role must be one of {APERTURE_ROLES}")
    return problems


def parse_range(text) -> tuple:
    """'lo-hi' (also 'lo,hi', 'lo to hi', or a 2-sequence) -> (lo, hi) floats, supporting negative lo."""
    if isinstance(text, (list, tuple)):
        lo, hi = float(text[0]), float(text[1])
    else:
        s = str(text).strip().replace(" to ", ",").replace(";", ",")
        if "," in s:
            lo, hi = [float(t) for t in s.split(",")[:2]]
        else:
            # split on a '-' that is not at position 0 and not preceded by 'e'/'E' (exponent)
            idx = None
            for i in range(1, len(s)):
                if s[i] == "-" and s[i - 1] not in "eE":
                    idx = i
                    break
            if idx is None:
                raise ValueError(f"cannot parse range {text!r}")
            lo, hi = float(s[:idx]), float(s[idx + 1:])
    if lo > hi:
        lo, hi = hi, lo
    return lo, hi


def parse_windows(text: str) -> List[tuple]:
    """Multi-line 'lo-hi' text -> list of (lo, hi)."""
    out = []
    for line in str(text or "").replace(";", "\n").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        out.append(parse_range(line))
    return out


def parse_float_list(text: str) -> List[float]:
    return [float(t) for t in str(text or "").replace(";", ",").split(",") if t.strip()]


def parse_name_list(text: str) -> List[str]:
    return [t.strip() for t in str(text or "").replace(";", ",").split(",") if t.strip()]


def load_config(path: str) -> Dict[str, Any]:
    """Load a JSON config, filling anything missing with defaults and coercing types."""
    with open(path) as fh:
        user = json.load(fh)
    cfg = defaults()
    base_dir = os.path.dirname(os.path.abspath(path))
    unknown = [k for k in user if k not in cfg and k != "apertures"]
    for k, v in user.items():
        if k == "apertures":
            cfg["apertures"] = [normalise_aperture(a) for a in v]
        elif k in FIELD_BY_KEY:
            cfg[k] = coerce(FIELD_BY_KEY[k], v)
        elif k == "region_file":
            cfg[k] = v
    # resolve relative paths against the config file's folder
    for key in ("blue_cube", "red_cube", "blue_var_cube", "red_var_cube", "region_file",
                "custom_template_path"):
        if cfg.get(key) and not os.path.isabs(os.path.expanduser(cfg[key])):
            cand = os.path.join(base_dir, cfg[key])
            cfg[key] = cand if os.path.exists(cand) else cfg[key]
    if unknown:
        cfg["_unknown_keys"] = unknown
    return cfg


def save_config(cfg: Dict[str, Any], path: str) -> None:
    out = {k: v for k, v in cfg.items() if not k.startswith("_")}
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)


def normalise_aperture(a: Dict[str, Any]) -> Dict[str, Any]:
    """Fill defaults and coerce numeric aperture fields."""
    out = dict(APERTURE_DEFAULTS)
    out.update(a)
    for k in ("x", "y", "a", "b", "angle"):
        out[k] = float(out[k])
    for k in ("z_known", "z_min", "z_max", "ra", "dec", "a_as", "b_as"):
        v = out.get(k, "")
        out[k] = "" if v in ("", None) else float(v)
    out["id"] = str(out["id"]).strip()
    out["group"] = str(out.get("group", "") or "").strip()
    return out
