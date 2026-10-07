"""
pipeline_qc.py — fixes for the defects found in the DESJ0206-0114 audit
(F140W<->F200LP ~65 mas misalignment + cosmic-ray-ridden F200LP mosaic) plus
two found later triaging real batch-processing failures (4b, 2026-10-01).

Each section below is self-contained and maps to one defect. Imported by
hst_reduction.py; nothing here runs on import. See README.md's "Astrometry &
CR-rejection fixes" section for the full writeup of each one.

Defect 1  — astrometric reference preflight + correction
Defect 2  — DQ CR-flag survival assertion
Defect 4  — tweakreg shiftfile NaN check
Defect 4b — tweakreg shiftfile fit-quality check (a fit that converges onto
            noise in a sparse field, rather than failing outright)
Defect 5  — N=2 cosmic-ray rejection (L.A.Cosmic + source protection +
            pairwise output-frame rejection)
Defect 6  — tweakreg skysigma/count-rate mismatch

(Defect 3 — tweakreg dqbits/use_sharp_round/conv_width — is fixed in place
at the TweakReg call sites in hst_reduction.py / config.py; there was no
reusable logic to factor out here.)
"""

from pathlib import Path
import shutil

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS


# ── Defect 1: astrometric reference preflight + correction ─────────────────
#
# WCSNAME alone is not a reliable signal: a "-GSC240"/"-GSC242" suffix can
# mean either "genuinely fit to that catalog" or "distortion-only solution,
# catalog name inherited from the default reference file" -- WCSNAME cannot
# tell the two apart. drizzlepac/stwcs *does* write an unambiguous keyword
# for exactly this purpose: WCSTYPE, e.g.
#   'undistorted a posteriori solution relatively aligned to GSC242'  (real fit)
#   'undistorted a priori solution based on GSC240'                   (no fit)
# Verified against this target's own files: F200LP flc SCI headers carry
# WCSTYPE="...a posteriori solution..." with RMS_RA/RMS_DEC/NMATCHES
# populated; F140W flt SCI headers carry WCSTYPE="...a priori solution..."
# with no HDRLET carrying a FIT_* WCS at all. That's the actual signal to
# key off, not a hardcoded "UV is always better" assumption -- a different
# visit could go the other way if the archive pipeline achieved a real fit
# for IR but not UV.

class AmbiguousAstrometricReference(RuntimeError):
    """Raised when the preflight check cannot safely pick a reference filter."""


def classify_wcs(fits_path, ext=('SCI', 1)):
    """Read the active WCS quality metadata from one exposure's SCI header."""
    with fits.open(fits_path) as hdul:
        hdr = hdul[ext].header
        wcsname = hdr.get('WCSNAME', '')
        wcstype = hdr.get('WCSTYPE', '')
        rms_ra = hdr.get('RMS_RA')
        rms_dec = hdr.get('RMS_DEC')
        nmatch = hdr.get('NMATCHES')
    return {
        'file': str(fits_path),
        'wcsname': wcsname,
        'wcstype': wcstype,
        'is_a_posteriori': 'a posteriori' in str(wcstype).lower(),
        'rms_ra_mas': rms_ra,
        'rms_dec_mas': rms_dec,
        'nmatches': nmatch,
    }


def select_astrometric_reference(ir_exposures, uv_exposures):
    """
    Decide which filter (IR or UV) carries the trustworthy absolute
    astrometry for one target, from the raw exposures' own WCSTYPE.

    Returns ('IR'|'UV', report_dict). Raises AmbiguousAstrometricReference
    if both filters (or neither) show an a posteriori fit -- refusing and
    telling the user beats silently guessing wrong, per the defect this is
    fixing (a wrong guess here is exactly what caused the 65 mas offset).
    """
    ir_checks = [classify_wcs(p) for p in ir_exposures]
    uv_checks = [classify_wcs(p) for p in uv_exposures]
    ir_fit = any(c['is_a_posteriori'] for c in ir_checks)
    uv_fit = any(c['is_a_posteriori'] for c in uv_checks)
    report = {'ir': ir_checks, 'uv': uv_checks}

    if ir_fit and not uv_fit:
        return 'IR', report
    if uv_fit and not ir_fit:
        return 'UV', report
    if ir_fit and uv_fit:
        raise AmbiguousAstrometricReference(
            "Both IR and UV exposures carry an a posteriori WCS fit -- "
            "pick a reference manually by comparing RMS_RA/RMS_DEC/NMATCHES.\n"
            f"IR: {ir_checks}\nUV: {uv_checks}")
    raise AmbiguousAstrometricReference(
        "Neither IR nor UV exposures carry an a posteriori WCS fit -- "
        "absolute astrometry for BOTH filters is untrustworthy here; "
        "cross-filter alignment cannot be resolved automatically.\n"
        f"IR: {ir_checks}\nUV: {uv_checks}")


def measure_filter_offset(reference_sci_path, moving_sci_path,
                          search_radii_arcsec=(0.2, 0.5, 1.0, 2.0, 5.0),
                          min_matches=8, detect_nsigma=8.0):
    """
    Measure the sky-plane offset of `moving`'s astrometry relative to
    `reference`'s via matched-source centroids: detect compact sources
    independently in each filter's *native* mosaic (own pixel scale, own
    PSF, no reprojection involved), cross-match by sky position, and take
    the sigma-clipped median of the matched pairs' (ref - moving) vector.

    Returns (dra_arcsec, ddec_arcsec) -- the vector to ADD to `moving`'s
    CRVAL to bring it into agreement with `reference`.

    Why not whole-frame cross-correlation: an earlier version of this
    function reprojected `moving` onto `reference`'s grid and ran
    skimage.registration.phase_cross_correlation on the full frames. On
    this dataset that measured a ~19 ARCSEC shift -- three orders of
    magnitude too large for what this defect actually is (~65 mas) -- and
    applying it moved F140W's field center *away* from the true target
    position, not toward it. The self-verification check in
    reconcile_filter_astrometry() didn't catch it because both the "before"
    and "after" correlation measurements were locking onto the same
    spurious signal (CR spikes, sky-background mismatch, and/or the sharp
    edge of the reprojected footprint dominate a sparse, CR-heavy field's
    FFT correlation far more than the ~30-60 real matchable sources do --
    exactly the per-exposure whole-frame matching failure mode the audit
    warned about). Matched-source centroids don't have that failure mode:
    with a tight (0.2") initial search radius, ~55 of ~250 F140W-detected
    sources already matched F200LP sources on this target, which is only
    possible if the true offset is already well under an arcsec -- the
    search radius grows only if that tight radius doesn't find enough
    matches to be statistically robust.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u
    from astropy.stats import sigma_clip

    sky_ref, _ = _detect_sources_sky(reference_sci_path, nsigma=detect_nsigma)
    sky_mov, _ = _detect_sources_sky(moving_sci_path, nsigma=detect_nsigma)
    if len(sky_ref) == 0 or len(sky_mov) == 0:
        raise RuntimeError(
            f"No sources detected in {reference_sci_path} and/or "
            f"{moving_sci_path} -- cannot measure an astrometric offset.")

    idx, sep2d, _ = sky_mov.match_to_catalog_sky(sky_ref)
    for radius in search_radii_arcsec:
        matched = sep2d < (radius * u.arcsec)
        if matched.sum() >= min_matches:
            break
    else:
        raise RuntimeError(
            f"Fewer than {min_matches} matched sources between "
            f"{reference_sci_path} and {moving_sci_path} even at "
            f"{search_radii_arcsec[-1]}\" search radius "
            f"(best: {matched.sum()}) -- cannot measure a robust offset.")

    m_mov, m_ref = sky_mov[matched], sky_ref[idx[matched]]
    dra = (m_ref.ra - m_mov.ra).to(u.arcsec) * np.cos(m_mov.dec.to(u.rad).value)
    ddec = (m_ref.dec - m_mov.dec).to(u.arcsec)
    dra_clipped = sigma_clip(dra.value, sigma=3, maxiters=5)
    ddec_clipped = sigma_clip(ddec.value, sigma=3, maxiters=5)

    n_used = int(dra_clipped.count())
    if n_used < min_matches:
        raise RuntimeError(
            f"Only {n_used} matches survived sigma-clipping between "
            f"{reference_sci_path} and {moving_sci_path} -- too few to "
            "trust a robust median offset.")

    dra_arcsec = float(np.ma.median(dra_clipped))
    ddec_arcsec = float(np.ma.median(ddec_clipped))
    print(f"    measure_filter_offset: {n_used} matched sources "
          f"(radius={radius}\"), scatter=({np.ma.std(dra_clipped)*1000:.1f}, "
          f"{np.ma.std(ddec_clipped)*1000:.1f}) mas")
    return dra_arcsec, ddec_arcsec


def _detect_sources_sky(sci_path, fwhm=2.0, nsigma=8.0):
    """DAOStarFinder source list for one mosaic, in sky coordinates, with
    basic shape filtering to reject cosmic-ray-like detections (the same
    concern as Defect 3's use_sharp_round, applied here to source
    detection rather than TweakReg)."""
    from astropy.stats import sigma_clipped_stats
    from photutils.detection import DAOStarFinder
    from astropy.coordinates import SkyCoord

    with fits.open(sci_path) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        wcs = WCS(hdu.header)
        data = np.asarray(hdu.data, dtype=float)

    _mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    finder = DAOStarFinder(fwhm=fwhm, threshold=nsigma * std,
                           sharplo=0.3, sharphi=0.9, roundlo=-0.5, roundhi=0.5)
    sources = finder(data - median)
    if sources is None or len(sources) == 0:
        return SkyCoord([], [], unit='deg'), np.empty(0)
    xy = np.column_stack([sources['xcentroid'], sources['ycentroid']])
    sky = wcs.pixel_to_world(xy[:, 0], xy[:, 1])
    return sky, np.asarray(sources['flux'])


def measure_filter_offset_xcorr(reference_sci_path, moving_sci_path):
    """
    Independent cross-check for measure_filter_offset: whole-frame phase
    cross-correlation, but with background subtraction and outlier
    clipping before correlating (the earlier, unguarded version of this
    measurement produced a wildly wrong ~19 arcsec answer on this dataset
    -- see measure_filter_offset's docstring). This is NOT used to drive
    any header correction; it exists only so the acceptance test can
    report agreement between two independent methods, per the audit's
    explicit request. Do not trust it in isolation on a sparse/CR-heavy
    field -- treat disagreement with measure_filter_offset as a real
    warning, not noise to average away.

    Sign convention, verified by simulation (see git history / audit
    notes): skimage.registration.phase_cross_correlation(ref, moving)
    returns the pixel shift that, applied to `moving`, registers it onto
    `ref`. That shift, converted to a sky vector via `reference`'s own WCS
    Jacobian at the field center, is the correction to add to `moving`'s
    CRVAL -- *not* its negative.
    """
    from reproject import reproject_interp
    from skimage.registration import phase_cross_correlation
    from astropy.stats import sigma_clipped_stats

    with fits.open(reference_sci_path) as href, fits.open(moving_sci_path) as hmov:
        ref_hdu = next(h for h in href if h.data is not None and h.data.ndim == 2)
        ref_wcs = WCS(ref_hdu.header)
        ref_data = np.asarray(ref_hdu.data, dtype=float)
        mov_hdu = next(h for h in hmov if h.data is not None and h.data.ndim == 2)
        mov_reproj, footprint = reproject_interp(mov_hdu, ref_hdu.header)

    good = np.isfinite(ref_data) & np.isfinite(mov_reproj) & (footprint > 0)
    if good.sum() < 1000:
        raise RuntimeError(
            f"Not enough overlap between {reference_sci_path} and "
            f"{moving_sci_path} to measure an astrometric offset "
            f"({good.sum()} finite overlapping pixels).")

    _mean_r, med_r, std_r = sigma_clipped_stats(ref_data[good], sigma=3.0, maxiters=5)
    _mean_m, med_m, std_m = sigma_clipped_stats(mov_reproj[good], sigma=3.0, maxiters=5)
    a = np.where(good, np.clip(ref_data - med_r, -8 * std_r, 8 * std_r), 0.0)
    b = np.where(good, np.clip(mov_reproj - med_m, -8 * std_m, 8 * std_m), 0.0)

    shift_px, _error, _diffphase = phase_cross_correlation(a, b, upsample_factor=20)
    dy_px, dx_px = shift_px

    ny, nx = ref_data.shape
    x0, y0 = nx / 2.0, ny / 2.0
    p0 = ref_wcs.pixel_to_world(x0, y0)
    p1 = ref_wcs.pixel_to_world(x0 + dx_px, y0 + dy_px)
    dra_arcsec = (p1.ra - p0.ra).wrap_at('180d').arcsec * np.cos(np.radians(p0.dec.deg))
    ddec_arcsec = (p1.dec - p0.dec).arcsec
    return float(dra_arcsec), float(ddec_arcsec)


def _apply_crval_shift(sci_path, dra_arcsec, ddec_arcsec):
    with fits.open(sci_path, mode='update') as hdul:
        for hdu in hdul:
            if hdu.data is None or 'CRVAL1' not in hdu.header:
                continue
            dec0 = hdu.header['CRVAL2']
            hdu.header['CRVAL1'] += dra_arcsec / 3600.0 / np.cos(np.radians(dec0))
            hdu.header['CRVAL2'] += ddec_arcsec / 3600.0
            hdu.header['HISTORY'] = (
                f'pipeline_qc: astrometry corrected to reference filter, '
                f'dRA={dra_arcsec * 1000:.1f}mas dDec={ddec_arcsec * 1000:.1f}mas')
        hdul.flush()


def reconcile_filter_astrometry(reference_sci_path, moving_sci_path):
    """
    Correct `moving`'s absolute WCS in place to agree with `reference`.

    This is the actual Defect-1 fix, and it has to touch the *native*
    drizzled science file -- not a display-only reprojected copy -- because
    hst_products.py's cutout builder reads each filter's own _sci.fits
    header directly (build_cutout_h5 -> WCS(header)) and never looks at the
    run_step3 *_scaled_L3 outputs. Correcting the header here means the
    correction is inherited automatically by (a) the lenstronomy cutouts,
    (b) run_step3's reprojection (which uses this same file as hdu_ref),
    and (c) everything downstream of those.

    Self-verifying: re-measures the residual offset after writing the
    correction and raises (reverting from a backup) if it did not shrink.
    This exists specifically to catch a sign error in measure_filter_offset
    -- shipping a WCS header that *claims* to be corrected but actually
    doubled the error would be worse than doing nothing, so this refuses
    to trust an unverified numerical correction.
    """
    reference_sci_path = Path(reference_sci_path)
    moving_sci_path = Path(moving_sci_path)

    before = measure_filter_offset(reference_sci_path, moving_sci_path)
    before_mag = float(np.hypot(*before))

    backup = moving_sci_path.with_suffix(moving_sci_path.suffix + '.pre_astrometry_fix.bak')
    if not backup.exists():
        shutil.copy2(moving_sci_path, backup)

    _apply_crval_shift(moving_sci_path, *before)

    after = measure_filter_offset(reference_sci_path, moving_sci_path)
    after_mag = float(np.hypot(*after))

    if after_mag >= before_mag:
        shutil.copy2(backup, moving_sci_path)
        raise RuntimeError(
            f"Astrometric correction for {moving_sci_path} did not improve "
            f"agreement with {reference_sci_path} "
            f"(before={before_mag * 1000:.1f} mas, after={after_mag * 1000:.1f} mas) "
            "-- reverted from backup. The measured offset or its sign is "
            "probably wrong for this pair; do not trust either header.")

    print(f"  [Defect 1 fix] {moving_sci_path.name}: corrected by "
          f"dRA={before[0] * 1000:+.1f} mas, dDec={before[1] * 1000:+.1f} mas "
          f"(residual vs {reference_sci_path.name}: {after_mag * 1000:.2f} mas)")
    return before, after


# ── Defect 2: DQ CR-flag survival assertion ─────────────────────────────────
#
# resetbits in this codebase is already explicit everywhere (never left at
# drizzlepac's "4096" default in a spot where that would erase real CR
# flags) -- see the audit. This assertion is the safety net the audit asked
# for regardless: a cheap, loud check that bit-4096 flags actually survived
# a drizzle call, so a future edit that reintroduces the erasure bug fails
# a run instead of silently producing an unmasked mosaic.

def assert_dq_flag_fraction(fits_paths, bit=4096, min_frac=0.0, label=""):
    """Raise if none of `fits_paths` still carries `bit` in its DQ arrays."""
    total = 0
    flagged = 0
    for p in fits_paths:
        with fits.open(p) as hdul:
            for hdu in hdul:
                if hdu.name == 'DQ' and hdu.data is not None:
                    dq = hdu.data
                    total += dq.size
                    flagged += int(np.sum((dq & bit) != 0))
    frac = (flagged / total) if total else 0.0
    if frac <= min_frac:
        raise AssertionError(
            f"[Defect 2 check{' — ' + label if label else ''}] "
            f"DQ bit {bit} fraction is {frac * 100:.4f}% (<= {min_frac * 100:.4f}%) "
            f"across {list(fits_paths)} -- CR flags appear to have been erased "
            "by this drizzle call (resetbits cleared the wrong bit, or a prior "
            "step never wrote them). Refusing to continue with an unmasked stack."
        )
    return frac


# ── Defect 4: tweakreg shiftfile NaN check ──────────────────────────────────
#
# Confirmed against the installed drizzlepac (imgclasses.py:performFit): when
# an exposure has too few matches (self.goodmatch=False, i.e. below minobj),
# drizzlepac sets fit['offset']=[nan,nan], fit['rot']=nan, fit['scale']=[nan]
# but still writes that image's row to the shiftfile, and with
# updatehdr=True it *still* stamps the new WCSNAME into that exposure's
# header regardless. So the shiftfile is the only place the failure is
# visible after the fact -- parse it and refuse to continue rather than
# ship a header that lies about the data.

def check_shiftfile_for_nan(shiftfile_path):
    path = Path(shiftfile_path)
    if not path.exists():
        raise FileNotFoundError(f"Expected TweakReg shiftfile not found: {path}")

    bad_rows = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            # columns per drizzlepac.imgclasses.Image.get_shiftfile_row:
            # name, xshift, yshift, rot, scale, xrms, yrms
            if any(p.lower() == 'nan' for p in parts[1:]):
                bad_rows.append(line)

    if bad_rows:
        raise RuntimeError(
            "[Defect 4 check] TweakReg produced NaN shift/rotation/scale for "
            f"the following exposure(s) in {path} (fit failed, most likely "
            "below minobj):\n  " + "\n  ".join(bad_rows) +
            "\nWith updatehdr=True, drizzlepac stamps the 'aligned' WCSNAME "
            "onto these files regardless of the failed fit -- refusing to "
            "continue rather than ship headers that misrepresent the data. "
            "Re-run with a lower threshold/minobj or inspect the exposure."
        )


# ── Defect 4b: tweakreg shiftfile fit-quality check ─────────────────────────
#
# check_shiftfile_for_nan above only catches the case where drizzlepac gives
# up entirely (NaN row). It does NOT catch a fit that converged onto junk,
# which is the more dangerous outcome because nothing downstream complains.
#
# Measured case (2026-10-01, AGEL103255+751854A 15867/F140W, from the
# 2026-09-17 run's surviving temp/shiftir_flt.txt):
#
#   ie5060bcq_crclean.fits  -0.664 1.384  rot=359.960595  scale=1.001509
#                                         xrms=4.365  yrms=4.724
#
# Those are three exposures of a single visit, dithered by POSTARG alone, so
# the true relative solution is a sub-pixel shift with rot=0 and scale=1 to
# ~1e-5. A 4.4 px residual with a 0.04 deg rotation and a 1.5e-3 scale change
# means xyxymatch paired noise peaks, not sources -- a direct consequence of
# IR_DRIZZLE_DEFAULTS['tolerance']=10.0 in a sparse field.
#
# Deliberately a warning, not a raise: this check was added after 193 targets
# were already reduced, and promoting it to an exception would retroactively
# fail re-runs of targets whose shipped products are fine. Promote it (set
# raise_on_fail=True at the call site) for the next full re-reduction.

# Two exposures of the same visit should align to well under a pixel; 1.0 px
# is a deliberately generous line between "real fit" and "matched noise".
SHIFTFILE_MAX_RMS_PX = 1.0
SHIFTFILE_MAX_ROT_DEG = 0.1
SHIFTFILE_MAX_SCALE_DEV = 1e-3


def check_shiftfile_fit_quality(shiftfile_path, raise_on_fail=False):
    """Flag TweakReg rows that converged but onto an implausible solution.

    Returns the list of suspect rows (empty if all rows look sane).
    """
    path = Path(shiftfile_path)
    if not path.exists():
        raise FileNotFoundError(f"Expected TweakReg shiftfile not found: {path}")

    suspect = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            try:
                xsh, ysh, rot, scale, xrms, yrms = (float(v) for v in parts[1:7])
            except (ValueError, IndexError):
                continue
            # The reference image's own row is an exact identity by construction.
            if (xsh, ysh, rot, scale) == (0.0, 0.0, 0.0, 1.0):
                continue
            # rot is reported in [0, 360); fold it to a signed deviation from 0.
            rot_dev = min(abs(rot), abs(rot - 360.0))
            reasons = []
            if max(xrms, yrms) > SHIFTFILE_MAX_RMS_PX:
                reasons.append(f'rms={max(xrms, yrms):.3f}px > {SHIFTFILE_MAX_RMS_PX}')
            if rot_dev > SHIFTFILE_MAX_ROT_DEG:
                reasons.append(f'rot={rot_dev:.4f}deg > {SHIFTFILE_MAX_ROT_DEG}')
            if abs(scale - 1.0) > SHIFTFILE_MAX_SCALE_DEV:
                reasons.append(f'|scale-1|={abs(scale - 1.0):.2e} > {SHIFTFILE_MAX_SCALE_DEV:.0e}')
            if reasons:
                suspect.append((Path(parts[0]).name, '; '.join(reasons), line))

    if suspect:
        msg = (
            "[Defect 4b check] TweakReg converged onto an implausible fit for "
            f"the following exposure(s) in {path}:\n  "
            + "\n  ".join(f'{name}: {why}' for name, why, _ in suspect)
            + "\nFor same-visit exposures expect sub-pixel rms, rot~0 and "
              "scale~1. A large residual usually means 'tolerance' is too "
              "loose for how sparse this field is -- lower 'tolerance' (and "
              "'tweak_threshold' to find more real sources) for this target "
              "via TWEAK_PARAM_OVERRIDES in config.py."
        )
        if raise_on_fail:
            raise RuntimeError(msg)
        print(f'[WARN] {msg}')
    return suspect


# ── Defect 5: N=2 cosmic-ray rejection ──────────────────────────────────────
#
# driz_cr's own median-based rejection degenerates to an average at N=2 and
# cannot reliably separate a CR from an undersampled source core. This
# section implements the four-part fix from the audit:
#   (a) L.A.Cosmic in PSF-convolution mode (astroscrappy)
#   (b) explicit source protection via repeat detection between the two
#       exposures
#   (c) combine_type='minmed' (already drizzlepac's default here, pinned
#       explicitly at the call site in hst_reduction.py)
#   (d) pairwise rejection in the output frame, mapped back to detector
#       coordinates through the full distortion model

def lacosmic_flag(data, gain=1.0, readnoise=3.1, sigclip=4.0, objlim=3.0,
                   psffwhm=1.9, inmask=None):
    """
    L.A.Cosmic CR mask for one exposure, in PSF-convolution mode.

    fsmode='median' (astroscrappy's default) compares each pixel to a
    median-filtered version of itself and cannot tell a sharp CR from an
    undersampled stellar core at F200LP's ~1.9 px FWHM -- fsmode='convolve'
    with an explicit Gaussian PSF model avoids that failure mode (see the
    audit: median mode destroyed 14-29% of real source cores in testing on
    this dataset).
    """
    import astroscrappy
    mask, _clean = astroscrappy.detect_cosmics(
        np.asarray(data, dtype=float), inmask=inmask,
        gain=gain, readnoise=readnoise,
        sigclip=sigclip, objlim=objlim,
        fsmode='convolve', psfmodel='gauss', psffwhm=psffwhm,
        cleantype='meanmask', verbose=False,
    )
    return mask


def detect_source_positions(data, fwhm=1.9, nsigma=5.0):
    """
    DAOStarFinder source list (x, y, flux) for repeat detection.

    sharplo/sharphi/roundlo/roundhi added after a real run on this target
    showed why they matter: without shape filtering, this found 636
    "confirmed" repeat-detected sources between the two F200LP exposures --
    far more than the ~47-54 real sources expected for this field. The
    excess is chance coincidences between cosmic-ray-driven spurious
    DAOStarFinder detections in each exposure (there are thousands of CR
    pixels pre-masking; even a small per-pixel coincidence probability
    across a multi-megapixel frame adds up), i.e. the same "DQ/shape-blind
    detection picks up CRs as sources" failure mode as Defect 3, just in
    this module's own source finder rather than TweakReg's. These cuts
    match what worked for the Defect 1 offset measurement's source finder.
    """
    from astropy.stats import sigma_clipped_stats
    from photutils.detection import DAOStarFinder

    mean, median, std = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
    finder = DAOStarFinder(fwhm=fwhm, threshold=nsigma * std,
                           sharplo=0.3, sharphi=0.9, roundlo=-0.5, roundhi=0.5)
    sources = finder(data - median)
    if sources is None or len(sources) == 0:
        return np.empty((0, 2)), np.empty((0,))
    xy = np.column_stack([sources['xcentroid'], sources['ycentroid']])
    return xy, np.asarray(sources['flux'])


def build_source_protection_mask(shape, xy_pixels, radius_px=4):
    """Boolean mask, True within `radius_px` of any (x, y) in `xy_pixels`."""
    mask = np.zeros(shape, dtype=bool)
    ny, nx = shape
    r2 = radius_px ** 2
    for x, y in xy_pixels:
        x0, x1 = max(0, int(x - radius_px)), min(nx, int(x + radius_px) + 1)
        y0, y1 = max(0, int(y - radius_px)), min(ny, int(y + radius_px) + 1)
        if x1 <= x0 or y1 <= y0:
            continue
        yy, xx = np.mgrid[y0:y1, x0:x1]
        mask[y0:y1, x0:x1] |= ((xx - x) ** 2 + (yy - y) ** 2) <= r2
    return mask


def repeat_detected_protection_mask(exp1_drz_path, exp2_drz_path, match_radius_arcsec=0.15,
                                     fwhm=3.0, nsigma=5.0):
    """
    Cross-match sources independently detected in two single-exposure
    drizzles (same output grid) and return their matched sky positions.

    Anything detected at the same sky position in both single-exposure
    drizzles cannot be a cosmic ray (a CR in one exposure has ~zero chance
    of landing on the same sky pixel in the other, independently-dithered
    exposure) -- this is a hard guarantee independent of L.A.Cosmic tuning,
    which is what lets step (a) run aggressively without risking real flux.
    """
    from astropy.coordinates import SkyCoord
    from astropy import units as u

    with fits.open(exp1_drz_path) as h1, fits.open(exp2_drz_path) as h2:
        d1 = next(h for h in h1 if h.data is not None and h.data.ndim == 2)
        d2 = next(h for h in h2 if h.data is not None and h.data.ndim == 2)
        wcs1, wcs2 = WCS(d1.header), WCS(d2.header)
        xy1, _ = detect_source_positions(np.asarray(d1.data, dtype=float), fwhm=fwhm, nsigma=nsigma)
        xy2, _ = detect_source_positions(np.asarray(d2.data, dtype=float), fwhm=fwhm, nsigma=nsigma)

    if len(xy1) == 0 or len(xy2) == 0:
        return SkyCoord([], [], unit='deg')

    sky1 = wcs1.pixel_to_world(xy1[:, 0], xy1[:, 1])
    sky2 = wcs2.pixel_to_world(xy2[:, 0], xy2[:, 1])
    idx, sep2d, _ = sky1.match_to_catalog_sky(sky2)
    matched = sep2d < (match_radius_arcsec * u.arcsec)
    return sky1[matched]


def sky_to_detector_pixels(sky_coords, exposure_sci_path):
    """
    Map sky positions to detector (x, y) pixel coordinates of one raw
    exposure's SCI extension, through its full distortion model (SIP +
    D2IM), via astropy WCS all_world2pix on the exposure's own header WCS
    (which for a calibrated flt/flc already includes the SIP + D2IM/NPOL
    corrections drizzlepac wrote into it).
    """
    with fits.open(exposure_sci_path) as hdul:
        wcs = WCS(hdul['SCI', 1].header, hdul)
    x, y = wcs.world_to_pixel(sky_coords)
    return np.column_stack([x, y])


def apply_custom_cr_dq(fits_path, ext_ver, cr_mask, protect_mask=None, bit=8192):
    """
    Write `cr_mask` (True = cosmic ray) into DQ bit `bit` of one SCI/DQ
    extension pair, clearing any protected pixels first. Uses bit 8192
    (not 4096) deliberately -- see Defect 2 -- so this never collides with
    archive or driz_cr flags, and downstream AstroDrizzle calls only need
    resetbits=0 to keep these flags alive to the final drizzle.
    """
    mask = cr_mask.copy()
    if protect_mask is not None:
        mask &= ~protect_mask
    with fits.open(fits_path, mode='update') as hdul:
        dq = hdul['DQ', ext_ver].data
        dq[mask] = dq[mask] | bit
        hdul.flush()
    return int(mask.sum())


def pairwise_reject_output_frame(single_drz_paths, exposure_sci_paths, threshold_k=4.0,
                                 dilation_px=1, protect_sky=None, protect_radius_px=4):
    """
    Compare the two single-exposure drizzles (same output grid) directly
    and flag pixels where one is a large positive excursion relative to the
    other -- cosmic rays are always positive excursions, so a simple
    two-sided z-score-style comparison catches residual CRs that survive
    (a)-(c). Threshold k=4.0 with single-pixel dilation is what the audit's
    reference implementation found worked (k=2.5 flagged ~25% of all pixels
    -- far too aggressive).

    `protect_sky` (SkyCoord, optional): confirmed real sources from the same
    repeat-detection catalog step (b) uses. A real run on DESJ0206-0114
    without this showed why it's needed: 15 of 20 confirmed sources lost
    flux below the CR-free floor even though steps (a)-(c) protected them,
    because this step -- comparing the two exposures directly in the output
    frame -- has no way to know which pixels are protected unless told.
    Bright, dithered sources naturally show larger per-pixel differences
    between two exposures near their PSF core/wings (sub-pixel registration
    residuals, not cosmic rays), and this step was flagging exactly that.

    Flagged output-frame pixels are mapped back to each exposure's detector
    coordinates via drizzlepac.pixtopix.tran (direction='backward'), which
    composes the full IDCTAB/SIP + D2IM distortion model rather than a
    linear WCS approximation, and the corresponding DQ bit 8192 pixels are
    set so a final re-drizzle picks them up.

    Returns the number of detector pixels flagged, per exposure.
    """
    from scipy.ndimage import binary_dilation
    from drizzlepac import pixtopix

    assert len(single_drz_paths) == 2 and len(exposure_sci_paths) == 2

    with fits.open(single_drz_paths[0]) as h0, fits.open(single_drz_paths[1]) as h1:
        d0 = next(h for h in h0 if h.data is not None and h.data.ndim == 2)
        d1 = next(h for h in h1 if h.data is not None and h.data.ndim == 2)
        a = np.asarray(d0.data, dtype=float)
        b = np.asarray(d1.data, dtype=float)
        out_wcs = WCS(d0.header)

    diff = a - b
    sigma = np.nanstd(diff[np.isfinite(diff) & (diff != 0)])
    if not np.isfinite(sigma) or sigma == 0:
        return [0, 0]

    flag_a = diff > threshold_k * sigma   # a is a positive excursion vs b -> CR in exposure a
    flag_b = diff < -threshold_k * sigma  # b is a positive excursion vs a -> CR in exposure b
    if dilation_px:
        flag_a = binary_dilation(flag_a, iterations=dilation_px)
        flag_b = binary_dilation(flag_b, iterations=dilation_px)

    if protect_sky is not None and len(protect_sky):
        out_xy = np.column_stack(out_wcs.world_to_pixel(protect_sky))
        protect_mask = build_source_protection_mask(a.shape, out_xy, radius_px=protect_radius_px)
        flag_a &= ~protect_mask
        flag_b &= ~protect_mask

    n_flagged = []
    for flag_mask, single_path, exp_path in (
        (flag_a, single_drz_paths[0], exposure_sci_paths[0]),
        (flag_b, single_drz_paths[1], exposure_sci_paths[1]),
    ):
        ys, xs = np.nonzero(flag_mask)
        if len(xs) == 0:
            n_flagged.append(0)
            continue
        det_x, det_y = pixtopix.tran(
            str(single_path), str(exp_path), direction='backward',
            x=xs.astype(float), y=ys.astype(float), verbose=False)
        det_x = np.atleast_1d(det_x)
        det_y = np.atleast_1d(det_y)

        with fits.open(exp_path, mode='update') as hdul:
            n_sci = sum(1 for h in hdul if h.name == 'SCI')
            for ext_ver in range(1, n_sci + 1):
                sci_hdr = hdul['SCI', ext_ver].header
                naxis1, naxis2 = sci_hdr['NAXIS1'], sci_hdr['NAXIS2']
                in_chip = (det_x >= 0) & (det_x < naxis1) & (det_y >= 0) & (det_y < naxis2)
                if not np.any(in_chip):
                    continue
                dq = hdul['DQ', ext_ver].data
                ix = np.clip(det_x[in_chip].round().astype(int), 0, naxis1 - 1)
                iy = np.clip(det_y[in_chip].round().astype(int), 0, naxis2 - 1)
                dq[iy, ix] = dq[iy, ix] | 8192
            hdul.flush()
        n_flagged.append(len(xs))

    return n_flagged


# ── Defect 6: tweakreg skysigma/count-rate mismatch ─────────────────────────
#
# drizzlepac's automatic source-detection sigma (computesig=True, the
# default whenever imagefindcfg doesn't set skysigma) estimates background
# noise as sigma = sqrt(2 * mode) -- a Poisson-noise formula that assumes
# the image is in raw COUNTS. Every image this pipeline hands to TweakReg
# is BUNIT='ELECTRONS/S' (count RATE, post-calibration), so that estimate
# is inflated by roughly sqrt(exposure time) -- e.g. for a 399s exposure,
# sqrt(399) ~= 20x. Confirmed directly against a real batch failure
# (AGEL103255+751854A, 15867/F140W): drizzlepac's auto-sigma was 1.265
# against an empirical (sigma-clipped, DQ-masked) pixel sigma of 0.073 --
# a 17x inflation, matching the sqrt(exptime) prediction. The effect is a
# real detection threshold roughly (inflation)x stricter than the
# 'threshold' parameter asks for, silently discarding genuine sources on
# every target -- it only shows up as an outright failure (zero sources
# found, "Fewer than two images available for alignment") in fields with no
# companion bright enough to survive the inflated cut anyway, but degrades
# alignment quality (fewer real matches to fit against) everywhere else.
#
# Fix: measure sigma ourselves from the actual pixel data (sigma-clipped,
# DQ-masked, same convention TweakReg itself would want) and pass it via
# skysigma with computesig=False, bypassing the flawed auto-estimate.

def compute_skysigma(fits_paths, dq_mask_bits=4096, sigma=3.0):
    """Empirical, DQ-masked, sigma-clipped background sigma for a set of
    crclean.fits files, to pass as TweakReg's imagefindcfg['skysigma'] with
    computesig=False. One representative value (median across every SCI
    extension of every input file) is returned since TweakReg's skysigma
    is a single shared value across all images in one call; the per-image
    values are typically within ~15% of each other for exposures of the
    same target/filter, so this is a reasonable shared estimate.
    """
    sigmas = []
    for path in fits_paths:
        with fits.open(path) as hdul:
            n_sci = sum(1 for h in hdul if h.name == 'SCI')
            for ext_ver in range(1, n_sci + 1):
                data = hdul['SCI', ext_ver].data
                dq = hdul['DQ', ext_ver].data if ('DQ', ext_ver) in hdul else None
                mask = (dq & dq_mask_bits) != 0 if dq is not None else None
                _, _, std = sigma_clipped_stats(data, mask=mask, sigma=sigma)
                sigmas.append(std)
    return float(np.median(sigmas))
