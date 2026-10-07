"""
hst_reduction.py — HST data reduction pipeline

Runs four steps in sequence for the proposal defined in config.py:
  1. Download  : fetch FLC/FLT files from MAST
  2. Drizzle   : align and drizzle IR (F140W) and UV (all other) exposures
  3. Reproject : resample UV frames onto the IR pixel grid, rename to *_L1.fits
  4. Green     : build synthetic green image from IR + UV for RGB composites

Usage:
    python hst_reduction.py [--steps 1234] [--target AGEL...] [--interactive]

    --steps        string of step numbers to run (default: all, e.g. '234' skips download)
    --target       run only one target (useful for testing / re-running a single object)
    --interactive, -i
                   review/override AstroDrizzle+TweakReg parameters before step 2 runs,
                   and confirm (or skip) each target's drizzle interactively
"""

import argparse
import gc
import glob
import os
import shutil
import stat
import time
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp

from config import (
    ACTIVE_CAMERA,
    ACTIVE_FILTERS,
    ACTIVE_PROPOSAL_ID,
    ACTIVE_TARGETS_CSV,
    IR_DRIZZLE_DEFAULTS,
    LACOSMIC_DEFAULTS,
    MAIN_DIR,
    PAIRWISE_CR_DEFAULTS,
    SOURCE_PROTECTION_DEFAULTS,
    TWEAK_PARAM_OVERRIDES,
    UV_DRIZZLE_DEFAULTS,
)
import pipeline_qc

# ── Interactive helpers ─────────────────────────────────────────────────────────

def _prompt_yes_no(prompt: str, default: bool = True) -> bool:
    suffix = '[Y/n]' if default else '[y/N]'
    ans = input(f"{prompt} {suffix}: ").strip().lower()
    if ans == '':
        return default
    return ans in ('y', 'yes')


def _prompt_download_dir(default_dir: Path) -> Path:
    """Confirm (or override) the base directory downloaded data will be stored under."""
    print(f"\nData will be downloaded to: {default_dir}")
    raw = input("Press Enter to accept, or type a different path: ").strip()
    if raw == '':
        return default_dir
    return Path(raw).expanduser()


def _prompt_edit_params(label: str, defaults: dict) -> dict:
    """Show each default value and let the user override it (Enter keeps default)."""
    print(f"\n--- {label} drizzle/tweakreg parameters (Enter to keep default) ---")
    params = dict(defaults)
    for key, val in defaults.items():
        raw = input(f"  {key} [{val}]: ").strip()
        if raw == '':
            continue
        if isinstance(val, bool):
            params[key] = raw.lower() in ('1', 'true', 'y', 'yes')
        elif isinstance(val, int) and not isinstance(val, bool):
            try:
                params[key] = int(raw)
            except ValueError:
                params[key] = raw
        elif isinstance(val, float):
            try:
                params[key] = float(raw)
            except ValueError:
                params[key] = raw
        else:
            params[key] = raw
    return params

# ── Utility helpers ────────────────────────────────────────────────────────────

def _fix_perms_and_retry(func, path, exc_info):
    try:
        os.chmod(path, stat.S_IWRITE | stat.S_IREAD | stat.S_IEXEC)
    except Exception:
        pass
    try:
        func(path)
    except Exception:
        pass


def ensure_cwd_not_inside(target: Path):
    try:
        if Path.cwd().resolve().is_relative_to(target.resolve()):
            os.chdir(target.parent.as_posix())
    except AttributeError:
        cwd = Path.cwd().resolve()
        tgt = target.resolve()
        try:
            cwd.relative_to(tgt)
            os.chdir(tgt.parent.as_posix())
        except ValueError:
            pass


def clean_mac_junk_recursive(root: Path, delete_cwd_ancillary: bool = True,
                              ancillary_exts=None) -> int:
    removed = 0
    for patt in ("**/._*", "**/.DS_Store"):
        for p in root.rglob(patt):
            try:
                p.unlink()
                removed += 1
            except FileNotFoundError:
                pass
    if ancillary_exts is None:
        ancillary_exts = ['coo', 'png', 'fits', 'match', 'log', 'list']
    if delete_cwd_ancillary:
        cwd = Path(os.getcwd()).resolve()
        try:
            same_dir = (cwd == root.resolve())
        except Exception:
            same_dir = False
        if not same_dir:
            for ext in ancillary_exts:
                for f in glob.glob(f"./*.{ext}"):
                    try:
                        os.remove(f)
                        removed += 1
                    except FileNotFoundError:
                        pass
    return removed


def safe_rmtree(p: Path, retries: int = 3, delay: float = 0.25):
    p = Path(p)
    try:
        if Path.cwd().resolve().is_relative_to(p.resolve()):
            raise RuntimeError(f"CWD is inside {p}; refusing to remove")
    except AttributeError:
        cwd = Path.cwd().resolve()
        try:
            cwd.relative_to(p.resolve())
            raise RuntimeError(f"CWD is inside {p}; refusing to remove")
        except ValueError:
            pass
    for i in range(retries + 1):
        try:
            shutil.rmtree(p, onerror=_fix_perms_and_retry)
            return
        except FileNotFoundError:
            return
        except OSError:
            if i == retries:
                raise
            time.sleep(delay)


def list_fits(dirpath: Path, pattern: str):
    return [p.as_posix() for p in dirpath.glob(pattern) if not p.name.startswith('._')]


def derive_rootname(fpath: str):
    try:
        return fits.getval(fpath, 'ROOTNAME')
    except Exception:
        return Path(fpath).name.split('_')[0]


def _float_or_none(v):
    try:
        if v is None:
            return None
        s = str(v).strip()
        if s == "" or s.lower() in {"none", "nan", "null"}:
            return None
        return float(s)
    except Exception:
        return None


def _resolve_z(row, src_cols, def_cols):
    """Return (z_deflector, z_source) from a CSV row, trying columns in order."""
    z_src = next((_float_or_none(row[c]) for c in src_cols if c in row and _float_or_none(row[c]) is not None), None)
    z_def = next((_float_or_none(row[c]) for c in def_cols if c in row and _float_or_none(row[c]) is not None), None)
    return z_def, z_src


SRC_COLS = ['z_spec_SRC_Spectral_Observations_Tally', 'z_source (from DR2 Redshifts)', 'z_spec_SRC']
DEF_COLS = ['z_spec_DE_Spectral_Observations_Tally',  'z_deflector (from DR2 Redshifts)',  'z_spec_DE']

# ── Step 1: Download ───────────────────────────────────────────────────────────

def download(lens_name: str, target_name: str, band: str, proposal_id: str, base_dir: Path = MAIN_DIR):
    """Download FLC/FLT files for one target from MAST."""
    from astroquery.mast import Observations  # import here so the rest of the script works without astroquery

    download_dir = base_dir / f'{lens_name}/HST/{proposal_id}_{band}/raw_data'
    download_dir.mkdir(parents=True, exist_ok=True)

    obs_table = Observations.query_criteria(
        proposal_id=proposal_id, target_name=target_name, filters=band)
    print(obs_table)

    file_format = ['FLT'] if band == 'F140W' else ['FLC']
    download_tab = Observations.download_products(
        obs_table['obsid'], mrp_only=False,
        download_dir=str(download_dir),
        productSubGroupDescription=file_format)

    science_files = glob.glob(
        str(download_dir / 'mastDownload' / 'HST' / '*' / '*fits'))
    for im in science_files:
        shutil.copy(im, str(download_dir))

    shutil.rmtree(str(download_dir / 'mastDownload'))
    for fp in glob.glob(str(download_dir / 'hst*')):
        try:
            os.remove(fp)
        except Exception as e:
            print(f"Could not delete {fp}: {e}")

    print(f"  Saved to: {download_dir.resolve()}")


def run_step1(targets, filters, proposal_id):
    print("\n=== STEP 1: Download ===")
    base_dir = _prompt_download_dir(MAIN_DIR)
    for row in targets:
        agel_name, mast_name = row[0], row[1]
        if agel_name.endswith('B'):
            continue
        for band in filters:
            print(f"  Downloading {agel_name} / {band}")
            download(agel_name, mast_name, band, proposal_id, base_dir=base_dir)

# ── Step 2: Drizzle ────────────────────────────────────────────────────────────

def ir_drizzler(object_name: str, propno: str, band: str, camera: str, params: dict = None):
    """Drizzle IR (WFC3/IR FLT) exposures."""
    from drizzlepac import astrodrizzle, tweakreg, tweakback

    p = dict(IR_DRIZZLE_DEFAULTS)
    if params:
        p.update(params)
    # Per-target TweakReg overrides win over both the defaults and the run-wide
    # (interactive) params -- they exist because a specific target's field is
    # too sparse for the shared values. See TWEAK_PARAM_OVERRIDES in config.py.
    target_override = TWEAK_PARAM_OVERRIDES.get((str(propno), object_name))
    if target_override:
        print(f"[INFO] Applying per-target TweakReg override for {object_name} "
              f"({propno}): {target_override}")
        p.update(target_override)

    raw_data_dir    = MAIN_DIR / f'{object_name}/HST/{propno}_{band}/raw_data'
    output_data_dir = MAIN_DIR / f'{object_name}/HST/{propno}_{band}'
    output_data_dir.mkdir(exist_ok=True)

    temp_dir = output_data_dir / 'temp'
    safe_rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    for src in raw_data_dir.glob('*.fits'):
        shutil.copy(src.as_posix(), temp_dir.as_posix())

    flt_files = list_fits(temp_dir, '*flt.fits')
    if not flt_files:
        raise FileNotFoundError(f"No *flt.fits in {temp_dir}")

    output_prefix = (output_data_dir / f'{object_name}_{band}_{camera}').as_posix()

    if len(flt_files) < 2:
        # A lone exposure has nothing to relatively align to -- no second
        # frame, and no external reference catalog is configured (refcat=
        # None). Measured directly against this exact failure
        # (AGEL001030-431515A, 16773/F140W, 1 raw exposure): TweakReg does
        # not raise or write a NaN row in this situation, it prints
        # "Fewer than two images available for alignment. Quitting..." and
        # exits WITHOUT writing a shiftfile at all, which made
        # check_shiftfile_for_nan's "missing file" branch fire on every
        # single-exposure target -- a structural certainty, not a fit
        # failure worth flagging. Skipping the call entirely leaves the
        # exposure on its native/archival WCS, which is what
        # reconcile_target_astrometry (run later, in run_step2) already
        # expects and corrects against the other filter's solution anyway.
        print(f"[WARN] Only {len(flt_files)} FLT file(s); skipping initial CR rejection "
              "and TweakReg (nothing to align a lone exposure against).")
    else:
        astrodrizzle.AstroDrizzle(
            flt_files, output=output_prefix,
            resetbits=4096, driz_cr_corr=True, preserve=False,
            final_wht_type=p['final_wht_type'], final_pixfrac=p['final_pixfrac'])

        shiftfile = (temp_dir / 'shiftir_flt.txt').as_posix()
        # See pipeline_qc.compute_skysigma (Defect 6): drizzlepac's own
        # auto-sigma assumes raw counts, but these images are count-rate
        # (ELECTRONS/S) -- computesig=False + an empirically-measured
        # skysigma here bypasses that ~sqrt(exptime) threshold inflation.
        skysigma = pipeline_qc.compute_skysigma(
            list_fits(temp_dir, '*crclean.fits'), dq_mask_bits=4096)
        tweakreg.TweakReg(
            (temp_dir / '*crclean.fits').as_posix(),
            # use_sharp_round=True added for the same reason as
            # above (belt-and-suspenders with the crclean images, which are
            # CR-interpolated but not CR-masked from source detection).
            imagefindcfg={'threshold': p['tweak_threshold'], 'conv_width': p['tweak_conv_width'],
                          'dqbits': ~4096, 'use_sharp_round': True,
                          'computesig': False, 'skysigma': skysigma},
            # Previously omitted (relying on drizzlepac's own refimagefindpars
            # defaults) -- but that default also has computesig=True, so the
            # reference image was just as affected by Defect 6 as every other
            # input image (confirmed: the original failure showed 0 objects
            # found in EVERY frame, including whichever was chosen as
            # reference). conv_width/dqbits/use_sharp_round mirror the
            # imagefindcfg values, matching the established pattern in the UV
            # branches below; 'threshold' is deliberately left at its own
            # default there too, same as those branches.
            refimagefindcfg={'conv_width': p['tweak_conv_width'], 'dqbits': ~4096,
                             'use_sharp_round': True,
                             'computesig': False, 'skysigma': skysigma},
            expand_refcat=True, enforce_user_order=False, shiftfile=True,
            outshifts=shiftfile,
            searchrad=p['tweak_searchrad'], ylimit=p['tweak_ylimit'],
            use2dhist=p['tweak_use2dhist'], xoffset=p['tweak_xoffset'],
            yoffset=p['tweak_yoffset'], tolerance=p['tolerance'],
            separation=p['tweak_separation'], minobj=p['tweak_minobj'],
            updatehdr=True,
            reusename=True, wcsname='IR_FLT', interactive=False)
        pipeline_qc.check_shiftfile_for_nan(shiftfile)
        pipeline_qc.check_shiftfile_fit_quality(shiftfile)

        for flt in flt_files:
            flid = derive_rootname(flt)
            crclean = temp_dir / f'{flid}_crclean.fits'
            if not crclean.exists():
                print(f"[WARN] Missing {crclean.name}; skipping tweakback.")
                continue
            tweakback.tweakback(crclean.as_posix(), input=flt, wcsname='IR_FLT')

    astrodrizzle.AstroDrizzle(
        flt_files, output=output_prefix,
        resetbits=8192, driz_cr_corr=False, preserve=False,
        final_pixfrac=p['final_pixfrac'], final_wcs=True, final_scale=p['pixel_size'])
    if len(flt_files) >= 2:
        pipeline_qc.assert_dq_flag_fraction(
            flt_files, bit=4096, min_frac=0.0,
            label=f'{object_name} {band}: pre final IR drizzle')

    gc.collect()
    clean_mac_junk_recursive(temp_dir)
    ensure_cwd_not_inside(temp_dir)
    time.sleep(0.1)
    safe_rmtree(temp_dir)
    print(f"[CLEAN] Removed temp: {temp_dir}")


def uv_drizzler(object_name: str, propno: str, band: str, camera: str, params: dict = None):
    """Drizzle UV/UVIS (ACS/WFC3 FLC) exposures."""
    from drizzlepac import astrodrizzle, tweakreg, tweakback

    p = dict(UV_DRIZZLE_DEFAULTS)
    if params:
        p.update(params)

    raw_data_dir    = MAIN_DIR / f'{object_name}/HST/{propno}_{band}/raw_data'
    output_data_dir = MAIN_DIR / f'{object_name}/HST/{propno}_{band}'
    output_data_dir.mkdir(exist_ok=True)

    temp_dir = output_data_dir / 'temp'
    safe_rmtree(temp_dir)
    temp_dir.mkdir(parents=True, exist_ok=True)

    for src in raw_data_dir.glob('*.fits'):
        shutil.copy(src.as_posix(), temp_dir.as_posix())

    flc_files = list_fits(temp_dir, '*flc.fits')
    if not flc_files:
        raise FileNotFoundError(f"No *flc.fits in {temp_dir}")

    output_prefix = (output_data_dir / f'{object_name}_{band}_{camera}').as_posix()

    if len(flc_files) < 2:
        print(f"[INFO] Single exposure — running AstroDrizzle without tweakreg.")
        # resetbits=8192 leaves archive CR flags (bit 4096) alone;
        # there's no second exposure for driz_cr to reject against anyway.
        astrodrizzle.AstroDrizzle(
            flc_files, output=output_prefix, preserve=False,
            resetbits=8192, driz_cr_scale=p['driz_cr_scale'], driz_cr_grow=p['driz_cr_grow'],
            final_wht_type=p['final_wht_type'], final_pixfrac=p['final_pixfrac_final'],
            final_wcs=True, final_scale=p['pixel_size'])
    elif len(flc_files) == 2:
        # N=2 median-based driz_cr rejection degenerates to a
        # simple average and cannot reliably separate a cosmic ray from an
        # undersampled point source at this dither (~2.5 UVIS px). See
        # pipeline_qc.py and _uv_drizzle_n2 below for the full L.A.Cosmic +
        # source-protection + pairwise-rejection stack this case needs.
        _uv_drizzle_n2(temp_dir, flc_files, output_prefix, p)
    else:
        # resetbits=8192 (not the drizzlepac default of 4096)
        # deliberately leaves archive CR flags (bit 4096) untouched through
        # this call, so driz_cr's freshly-derived flags (also bit 4096,
        # crbit default) are added on top rather than replacing them.
        astrodrizzle.AstroDrizzle(
            flc_files, output=output_prefix,
            resetbits=8192, driz_cr_corr=True, preserve=False,
            driz_cr_scale=p['driz_cr_scale'], driz_cr_grow=p['driz_cr_grow'],
            final_wht_type=p['final_wht_type'], final_pixfrac=p['final_pixfrac_initial'])

        shiftfile = (temp_dir / 'shiftuv_flt.txt').as_posix()
        # See pipeline_qc.compute_skysigma (Defect 6): bypass drizzlepac's
        # count-rate-blind auto-sigma the same way as the IR branch above.
        skysigma = pipeline_qc.compute_skysigma(
            list_fits(temp_dir, '*crclean.fits'), dq_mask_bits=4096)
        tweakreg.TweakReg(
            (temp_dir / '*crclean.fits').as_posix(),
            enforce_user_order=False,
            # dqbits ~4096, the bit driz_cr and the archive actually use, and added
            # use_sharp_round=True. refimagefindcfg now also carries both,
            # so the reference-image source catalog gets the same DQ
            # filtering as every other exposure (it previously had none).
            imagefindcfg={'threshold': p['tweak_threshold'], 'conv_width': p['tweak_conv_width'],
                          'dqbits': ~4096, 'use_sharp_round': True,
                          'computesig': False, 'skysigma': skysigma},
            refimagefindcfg={'conv_width': p['tweak_refconv_width'],
                             'dqbits': ~4096, 'use_sharp_round': True,
                             'computesig': False, 'skysigma': skysigma},
            shiftfile=True,
            outshifts=shiftfile,
            searchrad=p['tweak_searchrad'], ylimit=p['tweak_ylimit'],
            use2dhist=p['tweak_use2dhist'], xoffset=p['tweak_xoffset'],
            yoffset=p['tweak_yoffset'], tolerance=p['tolerance'],
            separation=p['tweak_separation'], updatehdr=True,
            wcsname='UVIS_FLC', reusename=True, interactive=False)
        pipeline_qc.check_shiftfile_for_nan(shiftfile)

        for flc in flc_files:
            flid = derive_rootname(flc)
            crclean = temp_dir / f'{flid}_crclean.fits'
            tweakback.tweakback(crclean.as_posix(), input=flc, wcsname='UVIS_FLC')

        astrodrizzle.AstroDrizzle(
            flc_files, output=output_prefix,
            resetbits=8192, driz_cr_corr=False, preserve=False,
            driz_cr_scale=p['driz_cr_scale'], driz_cr_grow=p['driz_cr_grow'],
            final_pixfrac=p['final_pixfrac_final'], final_wcs=True, final_scale=p['pixel_size'])
        pipeline_qc.assert_dq_flag_fraction(
            flc_files, bit=4096, min_frac=0.0,
            label=f'{object_name} {band}: pre final UV drizzle')

    gc.collect()
    clean_mac_junk_recursive(temp_dir)
    ensure_cwd_not_inside(temp_dir)
    time.sleep(0.1)
    safe_rmtree(temp_dir)
    print(f"[CLEAN] Removed temp: {temp_dir}")


def _uv_drizzle_n2(temp_dir: Path, flc_files: list, output_prefix: str, p: dict):
    """
    Defect 5: F200LP-style N=2 drizzle with the full CR-rejection stack the
    audit called for, in order:
      (a) L.A.Cosmic in PSF-convolution mode, per exposure
      (b) source protection via repeat detection between the two exposures
          (un-flags anything seen in both -- can't be a cosmic ray)
      (c) AstroDrizzle driz_cr with combine_type='minmed' (pinned explicitly;
          it happens to already be drizzlepac's default, but pin it so a
          future drizzlepac upgrade can't silently change it)
      (d) pairwise rejection comparing the two single-exposure drizzles in
          the output frame, mapped back to detector coordinates and
          re-drizzled

    All custom CR flags go into DQ bit 8192, never 4096 (Defect 2) -- keeps
    them out of the way of archive/driz_cr flags and means every call below
    can safely use resetbits=0 without erasing anything.

    This is new code, not a mechanical fix to an existing call site like the
    rest of this audit -- it has not been validated end-to-end against a
    real drizzled mosaic. Treat the first run on DESJ0206-0114 as validating
    this function, not just the target.
    """
    from drizzlepac import astrodrizzle, tweakreg, tweakback

    la = dict(LACOSMIC_DEFAULTS)
    sp = dict(SOURCE_PROTECTION_DEFAULTS)
    pw = dict(PAIRWISE_CR_DEFAULTS)

    # --- (a)+(b): preliminary single-exposure drizzles for repeat detection
    # Built before any CR masking so the source catalog isn't biased by
    # asymmetric CR cleaning between the two exposures. Sky-based matching
    # (not pixel-grid matching) is used, so these don't need a shared grid.
    prelim_singles = []
    for flc in flc_files:
        prelim_prefix = str(Path(flc).with_name(f'{Path(flc).stem}_prelim'))
        astrodrizzle.AstroDrizzle(
            [flc], output=prelim_prefix, resetbits=0, driz_cr_corr=False,
            preserve=False, build=False, final_wht_type=p['final_wht_type'], final_wcs=True)
        match = glob.glob(prelim_prefix + '_dr?_sci.fits')
        if not match:
            raise FileNotFoundError(f"Expected preliminary single drizzle not found: {prelim_prefix}_dr?_sci.fits")
        prelim_singles.append(match[0])

    protect_sky = pipeline_qc.repeat_detected_protection_mask(
        prelim_singles[0], prelim_singles[1],
        match_radius_arcsec=sp['match_radius_arcsec'],
        fwhm=sp['detect_fwhm_px'], nsigma=sp['detect_nsigma'])
    print(f"  [Defect 5b] {len(protect_sky)} source(s) confirmed by repeat detection; "
          "protecting from CR flagging.")

    for flc in flc_files:
        with fits.open(flc) as hdul:
            n_sci = sum(1 for h in hdul if h.name == 'SCI')
            hdr0 = hdul[0].header
        readnoise_vals = [hdr0.get(k) for k in ('READNSEA', 'READNSEB', 'READNSEC', 'READNSED')
                          if hdr0.get(k) is not None]
        readnoise = float(np.mean(readnoise_vals)) if readnoise_vals else 3.1
        for ext_ver in range(1, n_sci + 1):
            with fits.open(flc) as hdul:
                data = np.asarray(hdul['SCI', ext_ver].data, dtype=float)
            cr_mask = pipeline_qc.lacosmic_flag(data, gain=1.0, readnoise=readnoise, **la)
            protect_mask = None
            if len(protect_sky):
                det_xy = pipeline_qc.sky_to_detector_pixels(protect_sky, flc)
                # sky_to_detector_pixels uses the SCI,1 WCS; for SCI,2 (chip
                # 2 of UVIS) re-resolve against that extension specifically.
                if ext_ver != 1:
                    with fits.open(flc) as hdul:
                        wcs_ext = WCS(hdul['SCI', ext_ver].header, hdul)
                    x, y = wcs_ext.world_to_pixel(protect_sky)
                    det_xy = np.column_stack([x, y])
                protect_mask = pipeline_qc.build_source_protection_mask(
                    data.shape, det_xy, radius_px=sp['protect_radius_px'])
            n = pipeline_qc.apply_custom_cr_dq(flc, ext_ver, cr_mask, protect_mask=protect_mask, bit=8192)
            print(f"    {Path(flc).name}[SCI,{ext_ver}]: {n} px flagged CR (bit 8192, protection applied)")

    # prelim_singles / their _wht/_context siblings live in temp_dir and are
    # removed along with everything else by uv_drizzler's safe_rmtree(temp_dir).

    # --- (c): driz_cr pass, combine_type pinned to minmed, protecting our
    # new bit-8192 flags AND the archive's bit-4096 flags (resetbits=0).
    astrodrizzle.AstroDrizzle(
        flc_files, output=output_prefix, preserve=False,
        resetbits=0, driz_cr_corr=True, combine_type='minmed',
        driz_cr_scale=p['driz_cr_scale'], driz_cr_grow=p['driz_cr_grow'],
        final_wht_type=p['final_wht_type'], final_pixfrac=p['final_pixfrac_initial'])
    pipeline_qc.assert_dq_flag_fraction(
        flc_files, bit=8192, min_frac=0.0, label='post L.A.Cosmic pass (bit 8192)')

    shiftfile = (temp_dir / 'shiftuv_flt.txt').as_posix()
    # See pipeline_qc.compute_skysigma (Defect 6): same count-rate/auto-sigma
    # fix as the other TweakReg call sites; dq_mask_bits matches this
    # branch's own dqbits (archive CR bit 4096 + our L.A.Cosmic bit 8192).
    skysigma = pipeline_qc.compute_skysigma(
        list_fits(temp_dir, '*crclean.fits'), dq_mask_bits=(4096 | 8192))
    tweakreg.TweakReg(
        (temp_dir / '*crclean.fits').as_posix(),
        enforce_user_order=False,
        imagefindcfg={'threshold': p['tweak_threshold'], 'conv_width': p['tweak_conv_width'],
                      'dqbits': ~(4096 | 8192), 'use_sharp_round': True,
                      'computesig': False, 'skysigma': skysigma},
        refimagefindcfg={'conv_width': p['tweak_refconv_width'],
                         'dqbits': ~(4096 | 8192), 'use_sharp_round': True,
                         'computesig': False, 'skysigma': skysigma},
        shiftfile=True, outshifts=shiftfile,
        searchrad=p['tweak_searchrad'], ylimit=p['tweak_ylimit'],
        use2dhist=p['tweak_use2dhist'], xoffset=p['tweak_xoffset'],
        yoffset=p['tweak_yoffset'], tolerance=p['tolerance'],
        separation=p['tweak_separation'], updatehdr=True,
        wcsname='UVIS_FLC', reusename=True, interactive=False)
    pipeline_qc.check_shiftfile_for_nan(shiftfile)

    for flc in flc_files:
        flid = derive_rootname(flc)
        crclean = temp_dir / f'{flid}_crclean.fits'
        tweakback.tweakback(crclean.as_posix(), input=flc, wcsname='UVIS_FLC')

    astrodrizzle.AstroDrizzle(
        flc_files, output=output_prefix, preserve=False,
        resetbits=0, driz_cr_corr=False, combine_type='minmed',
        final_pixfrac=p['final_pixfrac_final'], final_wcs=True, final_scale=p['pixel_size'])
    pipeline_qc.assert_dq_flag_fraction(
        flc_files, bit=8192, min_frac=0.0, label='pre pairwise-rejection combined UV drizzle')

    # --- (d): pairwise rejection in the output frame ------------------------
    combined_sci = glob.glob(output_prefix + '_dr?_sci.fits')
    if not combined_sci:
        raise FileNotFoundError(f"Expected combined drizzle not found: {output_prefix}_dr?_sci.fits")
    combined_sci = combined_sci[0]

    single_paths = []
    for flc in flc_files:
        single_prefix = str(Path(flc).with_name(f'{Path(flc).stem}_single'))
        astrodrizzle.AstroDrizzle(
            [flc], output=single_prefix, resetbits=0, driz_cr_corr=False,
            preserve=False, build=False, final_wht_type=p['final_wht_type'],
            final_pixfrac=p['final_pixfrac_final'],
            final_wcs=True, final_refimage=combined_sci)
        match = glob.glob(single_prefix + '_dr?_sci.fits')
        if not match:
            raise FileNotFoundError(f"Expected single-exposure drizzle not found: {single_prefix}_dr?_sci.fits")
        single_paths.append(match[0])

    n_flagged = pipeline_qc.pairwise_reject_output_frame(
        single_paths, flc_files,
        threshold_k=pw['threshold_k'], dilation_px=pw['dilation_px'],
        # Reuse step (b)'s confirmed-source catalog: without this, a real
        # run on DESJ0206-0114 showed this step clipping real flux from
        # sources that (a)-(c) had correctly protected -- see
        # pairwise_reject_output_frame's docstring.
        protect_sky=protect_sky, protect_radius_px=sp['protect_radius_px'])
    print(f"  [Defect 5d] pairwise output-frame rejection flagged "
          f"{n_flagged} detector pixel(s) per exposure (bit 8192)")

    # --- final re-drizzle, now including the pairwise flags -----------------
    astrodrizzle.AstroDrizzle(
        flc_files, output=output_prefix, preserve=False,
        resetbits=0, driz_cr_corr=False, combine_type='minmed',
        final_pixfrac=p['final_pixfrac_final'], final_wcs=True, final_scale=p['pixel_size'])
    pipeline_qc.assert_dq_flag_fraction(
        flc_files, bit=8192, min_frac=0.0, label='final UV drizzle (post pairwise rejection)')


def _drz_suffix(band: str) -> str:
    return 'drz' if band == 'F140W' else 'drc'


def reconcile_target_astrometry(target: str, proposal_id: str, camera: str, filters: list):
    """
    Defect 1: after both filters for a target have been drizzled, check
    which one carries the trustworthy (a posteriori) absolute astrometric
    solution and correct the other one's WCS to match.

    This has to run against the *raw archive* flt/flc files for the
    WCSTYPE check (pipeline_qc.classify_wcs), since that's where the
    metadata drizzlepac's headerlets actually carry, and against the
    *drizzled native* _sci.fits for the correction (see
    pipeline_qc.reconcile_filter_astrometry's docstring for why it can't be
    the run_step3 display products instead).

    Only runs when this call actually processed an IR (F140W) filter and
    exactly one other (UV) filter for this target -- for any other
    combination there's nothing well-defined to reconcile automatically.
    """
    if 'F140W' not in filters:
        return
    uv_filters = [f for f in filters if f != 'F140W']
    if len(uv_filters) != 1:
        return
    uv_band = uv_filters[0]

    ir_raw = list_fits(MAIN_DIR / f'{target}/HST/{proposal_id}_F140W/raw_data', '*flt.fits')
    uv_raw = list_fits(MAIN_DIR / f'{target}/HST/{proposal_id}_{uv_band}/raw_data', '*flc.fits')
    if not ir_raw or not uv_raw:
        print(f"  [Defect 1] {target}: raw exposures not found for both filters — skipping preflight.")
        return

    try:
        reference, report = pipeline_qc.select_astrometric_reference(ir_raw, uv_raw)
    except pipeline_qc.AmbiguousAstrometricReference as exc:
        print(f"  [Defect 1] {target}: cannot automatically pick an astrometric "
              f"reference — leaving both filters' native WCS untouched.\n{exc}")
        return

    ir_sci = MAIN_DIR / f'{target}/HST/{proposal_id}_F140W/{target}_F140W_{camera}_{_drz_suffix("F140W")}_sci.fits'
    uv_sci = MAIN_DIR / f'{target}/HST/{proposal_id}_{uv_band}/{target}_{uv_band}_{camera}_{_drz_suffix(uv_band)}_sci.fits'
    if not ir_sci.exists() or not uv_sci.exists():
        print(f"  [Defect 1] {target}: drizzled science file(s) not found — skipping reconciliation.")
        return

    print(f"  [Defect 1] {target}: astrometric reference = {reference} "
          f"({report['ir' if reference == 'IR' else 'uv'][0]['wcstype']!r})")

    try:
        if reference == 'UV':
            pipeline_qc.reconcile_filter_astrometry(reference_sci_path=uv_sci, moving_sci_path=ir_sci)
        else:
            pipeline_qc.reconcile_filter_astrometry(reference_sci_path=ir_sci, moving_sci_path=uv_sci)
    except RuntimeError as exc:
        print(f"  [Defect 1] {target}: astrometric correction FAILED and was reverted — "
              f"filters remain misaligned, needs manual review.\n{exc}")


def run_step2(targets, filters, proposal_id, camera, interactive=False):
    print("\n=== STEP 2: Drizzle ===")

    ir_params = dict(IR_DRIZZLE_DEFAULTS)
    uv_params = dict(UV_DRIZZLE_DEFAULTS)

    if interactive:
        print("\nThese AstroDrizzle/TweakReg settings apply to every target in this run.")
        print("(Defaults live in config.py — edit IR_DRIZZLE_DEFAULTS / UV_DRIZZLE_DEFAULTS")
        print("to change them permanently.)")
        if 'F140W' in filters and _prompt_yes_no("Review IR (F140W) drizzle parameters?", default=False):
            ir_params = _prompt_edit_params('IR (F140W)', ir_params)
        if any(b != 'F140W' for b in filters) and _prompt_yes_no("Review UV drizzle parameters?", default=False):
            uv_params = _prompt_edit_params('UV', uv_params)

    skip_confirm = not interactive
    for row in targets:
        target = row[0]
        if target.endswith('B'):
            continue
        for band in filters:
            print(f"\n  Processing {target} {band} {camera}")

            if interactive and not skip_confirm:
                ans = input("    Proceed? [Y]es / [n]o (skip) / [a]ll (stop asking) / [q]uit: ").strip().lower()
                if ans == 'q':
                    print("  Stopping step 2 at user request.")
                    return
                if ans == 'n':
                    print(f"  Skipping {target} {band}.")
                    continue
                if ans == 'a':
                    skip_confirm = True

            if band == 'F140W':
                ir_drizzler(target, proposal_id, band, camera, params=ir_params)
            else:
                uv_drizzler(target, proposal_id, band, camera, params=uv_params)

        # Defect 1: only meaningful once every filter processed in this run
        # for this target has actually finished drizzling.
        reconcile_target_astrometry(target, proposal_id, camera, filters)

# ── Step 3: Reproject ──────────────────────────────────────────────────────────

def run_step3(target_filter=None):
    """
    Reproject all UV frames (F200LP, F606W) onto the pixel grid of whichever
    F140W image is available, then rename all images to the *_img_L1.fits scheme.

    Defect 1 note: this always uses F140W's native header as the resampling
    grid, which looks backwards given that F200LP more often carries the
    trustworthy (a posteriori) astrometric solution for GO-16773. It's safe
    anyway: run_step2's reconcile_target_astrometry() already corrects
    F140W's *native* _drz_sci.fits header (the file this function reads as
    ref_fits1/ref_fits2) to agree with whichever filter's preflight check
    picked as reference, before this step ever runs. Reprojecting onto an
    already-corrected F140W header does not reintroduce the misalignment.
    If step 2 was skipped for a target (e.g. `--steps 3` on its own) this
    reprojection inherits whatever astrometry the F140W file already has --
    run the preflight check manually first if you're not sure it was fixed.

    KNOWN LIMITATION: the three (proposal_id, filter) combinations below are
    hardcoded, not read from config.py's ALL_PROPOSALS (nothing reads that
    list -- see its comment in config.py). A fourth proposal or a filter
    other than F140W/F200LP/F606W is silently never reprojected -- no error,
    it just never enters the glob patterns below. Add the new combination
    here (and in run_step4 and hst_products.py's cam_map / 3-colour block)
    if you add one. See README.md's "Known limitations" section.
    """
    print("\n=== STEP 3: Reproject UV → IR pixel scale ===")

    rename = True  # set False to skip the _img_L1.fits rename step

    for target_dir in sorted(p for p in MAIN_DIR.iterdir() if p.is_dir()):
        target = target_dir.name
        if target_filter and target != target_filter:
            continue

        ref_fits1 = list(target_dir.glob('HST/16773_F140W/*_drz_sci.fits'))
        ref_fits2 = list(target_dir.glob('HST/15867_F140W/*_drz_sci.fits'))
        fits_file1 = list(target_dir.glob('HST/16773_F200LP/*_drc_sci.fits'))
        fits_file2 = list(target_dir.glob('HST/17307_F606W/*_drc_sci.fits'))

        # Reproject F200LP onto F140W grid
        if fits_file1:
            ref = ref_fits1 or ref_fits2
            if ref:
                hdu_ref = fits.open(ref[0])[0]
                hdu_uv  = fits.open(fits_file1[0])[0]
                array, _ = reproject_interp(hdu_uv, hdu_ref.header)
                out = target_dir / f'HST/16773_F200LP/{target}_F200LP_WFC3_drc_img_scaled_L3.fits'
                fits.writeto(str(out), array, hdu_ref.header, overwrite=True)
                lbl = '16773' if ref_fits1 else '15867'
                print(f"  {target}: F200LP → {lbl}_F140W reprojected")

        # Reproject F606W onto F140W grid
        if fits_file2:
            ref = ref_fits1 or ref_fits2
            if ref:
                hdu_ref = fits.open(ref[0])[0]
                hdu_uv  = fits.open(fits_file2[0])[0]
                array, _ = reproject_interp(hdu_uv, hdu_ref.header)
                out = target_dir / f'HST/17307_F606W/{target}_F606W_ACS_drc_img_scaled_L3.fits'
                fits.writeto(str(out), array, hdu_ref.header, overwrite=True)
                print(f"  {target}: F606W → F140W reprojected")

        if rename:
            for src, dst_name in [
                (ref_fits1,   f'HST/16773_F140W/{target}_F140W_WFC3_drz_img_L1.fits'),
                (ref_fits2,   f'HST/15867_F140W/{target}_F140W_WFC3_drz_img_L1.fits'),
                (fits_file1,  f'HST/16773_F200LP/{target}_F200LP_WFC3_drc_img_L1.fits'),
                (fits_file2,  f'HST/17307_F606W/{target}_F606W_ACS_drc_img_L1.fits'),
            ]:
                if src:
                    hdu = fits.open(src[0])[0]
                    fits.writeto(str(target_dir / dst_name), hdu.data, hdu.header, overwrite=True)

# ── Step 4: Green image ────────────────────────────────────────────────────────

def run_step4(target_filter=None):
    """
    Create a synthetic green channel: green = sqrt(IR² + UV²).
    Requires the scaled_L3.fits files produced by step 3.

    Same KNOWN LIMITATION as run_step3 above: '16773_F140W'/'15867_F140W'/
    '16773_F200LP'/'17307_F606W' are hardcoded glob patterns below, not
    driven by config.py's ALL_PROPOSALS.
    """
    print("\n=== STEP 4: Green image ===")

    for target_dir in sorted(p for p in MAIN_DIR.iterdir() if p.is_dir()):
        target = target_dir.name
        if target_filter and target != target_filter:
            continue
        if target.endswith('B'):
            continue

        # Locate IR reference
        ir_files = (list(target_dir.glob('HST/16773_F140W/*drz_sci.fits')) or
                    list(target_dir.glob('HST/15867_F140W/*drz_sci.fits')))
        if not ir_files:
            print(f"  {target}: no F140W found — skipping")
            continue

        # Locate UV scaled file (prefer F200LP, fall back to F606W)
        uv_file = (list(target_dir.glob('HST/16773_F200LP/*_scaled_L3.fits')) or
                   list(target_dir.glob('HST/17307_F606W/*_scaled_L3.fits')))
        if not uv_file:
            print(f"  {target}: no scaled UV file found — skipping")
            continue

        if list(target_dir.glob('HST/16773_F200LP/*_scaled_L3.fits')):
            outfile = target_dir / f'HST/16773_F200LP/{target}_F200LP_WFC3_green_img_scaled_L3.fits'
        else:
            outfile = target_dir / f'HST/17307_F606W/{target}_F606W_ACS_green_img_scaled_L3.fits'

        data_ir = fits.open(ir_files[0])[0].data
        data_uv = fits.open(uv_file[0])[0].data
        green   = np.sqrt(data_ir**2 + data_uv**2)

        hdu = fits.PrimaryHDU(green)
        hdu.header.update(WCS(fits.open(ir_files[0])[0].header).to_header())
        fits.HDUList([hdu]).writeto(str(outfile), overwrite=True)
        print(f"  {target}: green image written")

# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='HST reduction pipeline')
    parser.add_argument('--steps', default='1234',
                        help='Steps to run, e.g. "234" to skip download (default: 1234)')
    parser.add_argument('--target', default=None,
                        help='Run only this target (AGEL name)')
    parser.add_argument('--interactive', '-i', action='store_true',
                        help='Prompt to review/override AstroDrizzle+TweakReg parameters '
                             'and confirm each target before drizzling (step 2)')
    args = parser.parse_args()

    steps = set(args.steps)

    # Load target list from CSV
    if not ACTIVE_TARGETS_CSV.exists():
        raise SystemExit(
            f"Error: targets CSV not found at {ACTIVE_TARGETS_CSV}.\n"
            f"Check ACTIVE_PROPOSAL_ID in config.py (currently {ACTIVE_PROPOSAL_ID!r}) "
            f"and that {ACTIVE_PROPOSAL_ID}_targets.csv exists inside MAIN_DIR "
            f"({MAIN_DIR}). Copy targets_template.csv to create one.")
    df = pd.read_csv(ACTIVE_TARGETS_CSV)
    targets = []  # list of (agel_name, mast_name, ra, dec, z_def, z_src)
    for _, row in df.iterrows():
        agel = row['objname']
        mast = row['catalogue_objname']
        if args.target and agel != args.target:
            continue
        z_def, z_src = _resolve_z(row, SRC_COLS, DEF_COLS)
        targets.append((agel, mast, row.get('RAJ2000'), row.get('DECJ2000'), z_def, z_src))

    if args.target and not targets:
        raise SystemExit(
            f"--target {args.target!r} not found in 'objname' column of {ACTIVE_TARGETS_CSV}. "
            f"Check for typos or a mismatched ACTIVE_PROPOSAL_ID in config.py.")

    print(f"Loaded {len(targets)} target(s) from {ACTIVE_TARGETS_CSV}")
    print(f"Proposal {ACTIVE_PROPOSAL_ID}  |  filters {ACTIVE_FILTERS}  |  camera {ACTIVE_CAMERA}")
    print(f"Running steps: {', '.join(sorted(steps))}")

    if args.interactive:
        if not _prompt_yes_no("\nProceed with this run?", default=True):
            print("Aborted.")
            return

    if '1' in steps:
        run_step1(targets, ACTIVE_FILTERS, ACTIVE_PROPOSAL_ID)
    if '2' in steps:
        run_step2(targets, ACTIVE_FILTERS, ACTIVE_PROPOSAL_ID, ACTIVE_CAMERA, interactive=args.interactive)
    if '3' in steps:
        run_step3(target_filter=args.target)
    if '4' in steps:
        run_step4(target_filter=args.target)

    print("\nReduction pipeline complete.")


if __name__ == '__main__':
    main()
