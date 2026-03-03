"""
hst_reduction.py — HST data reduction pipeline

Runs four steps in sequence for the proposal defined in config.py:
  1. Download  : fetch FLC/FLT files from MAST
  2. Drizzle   : align and drizzle IR (F140W) and UV (all other) exposures
  3. Reproject : resample UV frames onto the IR pixel grid, rename to *_L1.fits
  4. Green     : build synthetic green image from IR + UV for RGB composites

Usage:
    python hst_reduction.py [--steps 1234] [--target AGEL...]

    --steps   string of step numbers to run (default: all, e.g. '234' skips download)
    --target  run only one target (useful for testing / re-running a single object)
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
    MAIN_DIR,
)

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

def download(lens_name: str, target_name: str, band: str, proposal_id: str):
    """Download FLC/FLT files for one target from MAST."""
    from astroquery.mast import Observations  # import here so the rest of the script works without astroquery

    download_dir = MAIN_DIR / f'{lens_name}/HST/{proposal_id}_{band}/raw_data'
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


def run_step1(targets, filters, proposal_id):
    print("\n=== STEP 1: Download ===")
    for agel_name, mast_name in targets:
        if agel_name.endswith('B'):
            continue
        for band in filters:
            print(f"  Downloading {agel_name} / {band}")
            download(agel_name, mast_name, band, proposal_id)

# ── Step 2: Drizzle ────────────────────────────────────────────────────────────

def ir_drizzler(object_name: str, propno: str, band: str, camera: str, pixel_size: float):
    """Drizzle IR (WFC3/IR FLT) exposures."""
    from drizzlepac import astrodrizzle, tweakreg, tweakback

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
        print(f"[WARN] Only {len(flt_files)} FLT file(s); skipping initial CR rejection.")
        tweakreg.TweakReg(
            (temp_dir / '*flt.fits').as_posix(),
            imagefindcfg={'threshold': 50, 'conv_width': 2.5},
            expand_refcat=True, enforce_user_order=False, shiftfile=True,
            outshifts=(temp_dir / 'shiftir_flt.txt').as_posix(),
            searchrad=2.0, ylimit=0.3, updatehdr=False,
            reusename=True, wcsname='IR_FLT', interactive=False)
    else:
        astrodrizzle.AstroDrizzle(
            flt_files, output=output_prefix,
            resetbits=4096, driz_cr_corr=True,
            final_wht_type='EXP', final_pixfrac=1.0)

        tweakreg.TweakReg(
            (temp_dir / '*crclean.fits').as_posix(),
            imagefindcfg={'threshold': 5, 'conv_width': 3.5, 'dqbits': ~4096},
            expand_refcat=True, enforce_user_order=False, shiftfile=True,
            outshifts=(temp_dir / 'shiftir_flt.txt').as_posix(),
            searchrad=2.0, ylimit=0.3, updatehdr=True,
            reusename=True, wcsname='IR_FLT', interactive=False)

        for flt in flt_files:
            flid = derive_rootname(flt)
            crclean = temp_dir / f'{flid}_crclean.fits'
            if not crclean.exists():
                print(f"[WARN] Missing {crclean.name}; skipping tweakback.")
                continue
            tweakback.tweakback(crclean.as_posix(), input=flt, wcsname='IR_FLT')

    astrodrizzle.AstroDrizzle(
        flt_files, output=output_prefix,
        resetbits=8192, driz_cr_corr=False,
        final_pixfrac=1.0, final_wcs=True, final_scale=pixel_size)

    gc.collect()
    clean_mac_junk_recursive(temp_dir)
    ensure_cwd_not_inside(temp_dir)
    time.sleep(0.1)
    safe_rmtree(temp_dir)
    print(f"[CLEAN] Removed temp: {temp_dir}")


def uv_drizzler(object_name: str, propno: str, band: str, camera: str, pixel_size: float):
    """Drizzle UV/UVIS (ACS/WFC3 FLC) exposures."""
    from drizzlepac import astrodrizzle, tweakreg, tweakback

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
        astrodrizzle.AstroDrizzle(
            flc_files, output=output_prefix,
            resetbits=8192, driz_cr_scale='0.9 0.6', driz_cr_grow=3,
            final_wht_type='EXP', final_pixfrac=1.0,
            final_wcs=True, final_scale=pixel_size)
    else:
        astrodrizzle.AstroDrizzle(
            flc_files, output=output_prefix,
            resetbits=8192, driz_cr_corr=True,
            driz_cr_scale='0.9 0.6', driz_cr_grow=3,
            final_wht_type='EXP', final_pixfrac=0.8)

        tweakreg.TweakReg(
            (temp_dir / '*crclean.fits').as_posix(),
            enforce_user_order=False,
            imagefindcfg={'threshold': 20, 'conv_width': 3.5, 'dqbits': ~8192},
            refimagefindcfg={'conv_width': 2.5},
            shiftfile=True,
            outshifts=(temp_dir / 'shiftuv_flt.txt').as_posix(),
            searchrad=5.0, ylimit=0.6, updatehdr=True,
            wcsname='UVIS_FLC', reusename=True, interactive=False)

        for flc in flc_files:
            flid = derive_rootname(flc)
            crclean = temp_dir / f'{flid}_crclean.fits'
            tweakback.tweakback(crclean.as_posix(), input=flc, wcsname='UVIS_FLC')

        astrodrizzle.AstroDrizzle(
            flc_files, output=output_prefix,
            resetbits=8192, driz_cr_corr=False,
            driz_cr_scale='0.9 0.6', driz_cr_grow=3,
            final_pixfrac=1.0, final_wcs=True, final_scale=pixel_size)

    gc.collect()
    clean_mac_junk_recursive(temp_dir)
    ensure_cwd_not_inside(temp_dir)
    time.sleep(0.1)
    safe_rmtree(temp_dir)
    print(f"[CLEAN] Removed temp: {temp_dir}")


def run_step2(targets, filters, proposal_id, camera):
    print("\n=== STEP 2: Drizzle ===")
    for row in targets:
        target = row[0]
        if target.endswith('B'):
            continue
        for band in filters:
            print(f"\n  Processing {target} {band} {camera}")
            if band == 'F140W':
                ir_drizzler(target, proposal_id, band, camera, pixel_size=0.08)
            else:
                uv_drizzler(target, proposal_id, band, camera, pixel_size=0.05)

# ── Step 3: Reproject ──────────────────────────────────────────────────────────

def run_step3():
    """
    Reproject all UV frames (F200LP, F606W) onto the pixel grid of whichever
    F140W image is available, then rename all images to the *_img_L1.fits scheme.
    """
    print("\n=== STEP 3: Reproject UV → IR pixel scale ===")

    rename = True  # set False to skip the _img_L1.fits rename step

    for target_dir in sorted(p for p in MAIN_DIR.iterdir() if p.is_dir()):
        target = target_dir.name

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

def run_step4():
    """
    Create a synthetic green channel: green = sqrt(IR² + UV²).
    Requires the scaled_L3.fits files produced by step 3.
    """
    print("\n=== STEP 4: Green image ===")

    for target_dir in sorted(p for p in MAIN_DIR.iterdir() if p.is_dir()):
        target = target_dir.name
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
    args = parser.parse_args()

    steps = set(args.steps)

    # Load target list from CSV
    df = pd.read_csv(ACTIVE_TARGETS_CSV)
    targets = []  # list of (agel_name, mast_name, ra, dec, z_def, z_src)
    for _, row in df.iterrows():
        agel = row['objname']
        mast = row['catalogue_objname']
        if args.target and agel != args.target:
            continue
        z_def, z_src = _resolve_z(row, SRC_COLS, DEF_COLS)
        targets.append((agel, mast, row.get('RAJ2000'), row.get('DECJ2000'), z_def, z_src))

    print(f"Loaded {len(targets)} target(s) from {ACTIVE_TARGETS_CSV}")
    print(f"Running steps: {', '.join(sorted(steps))}")

    if '1' in steps:
        run_step1(targets, ACTIVE_FILTERS, ACTIVE_PROPOSAL_ID)
    if '2' in steps:
        run_step2(targets, ACTIVE_FILTERS, ACTIVE_PROPOSAL_ID, ACTIVE_CAMERA)
    if '3' in steps:
        run_step3()
    if '4' in steps:
        run_step4()

    print("\nReduction pipeline complete.")


if __name__ == '__main__':
    main()
