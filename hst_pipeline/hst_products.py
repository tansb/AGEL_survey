"""
hst_products.py — HST science-products pipeline

Runs three steps in sequence using the drizzled data produced by hst_reduction.py:
  1. Cutouts      : build lenstronomy-ready HDF5 + FITS cutout files (nb 3)
  2. Offset fitter: fit Sersic centroids to measure inter-filter astrometric offsets,
                    write ra_shift / dec_shift back into the cutout files (nb 7)
  3. Postage stamps: make grayscale / 2-colour / 3-colour PNG thumbnails (nb 6)

Usage:
    python hst_products.py [--steps 123] [--target AGEL...] [--preview]
                           [--ra DEG] [--dec DEG] [--z-src Z] [--z-def Z]
                           [--proposal-id ID]

    --steps       string of step numbers to run (default: all, e.g. '23' skips cutouts)
    --target      run only this target (AGEL name); uses the full CSV if omitted
    --ra / --dec  coordinates in decimal degrees — required when --target is given and
                  the target is not already in the CSV
    --z-src       source redshift (optional, stored as -1 if omitted)
    --z-def       deflector redshift (optional, stored as -1 if omitted)
    --proposal-id override ACTIVE_PROPOSAL_ID from config.py for this run
    --preview     show a contact-sheet grid of proposed cutout sizes before running;
                  prompts for confirmation before proceeding

Single-target examples:
    # Target already in CSV
    python hst_products.py --target AGEL110725+245943A

    # Target not yet in CSV — provide coordinates directly
    python hst_products.py --target AGEL110725+245943A --ra 166.855 --dec 24.995

    # Same, and append the new target to the CSV for future runs
    # (the script will prompt interactively)
"""

import argparse
import collections
import glob
import json
import math
import os
from pathlib import Path
from typing import Optional, Tuple

import h5py
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clipped_stats
from astropy.visualization import simple_norm
from astropy.wcs import WCS

from config import (
    ACTIVE_CAMERA,
    ACTIVE_FILTERS,
    ACTIVE_PROPOSAL_ID,
    ACTIVE_TARGETS_CSV,
    DATA_DIR,
    DEFAULT_CUTOUT_SIZE_ARCSEC,
    LENS_PROC_DIR,
    MAIN_DIR,
    OFFSET_NUM_DUP,
    PARENT_CATALOGUE_CSV,
    POSTAGE_STAMP_CONFIG_CSV,
    POSTAGE_STAMP_DIR,
    PROPOSAL_CSVS,
)

# ── Shared helpers ─────────────────────────────────────────────────────────────

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


SRC_COLS = ['z_spec_SRC_Spectral_Observations_Tally', 'z_source (from DR2 Redshifts)', 'z_spec_SRC']
DEF_COLS = ['z_spec_DE_Spectral_Observations_Tally',  'z_deflector (from DR2 Redshifts)',  'z_spec_DE']


def _resolve_z(row, src_cols=SRC_COLS, def_cols=DEF_COLS):
    z_src = next((_float_or_none(row[c]) for c in src_cols if c in row and _float_or_none(row[c]) is not None), None)
    z_def = next((_float_or_none(row[c]) for c in def_cols if c in row and _float_or_none(row[c]) is not None), None)
    return z_def, z_src


def _load_hdf5_2d(path, dataset=None):
    with h5py.File(path, "r") as f:
        if dataset:
            arr = f[dataset][...]
            if arr.ndim != 2:
                raise ValueError(f"Dataset '{dataset}' is {arr.ndim}D, expected 2D")
            return np.asarray(arr, dtype="float64")
        for k in ("kernel", "psf", "psf_kernel", "data", "image"):
            if k in f and isinstance(f[k], h5py.Dataset) and f[k].ndim == 2:
                return np.asarray(f[k][...], dtype="float64")
        best = None; best_n = -1
        def walk(g):
            nonlocal best, best_n
            for name, obj in g.items():
                if isinstance(obj, h5py.Dataset) and obj.ndim == 2:
                    n = int(np.prod(obj.shape))
                    if n > best_n:
                        best, best_n = obj, n
                elif isinstance(obj, h5py.Group):
                    walk(obj)
        walk(f)
        if best is None:
            raise ValueError(f"No 2D dataset found in {path}")
        return np.asarray(best[...], dtype="float64")


def _load_fits_2d(path):
    with fits.open(path, memmap=False) as hdul:
        for h in hdul:
            if h.data is None:
                continue
            arr = np.asarray(h.data)
            if arr.ndim == 2:
                return arr.astype("float64")
    raise ValueError(f"{path} had no 2D image HDU")


def _load_psf_any(psf_path, dataset=None, normalize=False, ensure_odd=False, center_peak=False):
    if psf_path is None:
        return None
    p = str(psf_path).lower()
    if p.endswith((".fits", ".fit", ".fz")):
        psf = _load_fits_2d(psf_path)
    elif p.endswith((".h5", ".hdf5")):
        psf = _load_hdf5_2d(psf_path, dataset=dataset)
    else:
        raise ValueError(f"Unsupported PSF format: {psf_path}")
    if ensure_odd:
        ny, nx = psf.shape
        if ny % 2 == 0 or nx % 2 == 0:
            psf = np.pad(psf, ((0, ny % 2 == 0), (0, nx % 2 == 0)), mode="constant")
    if center_peak:
        try:
            from scipy.ndimage import fourier_shift
            yy, xx = np.unravel_index(np.nanargmax(psf), psf.shape)
            cy, cx = (np.array(psf.shape) - 1) / 2.0
            shift = (cy - yy, cx - xx)
            psf = np.fft.ifftn(fourier_shift(np.fft.fftn(psf), shift)).real
        except Exception:
            pass
    if normalize:
        s = np.nansum(psf)
        if not np.isfinite(s) or s == 0:
            raise ValueError("PSF sum is zero/NaN; cannot normalize.")
        psf = psf / s
    return psf


def _load_weight_any(weight_path, dataset=None):
    if weight_path is None:
        return None
    p = str(weight_path).lower()
    if p.endswith((".fits", ".fit", ".fz")):
        return _load_fits_2d(weight_path)
    elif p.endswith((".h5", ".hdf5")):
        return _load_hdf5_2d(weight_path, dataset=dataset)
    else:
        raise ValueError(f"Unsupported weight format: {weight_path}")


def _read_cd_or_pc(header, force_pc=False):
    def has_cd(h): return all(k in h for k in ("CD1_1","CD1_2","CD2_1","CD2_2"))
    def has_pc(h): return all(k in h for k in ("PC1_1","PC1_2","PC2_1","PC2_2")) and ("CDELT1" in h and "CDELT2" in h)
    if (not force_pc) and has_cd(header):
        return float(header["CD1_1"]), float(header["CD1_2"]), float(header["CD2_1"]), float(header["CD2_2"])
    elif has_pc(header):
        c1, c2 = float(header["CDELT1"]), float(header["CDELT2"])
        return (c1*float(header["PC1_1"]), c1*float(header["PC1_2"]),
                c2*float(header["PC2_1"]), c2*float(header["PC2_2"]))
    else:
        raise ValueError("No CD or PC*CDELT found in header.")


def _build_transform(CD1_1, CD1_2, CD2_1, CD2_2):
    return np.array([[CD1_1, CD1_2],[CD2_1, CD2_2]], float) * 3600.0


def _centered_frame_ra_dec(nx, ny, A):
    dra, ddec = A.dot(np.array([int(nx/2), int(ny/2)], float))
    return -float(dra), -float(ddec)


def _estimate_background_rms(data, sigma=3.0, maxiters=5):
    _, _, std = sigma_clipped_stats(data, sigma=sigma, maxiters=maxiters, grow=False)
    return float(std)


def _resolve_center(header, wcs, center_ra=None, center_dec=None, center_x=None, center_y=None):
    if center_x is not None and center_y is not None:
        return float(center_x), float(center_y)
    if center_ra is not None and center_dec is not None:
        x, y = wcs.world_to_pixel_values(center_ra, center_dec)
        return float(x), float(y)
    nx, ny = int(header.get("NAXIS1")), int(header.get("NAXIS2"))
    return float(nx/2.0), float(ny/2.0)


def _size_from_arcsec(size_arcsec, A):
    eff = float(np.sqrt(np.linalg.det(A.T @ A)) / np.sqrt(2.0))
    if not np.isfinite(eff) or eff <= 0:
        eff = float(np.median(np.abs(A)))
    return max(8, int(np.round(size_arcsec / eff)))

# ── Cutout preview ────────────────────────────────────────────────────────────

def preview_cutouts(targets_df, proposal_id, filters, camera, cutout_size_arcsec=None):
    """
    Display a contact-sheet grid showing the proposed cutout box (red) overlaid
    on a thumbnail of the drizzled science image for every target.

    Use this to verify that DEFAULT_CUTOUT_SIZE_ARCSEC is appropriate before
    committing to a full batch run. Adjust the value in config.py if any targets
    look too tight or too generous, then re-run the preview to confirm.
    """
    if cutout_size_arcsec is None:
        cutout_size_arcsec = DEFAULT_CUTOUT_SIZE_ARCSEC

    # Collect targets that have a science file on disk
    entries = []
    for _, row in targets_df.iterrows():
        object_name = row['objname']
        if object_name.endswith('B'):
            continue
        for band in filters:
            suffix = 'drz' if band == 'F140W' else 'drc'
            fits_path = MAIN_DIR / f'{object_name}/HST/{proposal_id}_{band}/{object_name}_{band}_{camera}_{suffix}_sci.fits'
            if fits_path.exists():
                entries.append((object_name, band, float(row.get('RAJ2000', 0)),
                                float(row.get('DECJ2000', 0)), str(fits_path)))
                break  # one band per target is enough for the preview
        else:
            print(f"  [preview] No science FITS found for {object_name} — skipped")

    if not entries:
        print("No science FITS files found. Run the drizzle step first.")
        return

    n = len(entries)
    ncols = min(4, n)
    nrows = math.ceil(n / ncols)

    fig, axes = plt.subplots(nrows, ncols, figsize=(4.5 * ncols, 4.5 * nrows),
                             squeeze=False)

    for idx, (object_name, band, ra_deg, dec_deg, fits_path) in enumerate(entries):
        ax = axes[idx // ncols][idx % ncols]
        pix_scale = 0.08 if band == 'F140W' else 0.05
        box_pix   = cutout_size_arcsec / pix_scale

        try:
            with fits.open(fits_path) as hdul:
                hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
                w   = WCS(hdu.header)
                cx, cy = w.world_to_pixel(SkyCoord(ra_deg * u.deg, dec_deg * u.deg))
                cx, cy = float(cx), float(cy)

                # Extract a region 1.6× the cutout size to give context around the box
                pad    = int(box_pix / 2 * 1.6)
                ny, nx = hdu.data.shape
                x0, x1 = max(0, int(cx) - pad), min(nx, int(cx) + pad)
                y0, y1 = max(0, int(cy) - pad), min(ny, int(cy) + pad)
                region = np.array(hdu.data[y0:y1, x0:x1], dtype=float)

            norm = simple_norm(region, stretch='sqrt', percent=99.5)
            ax.imshow(region, norm=norm, cmap='gray', origin='lower')

            # Draw the proposed cutout box in the thumbnail's coordinate space
            half  = box_pix / 2
            rect_x = (cx - x0) - half
            rect_y = (cy - y0) - half
            ax.add_patch(patches.Rectangle(
                (rect_x, rect_y), box_pix, box_pix,
                linewidth=1.5, edgecolor='red', facecolor='none'))
            ax.plot(cx - x0, cy - y0, 'r+', markersize=10, markeredgewidth=1.5)
            ax.set_title(f'{object_name}\n{band}  |  {cutout_size_arcsec:.0f}"',
                         fontsize=7, pad=3)

        except Exception as exc:
            ax.text(0.5, 0.5, f'{object_name}\nError loading:\n{exc}',
                    transform=ax.transAxes, ha='center', va='center',
                    fontsize=6, color='red')

        ax.axis('off')

    # Hide any unused grid cells
    for idx in range(n, nrows * ncols):
        axes[idx // ncols][idx % ncols].set_visible(False)

    fig.suptitle(
        f'Cutout preview  —  proposed size: {cutout_size_arcsec:.1f}"  (red box)\n'
        'Adjust DEFAULT_CUTOUT_SIZE_ARCSEC in config.py if needed, then re-run --preview to confirm.',
        fontsize=10)
    plt.tight_layout()
    plt.show()


# ── Step 1: Build cutouts ──────────────────────────────────────────────────────

def build_cutout_h5(
    fits_path, output_h5,
    size_pixels=None, size_arcsec=None,
    center_ra=None, center_dec=None, center_x=None, center_y=None,
    psf_path=None, weight_path=None,
    psf_normalize=False, psf_ensure_odd=False, psf_center_peak=False,
    background_rms=None, exposure_time=None,
    force_pc=False, overwrite=False,
    save_fits_with_header=False,
    z_src=None, z_def=None,
):
    if size_pixels is None and size_arcsec is None:
        raise ValueError("Specify either size_pixels or size_arcsec.")

    with fits.open(fits_path, memmap=False) as hdul:
        hdu = next(h for h in hdul if (h.data is not None and h.data.ndim == 2))
        header = hdu.header
        img = np.asarray(hdu.data).astype("float64")
        w = WCS(header)
        CD1_1, CD1_2, CD2_1, CD2_2 = _read_cd_or_pc(header, force_pc=force_pc)
        A = _build_transform(CD1_1, CD1_2, CD2_1, CD2_2)
        cx, cy = _resolve_center(header, w, center_ra, center_dec, center_x, center_y)
        if size_pixels is None:
            size_pixels = _size_from_arcsec(size_arcsec, A)
        cut = Cutout2D(img, (cx, cy), size=(int(size_pixels), int(size_pixels)), wcs=w, mode="trim", copy=True)
        data = cut.data.astype("float64")
        ny, nx = data.shape
        ra0, dec0 = _centered_frame_ra_dec(nx, ny, A)
        bkg = float(background_rms) if background_rms is not None else _estimate_background_rms(data)
        exptime = float(exposure_time) if exposure_time is not None else float(header.get("EXPTIME", 1.0))
        kernel_psf = _load_psf_any(psf_path, normalize=psf_normalize, ensure_odd=psf_ensure_odd, center_peak=psf_center_peak)
        weight_map = _load_weight_any(weight_path)

    out_path = Path(output_h5)
    if (not overwrite) and out_path.exists():
        raise FileExistsError(f"{output_h5} exists. Set overwrite=True.")
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with h5py.File(out_path, "w") as f:
        f.create_dataset("image_data",        data=data,                             dtype="float64")
        f.create_dataset("background_rms",    data=np.array(bkg,     dtype="float64"))
        f.create_dataset("exposure_time",     data=np.array(exptime,  dtype="float64"))
        f.create_dataset("ra_at_xy_0",        data=np.array(ra0,      dtype="float64"))
        f.create_dataset("dec_at_xy_0",       data=np.array(dec0,     dtype="float64"))
        f.create_dataset("transform_pix2angle", data=A,               dtype="float64")
        f.create_dataset("z_spec_SRC",        data=np.array(z_src if z_src is not None else -1.0, dtype="float64"))
        f.create_dataset("z_spec_DE",         data=np.array(z_def if z_def is not None else -1.0, dtype="float64"))
        if kernel_psf is not None:
            f.create_dataset("kernel_psf",    data=kernel_psf,        dtype="float64")
        if weight_map is not None:
            f.create_dataset("weight_map",    data=weight_map,        dtype="float64")
        meta = {
            "source_fits": os.path.abspath(fits_path),
            "size_pixels": int(size_pixels),
            "CD_matrix_deg_per_pix": [CD1_1, CD1_2, CD2_1, CD2_2],
            "transform_pix2angle_arcsec": A.tolist(),
            "background_rms_method": "sigma_clipped_stats" if background_rms is None else "user_provided",
            "exptime_source": "FITS[EXPTIME]" if exposure_time is None else "user_provided",
            "psf_attached": kernel_psf is not None,
            "weight_attached": weight_map is not None,
            "z_spec_SRC": float(z_src) if z_src is not None else None,
            "z_spec_DE":  float(z_def) if z_def is not None else None,
        }
        f.create_dataset("meta", data=json.dumps(meta))

    if save_fits_with_header:
        with fits.open(fits_path, memmap=False) as hduL1:
            h2 = next(h for h in hduL1 if (h.data is not None and h.data.ndim == 2))
            w2 = WCS(h2.header)
            cut2 = Cutout2D(h2.data, (cx, cy), size=(int(size_pixels), int(size_pixels)), wcs=w2, mode="trim", copy=True)
        hdr = cut2.wcs.to_header()
        hdr["SIZE_PIX"] = (int(size_pixels), "cutout size in pixels")
        hdr["CEN_RA"]   = (center_ra,  "cutout center RA (deg)")
        hdr["CEN_DEC"]  = (center_dec, "cutout center Dec (deg)")
        hdr["CD1_1"]    = (CD1_1, "CD matrix [0,0] (deg/pix)")
        hdr["CD1_2"]    = (CD1_2, "CD matrix [0,1] (deg/pix)")
        hdr["CD2_1"]    = (CD2_1, "CD matrix [1,0] (deg/pix)")
        hdr["CD2_2"]    = (CD2_2, "CD matrix [1,1] (deg/pix)")
        hdr["LNS_RA0"]  = (ra0,  "lenstronomy ra_at_xy_0 (arcsec)")
        hdr["LNS_DE0"]  = (dec0, "lenstronomy dec_at_xy_0 (arcsec)")
        hdr["EXPTIME"]  = (exptime, "exposure time (s)")
        hdr["BKG_RMS"]  = (float(bkg), "background RMS (image units)")
        if z_src is not None:
            hdr["Z_SRC"] = (float(z_src), "source z_spec")
        if z_def is not None:
            hdr["Z_DEF"] = (float(z_def), "deflector z_spec")
        fits_name = str(out_path.with_suffix(".fits"))
        hdu_list = [fits.PrimaryHDU(data=data, header=hdr),
                    fits.ImageHDU(data=A.astype("float64"), name="PIX2ANG")]
        if kernel_psf is not None:
            hdu_list.append(fits.ImageHDU(data=kernel_psf.astype("float64"), name="PSF"))
        if weight_map is not None:
            hdu_list.append(fits.ImageHDU(data=weight_map.astype("float64"), name="WHT"))
        fits.HDUList(hdu_list).writeto(fits_name, overwrite=True)

    return str(out_path)


def run_step1(targets, proposal_id, filters, camera):
    print("\n=== STEP 1: Build cutouts ===")
    for row in targets:
        object_name = row['objname']
        if object_name.endswith('B'):
            continue

        ra_deg  = row.get('RAJ2000')
        dec_deg = row.get('DECJ2000')
        z_def, z_src = _resolve_z(row)

        for band in filters:
            pix_scale = 0.08 if band == 'F140W' else 0.05
            suffix    = 'drz' if band == 'F140W' else 'drc'
            fits_path = f"{MAIN_DIR}/{object_name}/HST/{proposal_id}_{band}/{object_name}_{band}_{camera}_{suffix}_sci.fits"
            if not Path(fits_path).exists():
                print(f"  {object_name} {band}: science FITS not found — skipping")
                continue

            output_h5 = f"{MAIN_DIR}/{object_name}/HST/{proposal_id}_{band}/{object_name}_{band}_{camera}_cutout_L3.h5"
            psf_path  = LENS_PROC_DIR / f"lens_processing/psf_model_{band}.h5"
            weight_path = f"{MAIN_DIR}/{object_name}/HST/{proposal_id}_{band}/{object_name}_{band}_{camera}_{suffix}_wht.fits"

            sc = SkyCoord(ra_deg * u.deg, dec_deg * u.deg)
            size_pixels = DEFAULT_CUTOUT_SIZE_ARCSEC / pix_scale

            print(f"  {object_name} {band}: building cutout ({DEFAULT_CUTOUT_SIZE_ARCSEC}\")")
            build_cutout_h5(
                fits_path, output_h5,
                size_pixels=size_pixels, size_arcsec=DEFAULT_CUTOUT_SIZE_ARCSEC,
                center_ra=sc.ra.deg, center_dec=sc.dec.deg,
                psf_path=str(psf_path) if psf_path.exists() else None,
                weight_path=weight_path if Path(weight_path).exists() else None,
                overwrite=True, save_fits_with_header=True,
                z_src=z_src, z_def=z_def,
            )

# ── Step 2: Offset fitter ──────────────────────────────────────────────────────

def _import_data(object_name, propno, band, camera):
    data_file = DATA_DIR / f"{object_name}/HST/{propno}_{band}/{object_name}_{band}_{camera}_cutout_L3.h5"
    if not data_file.exists():
        return None, None, None

    with h5py.File(data_file, 'r') as f:
        kwargs_data = {
            'image_data':        f['image_data'][()],
            'background_rms':    f['background_rms'][()],
            'noise_map':         None,
            'exposure_time':     f['exposure_time'][()],
            'ra_at_xy_0':        f['ra_at_xy_0'][()],
            'dec_at_xy_0':       f['dec_at_xy_0'][()],
            'transform_pix2angle': f['transform_pix2angle'][()],
        }

    psf_path = LENS_PROC_DIR / f"lens_processing/psf_model_{band}.h5"
    if psf_path.exists():
        with h5py.File(psf_path, 'r') as f:
            kernel = f['kernel_point_source'][()]
        kwargs_psf = {
            'psf_type': "PIXEL",
            'kernel_point_source': kernel,
            'kernel_point_source_init': kernel,
        }
    else:
        print(f"  [WARN] No PSF found for {band}")
        kwargs_psf = None

    return kwargs_data, kwargs_psf, None


def _process_lens(kwargs_data, kwargs_psf, kwargs_mask, radius, x_offset, y_offset):
    from lenstronomy.Data.coord_transforms import Coordinates
    from lenstronomy.Util import mask_util
    import lenstronomy.Util.util as util
    from lenstronomy.Util.util import array2image

    ra_at_xy_0      = kwargs_data['ra_at_xy_0']
    dec_at_xy_0     = kwargs_data['dec_at_xy_0']
    transform_pix2angle = kwargs_data['transform_pix2angle']
    coords = Coordinates(transform_pix2angle, ra_at_xy_0, dec_at_xy_0)
    data_cutout = kwargs_data['image_data']
    numPix = len(data_cutout)
    lens_center_ra, lens_center_dec = coords.map_pix2coord(numPix/2, numPix/2)

    x_coords, y_coords = coords.coordinate_grid(numPix, numPix)
    mask_outer = mask_util.mask_center_2d(
        lens_center_ra + x_offset, lens_center_dec + y_offset, radius,
        util.image2array(x_coords), util.image2array(y_coords))
    mask = 1 - mask_outer

    if kwargs_mask is not None:
        provided = np.abs(np.asarray(kwargs_mask) - 1)
        mask *= provided.flatten()
    mask = np.clip(mask, 0, 1)
    return array2image(mask)


def _fit_perturber(object_name, band, camera, kwargs_data, kwargs_psf, mask,
                   x_offset, y_offset, bound=0.5):
    import seaborn as sns
    from lenstronomy.Workflow.fitting_sequence import FittingSequence
    from lenstronomy.Plots.model_plot import ModelPlot

    kwargs_numerics = {'supersampling_factor': 1}
    kwargs_likelihood = {'check_bounds': True, 'image_likelihood_mask_list': [mask]}
    multi_band_list = [[kwargs_data, kwargs_psf, kwargs_numerics]]
    kwargs_data_joint = {'multi_band_list': multi_band_list, 'multi_band_type': 'multi-linear'}
    kwargs_model = {'lens_light_model_list': ['SERSIC_ELLIPSE']}

    fixed_ll = [{'n_sersic': 4.}]
    init_ll  = [{'R_sersic': .1, 'n_sersic': 4, 'e1': 0, 'e2': 0,
                 'center_x': x_offset, 'center_y': y_offset}]
    sigma_ll = [{'n_sersic': 0.5, 'R_sersic': 0.2, 'e1': 0.1, 'e2': 0.1,
                 'center_x': 0.5, 'center_y': 0.5}]
    lower_ll = [{'e1': -0.5, 'e2': -0.5, 'R_sersic': 0.01, 'n_sersic': 0.5,
                 'center_x': x_offset - bound, 'center_y': y_offset - bound}]
    upper_ll = [{'e1': 0.5, 'e2': 0.5, 'R_sersic': 10, 'n_sersic': 8,
                 'center_x': x_offset + bound, 'center_y': y_offset + bound}]

    kwargs_params = {'lens_light_model': [init_ll, sigma_ll, fixed_ll, lower_ll, upper_ll]}
    fitting_seq = FittingSequence(kwargs_data_joint, kwargs_model, {}, kwargs_likelihood, kwargs_params)
    fitting_seq.fit_sequence([['PSO', {'sigma_scale': 1., 'n_particles': 100, 'n_iterations': 50}]])
    result = fitting_seq.best_fit()

    ll = result['kwargs_lens_light'][0]
    return ll['n_sersic'], ll['R_sersic'], ll['center_x'], ll['center_y']


def _offset_fitter(object_name, propid_list, filter_list, camera_list, num_dup=4):
    print(f"  Fitting offsets: {object_name}")
    data_by_band = {}
    for i, (propid, band, cam) in enumerate(zip(propid_list, filter_list, camera_list)):
        kwargs_data, kwargs_psf, kwargs_mask = _import_data(object_name, propid, band, cam)
        if kwargs_data is None:
            print(f"    Skipping {band} {cam} — no cutout file.")
            return None

        mask = _process_lens(kwargs_data, kwargs_psf, kwargs_mask, 1.0, 0, 0)
        _, _, cx, cy = _fit_perturber(object_name, band, cam, kwargs_data, kwargs_psf, mask, 0, 0, bound=1.)

        mask = _process_lens(kwargs_data, kwargs_psf, kwargs_mask, 0.6, cx, cy)
        xs, ys = [], []
        for _ in range(num_dup):
            _, _, cx, cy = _fit_perturber(object_name, band, cam, kwargs_data, kwargs_psf, mask, cx, cy, bound=1.)
            xs.append(cx); ys.append(cy)
        data_by_band[i] = (np.mean(xs), np.mean(ys))

    for i, (propid, band, cam) in enumerate(zip(propid_list, filter_list, camera_list)):
        avg_x, avg_y = data_by_band[i]
        print(f"    {band} offset: ra={-avg_x:.4f}  dec={-avg_y:.4f} arcsec")

        h5_path = DATA_DIR / f"{object_name}/HST/{propid}_{band}/{object_name}_{band}_{cam}_cutout_L3.h5"
        with h5py.File(h5_path, 'a') as f:
            for key in ('ra_shift', 'dec_shift'):
                if key in f:
                    del f[key]
            f.create_dataset('ra_shift',  data=-avg_x)
            f.create_dataset('dec_shift', data=-avg_y)

        fits_path = DATA_DIR / f"{object_name}/HST/{propid}_{band}/{object_name}_{band}_{cam}_cutout_L3.fits"
        if fits_path.exists():
            with fits.open(fits_path, mode="update") as hdul:
                hdul[0].header['RA_SHFT'] = (-avg_x, 'ra shift applied')
                hdul[0].header['DE_SHFT'] = (-avg_y, 'dec shift applied')
                hdul.flush()


def run_step2(targets, proposal_id, filters, camera):
    print("\n=== STEP 2: Offset fitter ===")
    for _, row in targets.iterrows():
        target = row['objname']
        if target.endswith('B'):
            continue
        propid_list = [proposal_id] * len(filters)
        cam_list    = [camera] * len(filters)
        _offset_fitter(target, propid_list, filters, cam_list, num_dup=OFFSET_NUM_DUP)

# ── Step 3: Postage stamps ─────────────────────────────────────────────────────

label_size    = 20
tick_size     = 16
cbar_tick_size = 12
plt.rc('font', family='serif')


def _make_grayscale_postage(img1, output_dir, object_name, ra, dec,
                             zsrc1, zdef, band, zsrc2=None,
                             thumb_size_arcsec=20, scalebar_size_arcsec=1, clims_set=1):
    import aplpy
    fig = plt.figure(figsize=(12, 6))
    f  = aplpy.FITSFigure(img1, figure=fig, subplot=(1, 2, 1), north=True)
    f2 = aplpy.FITSFigure(img1, figure=fig, subplot=(1, 2, 2), north=True)

    if clims_set == 2:
        f.show_grayscale(invert=True, vmin=-0.2, vmax=2)
        f2.show_grayscale(invert=True, vmin=-0.1, vmax=0.4)
    else:
        f.show_grayscale(invert=True, vmin=-0.05, vmax=0.2)
        f2.show_grayscale(invert=True, vmin=-0.05, vmax=0.1)

    img_size = (thumb_size_arcsec * u.arcsec).to(u.deg).value
    for panel in (f, f2):
        panel.recenter(ra, dec, width=img_size, height=img_size)
        panel.axis_labels.hide(); panel.tick_labels.hide(); panel.ticks.hide()
        panel.add_scalebar(scalebar_size_arcsec * u.arcsec, label=f'{scalebar_size_arcsec}"', corner='bottom right')
        panel.scalebar.set_font(size=tick_size, weight='medium')
        panel.scalebar.set(color='k', linewidth=2)
        panel.ax.text(0.02, 0.93, band, color='k',
                      bbox=dict(facecolor='w', alpha=1), size=cbar_tick_size,
                      ha='left', va='bottom', transform=panel.ax.transAxes)

    title = (f'{object_name}     z_def = {zdef}     z_src1 = {zsrc1}     z_src2 = {zsrc2}'
             if zsrc2 else f'{object_name}    z_def = {zdef}     z_src = {zsrc1}')
    fig.suptitle(title, size=label_size - (3 if zsrc2 else 0))
    plt.subplots_adjust(wspace=0.05)
    plt.savefig(f'{output_dir}/{object_name}_{band}_{thumb_size_arcsec}arcs_image_L3.png', bbox_inches='tight')
    plt.close(fig)


_CLIMS = {
    1:  dict(vmin_r=-0.1, vmax_r=0.8,  vmin_b=-0.01, vmax_b=0.1,  vmin_g=-0.1, vmax_g=1),
    2:  dict(vmin_r=-0.05, vmax_r=1.1, vmin_b=0.07,  vmax_b=0.7,  vmin_g=-0.05, vmax_g=1.5),
    3:  dict(vmin_r=0.4,  vmax_r=2.,   vmin_b=0.0,   vmax_b=0.2,  vmin_g=0.5,  vmax_g=2.5),
    32: dict(vmin_r=0.,   vmax_r=3.,   vmin_b=0.0,   vmax_b=0.15, vmin_g=0.,   vmax_g=5.),
    4:  dict(vmin_r=-0.6, vmax_r=6.,   vmin_b=-0.01, vmax_b=0.45, vmin_g=-0.9, vmax_g=9.),
    5:  dict(vmin_r=-0.1, vmax_r=1.5,  vmin_b=0,     vmax_b=0.1,  vmin_g=-0.1, vmax_g=2.),
    6:  dict(vmin_r=-0.1, vmax_r=1.5,  vmin_b=0,     vmax_b=0.2,  vmin_g=-0.1, vmax_g=2.),
    7:  dict(vmin_r=0.0,  vmax_r=3.0,  vmin_b=0,     vmax_b=0.25, vmin_g=0.0,  vmax_g=4.),
    8:  dict(vmin_r=0.0,  vmax_r=3.0,  vmin_b=0,     vmax_b=0.3,  vmin_g=0.0,  vmax_g=4.),
    9:  dict(vmin_r=-0.7, vmax_r=5.0,  vmin_b=-0.07, vmax_b=0.5,  vmin_g=-0.7, vmax_g=7.),
    10: dict(vmin_r=-0.03, vmax_r=1.0, vmin_b=-0.1,  vmax_b=1.0,  vmin_g=0.0,  vmax_g=1.3),
    11: dict(vmin_r=-0.1, vmax_r=2.0,  vmin_b=-0.01, vmax_b=0.2,  vmin_g=-0.1, vmax_g=2.5),
    12: dict(vmin_r=-0.1, vmax_r=1.5,  vmin_b=-0.02, vmax_b=0.15, vmin_g=-0.15, vmax_g=2.0),
    13: dict(vmin_r=0.,   vmax_r=1.,   vmin_b=0.,    vmax_b=0.6,  vmin_g=0.,   vmax_g=1.5),
    14: dict(vmin_r=-0.05, vmax_r=0.6, vmin_b=0.,    vmax_b=0.4,  vmin_g=-0.07, vmax_g=0.8),
    15: dict(vmin_r=0.,   vmax_r=0.7,  vmin_b=0.,    vmax_b=0.8,  vmin_g=0.,   vmax_g=0.9),
}


def _make_color_postage(red, green, blue, input_dir, output_dir,
                        target, ra, dec, filt_labels,
                        zsrc1, zdef, zsrc2=None,
                        thumb_size_arcsec=20, scalebar_size_arcsec=1, clims_set=1):
    import aplpy
    clims = _CLIMS.get(clims_set, _CLIMS[1])
    rgb_cube = f'{input_dir}/rgb_cube_L3.fits'
    rgb_img  = f'{input_dir}/{target}_{"_".join(filt_labels)}_RGB_image_L3.png'

    aplpy.make_rgb_cube([red, green, blue], rgb_cube, north=True)
    aplpy.make_rgb_image(rgb_cube, rgb_img, **clims)

    f = aplpy.FITSFigure(f'{input_dir}/rgb_cube_L3_2d.fits')
    f.show_rgb(rgb_img)
    f.axis_labels.hide(); f.tick_labels.hide(); f.ticks.hide()

    img_size = (thumb_size_arcsec * u.arcsec).to(u.deg).value
    f.recenter(ra, dec, width=img_size, height=img_size)
    f.add_scalebar(scalebar_size_arcsec * u.arcsec, label=f'{scalebar_size_arcsec}"', corner='bottom right')
    f.scalebar.set_font(size=tick_size, weight='medium')
    f.scalebar.set(color='white', linewidth=2)

    title = (f'{target}     z_def = {zdef}     z_src1 = {zsrc1}     z_src2 = {zsrc2}'
             if zsrc2 else f'{target}    z_def = {zdef}     z_src = {zsrc1}')
    f.set_title(title, size=tick_size)

    colours = ['r', 'g', 'b'] if len(filt_labels) == 3 else ['r', 'b']
    y_positions = [0.95, 0.90, 0.85]
    for label, colour, ypos in zip(filt_labels, colours, y_positions):
        f.ax.text(0.02, ypos, label, color=colour,
                  bbox=dict(facecolor='w', alpha=1), size=cbar_tick_size,
                  ha='left', va='bottom', transform=f.ax.transAxes)

    out_name = f'{output_dir}/{target}_{"_".join(filt_labels)}_{thumb_size_arcsec}arcs_image_L3.png'
    f.save(out_name, dpi=1200)
    plt.close('all')


def run_step3(target_filter=None, extra_target_info=None):
    """
    target_filter     : if set, only this target is processed (respects --target flag)
    extra_target_info : dict with keys ra/dec/z_src/z_def — used when a target was
                        supplied via CLI and may not be in any proposal CSV
    """
    print("\n=== STEP 3: Postage stamps ===")
    if target_filter:
        print(f"    (single-target mode: {target_filter})")

    # Load per-target config CSV
    ps_cfg = pd.read_csv(POSTAGE_STAMP_CONFIG_CSV).set_index('objname')

    # Load parent catalogue for second-source redshifts
    parent_cat = pd.read_csv(PARENT_CATALOGUE_CSV) if PARENT_CATALOGUE_CSV.exists() else pd.DataFrame()

    # Load all proposal CSVs and build a lookup: objname → (ra, dec, z_def, z_src)
    target_info = {}
    for propid, csv_path in PROPOSAL_CSVS.items():
        if not csv_path.exists():
            continue
        df = pd.read_csv(csv_path)
        for _, row in df.iterrows():
            name = row['objname']
            if name not in target_info:
                z_def, z_src = _resolve_z(row)
                target_info[name] = {
                    'ra': row.get('RAJ2000'), 'dec': row.get('DECJ2000'),
                    'z_def': z_def, 'z_src': z_src,
                }

    # Inject CLI-supplied info for targets not covered by any proposal CSV
    if target_filter and extra_target_info and target_filter not in target_info:
        target_info[target_filter] = extra_target_info

    POSTAGE_STAMP_DIR.mkdir(parents=True, exist_ok=True)

    for target_dir in sorted(p for p in MAIN_DIR.iterdir() if p.is_dir()):
        target = target_dir.name

        # Single-target filter
        if target_filter and target != target_filter:
            continue

        # Check skip flag in config CSV
        if target in ps_cfg.index and ps_cfg.loc[target, 'skip']:
            continue
        if target.endswith('B'):
            continue

        info = target_info.get(target, {})
        ra   = info.get('ra')
        dec  = info.get('dec')
        if ra is None or dec is None:
            print(f"  {target}: no coordinates found — skipping")
            continue

        z_def = info.get('z_def')
        z_src = info.get('z_src')

        # Look up second-source z (from B entry in parent catalogue)
        z_src_2 = None
        if not parent_cat.empty:
            match = parent_cat[parent_cat['objname'] == f'{target[:-1]}B']
            if len(match):
                z_src_2 = _resolve_z(match.iloc[0])[1]

        # Per-target postage stamp settings
        cfg_row = ps_cfg.loc[target] if target in ps_cfg.index else None
        thumb_size      = int(cfg_row['thumb_size_arcsec']) if cfg_row is not None else 20
        clims_set       = int(cfg_row['clims_set'])         if cfg_row is not None else 1
        alt_run         = bool(cfg_row['alt_run'])           if cfg_row is not None else False
        scalebar_size   = 2

        # Discover which filters are present for this target
        filts_list  = []
        propids_list = []
        for subdir in sorted(d for d in (target_dir / 'HST').iterdir() if d.is_dir()):
            parts = subdir.name.split('_')
            if len(parts) < 2:
                continue
            propid, filt = parts[0], parts[1]
            if filt not in filts_list:
                filts_list.append(filt)
                propids_list.append(propid)

        if not filts_list:
            continue

        print(f"  {target}: {filts_list}")

        # Determine camera per filter
        cam_map = {'F140W': 'WFC3', 'F200LP': 'WFC3', 'F606W': 'ACS'}

        hst_dir = target_dir / 'HST'

        if len(filts_list) == 1:
            band = filts_list[0]
            cam  = cam_map.get(band, 'WFC3')
            propid = propids_list[0]
            suffix = 'drz' if band == 'F140W' else 'drc'
            imgs = (list((hst_dir / f'{propid}_{band}').glob(f'*{suffix}_img_L1.fits')) or
                    list((hst_dir / f'{propid}_{band}').glob(f'*{suffix}_sci.fits')))
            if not imgs:
                print(f"    No image file found for {target} {band}")
                continue
            _make_grayscale_postage(
                str(imgs[0]), str(POSTAGE_STAMP_DIR),
                target, ra, dec, zsrc1=z_src, zdef=z_def, band=band, zsrc2=z_src_2,
                thumb_size_arcsec=thumb_size, scalebar_size_arcsec=scalebar_size, clims_set=clims_set)

        elif len(filts_list) == 2:
            filt1, filt2 = filts_list[0], filts_list[1]
            cam1  = cam_map.get(filt1, 'WFC3')
            cam2  = cam_map.get(filt2, 'ACS' if filt2 == 'F606W' else 'WFC3')
            prop1, prop2 = propids_list[0], propids_list[1]

            if alt_run:
                # Files written by the standalone 6.5-Reproject_and_Rescale notebook
                red   = str(hst_dir / f'{prop1}_{filt1}/{target}_{filt1}_L3.fits')
                green = str(hst_dir / f'{prop1}_{filt1}/{target}_{filt1}_green_L3.fits')
                blue  = str(hst_dir / f'{prop2}_{filt2}/{target}_{filt2}_L3.fits')
            else:
                red   = str(hst_dir / f'{prop1}_{filt1}/{target}_{filt1}_{cam1}_drz_img_L1.fits')
                green = str(hst_dir / f'{prop2}_{filt2}/{target}_{filt2}_{cam2}_green_img_scaled_L3.fits')
                blue  = str(hst_dir / f'{prop2}_{filt2}/{target}_{filt2}_{cam2}_drc_img_scaled_L3.fits')

            _make_color_postage(
                red, green, blue, str(hst_dir), str(POSTAGE_STAMP_DIR),
                target, ra, dec, filt_labels=[filt1, filt2],
                zsrc1=z_src, zdef=z_def, zsrc2=z_src_2,
                thumb_size_arcsec=thumb_size, scalebar_size_arcsec=scalebar_size, clims_set=clims_set)

        elif len(filts_list) >= 3:
            # 3-colour: F140W=red, F606W=green (scaled), F200LP=blue (scaled)
            filt1, filt2, filt3 = 'F140W', 'F606W', 'F200LP'
            p1 = next((propids_list[i] for i, f in enumerate(filts_list) if f == filt1), None)
            p2 = next((propids_list[i] for i, f in enumerate(filts_list) if f == filt2), None)
            p3 = next((propids_list[i] for i, f in enumerate(filts_list) if f == filt3), None)
            red   = str(hst_dir / f'{p1}_{filt1}/{target}_{filt1}_WFC3_drz_img_L1.fits')
            green = str(hst_dir / f'{p2}_{filt2}/{target}_{filt2}_ACS_drc_img_scaled_L3.fits')
            blue  = str(hst_dir / f'{p3}_{filt3}/{target}_{filt3}_WFC3_drc_img_scaled_L3.fits')
            _make_color_postage(
                red, green, blue, str(hst_dir), str(POSTAGE_STAMP_DIR),
                target, ra, dec, filt_labels=[filt1, filt2, filt3],
                zsrc1=z_src, zdef=z_def, zsrc2=z_src_2,
                thumb_size_arcsec=thumb_size, scalebar_size_arcsec=scalebar_size, clims_set=clims_set)

# ── Entry point ────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='HST science-products pipeline')
    parser.add_argument('--steps', default='123',
                        help='Steps to run, e.g. "23" to skip cutout building (default: 123)')
    parser.add_argument('--target', default=None,
                        help='AGEL target name; uses the full CSV if omitted')
    parser.add_argument('--ra', type=float, default=None,
                        help='Target RA in decimal degrees (required if target not in CSV)')
    parser.add_argument('--dec', type=float, default=None,
                        help='Target Dec in decimal degrees (required if target not in CSV)')
    parser.add_argument('--z-src', type=float, default=None,
                        help='Source redshift (optional; stored as -1 if omitted)')
    parser.add_argument('--z-def', type=float, default=None,
                        help='Deflector redshift (optional; stored as -1 if omitted)')
    parser.add_argument('--proposal-id', default=None,
                        help='Override ACTIVE_PROPOSAL_ID from config.py for this run')
    parser.add_argument('--preview', action='store_true',
                        help='Show cutout-size preview grid before running; prompts to confirm')
    args = parser.parse_args()
    steps = set(args.steps)

    proposal_id = args.proposal_id or ACTIVE_PROPOSAL_ID

    # ── Build working dataframe and extra_target_info for step 3 ──────────────
    extra_target_info = None   # metadata for run_step3 if target is not in any proposal CSV

    if args.target:
        # Try to find the target in the active CSV first
        df_full = pd.read_csv(ACTIVE_TARGETS_CSV) if ACTIVE_TARGETS_CSV.exists() else pd.DataFrame()
        df = df_full[df_full['objname'] == args.target] if not df_full.empty else pd.DataFrame()

        if not df.empty:
            # Found in CSV — use that row; also pre-build extra_target_info as a
            # fallback for step 3 in case this proposal isn't in PROPOSAL_CSVS yet
            row = df.iloc[0]
            z_def, z_src = _resolve_z(row)
            extra_target_info = {
                'ra':    float(row.get('RAJ2000', 0)),
                'dec':   float(row.get('DECJ2000', 0)),
                'z_src': z_src,
                'z_def': z_def,
            }
            print(f"Found '{args.target}' in {ACTIVE_TARGETS_CSV}")

        else:
            # Not in CSV — coordinates must be supplied on the command line
            if args.ra is None or args.dec is None:
                print(f"Error: '{args.target}' was not found in {ACTIVE_TARGETS_CSV}.")
                print("Supply --ra and --dec to run without a CSV entry, or add the")
                print("target to the CSV manually and re-run.")
                return

            print(f"'{args.target}' not found in CSV — using coordinates from command line.")
            extra_target_info = {
                'ra':    args.ra,
                'dec':   args.dec,
                'z_src': args.z_src,
                'z_def': args.z_def,
            }

            # Build a minimal single-row dataframe for steps 1 and 2
            df = pd.DataFrame([{
                'objname':                              args.target,
                'catalogue_objname':                   args.target,
                'RAJ2000':                             args.ra,
                'DECJ2000':                            args.dec,
                'z_spec_SRC':                          args.z_src,
                'z_spec_DE':                           args.z_def,
                'z_spec_SRC_Spectral_Observations_Tally': float('nan'),
                'z_source (from DR2 Redshifts)':       float('nan'),
                'z_spec_DE_Spectral_Observations_Tally':  float('nan'),
                'z_deflector (from DR2 Redshifts)':    float('nan'),
            }])

            # Offer to append this target to the CSV for future runs
            if ACTIVE_TARGETS_CSV.exists():
                print(f"\nWould you like to append '{args.target}' to {ACTIVE_TARGETS_CSV}")
                print("so it is included in future full-CSV runs?")
                ans = input("Append? [y/N]: ").strip().lower()
                if ans == 'y':
                    updated = pd.concat([df_full, df], ignore_index=True)
                    updated.to_csv(ACTIVE_TARGETS_CSV, index=False)
                    print(f"  Appended to {ACTIVE_TARGETS_CSV}")
                    print("  Note: update 'catalogue_objname' in the CSV if you need")
                    print("  this target's MAST name for the download step.")

    else:
        # No --target: run the full CSV
        if not ACTIVE_TARGETS_CSV.exists():
            print(f"Error: CSV not found at {ACTIVE_TARGETS_CSV}.")
            print("Set ACTIVE_TARGETS_CSV in config.py, or use --target with --ra/--dec")
            print("for a single-target run without a CSV.")
            return
        df = pd.read_csv(ACTIVE_TARGETS_CSV)

    print(f"\nProcessing {len(df)} target(s)  |  proposal: {proposal_id}")
    print(f"Running steps: {', '.join(sorted(steps))}")

    # ── Optional cutout preview ───────────────────────────────────────────────
    if args.preview:
        print(f"\nShowing cutout preview (DEFAULT_CUTOUT_SIZE_ARCSEC = {DEFAULT_CUTOUT_SIZE_ARCSEC}\")...")
        preview_cutouts(df, proposal_id, ACTIVE_FILTERS, ACTIVE_CAMERA)
        ans = input("\nProceed with selected steps? [y/N]: ").strip().lower()
        if ans != 'y':
            print("Aborted.")
            return

    # ── Run selected steps ────────────────────────────────────────────────────
    if '1' in steps:
        run_step1(df.to_dict('records'), proposal_id, ACTIVE_FILTERS, ACTIVE_CAMERA)
    if '2' in steps:
        run_step2(df, proposal_id, ACTIVE_FILTERS, ACTIVE_CAMERA)
    if '3' in steps:
        run_step3(target_filter=args.target, extra_target_info=extra_target_info)

    print("\nProducts pipeline complete.")


if __name__ == '__main__':
    main()
