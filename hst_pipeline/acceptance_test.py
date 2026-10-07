"""
acceptance_test.py — measures the numbers from the DESJ0206-0114 audit
acceptance test against a target's drizzled output.

Deliberately standalone from hst_reduction.py / pipeline_qc.py's internal
CR-rejection machinery: the "intrinsically CR-free" floor for the aperture-
photometry check is built here directly from the RAW archive exposures
(no tweakreg, no L.A.Cosmic, no custom DQ flags) so this check does not
depend on any of that code being correct -- an independent check should not
share a bug with the thing it's checking.

Usage:
    python acceptance_test.py --main-dir local_test_run --target DESJ0206-0114 \
        --proposal 16773 --camera WFC3
"""
import argparse
import glob
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy import units as u

import pipeline_qc


def find_sci(main_dir, target, propno, band, camera):
    suffix = 'drz' if band == 'F140W' else 'drc'
    path = Path(main_dir) / target / 'HST' / f'{propno}_{band}' / f'{target}_{band}_{camera}_{suffix}_sci.fits'
    if not path.exists():
        raise FileNotFoundError(path)
    return path


def check_1_astrometric_residual(f140_sci, f200_sci):
    print("\n--- (1) F140W <-> F200LP residual: matched-source centroids ---")
    dra, ddec = pipeline_qc.measure_filter_offset(f200_sci, f140_sci)
    mag_mas = np.hypot(dra, ddec) * 1000
    print(f"  centroid method:        dRA={dra*1000:+.2f} mas  dDec={ddec*1000:+.2f} mas  "
          f"|offset|={mag_mas:.2f} mas")

    print("\n--- (2) Independent cross-check: whole-frame cross-correlation ---")
    try:
        dra_x, ddec_x = pipeline_qc.measure_filter_offset_xcorr(f200_sci, f140_sci)
        mag_x_mas = np.hypot(dra_x, ddec_x) * 1000
        print(f"  cross-correlation method: dRA={dra_x*1000:+.2f} mas  dDec={ddec_x*1000:+.2f} mas  "
              f"|offset|={mag_x_mas:.2f} mas")
        agreement_mas = np.hypot(dra - dra_x, ddec - ddec_x) * 1000
        print(f"  agreement between methods: {agreement_mas:.2f} mas "
              f"(acceptance target: < 5 mas)")
        if agreement_mas >= 5:
            print("  *** methods do NOT agree to < 5 mas on this field -- see note in "
                  "pipeline_qc.measure_filter_offset_xcorr's docstring: whole-frame "
                  "cross-correlation is known to be unreliable on this sparse, "
                  "CR-heavy field and should not be trusted over the source-matched "
                  "result. Reporting both rather than hiding the disagreement.")
    except Exception as exc:
        print(f"  cross-correlation method FAILED: {exc}")
        agreement_mas = None

    return mag_mas, agreement_mas


def check_2_cr_detections(f200_sci, label="after-fix"):
    print(f"\n--- (3) F200LP CR-like detections at 5sigma, sharpness>0.8 ({label}) ---")
    from photutils.detection import DAOStarFinder

    with fits.open(f200_sci) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        data = np.asarray(hdu.data, dtype=float)

    finite = np.isfinite(data)
    mean, median, std = sigma_clipped_stats(data[finite], sigma=3.0, maxiters=5)
    finder = DAOStarFinder(fwhm=1.9, threshold=5 * std, sharplo=0.8, sharphi=2.0)
    sources = finder(np.where(finite, data - median, 0.0))
    n = 0 if sources is None else len(sources)
    print(f"  {n} detections (acceptance target: few hundred, not >10,000)")
    return n


def check_3_pixel_sigma(f200_sci, f200_wht, label="after-fix"):
    print(f"\n--- (4) F200LP pixel sigma in full-depth region ({label}) ---")
    with fits.open(f200_sci) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        data = np.asarray(hdu.data, dtype=float)
    with fits.open(f200_wht) as hdul:
        hdu_w = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        wht = np.asarray(hdu_w.data, dtype=float)

    full_depth = wht >= 0.9 * np.nanmax(wht)
    region = data[full_depth & np.isfinite(data)]
    mean, median, std = sigma_clipped_stats(region, sigma=3.0, maxiters=5)
    print(f"  {full_depth.sum()} full-depth pixels; sigma={std:.4f} "
          f"(acceptance target: ~0.03, not ~0.23)")
    return std


def build_raw_single_drizzles(raw_flc_files, refimage, out_prefix):
    """Drizzle each raw archive exposure alone (no tweak, no CR rejection,
    no custom DQ) onto `refimage`'s grid -- the independent, unprocessed
    floor for the aperture-photometry check below."""
    from drizzlepac import astrodrizzle

    paths = []
    for i, flc in enumerate(raw_flc_files):
        prefix = f"{out_prefix}_raw{i}"
        astrodrizzle.AstroDrizzle(
            [str(flc)], output=prefix, resetbits=0, driz_cr_corr=False,
            build=False, final_wht_type='EXP', final_wcs=True, final_refimage=str(refimage))
        match = glob.glob(prefix + '_dr?_sci.fits')
        if not match:
            raise FileNotFoundError(f"{prefix}_dr?_sci.fits not produced")
        paths.append(match[0])
    return paths


def check_4_flux_floor(f200_sci_corrected, raw_single_paths, protect_sky):
    print("\n--- (5) Aperture photometry vs pairwise-minimum (CR-free) floor ---")
    from photutils.aperture import CircularAperture, aperture_photometry

    with fits.open(f200_sci_corrected) as hdul:
        hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        wcs = WCS(hdu.header)
        data = np.asarray(hdu.data, dtype=float)

    with fits.open(raw_single_paths[0]) as h0, fits.open(raw_single_paths[1]) as h1:
        d0 = next(h for h in h0 if h.data is not None and h.data.ndim == 2)
        d1 = next(h for h in h1 if h.data is not None and h.data.ndim == 2)
        floor = np.fmin(np.asarray(d0.data, dtype=float), np.asarray(d1.data, dtype=float))

    if len(protect_sky) == 0:
        print("  no confirmed sources available to check.")
        return []

    xy = np.column_stack(wcs.world_to_pixel(protect_sky))
    results = []
    for i, (x, y) in enumerate(xy[:20]):  # cap for runtime
        if not (0 <= x < data.shape[1] and 0 <= y < data.shape[0]):
            continue
        ap = CircularAperture((x, y), r=3.0)
        flux_final = float(aperture_photometry(data, ap)['aperture_sum'][0])
        flux_floor = float(aperture_photometry(floor, ap)['aperture_sum'][0])
        ok = flux_final >= flux_floor
        results.append((i, flux_final, flux_floor, ok))
    n_below = sum(1 for r in results if not r[3])
    print(f"  {len(results)} sources checked; {n_below} fell BELOW the CR-free floor "
          f"(acceptance target: 0 -- floor is a lower bound, corrected mosaic should "
          f"never fall below it)")
    for i, ff, fl, ok in results:
        flag = "OK" if ok else "*** BELOW FLOOR ***"
        print(f"    source {i}: final={ff:9.2f}  floor={fl:9.2f}  ratio={ff/fl if fl else float('nan'):6.3f}  {flag}")
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--main-dir', required=True)
    ap.add_argument('--target', required=True)
    ap.add_argument('--proposal', default='16773')
    ap.add_argument('--camera', default='WFC3')
    args = ap.parse_args()

    main_dir = Path(args.main_dir)
    f140_sci = find_sci(main_dir, args.target, args.proposal, 'F140W', args.camera)
    f200_sci = find_sci(main_dir, args.target, args.proposal, 'F200LP', args.camera)
    f200_wht = Path(str(f200_sci).replace('_sci.fits', '_wht.fits'))

    print(f"F140W: {f140_sci}")
    print(f"F200LP: {f200_sci}")

    mag_mas, agreement_mas = check_1_astrometric_residual(f140_sci, f200_sci)
    n_cr = check_2_cr_detections(f200_sci)
    sigma = check_3_pixel_sigma(f200_sci, f200_wht)

    raw_dir = main_dir / args.target / 'HST' / f'{args.proposal}_F200LP' / 'raw_data'
    raw_flc = sorted(raw_dir.glob('*flc.fits'))
    if len(raw_flc) == 2:
        protect_sky = pipeline_qc.repeat_detected_protection_mask(
            *[str(p) for p in _prelim_for_floor(raw_flc, f200_sci)],
            match_radius_arcsec=0.15, fwhm=3.0, nsigma=8.0)
        raw_singles = build_raw_single_drizzles(
            raw_flc, f200_sci, str(main_dir / args.target / 'HST' / f'{args.proposal}_F200LP' / 'floor_check'))
        check_4_flux_floor(f200_sci, raw_singles, protect_sky)
    else:
        print("\n--- (5) skipped: raw F200LP exposures not found ---")

    print("\n=== SUMMARY ===")
    print(f"F140W-F200LP residual (after correction): {mag_mas:.2f} mas "
          f"(acceptance target: consistent with zero at the few-mas level)")
    if agreement_mas is not None:
        print(f"Centroid vs cross-correlation agreement: {agreement_mas:.2f} mas "
              f"(acceptance target: < 5 mas)")
    print(f"F200LP CR-like detections (after fix): {n_cr} (acceptance target: few hundred)")
    print(f"F200LP full-depth pixel sigma (after fix): {sigma:.4f} (acceptance target: ~0.03)")


def _prelim_for_floor(raw_flc, refimage):
    """Quick single-exposure drizzles (native grid, no ref-matching needed
    since repeat_detected_protection_mask matches on sky, not pixel grid)
    reused just to build a source list for the flux-floor check."""
    from drizzlepac import astrodrizzle
    paths = []
    for i, flc in enumerate(raw_flc):
        prefix = str(Path(flc).with_name(f'{Path(flc).stem}_floorprelim'))
        astrodrizzle.AstroDrizzle(
            [str(flc)], output=prefix, resetbits=0, driz_cr_corr=False,
            build=False, final_wht_type='EXP', final_wcs=True)
        match = glob.glob(prefix + '_dr?_sci.fits')
        paths.append(match[0])
    return paths


if __name__ == '__main__':
    main()
