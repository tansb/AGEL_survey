"""
config.py — shared configuration for hst_reduction.py and hst_products.py

Edit the values in this file before running either pipeline script.
"""

from pathlib import Path

# ── Directories ───────────────────────────────────────────────────────────────
# !! Update every path below before running anything !!

# Root directory where all HST data lives (per-target subdirectories inside)
MAIN_DIR = Path('/Volumes/AGEL/agel_hst/HST_workingdata')

# Directory where cutout HDF5/FITS products are written (usually same as MAIN_DIR)
DATA_DIR = Path('/Volumes/AGEL/agel_hst/HST_workingdata')

# Directory containing PSF models (psf_model_<band>.h5) and other supplementary files
LENS_PROC_DIR = Path('/Volumes/AGEL/agel_hst/lens_processing')

# Output directory for postage-stamp PNGs
POSTAGE_STAMP_DIR = Path('/Volumes/AGEL/agel_hst/HST_postage_stamps')

# Local raw-data cache batch_runner.py checks before falling back to a MAST
# download (batch_runner.stage_raw_data) — same directory layout as MAIN_DIR,
# i.e. {objname}/HST/{proposal_id}_{filter}/raw_data/*.fl[tc].fits.
HST_DATA_CACHE = Path('/Volumes/AGEL/agel_hst/HST_data')

# ── Active reduction run ───────────────────────────────────────────────────────
# Set these to match the proposal you are currently downloading / drizzling.
# These are used by hst_reduction.py (steps 1 and 2).

ACTIVE_PROPOSAL_ID = '17307'
ACTIVE_FILTERS     = ['F606W']   # list so multiple filters can be processed in one run
ACTIVE_CAMERA      = 'ACS'       # 'ACS' or 'WFC3'

# CSV used for downloading and drizzling.
# Must contain columns: objname (AGEL name), catalogue_objname (MAST name),
# RAJ2000, DECJ2000, and redshift columns (see notebooks for exact column names).
ACTIVE_TARGETS_CSV = MAIN_DIR / f'{ACTIVE_PROPOSAL_ID}_targets.csv'

# ── All known proposals ────────────────────────────────────────────────────────
# NOT actually read by any script (checked 2026-10 -- grep the codebase for
# ALL_PROPOSALS and this is the only hit). The reprojection (run_step3) and
# green-image (run_step4) steps in hst_reduction.py, and the postage-stamp
# colour-channel logic in hst_products.py, each have the proposal/filter/
# camera combinations below HARDCODED independently instead of reading this
# list -- see the "Known limitations" section in README.md. Keep this updated
# as documentation of what those three places assume, but editing it alone
# does NOT add a new proposal to those steps; you have to edit each of them.
# Format: (proposal_id, filter_name, camera)
ALL_PROPOSALS = [
    ('15867', 'F140W',  'WFC3'),
    ('16773', 'F140W',  'WFC3'),
    ('16773', 'F200LP', 'WFC3'),
    ('17307', 'F606W',  'ACS'),
]

# Per-proposal target CSVs (used by the postage-stamp step to look up RA/Dec/z).
PROPOSAL_CSVS = {
    '15867': MAIN_DIR / '15867_targets.csv',
    '16773': MAIN_DIR / '16773_targets.csv',
    '17307': MAIN_DIR / '17307_targets.csv',
}

# ── Postage-stamp config CSV ───────────────────────────────────────────────────
# Path to the CSV that controls per-target thumbnail and colour-scale settings.
# Generated alongside this file; edit it to customise individual targets.
POSTAGE_STAMP_CONFIG_CSV = Path(__file__).parent / 'postage_stamp_config.csv'

# ── Parent catalogue ───────────────────────────────────────────────────────────
# Full parent catalogue; used to look up second-source redshifts (z_src_2).
# Set to None to skip second-source redshift lookup in postage stamps.
PARENT_CATALOGUE_CSV = Path('/Volumes/AGEL/agel_hst/Parent-Catalogue-All_targets.csv')

# ── Cutout settings (used by hst_products.py step 1) ──────────────────────────
DEFAULT_CUTOUT_SIZE_ARCSEC = 25.0   # arcsec; per-target overrides go in targets CSV

# ── Offset-fitter settings (used by hst_products.py step 2) ───────────────────
OFFSET_NUM_DUP = 4   # number of PSO repetitions for the 0.6-arcsec fitting stage

# ── Drizzle defaults (used by hst_reduction.py step 2) ─────────────────────────
# These are the starting values shown when running with --interactive; they are
# used as-is when --interactive is not passed. Edit here to change the defaults
# for every run, or override per-run at the prompt.

IR_DRIZZLE_DEFAULTS = {
    'pixel_size':              0.08,   # final_scale, arcsec/pix
    'final_pixfrac':           1.0,
    'final_wht_type':          'EXP',
    'tweak_threshold':         5,      # multi-exposure imagefindcfg threshold
    # Audit defect 3: drizzlepac's imagefindpars conv_width default (3.5) is
    # tuned for well-sampled UVIS/ACS PSFs. WFC3/IR is undersampled at
    # ~1.2 px FWHM, so 3.5 over-smooths and biases/loses source detections
    # feeding TweakReg; 2.5 matches the IR PSF width.
    'tweak_conv_width':        2.5,
    'tweak_searchrad':         2.0,
    'tweak_ylimit':            0.3,
    # drizzlepac's default separation=0.0 applies no minimum-separation
    # filtering to detected sources before matching. Measured against a
    # real batch failure (AGEL144431+241843A, 16773/F140W, 3 exposures,
    # only 18-27 sources detected per sparse IR frame): with separation=0,
    # stsci.stimage's compiled xyxymatch occasionally produces more
    # candidate match pairs than its output buffer (sized to the smaller
    # image's source count) can hold, crashing with "Number of output
    # coordinates exceeded allocation" deep inside TweakReg's own fit step
    # (drizzlepac/imgclasses.py Image.match -> xyxymatch) -- not caught by
    # our own Defect-4 NaN check, since TweakReg never gets far enough to
    # write a shiftfile. Reproduced intermittently even re-running the
    # identical images (small nondeterminism in driz_cr's detected source
    # positions pushes this marginal, very-sparse-field case over the
    # edge some runs and not others). separation=0.5 (px) filters out
    # pathologically close duplicate detections before matching and
    # resolved every case checked, with no effect on real, well-separated
    # sources at HST resolution.
    'tweak_separation':        0.5,
    # Same 2dhist coarse-shift instability as UV_DRIZZLE_DEFAULTS below (see
    # its comment for the full explanation) -- confirmed on the IR side too
    # once Defect 6's sigma fix actually let TweakReg see real sources
    # (AGEL103255+751854A / AGEL142822+031800A, 15867/F140W: 67-156 sources
    # detected per frame, still "Not enough matches (< 15)" with
    # use2dhist's default coarse search). Same fix, same justification: IR
    # dithers are also always sub-arcsec here, so use2dhist=False with an
    # explicit (0, 0) guess and a tolerance-bounded direct match is safe
    # and resolved both cases cleanly.
    'tweak_use2dhist':         False,
    'tweak_xoffset':           0.0,
    'tweak_yoffset':           0.0,
    'tolerance':               10.0,
    # drizzlepac's own default (minobj=15) for the minimum number of matched
    # pairs a fit needs. Named here so TWEAK_PARAM_OVERRIDES can lower it per
    # target without touching the TweakReg call site.
    'tweak_minobj':            15,
}

# ── Per-target TweakReg overrides (used by hst_reduction.ir_drizzler) ──────────
#
# Keyed by (proposal_id, objname); the dict is merged over IR_DRIZZLE_DEFAULTS
# for that target only, so the 50 already-reduced 15867 targets keep the exact
# parameters they were reduced with.
#
# Why these two targets need their own values (measured 2026-10-01, by running
# TweakReg over their existing temp/*crclean.fits with a parameter sweep):
#
#   The shared 'tolerance': 10.0 above is the real culprit. It is the radius
#   (output pixels) inside which xyxymatch will accept a source pair, and at
#   10 px it accepts essentially anything -- in a sparse IR field that means
#   matching noise peaks. Both targets are sparse at the default threshold=5
#   (AGEL142822+031800A: only 13-38 sources per 1014x1014 frame;
#   AGEL103255+751854A: 61-95), so the match either fell under minobj=15 and
#   wrote a NaN shiftfile row (the observed failure), or "succeeded" on junk:
#   the real 2026-09-17 run left AGEL103255+751854A/ie5060bcq with
#   rot=359.96 deg, scale=1.0015 and xrms/yrms = 4.37/4.72 px -- a fit that
#   should never have been accepted for two exposures of the same visit.
#
#   threshold=3 brings the detection count to ~1200-1500 per frame, and
#   tolerance=1.0 then keeps only genuinely coincident pairs. Measured result:
#
#     target                 variant                       worst rms   rot      dscale
#     AGEL103255+751854A     thr5/tol10 (current default)  NaN (fail)
#     AGEL103255+751854A     thr3/tol1.0/minobj=8          0.543 px    0.038deg 4.6e-4
#     AGEL142822+031800A     thr5/tol10 (current default)  NaN (fail)
#     AGEL142822+031800A     thr3/tol1.0/minobj=8          0.520 px    0.034deg 5.0e-5
#
#   minobj=8 is needed only by AGEL103255+751854A (tol=1.0 leaves it 14 good
#   pairs, just under drizzlepac's 15); AGEL142822+031800A clears 15 on its
#   own but is given the same value for consistency. 8 pairs still
#   over-constrains a 4-parameter shift+rot+scale fit.
#
# NOTE: tolerance=10.0 is loose for every IR target, not just these two --
# tolerance=1.0 scored equal or better on every target swept, including the
# source-rich ones (AGEL142104+002219A 0.142 -> 0.133 px, AGEL144431+241843A
# 0.212 -> 0.186 px). It is left alone as the global default here only so the
# already-reduced targets stay reproducible; consider lowering it for the next
# full re-reduction.
TWEAK_PARAM_OVERRIDES = {
    ('15867', 'AGEL103255+751854A'): {
        'tweak_threshold': 3,
        'tolerance':       1.0,
        'tweak_minobj':    8,
    },
    ('15867', 'AGEL142822+031800A'): {
        'tweak_threshold': 3,
        'tolerance':       1.0,
        'tweak_minobj':    8,
    },
}

UV_DRIZZLE_DEFAULTS = {
    'pixel_size':          0.05,   # final_scale, arcsec/pix
    'final_pixfrac_initial': 0.8,  # pixfrac for the CR-rejection drizzle pass
    'final_pixfrac_final':   1.0,  # pixfrac for the final drizzle pass
    'final_wht_type':      'EXP',
    'driz_cr_scale':       '0.9 0.6',
    'driz_cr_grow':        3,
    'tweak_threshold':     20,
    # UVIS/ACS are well-sampled, so the drizzlepac default of 3.5 is
    # appropriate here (unlike IR above) -- left unchanged.
    'tweak_conv_width':    3.5,
    'tweak_refconv_width': 2.5,
    'tweak_searchrad':     5.0,
    'tweak_ylimit':        0.6,
    # TweakReg's use2dhist=True (drizzlepac default) coarse-aligns the two
    # exposures with a 2D-histogram peak search over +/-searchrad (in
    # ARCSECONDS, not pixels -- easy to misread). Measured against a batch
    # run of 18 targets: this survey's exposure pairs are only ever
    # dithered by a fraction of an arcsec (POSTARG deltas of ~0.1-0.25"
    # across every target checked, for CR rejection, not mosaicking), so
    # searchrad should matter very little -- but the histogram's peak
    # location turned out to depend on searchrad in a non-monotonic,
    # essentially arbitrary way for crowded/extended-source fields (one
    # target's fit succeeded at 0.5" and 3.0" but failed at 1.0" and 2.0").
    # No single searchrad value is safe for every target. A careful
    # sky-position cross-match (pipeline_qc's repeat-detection, used for
    # Defect 5b) confirmed 80-400+ genuinely well-matched sources in every
    # field checked, including the ones that failed -- the true offset is
    # small and real, only the histogram-based coarse search is unreliable.
    # use2dhist=False skips that step entirely and matches directly around
    # an explicit (xoffset, yoffset) guess using plain tolerance-bounded
    # nearest-neighbor matching (xyxymatch) -- no binning, so no bin-edge
    # sensitivity. Since the true offset is always small, (0, 0) is a safe
    # guess and this resolved every case checked (2dhist- and
    # searchrad-tuning failures alike) with a clean, consistent fit.
    'tweak_use2dhist':     False,
    'tweak_xoffset':       0.0,
    'tweak_yoffset':       0.0,
    # Bounds the direct xyxymatch search around (xoffset, yoffset) now that
    # there's no coarse pre-alignment step to narrow it first. Every true
    # offset measured (via a fit or via POSTARG) has been under ~3 px;
    # 10 px gives generous headroom for that while staying tight enough
    # in a crowded field that it shouldn't cross-match unrelated sources.
    'tolerance':           10.0,
    # See IR_DRIZZLE_DEFAULTS['tweak_separation'] -- same xyxymatch output-
    # buffer overflow is possible here in principle (same underlying
    # drizzlepac/stsci.stimage code path). Not yet observed on the UV side,
    # applied defensively for the same near-zero cost.
    'tweak_separation':    0.5,
}

# ── N=2 cosmic-ray rejection (used by hst_reduction.py's uv_drizzler when
# exactly 2 UV exposures are present; see audit defect 5) ─────────────────
# driz_cr's own median-based rejection degenerates to a simple average at
# N=2 and cannot reliably tell a cosmic ray from an undersampled point
# source. These parameters drive pipeline_qc.py's L.A.Cosmic + source-
# protection + pairwise-rejection pass that supplements (not replaces)
# driz_cr for that specific case.
LACOSMIC_DEFAULTS = {
    'sigclip':  4.0,    # safe to run aggressively once source protection is in place
    'objlim':   3.0,
    'psffwhm':  1.9,    # F200LP PSF width, px -- NOT astroscrappy's 2.5 default
}

SOURCE_PROTECTION_DEFAULTS = {
    'match_radius_arcsec': 0.15,   # repeat-detection cross-match tolerance
    'protect_radius_px':   4,      # disc radius un-flagged around confirmed sources
    'detect_fwhm_px':      3.0,    # DAOStarFinder FWHM on the *drizzled* (0.05"/pix) grid
    # 5.0 (with no shape filtering) measured 636 "confirmed" sources on the
    # real DESJ0206-0114 F200LP pair -- mostly chance coincidences between
    # cosmic-ray-driven spurious detections in each exposure, not real
    # matches (see pipeline_qc.detect_source_positions). 8.0 sigma plus the
    # sharpness/roundness cuts added there bring this down to ~150; still
    # over-protects a little in exchange for never risking a real source,
    # which is the safer direction to err in.
    'detect_nsigma':       8.0,
}

PAIRWISE_CR_DEFAULTS = {
    'threshold_k':  4.0,   # k=2.5 flagged ~25% of all pixels -- far too aggressive
    'dilation_px':  1,
}
