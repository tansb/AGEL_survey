"""
config.py — shared configuration for hst_reduction.py and hst_products.py

Edit the values in this file before running either pipeline script.
"""

from pathlib import Path

# ── Directories ───────────────────────────────────────────────────────────────
# !! Update every path below before running anything !!

# Root directory where all HST data lives (per-target subdirectories inside)
MAIN_DIR = Path('/path/to/your/HST_data')

# Directory where cutout HDF5/FITS products are written (usually same as MAIN_DIR)
DATA_DIR = Path('/path/to/your/HST_data')

# Directory containing PSF models (psf_model_<band>.h5) and other supplementary files
LENS_PROC_DIR = Path('/path/to/your/supplementary')

# Output directory for postage-stamp PNGs
POSTAGE_STAMP_DIR = Path('/path/to/your/postage_stamps')

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
# Used by the reprojection, green-image, cutout, and postage-stamp steps, which
# scan all proposal subdirectories rather than a single active proposal.
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
PARENT_CATALOGUE_CSV = Path('/path/to/your/parent_catalogue.csv')

# ── Cutout settings (used by hst_products.py step 1) ──────────────────────────
DEFAULT_CUTOUT_SIZE_ARCSEC = 25.0   # arcsec; per-target overrides go in targets CSV

# ── Offset-fitter settings (used by hst_products.py step 2) ───────────────────
OFFSET_NUM_DUP = 4   # number of PSO repetitions for the 0.6-arcsec fitting stage
