# Preface Notes from CW: 
Before you begin running code I *highly* suggest

(1) organizing your data such that you have your codes in your main working directory and then all the observing data organized as /AGEL_name/HST/proposalID_filter/ (e.g. if you have F140W and F200LP data for the target AGEL110725+245943A from the Glazebrook 16773 proposal, you would probably want to organize your data as /AGEL110725+245943A/HST/16773_F140W/ and /AGEL110725+245943A/HST/16773_F200LP/)

(2) creating a CSV file (I decided to go with one file per proposal) that contains the target name recorded in MAST, the matching AGEL target ID, the target RA and DEC coordinates, and any additional redshift info that the source may have. Helpful to start with the Parent Catalog in Airtable, and filter/extract the info you need from there. What I ended up doing was download the CSV of observations for a given proposal from the MAST search, cross-matched that to the AGEL parents catalog, and then just created a new CSV file with the relevant columns. Then I just add sources to this CSV as they come in for any on-going programs (e.g. 17307).

# HST Pipeline — Final Scripts

Two Python scripts are used for batch processing of HST imaging. A shared config file and a per-target CSV control all settings without touching the scripts themselves.

---

## Files

| File | Purpose |
|---|---|
| `config.py` | All paths and proposal settings — **edit this first** |
| `postage_stamp_config.csv` | Per-target thumbnail size, colour scale, skip/alt flags |
| `hst_reduction.py` | download → drizzle → reproject → green image |
| `hst_products.py` | cutouts → offset fitter → postage stamps |
| `6.5-Reproject_and_Rescale.ipynb` | Standalone manual alignment for problem targets (see below) |

---

## Installation

### Prerequisites

- [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)
- A MAST account (free) if you intend to use the download step — set up credentials with `astroquery` before first run

### 1. Clone the repository

```bash
git clone https://github.com/<your-username>/my-astro-tools.git
cd my-astro-tools/hst-pipeline
```

### 2. Create the conda environment

```bash
conda env create -f environment.yml
conda activate hst-pipeline
```

> **Note on drizzlepac:** `drizzlepac` is an STScI package with compiled C extensions. Installing via conda (`conda-forge`) is strongly recommended over pip, which may fail on some platforms. If the environment solve is slow, try [mamba](https://mamba.readthedocs.io) as a drop-in replacement: `mamba env create -f environment.yml`.

### 3. Configure paths

Open `config.py` and update every path under the `Directories` block to point to your local data storage (see [Setup — config.py](#setup--configpy) below).

### 4. Prepare your targets CSV

Copy `targets_template.csv`, rename it, to e.g. `{proposal_id}_targets.csv`, and place it inside `MAIN_DIR`, and fill in your target list. See [Required input files](#required-input-files) for the full column specification.

---

## Quick start

```bash
# 1. Edit config.py for your proposal
# 2. Run the reduction pipeline
python hst_reduction.py

# 3. (Optional) Run Reproject_and_Rescale notebook for any alt_run targets
# 4. Run the products pipeline
python hst_products.py
```

Both scripts accept `--steps` and `--target` flags for partial runs and single-target testing (see Usage section).

---

## Setup — config.py

Open `config.py` and set the following before running anything.

### Directories

```python
MAIN_DIR          # root data directory; per-target folders live inside
DATA_DIR          # where cutout HDF5/FITS files are written (usually same as MAIN_DIR)
LENS_PROC_DIR     # directory containing PSF models and supplementary files
POSTAGE_STAMP_DIR # output directory for postage-stamp PNGs
```

### Active proposal (reduction and cutout steps)

```python
ACTIVE_PROPOSAL_ID = '17307'   # proposal ID string matching MAST and your CSV filename
ACTIVE_FILTERS     = ['F606W'] # list — multiple filters processed in one run
ACTIVE_CAMERA      = 'ACS'     # 'ACS' or 'WFC3'
ACTIVE_TARGETS_CSV             # auto-set to MAIN_DIR/{ACTIVE_PROPOSAL_ID}_targets.csv
```

### All proposals (reprojection, green image, postage stamps)

```python
ALL_PROPOSALS   # list of (proposal_id, filter, camera) tuples — add entries as new data arrives
PROPOSAL_CSVS   # dict mapping proposal_id → path to that proposal's targets CSV
```

### Other settings

```python
DEFAULT_CUTOUT_SIZE_ARCSEC = 25.0   # cutout size for all targets
OFFSET_NUM_DUP             = 2      # PSO repetitions in the offset-fitter step
PARENT_CATALOGUE_CSV               # full AGEL parent catalogue (for second-source redshifts)
```

---

## Required input files

### 1. Targets CSV  (`{proposal_id}_targets.csv`)

One row per target. Required columns:

| Column | Description |
|---|---|
| `objname` | AGEL target name (e.g. `AGEL110725+245943A`) |
| `catalogue_objname` | Name as it appears in the MAST archive (used for download only) |
| `RAJ2000` | Right ascension in decimal degrees |
| `DECJ2000` | Declination in decimal degrees |
| `z_spec_SRC_Spectral_Observations_Tally` | Source redshift (first priority) |
| `z_source (from DR2 Redshifts)` | Source redshift (second priority) |
| `z_spec_SRC` | Source redshift (fallback) |
| `z_spec_DE_Spectral_Observations_Tally` | Deflector redshift (first priority) |
| `z_deflector (from DR2 Redshifts)` | Deflector redshift (second priority) |
| `z_spec_DE` | Deflector redshift (fallback) |

The redshift columns are tried in order; the first non-null value is used. Missing redshifts are stored as `-1.0` in the HDF5 files.

> **Note:** `catalogue_objname` only needs to match the MAST observation header exactly. It can differ from `objname`. Check the MAST portal if a download returns no results.

> **Note:** B entries (e.g. `AGEL110725+245943B`) are automatically skipped in all steps. They are used only to look up a second-source redshift for postage stamp labels.

### 2. PSF model files

Place HDF5 PSF models at:

```
LENS_PROC_DIR/lens_processing/psf_model_{band}.h5
```

Each file must contain a dataset named `kernel_point_source` (2D array). If a PSF file is not found, the cutout is still built without one and a warning is printed.

### 3. Data directory structure

The scripts expect data organised as:

```
MAIN_DIR/
└── {AGEL_name}/
    └── HST/
        └── {proposal_id}_{filter}/
            ├── raw_data/          ← FLC/FLT files (created by download step)
            ├── *_sci.fits         ← drizzled science image
            ├── *_wht.fits         ← drizzled weight image
            ├── *_img_L1.fits      ← renamed science image (created by reproject step)
            ├── *_scaled_L3.fits   ← UV reprojected onto IR grid
            ├── *_green_img_scaled_L3.fits  ← synthetic green image
            ├── *_cutout_L3.h5     ← lenstronomy-ready cutout
            └── *_cutout_L3.fits   ← same cutout with full FITS header
```

Each step in the pipeline creates the files needed by the next step, so running them in order produces a complete directory tree automatically.

---

## hst_reduction.py

Covers data fetch, drizzle and align, reprojection, creation of median channel for RGBs

### Steps

| Step | Flag | What it does |
|---|---|---|
| 1 | `1` | Queries MAST and downloads FLC/FLT files |
| 2 | `2` | Aligns and drizzles exposures (IR at 0.08″/pix, UV at 0.05″/pix) |
| 3 | `3` | Reprojects UV frames onto the IR pixel grid; renames outputs to `*_img_L1.fits` |
| 4 | `4` | Creates synthetic green channel: `green = sqrt(IR² + UV²)` |

### Usage

```bash
# Run all four steps
python hst_reduction.py

# Run only drizzle + reproject (e.g. data already downloaded)
python hst_reduction.py --steps 23

# Test on one target before running the full catalogue
python hst_reduction.py --steps 234 --target AGEL110725+245943A
```

### Notes

- Steps 3 and 4 scan **all** subdirectories of `MAIN_DIR`, not just the active proposal, so they update the full dataset each time they run.
- Step 2 picks `pixel_size = 0.08` for F140W and `0.05` for all other filters automatically.
- AstroDrizzle and TweakReg parameters are set to sensible defaults but may need per-target adjustment for challenging data. Check the outputs of at least a handful of targets before running the full catalogue.

---

## Reproject_and_Rescale.ipynb (standalone)

This notebook is kept separate because it requires **manual parameter setting per target** and is only needed for a subset of objects where the WCS-based reprojection in step 3 is insufficiently accurate.

**When to use it:** targets listed with `alt_run = True` in `postage_stamp_config.csv`.

**What it does:** uses OpenCV template matching and `astroalign` to register the UV frame to the IR frame, then writes three output files per target:

```
HST/{proposal_id}_F140W/{target}_F140W_L3.fits
HST/{proposal_id}_F140W/{target}_F140W_green_L3.fits
HST/{proposal_id}_{uv_filter}/{target}_{uv_filter}_L3.fits
```

These names are exactly what `hst_products.py` step 3 looks for when `alt_run = True`, so no further changes are needed — just run the notebook for each affected target, then proceed with the products pipeline.

---

## hst_products.py

Covers cutout creation, offset alignment, and postage stamp creation.

### Steps

| Step | Flag | What it does |
|---|---|---|
| 1 | `1` | Builds lenstronomy-ready HDF5 + FITS cutouts centred on target RA/Dec |
| 2 | `2` | Fits Sersic centroids in each band; writes `ra_shift`/`dec_shift` into cutout files |
| 3 | `3` | Makes grayscale, 2-colour, or 3-colour PNG postage stamps |

### Usage

```bash
# Run all three steps
python hst_products.py

# Make postage stamps only (cutouts and offsets already done)
python hst_products.py --steps 3

# Test cutout building on one target
python hst_products.py --steps 1 --target AGEL110725+245943A
```

### Postage stamp config CSV

`postage_stamp_config.csv` controls per-target appearance. Edit this file rather than the script.

| Column | Type | Description |
|---|---|---|
| `objname` | string | AGEL target name |
| `thumb_size_arcsec` | int | Thumbnail diameter in arcseconds (default 20) |
| `clims_set` | int | Colour-scale preset 1–15 (default 1; see script for limits per preset) |
| `skip` | bool | Set `True` to exclude a target from all postage-stamp output |
| `alt_run` | bool | Set `True` for targets processed with the 6.5 notebook (uses different input filenames) |

To add a new target not in the CSV, the script falls back to `thumb_size = 20` and `clims_set = 1`. Add a row to the CSV only when you need non-default values.

### How the postage stamp type is chosen

- **1 filter found** → grayscale (two panels at different stretch)
- **2 filters found** → 2-colour RGB (IR = red, UV = blue, green = synthetic blend)
- **3+ filters found** → 3-colour RGB (F140W = red, F606W = green, F200LP = blue)

---

## Typical full-pipeline run (new proposal)

```bash
# 1. Add proposal to config.py: set ACTIVE_PROPOSAL_ID, ACTIVE_FILTERS, ACTIVE_CAMERA
#    Add entry to ALL_PROPOSALS and PROPOSAL_CSVS
#    Place {proposal_id}_targets.csv in MAIN_DIR

# 2. Download and reduce
python hst_reduction.py

# 3. Inspect drizzled outputs — check alignment and CR rejection for a sample of targets

# 4. For any alt_run targets, open and run 6.5-Reproject_and_Rescale.ipynb
#    then set alt_run = True in postage_stamp_config.csv for those targets

# 5. Build cutouts, fit offsets, make stamps
python hst_products.py
```
