## Preface Notes from CW: 
Before you begin running code I *highly* suggest
(1)FOLLOWING AGEL DMP: organize your data such that you have your codes in your main working directory and then all the observing data organized as /AGEL_name/HST/proposalID_filter/ (e.g. if you have F140W and F200LP data for the target AGEL110725+245943A from the Glazebrook 16773 proposal, you would probably want to organize your data as /AGEL110725+245943A/HST/16773_F140W/ and /AGEL110725+245943A/HST/16773_F200LP/)
(2) creating a CSV file (I decided to go with one file per proposal) that contains the target name recorded in MAST, the matching AGEL target ID, the target RA and DEC coordinates, and any additional redshift info that the source may have. Helpful to start with the Parent Catalog in Airtable, and filter/extract the info you need from there. What I ended up doing was download the CSV of observations for a given proposal from the MAST search, cross-matched that to the AGEL parents catalog, and then just created a new CSV file with the relevant columns. Then I just add sources to this CSV as they come in for any on-going programs (e.g. 17307).

# HST Pipeline — Final Scripts

Two Python scripts are used for batch processing of HST imaging. A shared config file and a per-target CSV control all settings without touching the scripts themselves.

---

## Files

| File | Purpose |
|---|---|
| `config.py` | All paths and proposal settings — **edit this first** |
| `batch_runner.py` | Recommended entry point for a whole catalogue: wraps `hst_reduction.py` + `hst_products.py` per target with resumable progress tracking — see [batch_runner.py](#batch_runnerpy) |
| `batch_processing_tracker.csv` | One row per `(proposal, objname)`; `status` is `pending` / `in_progress` / `complete` / `failed` (or a custom value like `bad_data` once you've triaged a failure — `batch_runner.py` only ever treats `pending` specially). Generated/updated by `batch_runner.py`; safe and expected to be committed to version control as a running log |
| `postage_stamp_config.csv` | Per-target thumbnail size, colour scale, skip/alt flags |
| `targets_template.csv` | Empty starting point for a new `{proposal_id}_targets.csv` |
| `hst_reduction.py` | download → drizzle → reproject → green image |
| `hst_products.py` | cutouts → offset fitter → postage stamps |
| `pipeline_qc.py` | astrometric preflight/correction, DQ/shiftfile QC checks, N=2 cosmic-ray rejection — see [Astrometry & CR-rejection fixes](#astrometry--cr-rejection-fixes) |
| `acceptance_test.py` | standalone validation script (astrometric residual, CR-detection count, pixel noise, flux-floor check) for one already-drizzled target |
| `preflight_report.py` | standalone report of `pipeline_qc.select_astrometric_reference()`'s verdict for every target in one proposal, without running any reduction — useful to scan a catalogue for `AMBIGUOUS` targets before a full batch run |
| `rerun_tolerance_fix.py` | Example of the "re-run a specific list of previously-failed targets" pattern (reuses `batch_runner.process_target`) — copy and adapt it rather than running it as-is; its target list is specific to one historical defect fix |
| `Reproject_and_Rescale.ipynb` | Standalone manual alignment for problem targets (see below) |

---

## Getting Started

A strict, numbered walkthrough from a fresh clone to a completed first target. Each step says what breaks if you skip it. For a deeper dive into any one piece, follow the link to its full section further down; for a list of things that will trip you up even after setup is correct, see [Known limitations & gotchas](#known-limitations--gotchas) at the end of this walkthrough.

**What this repo is, and isn't:** cloning it gets you code, `config.py`, and the tracking CSVs (`batch_processing_tracker.csv`, `postage_stamp_config.csv`) — `.gitignore` deliberately excludes every FITS/H5 file, so **no actual HST data comes with the repo.** You either point `config.py` at an existing data drive (e.g. a shared drive another team member already populated) or let the pipeline download fresh from MAST into a new one. Decide which before step 4.

### 1. Prerequisites

- [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html).
- **macOS or Linux.** `drizzlepac` (the core dependency — see step 2) is not published for native Windows on conda-forge; on Windows, use [WSL](https://learn.microsoft.com/en-us/windows/wsl/install) instead.
- **Disk space**, if downloading fresh rather than pointing at existing data: a single ACS or WFC3 exposure pair is commonly 150–350 MB, and the pipeline keeps raw, drizzled, and cutout products per target per filter — budget on the order of 1–2 GB per target per filter for a rough estimate, so a ~50–80 target catalogue in 2 filters can be a few hundred GB.
- A free [MAST account](https://auth.mast.stsci.edu) **only if** the proposal you're downloading is still inside its exclusive-access period (new HST GO data is normally proprietary for ~1 year). Already-public data (this includes both SNAP programs this pipeline was built for, 15867 and 16773) downloads anonymously — nothing to set up. If a download step returns zero results for a target you know was observed, this is the first thing to check (see [Known limitations](#known-limitations--gotchas) for how to tell, and what `download()` would need in order to authenticate — it doesn't do so automatically today).

### 2. Clone and create the environment

```bash
git clone https://github.com/<your-username>/my-astro-tools.git
cd my-astro-tools/hst-pipeline

conda env create -f environment.yml
conda activate hst-pipeline
```

> **Note on drizzlepac:** it's an STScI package with compiled C extensions. Installing via conda (`conda-forge`) is strongly recommended over the `requirements.txt` pip fallback, which may fail to build on some platforms. If the environment solve is slow, try [mamba](https://mamba.readthedocs.io) as a drop-in replacement: `mamba env create -f environment.yml`.

### 3. Verify the install

Before touching any real data, confirm every dependency actually imports — this catches a bad environment solve immediately instead of mid-batch-run on whichever package happens to be needed by the first target that reaches that code path:

```bash
python -c "import astropy, astroquery, drizzlepac, photutils, astroscrappy, reproject, h5py, lenstronomy, aplpy; print('OK')"
```

If this prints `OK`, the environment is sound. If drizzlepac fails to import, it very likely wasn't installed from conda-forge (see the note in step 2).

### 4. Configure `config.py`

Open `config.py` and update every path in the `Directories` block at the top to point at your actual data storage — see [Setup — config.py](#setup--configpy) below for what each one is for. This is the one file every script imports from, so it's also the one place a stale path silently breaks everything downstream.

### 5. Adding a new proposal? Edit `batch_runner.py` too

`config.py`'s paths are proposal-agnostic, but a few proposal-specific details live only in `batch_runner.py`'s `PROPOSAL_META` and `PROPOSAL_CSVS` dicts (filters, camera, CSV path) — add an entry there for any proposal ID not already present, or `batch_runner.py --proposal <id>` will reject it outright with "invalid choice". If you're only running on 15867/16773/17307, skip this step. In either case, read [Known limitations](#known-limitations--gotchas) below before processing a *filter* this pipeline hasn't seen before (F140W/F200LP/F606W) — several places assume exactly those three.

### 6. Prepare your targets CSV

Copy `targets_template.csv` to `{proposal_id}_targets.csv` (matching whatever you set `ACTIVE_PROPOSAL_ID` to in step 4) and place it inside `MAIN_DIR`. Fill in your target list — see [Required input files](#required-input-files) below for the exact column spec. `catalogue_objname` must match the target name as it appears in the MAST archive exactly, or the download step (`hst_reduction.py` step 1) silently returns zero results for that row.

### 7. (Optional) PSF models

Only needed for `hst_products.py` step 2 (offset-fitting) and for sharper cutouts in step 1 — the pipeline runs without them (step 1 just skips PSF deconvolution with a `[WARN] No PSF found`). See [Required input files](#required-input-files) for the exact path and HDF5 key expected.

### 8. First run

**Recommended, for a whole catalogue** — `batch_runner.py` wraps both scripts below per target with resumable progress tracking, which matters once you're past a handful of targets:

```bash
python batch_runner.py --proposal 17307 --seed          # once, to populate the tracker
python batch_runner.py --proposal 17307 --batch-size 10 # re-run to process the next 10
```

See [batch_runner.py](#batch_runnerpy) below for the full usage and how to handle a `failed` row.

**For a single target, or to debug one step at a time** — call the two scripts directly (no tracker, no resumability):

```bash
python hst_reduction.py --target AGEL110725+245943A   # download, drizzle, reproject, green image
python hst_products.py  --target AGEL110725+245943A   # cutouts, offset fit, postage stamp
```

Drop `--target ...` from either command to run the full CSV instead of one target. Both scripts also accept `--steps` for partial runs (e.g. `--steps 234` to skip the download step) — see their own sections below for every flag.

Either path, **run them in order and let each step finish before the next**: drizzle before reproject, reproject before the green image, and both reduction scripts entirely before `hst_products.py`. The pipeline does not check this for you — skipping ahead produces a deep, unhelpful crash inside `reproject`/`aplpy` (a missing-file error several stack frames down) rather than a clear "run step X first" message.

### 9. Check it worked

For one target, you should now have (inside `MAIN_DIR/{objname}/HST/{proposal_id}_{filter}/`): a `raw_data/` folder with the downloaded `*_flt.fits`/`*_flc.fits`, `*_drz_sci.fits`/`*_drc_sci.fits` (the drizzled mosaic), `*_cutout_L3.fits` and `*_cutout_L3.h5` (the lenstronomy-ready cutout), and a PNG thumbnail under `POSTAGE_STAMP_DIR`. If any of those are missing, work backwards through step 8's console output — every step prints a line per target, and a skip/failure says why.

From here, the detailed reference sections below cover each script's full flag set and the specific defects `pipeline_qc.py` works around.

---

## Known limitations & gotchas

Things that won't stop a correct setup from running, but will produce a confusing result or an unhelpful crash if you don't know about them going in.

### Hardcoded to exactly three filters

`hst_reduction.py`'s `run_step3` (reproject) and `run_step4` (green image), and `hst_products.py`'s postage-stamp colour-channel logic, each hardcode the specific combinations `15867/F140W`, `16773/F140W`, `16773/F200LP`, `17307/F606W` — independently of each other and of `config.py`'s `ALL_PROPOSALS`, which despite its name and comment is **not read by any script** (confirmed by grepping the codebase; it exists only as documentation of what the hardcoded spots below assume). Concretely:

- A target in a proposal/filter combination not in that list is silently never reprojected and never gets a green channel — no error, it just never matches the glob patterns inside `run_step3`/`run_step4`.
- `hst_products.py`'s `cam_map` defaults an unrecognised filter to camera `'WFC3'`, which is simply wrong if that filter is actually ACS (or something else).
- The ≥3-filter postage-stamp branch hardcodes red=F140W/green=F606W/blue=F200LP by literal filter name; a third-plus filter outside that exact set produces no 3-colour image at all.

Adding a fourth proposal or a filter besides F140W/F200LP/F606W means editing all of: `run_step3` and `run_step4` in `hst_reduction.py`, the `cam_map` and 3-colour block in `hst_products.py`, and `PROPOSAL_META`/`PROPOSAL_CSVS` in `batch_runner.py`. Each of those spots has an inline comment pointing back here.

### MAST downloads and proprietary data

`hst_reduction.py`'s `download()` queries MAST anonymously — correct and sufficient for already-public data (this covers both SNAP programs, 15867 and 16773). New HST GO/GAP observations are normally proprietary for about a year after they're taken; for those, an anonymous query returns zero results with no error explaining why (you'll just see an empty `obs_table` printed, then nothing downloaded). This codebase does not call `astroquery.mast.Observations.login()` anywhere. If you need proprietary data, you'd add that call yourself near the top of `download()` in `hst_reduction.py`, e.g.:

```python
from astroquery.mast import Observations
Observations.login(token='<your MAST API token>')  # from https://auth.mast.stsci.edu
```

### The download step always prompts, with or without `--interactive`

`run_step1` (download) calls `_prompt_download_dir()` unconditionally — it asks you to confirm or override the download directory every time, regardless of whether you passed `--interactive`. This is a blocking `input()` call: if you ever run `hst_reduction.py` with its stdin redirected or closed (cron, CI, a background process), step 1 will hang waiting for input that will never arrive. Run interactively, or skip step 1 (`--steps 234`) and stage raw data another way (see next point).

### `hst_reduction.py` never uses `HST_DATA_CACHE`

Only `batch_runner.py`'s `stage_raw_data()` checks `HST_DATA_CACHE` before falling back to MAST. Calling `hst_reduction.py` directly always re-downloads from MAST regardless of what's already cached — another reason `batch_runner.py` is the recommended entry point for anything beyond single-target debugging.

### `--preview` needs a display

`hst_products.py --preview` calls `matplotlib`'s `plt.show()` to pop up an interactive contact sheet. On a remote/headless machine with no GUI backend available, this will do nothing useful (or hang, depending on the backend matplotlib resolves to) rather than show you anything. It's a debugging convenience, not required for any pipeline step — just omit `--preview` on a headless machine.

### Step ordering is not enforced

Nothing checks that step 2 (drizzle) finished before step 3 (reproject) runs, that step 3 finished before step 4 (green image), or that `hst_reduction.py` finished before `hst_products.py` starts. Running out of order produces a deep, unhelpful crash several stack frames inside `reproject` or `aplpy` (a missing-file error, not a clear "run step X first") rather than failing cleanly at the point you actually made the mistake.

---

## Setup — config.py

Open `config.py` and set the following before running anything.

### Directories

```python
MAIN_DIR          # root data directory; per-target folders live inside
DATA_DIR          # where cutout HDF5/FITS files are written (usually same as MAIN_DIR)
LENS_PROC_DIR     # directory containing PSF models and supplementary files
POSTAGE_STAMP_DIR # output directory for postage-stamp PNGs
HST_DATA_CACHE    # batch_runner.py only: a pre-downloaded raw-data cache (same
                  # layout as MAIN_DIR) it checks before falling back to MAST
```

### Active proposal (reduction and cutout steps)

```python
ACTIVE_PROPOSAL_ID = '17307'   # proposal ID string matching MAST and your CSV filename
ACTIVE_FILTERS     = ['F606W'] # list — multiple filters processed in one run
ACTIVE_CAMERA      = 'ACS'     # 'ACS' or 'WFC3'
ACTIVE_TARGETS_CSV             # auto-set to MAIN_DIR/{ACTIVE_PROPOSAL_ID}_targets.csv
```

### All proposals

```python
ALL_PROPOSALS   # documentation only -- see "Known limitations" below; not actually
                # read by any script. Keep it updated anyway as a record of what
                # the hardcoded combinations elsewhere assume.
PROPOSAL_CSVS   # dict mapping proposal_id → path to that proposal's targets CSV;
                # this one IS used, by hst_products.py's postage-stamp step
```

### Other settings

```python
DEFAULT_CUTOUT_SIZE_ARCSEC = 25.0   # cutout size for all targets
OFFSET_NUM_DUP             = 4      # PSO repetitions in the offset-fitter step
PARENT_CATALOGUE_CSV               # full AGEL parent catalogue (for second-source redshifts)
```

### Drizzle / CR-rejection defaults

```python
IR_DRIZZLE_DEFAULTS         # AstroDrizzle + TweakReg params for F140W (step 2)
UV_DRIZZLE_DEFAULTS         # AstroDrizzle + TweakReg params for UV filters (step 2)
LACOSMIC_DEFAULTS           # astroscrappy params, used only when a UV filter has exactly 2 exposures
SOURCE_PROTECTION_DEFAULTS  # repeat-detection params for the N=2 CR path
PAIRWISE_CR_DEFAULTS        # output-frame pairwise-rejection params for the N=2 CR path
```

See [Astrometry & CR-rejection fixes](#astrometry--cr-rejection-fixes) for what the last three do and why they exist.

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

Place HDF5 PSF models directly inside `LENS_PROC_DIR` (not in a subfolder — note that `LENS_PROC_DIR` itself is typically already named `.../lens_processing` in `config.py`'s example path):

```
LENS_PROC_DIR/psf_model_{band}.h5
```

Each file must contain a dataset named `kernel_point_source` (2D array) for `hst_products.py` step 2 (offset-fitting) — step 1's cutout builder is more lenient and will also accept a dataset named `kernel`, `psf`, `psf_kernel`, `data`, or `image`, but `kernel_point_source` is the only name both steps agree on, so use it. If a PSF file is not found, the cutout is still built without one and a `[WARN] No PSF found` is printed — this is not fatal.

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

**`MAIN_DIR` should contain only per-target directories** (each with its own `HST/` subfolder as above) — several steps do `for target_dir in MAIN_DIR.iterdir(): ...` with no filtering beyond "is a directory", so a stray folder dropped into `MAIN_DIR` (a `.git/` if you mistakenly version-control it directly, a half-renamed target, etc.) is picked up and treated as a target. `hst_products.py`'s postage-stamp step now skips anything without an `HST/` subdirectory rather than crashing, but the reduction steps in `hst_reduction.py` don't have that guard.

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

# Interactively review/override drizzle parameters and confirm each target
python hst_reduction.py --steps 2 --interactive
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--steps` | `1234` | Steps to run, e.g. `"234"` to skip download |
| `--target` | *(full CSV)* | Run only this AGEL target name |
| `--interactive`, `-i` | `False` | Prompt to review AstroDrizzle/TweakReg parameters before step 2, and confirm/skip each target as it runs |

### Notes

- Steps 3 and 4 scan **all** subdirectories of `MAIN_DIR`, not just the active proposal, so they update the full dataset each time they run.
- Step 2 picks `pixel_size = 0.08` for F140W and `0.05` for all other filters automatically.
- AstroDrizzle and TweakReg parameters are set to sensible defaults (`IR_DRIZZLE_DEFAULTS` / `UV_DRIZZLE_DEFAULTS` in `config.py`) but may need per-target adjustment for challenging data. Check the outputs of at least a handful of targets before running the full catalogue.
- After both filters for a target finish drizzling, step 2 automatically runs the astrometric preflight/correction described below (`reconcile_target_astrometry`). This is silent when it succeeds (look for `[Defect 1 fix]` in the console output) and loud when it can't run automatically — see [Astrometry & CR-rejection fixes](#astrometry--cr-rejection-fixes).
- When a UV filter has exactly 2 exposures (as F200LP typically does for GO-16773), step 2 runs a heavier cosmic-ray-rejection path (`_uv_drizzle_n2` in `hst_reduction.py`) instead of the standard `driz_cr` pass — this takes noticeably longer than a 3+ exposure filter and prints `[Defect 5b]` / `[Defect 5d]` progress lines. This is expected, not a hang.

---

## Astrometry & CR-rejection fixes

A 2026 audit of a diagnostic reduction of DESJ0206-0114 (GO-16773) found the first five defects below, producing a ~65-170 mas F140W↔F200LP misalignment and a cosmic-ray-ridden F200LP mosaic (F200LP has only 2 exposures, where AstroDrizzle's own median-based CR rejection is unreliable). Defect 6 and 4b were found later, triaging real batch-processing failures. All are fixed in `pipeline_qc.py` and wired into `hst_reduction.py` / `config.py`; two of the original five are also fixed in `Reproject_and_Rescale.ipynb`. This section is the reference for what changed and why — all of these were silent failures (nothing errored, they just produced a slightly-wrong, CR-contaminated, or misaligned mosaic) except 4b, which produces the same silent-corruption risk as 4 but on a fit that *converges* onto garbage rather than failing outright.

| # | Defect | Fix |
|---|---|---|
| 1 | Reprojecting/aligning onto whichever filter's mosaic happens to be used as the display reference discards the only trustworthy absolute astrometry, if the two filters' independent WCS solutions disagree (one has a genuine *a posteriori* Gaia/GSC fit, the other doesn't) | `pipeline_qc.select_astrometric_reference()` reads the `WCSTYPE` keyword on the raw exposures (not a hardcoded "UV is always better" assumption) to pick a reference automatically, or raises `AmbiguousAstrometricReference` and refuses to proceed if it can't. `reconcile_target_astrometry()` then corrects the non-reference filter's native drizzled header in place via matched-source centroids, self-verifying that the correction actually reduced the residual before keeping it. |
| 2 | AstroDrizzle's `resetbits` defaults to `"4096"`, silently erasing cosmic-ray flags (archive's own + any written by a prior pass) if a later call uses that default | Already correct throughout this codebase (`resetbits=8192` used everywhere that must preserve bit 4096) — `pipeline_qc.assert_dq_flag_fraction()` adds a loud assertion after each drizzle call as a regression guard. |
| 3 | TweakReg's `imagefindpars` ignores the DQ array (`dqbits=""`) and source shape (`use_sharp_round=False`) by default, so cosmic rays get treated as candidate alignment sources | `dqbits` and `use_sharp_round=True` now set on every TweakReg call (the UV path previously excluded the wrong bit — `~8192`, which nothing sets — instead of `~4096`, where the real flags are); IR's `conv_width` fixed from 3.5 (tuned for well-sampled UVIS/ACS) to 2.5 (matches WFC3/IR's undersampled PSF). |
| 4 | `updatehdr=True` stamps the "aligned" WCSNAME onto an exposure even when its fit failed (too few matches) — drizzlepac sets `NaN` shifts internally but writes the header anyway | `pipeline_qc.check_shiftfile_for_nan()` parses the shiftfile after every TweakReg call and raises before that header ships. |
| 4b | A fit that *doesn't* fail can still be garbage: with the shared `tolerance=10.0` (loose, to tolerate real sub-pixel dithers) and a sparse source field, `xyxymatch` can pair noise peaks and "converge" with no NaNs — e.g. a real 2026-09-17 batch run matching two exposures of one visit at rot=359.96°, scale=1.0015, rms≈4.4 px, which should never have been accepted | `pipeline_qc.check_shiftfile_fit_quality()` runs after `check_shiftfile_for_nan()` on every IR TweakReg call and warns (does not yet raise — see its docstring) when rms/rotation/scale deviate from the sub-pixel, rot≈0, scale≈1 solution two exposures of one visit should produce. For a target this catches, lower `tweak_threshold`/`tolerance` (and usually `tweak_minobj`) for that target specifically via `TWEAK_PARAM_OVERRIDES` in `config.py` rather than loosening the shared defaults — see the comment above that dict for two worked examples. |
| 5 | N=2 median-based `driz_cr` rejection degenerates to an average and can't reliably separate a cosmic ray from an undersampled point source | New `_uv_drizzle_n2()` path in `hst_reduction.py`, used automatically whenever a UV filter has exactly 2 exposures: (a) L.A.Cosmic in PSF-convolution mode (`astroscrappy`), (b) source protection via repeat detection between the two exposures, (c) `combine_type='minmed'` pinned explicitly, (d) pairwise rejection comparing the two exposures directly in the output frame, mapped back to detector coordinates via `drizzlepac.pixtopix`. |
| 6 | TweakReg's `computesig=True` default estimates sky sigma assuming raw counts; these images are count-rate (electrons/s), so the auto-estimate is inflated by ~√exptime, which can suppress real source detections | `pipeline_qc.compute_skysigma()` computes an empirical, DQ-masked sigma from the actual `*crclean.fits` frames and every TweakReg call passes it explicitly with `computesig=False` (`imagefindcfg` and `refimagefindcfg` alike — the reference image needs the same fix, not just the inputs). |

New `config.py` blocks: `LACOSMIC_DEFAULTS`, `SOURCE_PROTECTION_DEFAULTS`, `PAIRWISE_CR_DEFAULTS` (all only used by the N=2 CR path); `TWEAK_PARAM_OVERRIDES` (Defect 4b, per-target IR TweakReg overrides keyed by `(proposal_id, objname)`).

Separately (not one of the audit's silent-corruption defects, just disk hygiene): every `AstroDrizzle()` call now passes `preserve=False`. Its default (`True`) backs up pristine copies of its inputs into `OrIg_files/` in the process's current working directory before modifying them — not under the target's own folder, under wherever you happened to launch the script from. Since every call here only ever touches `temp/` (a disposable copy of `raw_data/`, deleted at the end of a successful run), that backup is always redundant with the untouched original in `raw_data/`, and left unset it silently accumulates gigabytes of duplicate raw frames in your working directory, one `OrIg_files/` per target processed.

### Validating a target end-to-end

`acceptance_test.py` measures the numbers this audit cares about against an already-drizzled target:

```bash
python acceptance_test.py --main-dir /path/to/MAIN_DIR --target DESJ0206-0114 --proposal 16773 --camera WFC3
```

It reports, independently:
1. **Matched-source centroid** F140W↔F200LP residual (the primary, trusted measurement).
2. **Whole-frame cross-correlation** residual, as an independent cross-check — expect these two to disagree on sparse/CR-heavy fields; report the disagreement rather than picking whichever number looks better. An earlier version of the astrometric fix used *only* the cross-correlation method and it measured a ~19 arcsec "correction" on this exact target — three orders of magnitude too large — because a sparse, cosmic-ray-heavy field's whole-frame FFT correlation is dominated by CR/noise/edge structure rather than the handful of real sources. Matched-source centroids don't have that failure mode; keep the cross-correlation number as a sanity check, not the source of truth.
3. **F200LP CR-like detections** (DAOStarFinder, 5σ, sharpness > 0.8) in the full-depth mosaic.
4. **F200LP pixel σ** in the full-depth region.
5. **Aperture photometry vs. the per-pixel minimum of the two raw, unprocessed single-exposure drizzles** — a biased-low, CR-free flux floor. A real source should never fall below it.

Known limitation of check 5: it builds the flux floor from the *raw, untweaked* archive exposures, while the actual mosaic is built from *tweaked* (TweakReg-corrected) exposures. On DESJ0206-0114 this introduced its own ~100+ mas registration offset between the floor images and the real mosaic — enough to bias compact-source aperture photometry on its own, independent of any real CR-masking flux loss. Treat a "below floor" result from this check as inconclusive until the floor is rebuilt from the same tweaked exposures the real mosaic uses (not yet implemented — would need `_uv_drizzle_n2` to keep its post-tweakback FLC copies instead of deleting them with the rest of `temp/`).

### Interactive mode (`--interactive`)

Run step 2 with `--interactive` (or `-i`) to review drizzle settings without editing `config.py`:

1. At the start of step 2 you're asked whether to review the IR (F140W) and/or UV drizzle parameters. Answering yes walks through each value (pixel scale, `final_pixfrac`, TweakReg `threshold`/`conv_width`/`searchrad`/`ylimit`, `driz_cr_scale`, `driz_cr_grow`) with the current default shown — press Enter to keep it, or type a new value. These overrides apply to every target for the rest of the run.
2. Before each target's drizzle runs, you're prompted `[Y]es / [n]o (skip) / [a]ll (stop asking) / [q]uit` — useful for skipping a known-bad target or bailing out early without killing the process.

Overrides made this way are for that run only. To change the defaults permanently, edit `IR_DRIZZLE_DEFAULTS` / `UV_DRIZZLE_DEFAULTS` in `config.py`.

---

## Reproject_and_Rescale.ipynb (standalone)

This notebook is kept separate because it requires **manual parameter setting per target** and is only needed for a subset of objects where the WCS-based reprojection in step 3 is insufficiently accurate.

**When to use it:** targets listed with `alt_run = True` in `postage_stamp_config.csv`.

**What it does:** uses OpenCV template matching and `astroalign` to register one filter onto the other, then writes three output files per target:

```
HST/{proposal_id}_F140W/{target}_F140W_L3.fits
HST/{proposal_id}_F140W/{target}_F140W_green_L3.fits
HST/{proposal_id}_{uv_filter}/{target}_{uv_filter}_L3.fits
```

These names are exactly what `hst_products.py` step 3 looks for when `alt_run = True`, so no further changes are needed — just run the notebook for each affected target, then proceed with the products pipeline.

**Which filter gets warped:** a cell right before the `astroalign` registration step (section 2.4) runs the same astrometric preflight check as `pipeline_qc.select_astrometric_reference()` — it looks for the raw exposures under `raw_data/` and sets `IR_IS_REFERENCE` accordingly (defaulting to `False`, i.e. F200LP is trustworthy and F140W moves, if the raw exposures can't be found from this notebook's expected path). This used to be hardcoded to always warp F200LP onto F140W regardless of which filter's astrometry was actually good — see [Astrometry & CR-rejection fixes](#astrometry--cr-rejection-fixes). If you see `[Defect 1 preflight] ... defaulting to ... but VERIFY manually`, check `IR_IS_REFERENCE` makes sense for that target before trusting the output.

---

## hst_products.py

Covers cutout creation, offset alignment, and postage stamp creation.

### Steps

| Step | Flag | What it does |
|---|---|---|
| 1 | `1` | Builds lenstronomy-ready HDF5 + FITS cutouts centred on target RA/Dec |
| 2 | `2` | Fits Sersic centroids in each band; writes `ra_shift`/`dec_shift` into cutout files |
| 3 | `3` | Makes grayscale, 2-colour, or 3-colour PNG postage stamps; respects `--target` for single-target output |

### Usage

```bash
# ── Full CSV (default) ─────────────────────────────────────────────────
# Run all three steps across every target in ACTIVE_TARGETS_CSV
python hst_products.py

# ── Single target already in the CSV ──────────────────────────────────
# All three steps
python hst_products.py --target AGEL110725+245943A

# Postage stamps only
python hst_products.py --target AGEL110725+245943A --steps 3

# ── Single target NOT yet in the CSV ──────────────────────────────────
# --ra and --dec are required; the script will prompt whether to append
# the target to the CSV for future full-CSV runs
python hst_products.py --target AGEL999999+000000A --ra 12.345 --dec -45.678

# With optional redshifts
python hst_products.py --target AGEL999999+000000A \
    --ra 12.345 --dec -45.678 --z-src 1.83 --z-def 0.45

# ── Override proposal ID without editing config.py ────────────────────
python hst_products.py --target AGEL110725+245943A --proposal-id 16773

# ── Preview cutout sizes before committing to a run ───────────────────
python hst_products.py --preview
python hst_products.py --target AGEL110725+245943A --preview
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--steps` | `123` | Steps to run, e.g. `"3"` for postage stamps only |
| `--target` | *(full CSV)* | Run only this AGEL target name |
| `--ra` | — | Target RA in decimal degrees — required if target is not in the CSV |
| `--dec` | — | Target Dec in decimal degrees — required if target is not in the CSV |
| `--z-src` | `None` | Source redshift (stored as `-1` if omitted) |
| `--z-def` | `None` | Deflector redshift (stored as `-1` if omitted) |
| `--proposal-id` | config.py value | Override `ACTIVE_PROPOSAL_ID` for this run only |
| `--preview` | `False` | Show cutout-size contact sheet; prompts before continuing |

### Postage stamp config CSV

`postage_stamp_config.csv` controls per-target appearance. Edit this file rather than the script.

| Column | Type | Description |
|---|---|---|
| `objname` | string | AGEL target name |
| `thumb_size_arcsec` | int | Thumbnail diameter in arcseconds (default 20) |
| `clims_set` | int | Colour-scale preset 1–15 (default 1; see script for limits per preset) |
| `skip` | bool | Set `True` to exclude a target from all postage-stamp output |
| `alt_run` | bool | Set `True` for targets processed with the Reproject notebook (uses different input filenames) |

To add a new target not in the CSV, the script falls back to `thumb_size = 20` and `clims_set = 1`. Add a row to the CSV only when you need non-default values.

### How the postage stamp type is chosen

- **1 filter found** → grayscale (two panels at different stretch)
- **2 filters found** → 2-colour RGB (IR = red, UV = blue, green = synthetic blend)
- **3+ filters found** → 3-colour RGB (F140W = red, F606W = green, F200LP = blue)

---

## batch_runner.py

Drives a whole catalogue end to end — per target: stage raw data (reuse `HST_DATA_CACHE` if present, else download from MAST), `hst_reduction.py` steps 2-4 (drizzle/reproject/green), then `hst_products.py` steps 1+3 (cutouts + postage stamps; step 2, offset-fitting, is skipped by design for batch runs). Progress is tracked in `batch_processing_tracker.csv`, written after *every* target — not just at the end — so a crash, a `Ctrl-C`, or killing the process after target 40 of 80 loses nothing: the tracker already reflects targets 1-40, and the next invocation picks up at 41.

Proposal metadata (filters, camera) lives in `PROPOSAL_META` inside `batch_runner.py` itself, not `config.py` — add an entry there (and to `PROPOSAL_CSVS`) for a new proposal.

### Usage

```bash
# 1. First time for a proposal (and again whenever the proposal's CSV gains
#    targets, e.g. an ongoing GAP program): add any new targets to the
#    tracker as 'pending' rows. Safe to re-run — only adds rows that don't
#    already exist; never touches an existing row's status.
python batch_runner.py --proposal 17307 --seed

# 2. Process the next N pending targets (default 10). Re-run with the same
#    command to keep working through the queue.
python batch_runner.py --proposal 17307 --batch-size 10
python batch_runner.py --proposal 17307 --batch-size 10   # again, for the next 10
```

### CLI flags

| Flag | Default | Description |
|---|---|---|
| `--proposal` | *(required)* | Proposal ID — must be a key in `PROPOSAL_META` |
| `--seed` | `False` | Add new targets from the proposal's CSV to the tracker as `pending` rows, then exit without processing anything |
| `--batch-size` | `10` | Number of `pending` targets to process this run |

### Handling a failure

A target that raises anywhere in the pipeline is caught, marked `status=failed` with the exception's type and message in `notes`, and the run continues to the next target rather than aborting the whole batch. To retry it: either re-run `batch_runner.py` after manually setting that row's `status` back to `pending`, or write a small one-off script in the style of `rerun_tolerance_fix.py` that calls `batch_runner.process_target()` directly on a specific list.

Before retrying, read `notes` — it's the full exception message, and per [Astrometry & CR-rejection fixes](#astrometry--cr-rejection-fixes) and the diagnoses already recorded in `batch_processing_tracker.csv`, most recurring failures have a known cause and fix (sparse-field TweakReg matching → `TWEAK_PARAM_OVERRIDES`; a `NoOverlapError` → check the target's catalogue coordinates against `Parent-Catalogue-All_targets.csv` and the raw exposure headers' `RA_TARG`/`DEC_TARG` before assuming the data itself is bad). If a failure turns out to be genuinely bad data (lost exposure, wrong-pointing visit) rather than something fixable, set `status` to a distinct value like `bad_data` instead of `pending`/`failed` and record why in `notes` — `batch_runner.py` only ever looks for exactly `pending`, so anything else is simply left alone and visible in the tracker as triaged.

---

## Typical full-pipeline run (new proposal)

```bash
# 1. Add proposal to config.py: set ACTIVE_PROPOSAL_ID, ACTIVE_FILTERS, ACTIVE_CAMERA
#    Add entry to ALL_PROPOSALS and PROPOSAL_CSVS (config.py) and to
#    PROPOSAL_META and PROPOSAL_CSVS (batch_runner.py)
#    Place {proposal_id}_targets.csv in MAIN_DIR

# 2. Seed the tracker, then download and reduce every target
python batch_runner.py --proposal {proposal_id} --seed
python batch_runner.py --proposal {proposal_id} --batch-size 20   # repeat until "No pending targets left"

# 3. Inspect drizzled outputs — check alignment and CR rejection for a sample of targets

# 4. For any alt_run targets, open and run Reproject_and_Rescale.ipynb
#    then set alt_run = True in postage_stamp_config.csv for those targets
#    (batch_runner.py's cutout/postage-stamp steps already handle alt_run
#    targets automatically on their next run once that flag is set)
```

Or, to run the two scripts directly instead (no tracker, no resumability — mainly useful for a single target or for debugging one step):

```bash
# 1. Add proposal to config.py: set ACTIVE_PROPOSAL_ID, ACTIVE_FILTERS, ACTIVE_CAMERA
#    Add entry to ALL_PROPOSALS and PROPOSAL_CSVS
#    Place {proposal_id}_targets.csv in MAIN_DIR

# 2. Download and reduce
python hst_reduction.py

# 3. Inspect drizzled outputs — check alignment and CR rejection for a sample of targets

# 4. For any alt_run targets, open and run Reproject_and_Rescale.ipynb
#    then set alt_run = True in postage_stamp_config.csv for those targets

# 5. Build cutouts, fit offsets, make stamps
python hst_products.py
```
