# AGEL Survey Tools

A collection of Python pipelines for processing and analysing data from the
[AGEL](https://agel.readthedocs.io) (the ARC Gravitational-arc Eye Legacysurvey) gravitational-lens catalogue.
Two independent packages are included; each has its own environment and README.

---

## Packages

### 1. `agel_xmatch/` — Archive Cross-Match Pipeline

A command-line tool that cross-matches the AGEL lens catalogue against major
astronomical archives, downloads matched data products, and generates postage
stamp images.

**Key features:**
- Cross-matches against JWST (MAST), Chandra, XRISM, XMM-Newton, eROSITA, and
  Swift (HEASARC) in a single run
- Smart JWST filter selection — scores observations by calibration level,
  instrument priority, and wavelength spread to pick the best 1–3 filters per
  source for false-colour imaging
- Parallel FITS downloads for JWST; shell-script generation for X-ray archives
- Postage stamp generation (grayscale, 2-colour, or RGB) with WCS scale bar and
  redshift annotation; per-source FOV and stretch overrides
- Modular stage flags (`--skip-crossmatch`, `--skip-filter`, `--skip-download`,
  `--only-stamps`, `--clean`) so any stage can be re-run in isolation

**Entry point:** `agel_xmatch/archive_pipeline.py`

**Quick start:**
```bash
cd agel_xmatch
conda env create -f environment.yml && conda activate archive-pipeline
python archive_pipeline.py --archive jwst
python archive_pipeline.py --archive jwst chandra xmm
```

See [agel_xmatch/README.md](agel_xmatch/README.md) for the full configuration
reference, pipeline stage diagram, and output directory layout.

---

### 2. `hst_pipeline/` — HST Reduction & Products Pipeline

Two Python scripts that take raw HST FLC/FLT files through drizzling,
reprojection, lenstronomy-ready cutout creation, and final postage stamp
generation. A shared `config.py` and a per-target CSV control all settings
without touching the scripts themselves.

**Key features:**
- Downloads FLC/FLT files from MAST, aligns exposures with TweakReg, and
  drizzles with AstroDrizzle (IR at 0.08″/pix, UV at 0.05″/pix)
- Reprojects UV frames onto the IR pixel grid and builds a synthetic green
  channel (`green = sqrt(IR² + UV²)`) for RGB composites
- Produces lenstronomy-ready HDF5 + FITS cutouts with Sérsic centroid offsets
- Grayscale, 2-colour, or 3-colour postage stamps depending on filter
  availability; per-target thumbnail size and colour scale via a CSV config file
- `Reproject_and_Rescale.ipynb` for the subset of targets where WCS reprojection
  is insufficient — uses OpenCV template matching and `astroalign` for manual
  fine-tuning
- `--steps`, `--target`, `--preview`, and other CLI flags for partial runs and
  single-target testing

**Entry points:**
| Script | Role |
|---|---|
| `hst_pipeline/hst_reduction.py` | Download → drizzle → reproject → synthetic green |
| `hst_pipeline/hst_products.py` | Cutouts → offset fitting → postage stamps |

**Quick start:**
```bash
cd hst_pipeline
conda env create -f environment.yml && conda activate hst-pipeline
# Edit config.py to set your paths and active proposal
python hst_reduction.py
python hst_products.py
```

See [hst_pipeline/README.md](hst_pipeline/README.md) for setup instructions,
the full `config.py` reference, required input file formats, and a walkthrough
for adding a new proposal.

---

## Repository layout

```
AGEL_survey/
├── agel_xmatch/
│   ├── archive_pipeline.py   # single-script archive cross-match + stamp pipeline
│   ├── environment.yml
│   └── requirements.txt
└── hst_pipeline/
    ├── hst_reduction.py      # download, drizzle, reproject, green channel
    ├── hst_products.py       # cutouts, offset fitting, postage stamps
    ├── config.py             # all paths and proposal settings — edit first
    ├── postage_stamp_config.csv  # per-target thumbnail and colour-scale overrides
    ├── targets_template.csv  # template for your per-proposal targets CSV
    ├── Reproject_and_Rescale.ipynb  # manual alignment for problem targets
    ├── environment.yml
    └── requirements.txt
```

---

## Requirements

Both packages require **Python ≥ 3.9** and conda (recommended). The two
environments are independent — install them separately.

| Package | Core extra dependencies |
|---|---|
| `agel_xmatch` | `astropy`, `astroquery`, `reproject` |
| `hst_pipeline` | all of the above + `drizzlepac`, `lenstronomy`, `h5py`, `astroalign`, `opencv-python` |

> **drizzlepac note:** `drizzlepac` has compiled C extensions. Installing via
> conda (`conda-forge`) is strongly recommended over pip, which may fail on some
> platforms. If the environment solve is slow, try
> [mamba](https://mamba.readthedocs.io) as a drop-in replacement:
> `mamba env create -f environment.yml`.

A free [MAST account](https://auth.mast.stsci.edu/login) with `astroquery`
credentials configured is required for any step that downloads HST or JWST data.

---

## Data organisation (hst_pipeline)

The HST pipeline expects data laid out as:

```
MAIN_DIR/
└── {AGEL_name}/
    └── HST/
        └── {proposal_id}_{filter}/
            └── raw_data/    ← FLC/FLT files
```

See the preface note in [hst_pipeline/README.md](hst_pipeline/README.md) for
the recommended workflow for building this directory tree from a MAST
observation CSV and the AGEL parent catalogue.

---

## License

MIT
