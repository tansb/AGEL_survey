# AGEL Archive Pipeline

A command-line pipeline that cross-matches the AGEL gravitational lens catalog
against major astronomical archives, downloads the matched data, and generates
postage stamp images.

## Features

- **Multi-archive cross-matching** via `astroquery`: JWST (MAST), Chandra,
  XRISM, XMM-Newton, eROSITA, and Swift (HEASARC) — query one or several in a
  single run
- **Smart JWST filter selection**: scores observations by calibration level,
  instrument priority, and wavelength spread to pick the best 1–3 filters per
  source for false-colour imaging
- **Parallel FITS downloads** for JWST science files; shell-script generation
  for X-ray archives
- **Postage stamp generation**: grayscale (single filter), 2-colour, or RGB
  composites with WCS scalebar and redshift titles; per-source FOV and stretch
  overrides
- **`--clean` flag**: removes full field-of-view FITS files after stamps are
  saved, keeping only the cutout PNGs

## Supported Archives

| Key        | Service | Catalog       |
|------------|---------|---------------|
| `jwst`     | MAST    | —             |
| `chandra`  | HEASARC | `chanmaster`  |
| `xrism`    | HEASARC | `xrismmastr`  |
| `xmm`      | HEASARC | `xmmmaster`   |
| `erosita`  | HEASARC | `erassmastr`  |
| `erosmaster` | HEASARC | `erosmaster` |
| `swift`    | HEASARC | `swiftmastr`  |

## Installation

Create and activate the conda environment (recommended):

```bash
git clone https://github.com/<your-org>/agel-archive-pipeline.git
cd agel-archive-pipeline
conda env create -f environment.yml
conda activate archive-pipeline
```

Or install dependencies via pip only:

```bash
pip install -r requirements.txt
```

## Quick Start

```bash
# Cross-match against JWST only, download FITS, generate stamps
python archive_pipeline.py --archive jwst

# Cross-match against JWST and Chandra together
python archive_pipeline.py --archive jwst chandra

# Chandra + XRISM only
python archive_pipeline.py --archive chandra xrism

# Re-use existing cross-match CSVs (skip the slow MAST/HEASARC queries)
python archive_pipeline.py --archive jwst --skip-crossmatch

# Regenerate postage stamps from already-downloaded FITS files
python archive_pipeline.py --archive jwst --only-stamps

# Full run, then delete the large FoV FITS, keeping only the stamp PNGs
python archive_pipeline.py --archive jwst --clean
```

## Configuration

All user-editable settings live at the top of `archive_pipeline.py` in clearly
marked sections. Key variables:

| Variable | Description |
|---|---|
| `INPUT_CATALOG_PATH` | Path to your source catalog CSV (`ra`, `dec`, `objname` columns) |
| `PARENT_CATALOG_PATH` | Path to catalog with redshifts (`RAJ2000`, `DECJ2000`, `z_source`, `z_deflector`) |
| `MAIN_DIR` | Root output directory for all downloaded files and stamps |
| `SEARCH_RADIUS_ARCSEC` | Cross-match search radius in arcseconds (default 15) |
| `DEFAULT_ARCHIVES` | Archives used when `--archive` is not supplied on the CLI |
| `JWST_N_FILTERS` | Target number of distinct JWST filters per source (1–3) |
| `INSTRUMENT_PRIORITY` | Preferred JWST instruments in order: `['NIRCAM', 'MIRI', 'NIRISS']` |
| `CATALOG_SLICE` | Tuple `(start, end)` to process a row slice of the input catalog, or `None` for all |
| `TEST_LIMIT` | Limit downloads to N files for testing; `None` for all |
| `DEFAULT_THUMB_SIZE_ARCSEC` | Default postage stamp size in arcseconds |
| `CUSTOM_THUMB_SIZES` | Per-source FOV overrides `{source_id: arcsec}` |
| `CUSTOM_STRETCHES` | Per-source percentile stretch `{source_id: (pmin, pmax)}` |
| `CUSTOM_RGB_LIMITS` | Per-source per-channel vmin/vmax or percentile limits |

## Pipeline Stages

```
Stage 1 – Cross-match   →  per-archive CSV files in MAIN_DIR/{archive}/
Stage 2 – Filter        →  {archive}_filtered.csv (JWST: best imaging per source)
Stage 3 – Download      →  JWST FITS in MAIN_DIR/jwst/preview_files/
                           X-ray: download_{archive}_files.sh shell scripts
Stage 4 – Stamps        →  PNG images in MAIN_DIR/jwst/postage_stamps/
[--clean]               →  deletes FITS files, keeps PNGs
```

Each stage can be skipped independently:

```bash
python archive_pipeline.py --archive jwst \
    --skip-crossmatch \   # load existing CSVs from stage 1
    --skip-filter \       # load existing filtered CSVs from stage 2
    --skip-download       # assume FITS already on disk
```

## Output Structure

```
MAIN_DIR/
├── jwst/
│   ├── jwst_agel_archive_matches.csv   # raw cross-match results
│   ├── jwst_filtered.csv               # selected observations for download
│   ├── preview_files/                  # downloaded FITS files
│   │   ├── AGEL003428+022522A_jw…fits
│   │   └── source_file_mapping.csv
│   └── postage_stamps/
│       ├── AGEL003428+022522A_F115W_F200W_F444W_40arcs_RGB.png
│       └── individual_filters/         # per-filter grayscale (optional)
├── chandra/
│   ├── chandra_agel_archive_matches.csv
│   ├── chandra_filtered.csv
│   └── preview_files/
│       ├── download_chandra_files.sh   # run to fetch event files
│       └── chandra_obsids.txt
└── xrism/
    └── …
```

## Requirements

- Python ≥ 3.9
- See [environment.yml](environment.yml) (conda) or [requirements.txt](requirements.txt) (pip) for the full dependency list

## Legacy Notebooks

The `crossref.ipynb`, `jwst_chandra_retrieval_complete.ipynb`, and
`jwst_postage_stamps.ipynb` notebooks are retained for reference. All of their
functionality is consolidated in `archive_pipeline.py`.

## License

MIT
