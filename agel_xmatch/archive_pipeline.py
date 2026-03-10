#!/usr/bin/env python3
"""
AGEL Archive Cross-Match Pipeline
==================================
Combines the functionality of three notebooks into a single configurable script:
  1. crossref.ipynb                      - Cross-match catalog with archives
  2. jwst_chandra_retrieval_complete.ipynb - Filter observations and download data
  3. jwst_postage_stamps.ipynb           - Generate postage stamp images

Supported archives:
  - JWST     (via MAST / astroquery)
  - Chandra  (via HEASARC: chanmaster)
  - XRISM    (via HEASARC: xrismmastr)
  - XMM-Newton (via HEASARC: xmmmaster)
  - eROSITA  (via HEASARC: erassmastr, erosmaster)
  - Swift    (via HEASARC: swiftmastr)

Usage:
  python archive_pipeline.py --archive jwst
  python archive_pipeline.py --archive jwst chandra
  python archive_pipeline.py --archive chandra xrism
  python archive_pipeline.py --archive jwst --skip-crossmatch   # re-use saved CSVs
  python archive_pipeline.py --archive jwst --only-stamps       # regenerate stamps only
  python archive_pipeline.py --archive jwst --clean             # delete full-FoV FITS after stamps

Available archive names: jwst, chandra, xrism, xmm, erosita, erosmaster, swift
"""

# ============================================================
# ===                   CONFIGURATION                      ===
# ============================================================

# --- Input catalog (for cross-matching) ---
# Must contain columns for RA, Dec, and object name.
INPUT_CATALOG_PATH = "/Users/cwat/research/prev_survey_tabs/agel_coordsonly.csv"
INPUT_RA_COL        = "ra"
INPUT_DEC_COL       = "dec"
INPUT_NAME_COL      = "objname"

# Slice of the input catalog to process (Python slice notation).
# Set to None to process the entire catalog.
CATALOG_SLICE = None  # e.g. (3000, 4000) to process rows 3000–3999

# Only query sources whose name ends in 'A' (primary source in compound systems).
QUERY_PRIMARY_ONLY = True

# --- Parent catalog (for postage stamp coordinates and redshifts) ---
# May be the same file as INPUT_CATALOG_PATH if it contains RA/Dec + redshifts.
# Set to None to use the cross-match result coordinates instead.
PARENT_CATALOG_PATH     = "/Volumes/AGEL/Parent Catalogue-All_targets.csv"
PARENT_RA_COL           = "RAJ2000"
PARENT_DEC_COL          = "DECJ2000"
PARENT_NAME_COL         = "objname"
PARENT_Z_SOURCE_COL     = "z_source (from DR2 Redshifts)"
PARENT_Z_DEFLECTOR_COL  = "z_deflector (from DR2 Redshifts)"

# --- Output root directory ---
from pathlib import Path
MAIN_DIR = Path("/Volumes/AGEL/agel_xmatch/")

# --- Search radius ---
SEARCH_RADIUS_ARCSEC = 15   # arcseconds

# ============================================================
# ===              ARCHIVE SELECTION                       ===
# ============================================================
# Default archives to query when --archive is not supplied on the command line.
# Override at runtime with:  python archive_pipeline.py --archive jwst chandra
#
# Valid names: jwst, chandra, xrism, xmm, erosita, erosmaster, swift
DEFAULT_ARCHIVES = ["jwst", "chandra", "erosita", "erosmaster"]

# Runtime-resolved set – populated automatically in main(); do not edit.
ACTIVE_ARCHIVES: set = set()

# HEASARC catalog names and preferred columns for each X-ray archive.
# Columns that don't exist in a result are silently skipped.
HEASARC_ARCHIVE_CONFIGS = {
    "chandra": {
        "catalog":      "chanmaster",
        "display_name": "Chandra",
        "columns":      ["obsid", "name", "ra", "dec", "livetime", "detector", "grating", "type"],
        "ra_col":       "ra",
        "dec_col":      "dec",
        "obsid_col":    "obsid",
        "download_url_template":
            "https://cda.harvard.edu/csccli/retrieveFiles?obsid={obsid}&filetypes=evt2",
    },
    "xrism": {
        "catalog":      "xrismmastr",
        "display_name": "XRISM",
        "columns":      ["obsid", "name", "ra", "dec", "livetime", "detector"],
        "ra_col":       "ra",
        "dec_col":      "dec",
        "obsid_col":    "obsid",
        "download_url_template": None,   # Customise per your XRISM data access
    },
    "xmm": {
        "catalog":      "xmmmaster",
        "display_name": "XMM-Newton",
        "columns":      ["obsid", "name", "ra", "dec", "duration", "instruments"],
        "ra_col":       "ra",
        "dec_col":      "dec",
        "obsid_col":    "obsid",
        "download_url_template": None,
    },
    "erosita": {
        "catalog":      "erassmastr",
        "display_name": "eROSITA (eRASS)",
        "columns":      ["obsid", "name", "ra", "dec", "livetime"],
        "ra_col":       "ra",
        "dec_col":      "dec",
        "obsid_col":    "obsid",
        "download_url_template": None,
    },
    "erosmaster": {
        "catalog":      "erosmaster",
        "display_name": "eROSITA (Master)",
        "columns":      ["obsid", "name", "ra", "dec", "exposure"],
        "ra_col":       "ra",
        "dec_col":      "dec",
        "obsid_col":    "obsid",
        "download_url_template": None,
    },
    "swift": {
        "catalog":      "swiftmastr",
        "display_name": "Swift",
        "columns":      ["obsid", "name", "ra", "dec", "uvot_exp", "xrt_exp"],
        "ra_col":       "ra",
        "dec_col":      "dec",
        "obsid_col":    "obsid",
        "download_url_template": None,
    },
}

# ============================================================
# ===          FILTER / RETRIEVAL SETTINGS                 ===
# ============================================================

# JWST filter selection
JWST_N_FILTERS        = 3           # Target distinct filters per source (1–3)
INSTRUMENT_PRIORITY   = ["NIRCAM", "MIRI", "NIRISS"]

# Download
MAX_WORKERS  = 2     # Parallel JWST download threads
TEST_LIMIT   = None  # Set to e.g. 3 to test a few files; None for all

# Directories (created automatically)
JWST_DATA_DIR    = MAIN_DIR / "jwst" / "preview_files"
CHANDRA_DATA_DIR = MAIN_DIR / "chandra" / "preview_files"

# ============================================================
# ===         POSTAGE STAMP SETTINGS                       ===
# ============================================================

STAMP_OUTPUT_DIR          = MAIN_DIR / "jwst" / "postage_stamps"
TEMP_DIR                  = JWST_DATA_DIR / "temp"
DEFAULT_THUMB_SIZE_ARCSEC = 20
DEFAULT_SCALEBAR_SIZE_ARCSEC = 2
CREATE_INDIVIDUAL_GRAYSCALE  = False  # Also save per-filter grayscale stamps

# Per-target thumbnail size overrides  {source_id: arcsec}
CUSTOM_THUMB_SIZES = {
    "AGEL003428+022522A": 40, "AGEL003526-201545A": 80, "AGEL004325-203717A": 100,
    "AGEL010258-491619A": 40, "AGEL022434-000228A": 30, "AGEL022421+000421A": 40,
    "AGEL023953-013456A": 40, "AGEL024525-530145A": 25, "AGEL024803-033145A": 50,
    "AGEL032727-132623A": 60, "AGEL040559-491558A": 60, "AGEL045407-101322A": 30,
    "AGEL055747-411327A": 50, "AGEL065832-555637A": 80, "AGEL100434+411244A": 30,
    "AGEL121219+273353A": 30, "AGEL124817+074258A": 60, "AGEL172336+341158A": 40,
    "AGEL201108-522814A": 50, "AGEL211119-011423A": 40, "AGEL213105-401921A": 30,
    "AGEL224844-443151A": 80, "AGEL230822-021132A": 80, "AGEL232511-411125A": 30,
}

# Per-target stretch overrides  {source_id: (pmin, pmax)}
CUSTOM_STRETCHES = {
    # "AGEL003428+022522A": (2, 98),
}

DEFAULT_STRETCH = (1, 99)   # (pmin, pmax) percentile range

# Per-target per-channel RGB scaling  {source_id: {vmin_r, vmax_r, ...}}
CUSTOM_RGB_LIMITS = {
    "AGEL020613-011417A": {"vmin_r": 0,    "vmax_r": 1.9, "vmin_g": 0.05, "vmax_g": 2.2, "vmin_b": 0.02, "vmax_b": 2.0},
    "AGEL003428+022522A": {"vmin_r": 0,    "vmax_r": 1.9, "vmin_g": 0.05, "vmax_g": 2.2, "vmin_b": 0.02, "vmax_b": 2.0},
    "AGEL003526-201545A": {"vmin_r": 0.1,  "vmax_r": 1.7, "vmin_g": 0.01, "vmax_g": 3.3, "vmin_b": 0.1,  "vmax_b": 1.5},
    "AGEL010258-491619A": {"vmin_r": 0.13, "vmax_r": 0.88,"vmin_g": 0,    "vmax_g": 1.2, "vmin_b": 0.15, "vmax_b": 0.4},
    "AGEL021719-050759A": {"vmin_r": 0.2,  "vmax_r": 1.1, "vmin_g": 0.05, "vmax_g": 1.8, "vmin_b": 0.2,  "vmax_b": 0.5},
    "AGEL022421+000421A": {"vmin_r": 0.0,  "vmax_r": 1.0, "vmin_g": 0.0,  "vmax_g": 1.5, "vmin_b": 0.0,  "vmax_b": 1.5},
    "AGEL022434-000228A": {"vmin_r": 0.35, "vmax_r": 1.6, "vmin_g": 0.05, "vmax_g": 1.8, "vmin_b": 0.3,  "vmax_b": 0.9},
    "AGEL023953-013456A": {"vmin_r": 0.2,  "vmax_r": 3.0, "vmin_g": 0.3,  "vmax_g": 4.0, "vmin_b": -0.3, "vmax_b": 4.5},
    "AGEL024803-033145A": {"vmin_r": 0.0,  "vmax_r": 1.9, "vmin_g": 0.05, "vmax_g": 2.9, "vmin_b": 0.02, "vmax_b": 2.0},
    "AGEL025739-220927A": {"vmin_r": 0.1,  "vmax_r": 1.8, "vmin_g": 0.0,  "vmax_g": 3.6, "vmin_b": 0.1,  "vmax_b": 1.5},
    "AGEL142032+525822A": {"vmin_r": 0.3,  "vmax_r": 1.9, "vmin_g": 0.05, "vmax_g": 4.0, "vmin_b": 0.02, "vmax_b": 2.0},
    "AGEL142955+120236A": {"vmin_r": 0.1,  "vmax_r": 2.0, "vmin_g": 0.0,  "vmax_g": 4.0, "vmin_b": 0.2,  "vmax_b": 1.5},
    "AGEL172336+341158A": {"vmin_r": 0.0,  "vmax_r": 1.5, "vmin_g": 0.05, "vmax_g": 3.0, "vmin_b": 0.1,  "vmax_b": 1.5},
    "AGEL213105-401921A": {"vmin_r": 0.1,  "vmax_r": 3.0, "vmin_g": 0.05, "vmax_g": 5.0, "vmin_b": 0.1,  "vmax_b": 1.9},
    "AGEL230825-021214A": {"vmin_r": 0.1,  "vmax_r": 2.0, "vmin_g": 0.0,  "vmax_g": 4.5, "vmin_b": 0.1,  "vmax_b": 2.5},
}

# ============================================================
# ===                    IMPORTS                           ===
# ============================================================

import os
import sys
import argparse
import warnings
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import requests
import matplotlib
matplotlib.use("Agg")   # non-interactive backend for batch mode
import matplotlib.pyplot as plt
from matplotlib import patheffects

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.table import vstack, Table
from astropy.visualization import (ImageNormalize, AsymmetricPercentileInterval,
                                   PercentileInterval)
from astropy.wcs import WCS

from astroquery.mast import Observations
from astroquery.heasarc import Heasarc

try:
    from reproject import reproject_interp
    REPROJECT_AVAILABLE = True
except ImportError:
    REPROJECT_AVAILABLE = False
    print("WARNING: reproject not installed. RGB channel alignment will be skipped.\n"
          "Install with: pip install reproject")

warnings.filterwarnings("ignore")

# Plot styling (mirrors the postage stamps notebook)
label_size      = 20
cbar_tick_size  = 12
plt.rc("font", family="serif")


# ============================================================
# ===              HELPER UTILITIES                        ===
# ============================================================

def _active_heasarc_archives():
    """Return the list of non-JWST archive keys that are active this run."""
    return [k for k in HEASARC_ARCHIVE_CONFIGS if k in ACTIVE_ARCHIVES]


def save_table_to_csv(astropy_table, filename, output_dir):
    """Save an Astropy Table to CSV."""
    if astropy_table is None or len(astropy_table) == 0:
        print(f"  Skipping save: no rows for '{filename}'")
        return
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    full_path = output_dir / filename
    astropy_table.write(str(full_path), format="ascii.csv", overwrite=True)
    print(f"  Saved {len(astropy_table)} rows -> {full_path}")


# ============================================================
# ===     STAGE 1 : CROSS-MATCHING                        ===
# ============================================================

def _coord_string(coord):
    return (coord.ra.to_string(unit=u.hour, sep=":", precision=2)
            + " "
            + coord.dec.to_string(unit=u.deg, sep=":", precision=1))


def query_jwst(source_data, radius):
    """Query MAST for JWST observations around each source."""
    n = len(source_data)
    print(f"\n--- Querying JWST (MAST) with radius {radius.to_value(u.arcsec):.1f} arcsec "
          f"({n} sources) ---")
    tables = []
    t0 = time.time()
    for i, source in enumerate(source_data, 1):
        coord   = source["coord"]
        obj_id  = source["objname"]
        elapsed = time.time() - t0
        print(f"  [{i}/{n}] {obj_id}  (elapsed {elapsed:.0f}s)", end="")
        try:
            result = Observations.query_region(coord, radius=radius)
            jwst   = result[result["obs_collection"] == "JWST"]
            if len(jwst) > 0:
                print(f"  → {len(jwst)} JWST observation(s)")
                jwst["Input_Coord"] = _coord_string(coord)
                jwst["Source_ID"]   = obj_id
                tables.append(jwst)
            else:
                print("  → no matches")
        except Exception as exc:
            print(f"\n  ERROR querying MAST for {obj_id}: {exc}")
    if not tables:
        return None
    result = vstack(tables)
    result.meta = {}
    return result


def query_heasarc_archive(source_data, radius, archive_key):
    """Query a HEASARC archive for observations around each source."""
    cfg  = HEASARC_ARCHIVE_CONFIGS[archive_key]
    cat  = cfg["catalog"]
    name = cfg["display_name"]
    cols = cfg["columns"]
    n    = len(source_data)

    print(f"\n--- Querying {name} (HEASARC: {cat}) with radius {radius.to_value(u.arcsec):.1f} arcsec "
          f"({n} sources) ---")
    h = Heasarc()
    tables = []
    t0 = time.time()

    for i, source in enumerate(source_data, 1):
        coord  = source["coord"]
        obj_id = source["objname"]
        elapsed = time.time() - t0
        print(f"  [{i}/{n}] {obj_id}  (elapsed {elapsed:.0f}s)", end="")
        try:
            result = h.query_region(coord, catalog=cat, radius=radius)
            if result is not None and len(result) > 0:
                print(f"  → {len(result)} {name} observation(s)")
                result.meta = {}
                present_cols = [c for c in cols if c in result.colnames]
                clean = result[present_cols].copy()
                clean["Input_Coord"] = _coord_string(coord)
                clean["Source_ID"]   = obj_id
                tables.append(clean)
            else:
                print("  → no matches")
        except Exception as exc:
            print(f"\n  ERROR querying HEASARC ({cat}) for {obj_id}: {exc}")

    if not tables:
        return None
    return vstack(tables)


def run_crossmatch():
    """
    Stage 1 – Load the input catalog, query all enabled archives, save CSVs.

    Returns
    -------
    dict  {archive_key: astropy.table.Table or None}
    """
    t_stage = time.time()
    print("\n" + "="*70)
    print("STAGE 1 : CROSS-MATCHING")
    print("="*70)

    # --- load input catalog ---
    df = pd.read_csv(INPUT_CATALOG_PATH, encoding="utf-8-sig")
    if CATALOG_SLICE is not None:
        lo, hi = CATALOG_SLICE
        df = df.iloc[lo:hi]
    print(f"Loaded {len(df)} sources from {INPUT_CATALOG_PATH}")

    # --- build source list ---
    source_data = []
    for _, row in df.iterrows():
        obj_id = str(row[INPUT_NAME_COL])
        if QUERY_PRIMARY_ONLY and obj_id[-1] != "A":
            continue
        coord = SkyCoord(row[INPUT_RA_COL], row[INPUT_DEC_COL], unit=(u.deg, u.deg))
        source_data.append({"coord": coord, "objname": obj_id})

    print(f"Querying {len(source_data)} primary sources with radius {SEARCH_RADIUS_ARCSEC} arcsec")
    radius = SEARCH_RADIUS_ARCSEC * u.arcsec

    results = {}

    # --- JWST ---
    if "jwst" in ACTIVE_ARCHIVES:
        jwst_table = query_jwst(source_data, radius)
        results["jwst"] = jwst_table
        save_table_to_csv(jwst_table, "jwst_agel_archive_matches.csv",
                          MAIN_DIR / "jwst")

    # --- HEASARC archives ---
    for key in _active_heasarc_archives():
        tbl = query_heasarc_archive(source_data, radius, key)
        results[key] = tbl
        save_table_to_csv(tbl, f"{key}_agel_archive_matches.csv",
                          MAIN_DIR / key)

    # --- summary ---
    elapsed = time.time() - t_stage
    print(f"\n--- Cross-Match Summary (total {elapsed:.0f}s) ---")
    for key, tbl in results.items():
        n = len(tbl) if tbl is not None else 0
        print(f"  {key.upper():12s}: {n} matched observations")

    return results


def load_crossmatch_results():
    """
    Load previously saved cross-match CSVs (used when skipping Stage 1).

    Returns
    -------
    dict  {archive_key: pandas.DataFrame or None}
    """
    results = {}
    if "jwst" in ACTIVE_ARCHIVES:
        p = MAIN_DIR / "jwst" / "jwst_agel_archive_matches.csv"
        results["jwst"] = pd.read_csv(p) if p.exists() else None
    for key in _active_heasarc_archives():
        p = MAIN_DIR / key / f"{key}_agel_archive_matches.csv"
        results[key] = pd.read_csv(p) if p.exists() else None
    return results


# ============================================================
# ===     STAGE 2 : FILTERING                             ===
# ============================================================

# JWST filter wavelengths (microns) – used for diverse-filter selection
_FILTER_WAVELENGTHS = {
    "F070W": 0.70, "F090W": 0.90, "F115W": 1.15, "F140M": 1.40,
    "F150W": 1.50, "F150W2": 1.50, "F162M": 1.62, "F164N": 1.64,
    "F182M": 1.82, "F187N": 1.87, "F200W": 2.00, "F210M": 2.10, "F212N": 2.12,
    "F250M": 2.50, "F277W": 2.77, "F300M": 3.00, "F322W2": 3.22, "F323N": 3.23,
    "F335M": 3.35, "F356W": 3.56, "F360M": 3.60, "F405N": 4.05, "F410M": 4.10,
    "F430M": 4.30, "F444W": 4.44, "F460M": 4.60, "F466N": 4.66, "F470N": 4.70, "F480M": 4.80,
    "F560W": 5.60, "F770W": 7.70, "F1000W": 10.0, "F1065C": 10.65, "F1130W": 11.30,
    "F1140C": 11.40, "F1280W": 12.80, "F1500W": 15.0, "F1550C": 15.50, "F1800W": 18.0,
    "F2100W": 21.0, "F2300C": 23.0, "F2550W": 25.50,
    "F158M": 1.58, "F380M": 3.80,
}


def filter_best_jwst_for_thumbnails(df, n_filters=3, instrument_priority=None,
                                     prefer_same_instrument=True):
    """
    Select up to n_filters imaging observations per source, prioritising
    diverse wavelength coverage for false-colour image creation.

    Calibration level -1 is excluded (bad calibration).
    """
    if instrument_priority is None:
        instrument_priority = INSTRUMENT_PRIORITY

    # Keep only imaging products; exclude NIRSPEC and bad calibration
    df = df[df["dataproduct_type"].isin(["image", "cube"])].copy()
    df = df[~df["instrument_name"].str.contains("NIRSPEC", case=False, na=False)].copy()
    df = df[df["calib_level"] != -1].copy()

    if len(df) == 0:
        return pd.DataFrame()

    df["base_instrument"] = df["instrument_name"].str.split("/").str[0].str.upper()

    def extract_filter(fstr):
        if pd.isna(fstr):
            return None
        parts = str(fstr).replace(";", "-").replace("CLEAR-", "").replace("-CLEAR", "").split("-")
        for p in parts:
            if p.upper().startswith("F") and any(c.isdigit() for c in p):
                return p.upper()
        return str(fstr).upper()

    df["filter_name"] = df["filters"].apply(extract_filter)
    df["wavelength"]  = df["filter_name"].map(_FILTER_WAVELENGTHS)

    results = []
    for source_id in df["Source_ID"].unique():
        src = df[df["Source_ID"] == source_id].copy()

        # Score each observation
        src["score"] = 0.0
        src["score"] += src["calib_level"].map({3: 1000, 2: 500}).fillna(0)
        for idx, inst in enumerate(instrument_priority):
            src.loc[src["base_instrument"] == inst, "score"] += (len(instrument_priority) - idx) * 100
        src.loc[src["dataproduct_type"] == "image", "score"] += 50
        if src["t_exptime"].max() > 0:
            src["score"] += (src["t_exptime"] / src["t_exptime"].max()) * 20

        # Best observation per unique filter
        best_per_filter = (src.sort_values("score", ascending=False)
                             .groupby("filter_name", dropna=True)
                             .first()
                             .reset_index())

        if len(best_per_filter) == 0:
            continue

        best_per_filter = best_per_filter.sort_values("wavelength", ascending=True,
                                                       na_position="last")

        # Optionally restrict to a single instrument that has enough filters
        if prefer_same_instrument and len(best_per_filter) >= n_filters:
            for inst in instrument_priority:
                inst_df = best_per_filter[best_per_filter["base_instrument"] == inst]
                if len(inst_df) >= n_filters:
                    best_per_filter = inst_df
                    break

        n_avail = len(best_per_filter)
        if n_avail >= 3 and n_filters >= 3:
            indices = list(dict.fromkeys([0, n_avail // 2, n_avail - 1]))
            selected = best_per_filter.iloc[indices]
        elif n_avail >= 2 and n_filters >= 2:
            selected = best_per_filter.iloc[[0, n_avail - 1]]
        else:
            selected = best_per_filter.head(min(n_filters, n_avail))

        results.append(selected)

    if not results:
        return pd.DataFrame()

    result_df = pd.concat(results, ignore_index=True)
    result_df = result_df.drop(columns=[c for c in ["base_instrument", "filter_name",
                                                      "wavelength", "score"]
                                         if c in result_df.columns])
    return result_df


def run_filter(crossmatch_results):
    """
    Stage 2 – Apply observation filters and save filtered CSVs.

    Returns
    -------
    dict  {archive_key: pandas.DataFrame}
    """
    t_stage = time.time()
    print("\n" + "="*70)
    print("STAGE 2 : FILTERING")
    print("="*70)

    filtered = {}

    if "jwst" in ACTIVE_ARCHIVES and crossmatch_results.get("jwst") is not None:
        jwst_df = crossmatch_results["jwst"]
        if isinstance(jwst_df, Table):
            jwst_df = jwst_df.to_pandas()
        jwst_filtered = filter_best_jwst_for_thumbnails(
            jwst_df, n_filters=JWST_N_FILTERS, instrument_priority=INSTRUMENT_PRIORITY
        )
        print(f"\nJWST: {len(jwst_df)} -> {len(jwst_filtered)} observations "
              f"({jwst_filtered['Source_ID'].nunique()} sources, "
              f"{JWST_N_FILTERS}-filter selection)")
        filtered["jwst"]    = jwst_filtered
        filtered["jwst_raw"] = jwst_df
        jwst_filtered.to_csv(MAIN_DIR / "jwst" / "jwst_filtered.csv", index=False)

    for key in _active_heasarc_archives():
        tbl = crossmatch_results.get(key)
        if tbl is None:
            filtered[key] = None
            continue
        if isinstance(tbl, Table):
            tbl = tbl.to_pandas()
        # X-ray archives: keep all unique observations (no spatial filter needed)
        filtered[key] = tbl.drop_duplicates(
            subset=["Source_ID", HEASARC_ARCHIVE_CONFIGS[key]["obsid_col"]]
        )
        print(f"\n{key.upper()}: {len(tbl)} -> {len(filtered[key])} observations")
        filtered[key].to_csv(MAIN_DIR / key / f"{key}_filtered.csv", index=False)

    print(f"\nFiltering complete ({time.time() - t_stage:.0f}s)")
    return filtered


# ============================================================
# ===     STAGE 3 : DOWNLOADING                           ===
# ============================================================

def _download_file(url, output_path, timeout=120):
    """Download a single file; return (success, size_mb, error_msg)."""
    try:
        r = requests.get(url, timeout=timeout, stream=True)
        r.raise_for_status()
        downloaded = 0
        with open(output_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
                downloaded += len(chunk)
        return True, downloaded / (1024 * 1024), None
    except Exception as exc:
        return False, 0, str(exc)


def download_jwst_science_files(df, output_dir, max_workers=2, limit=None):
    """
    Download JWST science FITS files from MAST, skipping duplicates.
    Files are named  {Source_ID}_{original_filename}.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    existing_files     = list(output_dir.glob("*.fits"))
    existing_basenames = {f.name.split("_", 1)[1] if "_" in f.name else f.name
                          for f in existing_files}
    scheduled          = set()
    tasks              = []
    mapping            = []
    n_skipped_dup      = 0

    for _, row in df.iterrows():
        uri       = row.get("dataURL")
        if pd.isna(uri):
            continue
        source_id      = row["Source_ID"]
        orig_filename  = uri.split("/")[-1]
        new_filename   = f"{source_id}_{orig_filename}"
        out_path       = output_dir / new_filename

        if orig_filename in existing_basenames:
            existing_f = next((f for f in existing_files
                               if f.name.endswith(orig_filename)), None)
            if existing_f:
                mapping.append({"Source_ID": source_id,
                                 "base_filename": orig_filename,
                                 "actual_file": existing_f.name,
                                 "status": "existing"})
                continue

        if orig_filename in scheduled:
            n_skipped_dup += 1
            mapping.append({"Source_ID": source_id,
                             "base_filename": orig_filename,
                             "actual_file": f"(scheduled)_{orig_filename}",
                             "status": "duplicate_skipped"})
            continue

        if limit and len(tasks) >= limit:
            break

        scheduled.add(orig_filename)
        url = (f"https://mast.stsci.edu/api/v0.1/Download/file?uri={uri}"
               if uri.startswith("mast:") else uri)
        tasks.append((url, out_path, source_id, row.get("instrument_name", ""), orig_filename))
        mapping.append({"Source_ID": source_id,
                        "base_filename": orig_filename,
                        "actual_file": new_filename,
                        "status": "to_download"})

    if n_skipped_dup:
        print(f"  Skipped {n_skipped_dup} duplicate files (same mosaic, different source)")

    successful = failed = 0
    total_mb   = 0.0

    if tasks:
        print(f"  Downloading {len(tasks)} unique JWST files (workers={max_workers})…")
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            future_map = {ex.submit(_download_file, url, path): (url, path, src, inst, base)
                          for url, path, src, inst, base in tasks}
            for i, future in enumerate(as_completed(future_map), 1):
                url, path, src, inst, base = future_map[future]
                ok, mb, err = future.result()
                if ok:
                    successful += 1
                    total_mb   += mb
                    existing_basenames.add(base)
                    print(f"  [{i}/{len(tasks)}] OK  {path.name} ({mb:.1f} MB) [{src}]")
                else:
                    failed += 1
                    print(f"  [{i}/{len(tasks)}] ERR {path.name} – {err}")
    else:
        print("  All unique JWST files already downloaded.")

    map_path = output_dir / "source_file_mapping.csv"
    pd.DataFrame(mapping).to_csv(map_path, index=False)
    return {"successful": successful, "failed": failed, "total_mb": total_mb}


def create_heasarc_download_script(df, archive_key, output_dir):
    """
    Write a shell script to download data for a HEASARC archive.
    """
    cfg       = HEASARC_ARCHIVE_CONFIGS[archive_key]
    name      = cfg["display_name"]
    obsid_col = cfg["obsid_col"]
    url_tpl   = cfg["download_url_template"]

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    script_path   = output_dir / f"download_{archive_key}_files.sh"
    obsid_list    = output_dir / f"{archive_key}_obsids.txt"

    with open(script_path, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"# {name} data download script\n")
        f.write("# Generated by archive_pipeline.py\n\n")
        for _, row in df.iterrows():
            source_id = row["Source_ID"]
            obsid     = row[obsid_col]
            out_file  = f"{source_id}_{archive_key}_{obsid}.tar"
            f.write(f"# Source: {source_id} (ObsID: {obsid})\n")
            f.write(f"echo 'Downloading {source_id}...'\n")
            if url_tpl:
                url = url_tpl.format(obsid=obsid)
                f.write(f"wget -O {out_file} '{url}'\n")
            else:
                f.write(f"# TODO: customise download URL for {name} obsid={obsid}\n")
            f.write("\n")

    script_path.chmod(0o755)

    with open(obsid_list, "w") as f:
        f.write("Source_ID\tObsID\tCoordinates\n")
        for _, row in df.iterrows():
            f.write(f"{row['Source_ID']}\t{row[obsid_col]}\t{row.get('Input_Coord', '')}\n")

    print(f"  {name} download script -> {script_path}")
    print(f"  ObsID list            -> {obsid_list}")
    return script_path


def run_download(filtered):
    """
    Stage 3 – Download JWST FITS files and generate download scripts for
    HEASARC archives.
    """
    t_stage = time.time()
    print("\n" + "="*70)
    print("STAGE 3 : DOWNLOADING")
    print("="*70)

    if "jwst" in ACTIVE_ARCHIVES and filtered.get("jwst") is not None:
        jwst_filt = filtered["jwst"]
        if len(jwst_filt) == 0:
            print("No JWST observations to download.")
        else:
            print(f"\nJWST: downloading up to {TEST_LIMIT or 'all'} files…")
            df_to_dl = jwst_filt.head(TEST_LIMIT) if TEST_LIMIT else jwst_filt
            stats = download_jwst_science_files(
                df_to_dl, JWST_DATA_DIR,
                max_workers=MAX_WORKERS, limit=TEST_LIMIT
            )
            print(f"  Success: {stats['successful']}  "
                  f"Failed: {stats['failed']}  "
                  f"Total: {stats['total_mb']:.1f} MB")

    for key in _active_heasarc_archives():
        df = filtered.get(key)
        if df is None or len(df) == 0:
            continue
        cfg = HEASARC_ARCHIVE_CONFIGS[key]
        print(f"\n{cfg['display_name']}: creating download script for {len(df)} observations…")
        out_dir = MAIN_DIR / key / "preview_files"
        df_to_script = df.head(TEST_LIMIT) if TEST_LIMIT else df
        create_heasarc_download_script(df_to_script, key, out_dir)

    print(f"\nDownloads complete ({time.time() - t_stage:.0f}s)")


# ============================================================
# ===     STAGE 4 : POSTAGE STAMPS                        ===
# ============================================================

def get_thumb_size(source_id):
    if source_id in CUSTOM_THUMB_SIZES:
        return CUSTOM_THUMB_SIZES[source_id]
    return DEFAULT_THUMB_SIZE_ARCSEC


def get_stretch(source_id):
    return CUSTOM_STRETCHES.get(source_id, DEFAULT_STRETCH)


def get_rgb_limits(source_id):
    defaults = dict(vmin_r=None, vmax_r=None, vmin_g=None, vmax_g=None,
                    vmin_b=None, vmax_b=None,
                    pmin_r=None, pmax_r=None, pmin_g=None, pmax_g=None,
                    pmin_b=None, pmax_b=None)
    if source_id in CUSTOM_RGB_LIMITS:
        defaults.update(CUSTOM_RGB_LIMITS[source_id])
    return defaults


def get_pixel_scale_arcsec(wcs):
    """Robustly determine pixel scale from WCS (handles CD matrix and CDELT)."""
    try:
        cd = wcs.wcs.cd
        if cd is not None:
            ps = abs(np.sqrt(cd[0, 0]**2 + cd[1, 0]**2)) * 3600
            if 0.001 < ps < 100:
                return ps
    except AttributeError:
        pass
    try:
        cdelt = wcs.wcs.get_cdelt()
        if cdelt is not None and len(cdelt) >= 2 and abs(cdelt[0]) not in (1.0, 0.0):
            ps = abs(cdelt[0]) * 3600
            if 0.001 < ps < 100:
                return ps
    except Exception:
        pass
    try:
        from astropy.wcs.utils import proj_plane_pixel_scales
        ps = proj_plane_pixel_scales(wcs)[0] * 3600
        if 0.001 < ps < 100:
            return ps
    except Exception:
        pass
    return 0.031   # fallback: NIRCam short-wavelength default


def get_cutout(fits_path, ra, dec, size_arcsec, ext=None, verbose=False):
    """Extract a WCS-aware cutout from a FITS file."""
    with fits.open(fits_path, memmap=False) as hdul:
        if ext is None:
            sci_exts   = []
            image_exts = []
            for i, hdu in enumerate(hdul):
                name = (hdu.name or "").upper()
                if name.startswith("SCI"):
                    sci_exts.append(i)
                naxis = hdu.header.get("NAXIS", 0)
                if naxis == 2 and hdu.header.get("NAXIS1") and hdu.header.get("NAXIS2"):
                    image_exts.append(i)
            candidates = ([i for i in sci_exts if i in image_exts] or sci_exts
                          or image_exts or [1 if len(hdul) > 1 else 0])
            data = wcs_obj = None
            for cand in candidates:
                try:
                    data    = hdul[cand].data
                    wcs_obj = WCS(hdul[cand].header)
                    break
                except Exception:
                    continue
            if data is None:
                return None, None
        else:
            try:
                data    = hdul[ext].data
                wcs_obj = WCS(hdul[ext].header)
            except Exception as exc:
                print(f"  Error reading ext [{ext}] in {fits_path}: {exc}")
                return None, None

    position    = SkyCoord(ra, dec, unit="deg", frame="icrs")
    pixel_scale = get_pixel_scale_arcsec(wcs_obj)
    size_pixels = max(int(np.ceil(size_arcsec / pixel_scale)), 10)

    try:
        cutout = Cutout2D(data, position, size_pixels, wcs=wcs_obj)
        return cutout.data, cutout.wcs
    except Exception as exc:
        if verbose:
            print(f"  Warning: cutout failed: {exc}")
        return None, None


def get_source_info(source_id, parent_df):
    """
    Look up RA, Dec, and redshifts for a source from the parent catalog.
    Tries exact match, then match without trailing letter, then match with 'A'.
    """
    for query in [source_id,
                  source_id[:-1] if source_id[-1] in "ABCD" else None,
                  source_id + "A" if source_id[-1] not in "ABCD" else None]:
        if query is None:
            continue
        match = parent_df[parent_df[PARENT_NAME_COL] == query]
        if len(match) > 0:
            row = match.iloc[0]
            z_s = row.get(PARENT_Z_SOURCE_COL, None)
            z_d = row.get(PARENT_Z_DEFLECTOR_COL, None)
            return {
                "ra":          row[PARENT_RA_COL],
                "dec":         row[PARENT_DEC_COL],
                "z_source":    None if pd.isna(z_s) else z_s,
                "z_deflector": None if pd.isna(z_d) else z_d,
            }
    return None


def format_redshift_title(source_id, z_source=None, z_deflector=None):
    title  = source_id
    z_parts = []
    if z_deflector is not None:
        z_parts.append(f"$z_{{\\rm d}}$={z_deflector:.3f}")
    if z_source is not None:
        z_parts.append(f"$z_{{\\rm s}}$={z_source:.3f}")
    if z_parts:
        title += "\n" + ", ".join(z_parts)
    return title


def add_scalebar(ax, wcs, size_arcsec, color="black", fontsize=12, loc="lower right"):
    """Draw a physical scalebar on a WCS axis."""
    pixel_scale    = get_pixel_scale_arcsec(wcs)
    scalebar_pix   = size_arcsec / pixel_scale
    xlim, ylim     = ax.get_xlim(), ax.get_ylim()
    x_start = xlim[1] - 0.15 * (xlim[1] - xlim[0]) if "right" in loc else xlim[0] + 0.05 * (xlim[1] - xlim[0])
    y_pos   = ylim[0] + 0.08 * (ylim[1] - ylim[0]) if "lower" in loc else ylim[1] - 0.08 * (ylim[1] - ylim[0])
    x_end   = x_start - scalebar_pix
    ax.plot([x_end, x_start], [y_pos, y_pos], color=color, linewidth=2)
    ax.text((x_start + x_end) / 2,
            y_pos + 0.03 * (ylim[1] - ylim[0]),
            f'{size_arcsec}"', color=color, fontsize=fontsize,
            ha="center", va="bottom", fontweight="medium")


def get_filter_from_filename(filename):
    """Extract filter name from a JWST filename."""
    name_lower = filename.lower()
    filters = [
        "f070w", "f090w", "f115w", "f140m", "f150w2", "f150w",
        "f162m", "f182m", "f200w", "f210m", "f212n",
        "f250m", "f277w", "f300m", "f322w2", "f335m",
        "f356w", "f360m", "f410m", "f430m", "f444w", "f460m", "f480m",
        "f560w", "f770w", "f1000w", "f1065c", "f1130w", "f1140c",
        "f1280w", "f1500w", "f1550c", "f1800w", "f2100w", "f2300c", "f2550w",
    ]
    for filt in filters:
        if filt in name_lower:
            return filt.upper()
    return "UNKNOWN"


def get_wavelength_order(filter_name):
    return _FILTER_WAVELENGTHS.get(filter_name.upper(), 99)


def _make_grayscale_postage(fits_path, output_dir, source_id, ra, dec,
                             filter_name, thumb_size_arcsec, scalebar_size_arcsec=2,
                             pmin=None, pmax=None, show_dual=True,
                             z_source=None, z_deflector=None):
    """Save a grayscale postage stamp (single or dual-panel)."""
    if pmin is None or pmax is None:
        pmin, pmax = get_stretch(source_id)

    cutout_data, cutout_wcs = get_cutout(fits_path, ra, dec, thumb_size_arcsec)
    if cutout_data is None:
        return None

    valid = cutout_data[np.isfinite(cutout_data)]
    if len(valid) == 0:
        print(f"  Warning: no valid pixels in cutout for {source_id}")
        return None

    title = format_redshift_title(source_id, z_source, z_deflector)

    if show_dual:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6),
                                        subplot_kw={"projection": cutout_wcs})
        norm1 = ImageNormalize(cutout_data, interval=AsymmetricPercentileInterval(pmin, pmax))
        ax1.imshow(cutout_data, origin="lower", cmap="gray_r", norm=norm1)
        ax1.set_title(f"{filter_name} p({pmin},{pmax})", fontsize=cbar_tick_size)
        ax1.coords.frame.set_color("none")
        ax1.tick_params(axis="x", bottom=False, labelbottom=False)
        ax1.tick_params(axis="y", left=False, labelleft=False)
        add_scalebar(ax1, cutout_wcs, scalebar_size_arcsec)

        pmax_tight = pmin + (pmax - pmin) / 2
        norm2 = ImageNormalize(cutout_data,
                               interval=AsymmetricPercentileInterval(pmin, pmax_tight))
        ax2.imshow(cutout_data, origin="lower", cmap="gray_r", norm=norm2)
        ax2.set_title(f"{filter_name} p({pmin},{pmax_tight:.1f})", fontsize=cbar_tick_size)
        ax2.coords.frame.set_color("none")
        ax2.tick_params(axis="x", bottom=False, labelbottom=False)
        ax2.tick_params(axis="y", left=False, labelleft=False)
        add_scalebar(ax2, cutout_wcs, scalebar_size_arcsec)

        fig.suptitle(title, fontsize=label_size)
        plt.subplots_adjust(wspace=0.05)
    else:
        fig, ax = plt.subplots(1, 1, figsize=(8, 8),
                               subplot_kw={"projection": cutout_wcs})
        norm = ImageNormalize(cutout_data, interval=AsymmetricPercentileInterval(pmin, pmax))
        ax.imshow(cutout_data, origin="lower", cmap="gray_r", norm=norm)
        ax.set_title(title, fontsize=label_size)
        ax.coords.frame.set_color("none")
        ax.tick_params(axis="x", bottom=False, labelbottom=False)
        ax.tick_params(axis="y", left=False, labelleft=False)
        ax.text(0.02, 0.95, filter_name, color="k", fontsize=cbar_tick_size,
                ha="left", va="top", transform=ax.transAxes,
                bbox=dict(facecolor="w", alpha=0.8, edgecolor="none"))
        add_scalebar(ax, cutout_wcs, scalebar_size_arcsec)

    out = Path(output_dir) / f"{source_id}_{filter_name}_{thumb_size_arcsec}arcs.png"
    plt.savefig(out, bbox_inches="tight", dpi=150, facecolor="white")
    plt.close(fig)
    return out


def _reproject_to_reference(data, wcs_in, wcs_ref, shape_ref):
    """Reproject data to the reference WCS grid."""
    reprojected, _ = reproject_interp((data, wcs_in), wcs_ref,
                                       shape_out=shape_ref, order="bilinear")
    reprojected[~np.isfinite(reprojected)] = 0
    return reprojected


def _normalize_channel(data, vmin=None, vmax=None, pmin_ch=None, pmax_ch=None,
                        default_pmin=1, default_pmax=99):
    valid = data[np.isfinite(data)]
    if len(valid) == 0:
        return np.zeros_like(data)
    if vmin is None or vmax is None:
        p_lo  = pmin_ch if pmin_ch is not None else default_pmin
        p_hi  = pmax_ch if pmax_ch is not None else default_pmax
        vmin, vmax = np.percentile(valid, [p_lo, p_hi])
    if vmax == vmin:
        return np.zeros_like(data)
    out = np.clip((data - vmin) / (vmax - vmin), 0, 1)
    out[~np.isfinite(out)] = 0
    return out


def _make_rgb_postage(red_fits, green_fits, blue_fits, output_dir, source_id,
                       ra, dec, red_filter, green_filter, blue_filter,
                       thumb_size_arcsec, scalebar_size_arcsec=2,
                       pmin=1, pmax=99,
                       vmin_r=None, vmax_r=None, vmin_g=None, vmax_g=None,
                       vmin_b=None, vmax_b=None,
                       pmin_r=None, pmax_r=None, pmin_g=None, pmax_g=None,
                       pmin_b=None, pmax_b=None,
                       z_source=None, z_deflector=None):
    """Save an RGB composite postage stamp with WCS reprojection alignment."""
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    rd, rw = get_cutout(red_fits,   ra, dec, thumb_size_arcsec)
    gd, gw = get_cutout(green_fits, ra, dec, thumb_size_arcsec)
    bd, bw = get_cutout(blue_fits,  ra, dec, thumb_size_arcsec)

    if rd is None or gd is None or bd is None:
        print(f"  Error: cutout failed for {source_id}")
        return None

    # Save green channel FITS for reference
    if gw is not None:
        green_out = TEMP_DIR / f"{source_id}_{green_filter}_green_channel.fits"
        hdr = gw.to_header()
        hdr["OBJECT"] = source_id
        hdr["FILTER"] = green_filter
        fits.PrimaryHDU(data=gd, header=hdr).writeto(green_out, overwrite=True)

    # Align to blue (highest resolution) reference
    ref_wcs, ref_shape = bw, bd.shape
    needs_reproject    = (rd.shape != ref_shape or gd.shape != ref_shape)
    if not needs_reproject:
        needs_reproject = (abs(get_pixel_scale_arcsec(rw) - get_pixel_scale_arcsec(bw))
                           / get_pixel_scale_arcsec(bw)) > 0.01

    if needs_reproject and REPROJECT_AVAILABLE:
        print(f"  Reprojecting channels to align pixel scales (ref: {blue_filter} {ref_shape})")
        rd = _reproject_to_reference(rd, rw, ref_wcs, ref_shape)
        gd = _reproject_to_reference(gd, gw, ref_wcs, ref_shape)
    elif needs_reproject:
        print(f"  Warning: reproject not available; channels may be misaligned.")

    r = _normalize_channel(rd, vmin=vmin_r, vmax=vmax_r, pmin_ch=pmin_r, pmax_ch=pmax_r,
                            default_pmin=pmin, default_pmax=pmax)
    g = _normalize_channel(gd, vmin=vmin_g, vmax=vmax_g, pmin_ch=pmin_g, pmax_ch=pmax_g,
                            default_pmin=pmin, default_pmax=pmax)
    b = _normalize_channel(bd, vmin=vmin_b, vmax=vmax_b, pmin_ch=pmin_b, pmax_ch=pmax_b,
                            default_pmin=pmin, default_pmax=pmax)
    rgb   = np.dstack([r, g, b])
    title = format_redshift_title(source_id, z_source, z_deflector)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={"projection": ref_wcs})
    ax.imshow(rgb, origin="lower")
    ax.set_title(title, fontsize=label_size, color="white")
    ax.coords.frame.set_color("none")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.text(0.02, 0.98, red_filter,   color="red",   fontsize=cbar_tick_size,
            ha="left", va="top", transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    ax.text(0.02, 0.93, green_filter, color="green", fontsize=cbar_tick_size,
            ha="left", va="top", transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    ax.text(0.02, 0.88, blue_filter,  color="blue",  fontsize=cbar_tick_size,
            ha="left", va="top", transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    add_scalebar(ax, ref_wcs, scalebar_size_arcsec, color="white")

    out = (Path(output_dir)
           / f"{source_id}_{red_filter}_{green_filter}_{blue_filter}_{thumb_size_arcsec}arcs_RGB.png")
    plt.savefig(out, bbox_inches="tight", dpi=150, facecolor="black")
    plt.close(fig)
    return out


def _make_2color_postage(red_fits, blue_fits, output_dir, source_id,
                          ra, dec, red_filter, blue_filter,
                          thumb_size_arcsec, scalebar_size_arcsec=2,
                          pmin=1, pmax=99,
                          vmin_r=None, vmax_r=None, vmin_g=None, vmax_g=None,
                          vmin_b=None, vmax_b=None,
                          pmin_r=None, pmax_r=None, pmin_g=None, pmax_g=None,
                          pmin_b=None, pmax_b=None,
                          z_source=None, z_deflector=None):
    """Save a 2-colour composite (green = average of red + blue)."""
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    rd, rw = get_cutout(red_fits,  ra, dec, thumb_size_arcsec)
    bd, bw = get_cutout(blue_fits, ra, dec, thumb_size_arcsec)

    if rd is None or bd is None:
        print(f"  Error: cutout failed for {source_id}")
        return None

    ref_wcs, ref_shape = bw, bd.shape
    needs_reproject    = rd.shape != ref_shape
    if not needs_reproject:
        needs_reproject = (abs(get_pixel_scale_arcsec(rw) - get_pixel_scale_arcsec(bw))
                           / get_pixel_scale_arcsec(bw)) > 0.01

    if needs_reproject and REPROJECT_AVAILABLE:
        rd = _reproject_to_reference(rd, rw, ref_wcs, ref_shape)
    elif needs_reproject:
        print("  Warning: reproject not available; channels may be misaligned.")

    gd = (rd + bd) / 2.0

    # Save average green channel
    if ref_wcs is not None:
        green_out = TEMP_DIR / f"{source_id}_{red_filter}_{blue_filter}_green_avg.fits"
        hdr = ref_wcs.to_header()
        fits.PrimaryHDU(data=gd, header=hdr).writeto(green_out, overwrite=True)

    r = _normalize_channel(rd, vmin=vmin_r, vmax=vmax_r, pmin_ch=pmin_r, pmax_ch=pmax_r,
                            default_pmin=pmin, default_pmax=pmax)
    g = _normalize_channel(gd, vmin=vmin_g, vmax=vmax_g, pmin_ch=pmin_g, pmax_ch=pmax_g,
                            default_pmin=pmin, default_pmax=pmax)
    b = _normalize_channel(bd, vmin=vmin_b, vmax=vmax_b, pmin_ch=pmin_b, pmax_ch=pmax_b,
                            default_pmin=pmin, default_pmax=pmax)
    rgb   = np.dstack([r, g, b])
    title = format_redshift_title(source_id, z_source, z_deflector)

    fig, ax = plt.subplots(1, 1, figsize=(8, 8), subplot_kw={"projection": ref_wcs})
    ax.imshow(rgb, origin="lower")
    ax.set_title(title, fontsize=label_size, color="white")
    ax.coords.frame.set_color("none")
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", left=False, labelleft=False)
    ax.text(0.02, 0.98, red_filter,  color="red",  fontsize=cbar_tick_size,
            ha="left", va="top", transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    ax.text(0.02, 0.93, blue_filter, color="blue", fontsize=cbar_tick_size,
            ha="left", va="top", transform=ax.transAxes,
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none"))
    add_scalebar(ax, ref_wcs, scalebar_size_arcsec, color="white")

    out = (Path(output_dir)
           / f"{source_id}_{red_filter}_{blue_filter}_{thumb_size_arcsec}arcs_RGB.png")
    plt.savefig(out, bbox_inches="tight", dpi=150, facecolor="black")
    plt.close(fig)
    return out


def run_postage_stamps():
    """
    Stage 4 – Generate JWST postage stamp images from downloaded FITS files.

    Looks for FITS files in JWST_DATA_DIR, reads coordinates from
    PARENT_CATALOG_PATH, and writes PNG stamps to STAMP_OUTPUT_DIR.
    """
    t_stage = time.time()
    print("\n" + "="*70)
    print("STAGE 4 : POSTAGE STAMPS")
    print("="*70)

    # --- discover downloaded FITS files ---
    fits_files = sorted(JWST_DATA_DIR.glob("*.fits"))
    if not fits_files:
        print("No FITS files found in", JWST_DATA_DIR)
        return

    print(f"Found {len(fits_files)} FITS files in {JWST_DATA_DIR}")

    # Group by source ID and sort filters by wavelength
    files_by_source = defaultdict(list)
    for f in fits_files:
        source_id = f.name.split("_jw")[0]
        files_by_source[source_id].append(f)

    source_info = {}
    for source_id, files in files_by_source.items():
        filter_info = [(get_filter_from_filename(f.name), f) for f in files]
        filter_info.sort(key=lambda x: get_wavelength_order(x[0]))
        source_info[source_id] = filter_info

    # --- load parent catalog ---
    parent_df = None
    if PARENT_CATALOG_PATH and Path(PARENT_CATALOG_PATH).exists():
        parent_df = pd.read_csv(PARENT_CATALOG_PATH, encoding="utf-8-sig")
        parent_df.columns = [c.replace("\ufeff", "") for c in parent_df.columns]
        print(f"Loaded parent catalog: {len(parent_df)} entries")
    else:
        print("WARNING: Parent catalog not found. Using nominal coordinates from filenames.")

    # --- output directories ---
    STAMP_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)

    created_files   = []
    failed_sources  = []
    sorted_sources  = sorted(source_info.items())
    n_sources       = len(sorted_sources)

    for i, (source_id, filter_info) in enumerate(sorted_sources, 1):
        elapsed = time.time() - t_stage
        print(f"\n[{i}/{n_sources}] {source_id}  (elapsed {elapsed:.0f}s)…")

        # Resolve coordinates
        ra = dec = z_source = z_deflector = None
        if parent_df is not None:
            info = get_source_info(source_id, parent_df)
            if info:
                ra, dec = info["ra"], info["dec"]
                z_source, z_deflector = info["z_source"], info["z_deflector"]

        if ra is None:
            print(f"  Warning: no coordinates for {source_id} – skipping")
            failed_sources.append((source_id, "No coordinates"))
            continue

        thumb_size  = get_thumb_size(source_id)
        pmin, pmax  = get_stretch(source_id)
        rgb_limits  = get_rgb_limits(source_id)
        n_filters   = len(filter_info)

        print(f"  RA={ra:.5f}  Dec={dec:.5f}  thumb={thumb_size}\"  "
              f"p({pmin},{pmax})  n_filters={n_filters}")

        if n_filters == 1:
            flt, fp = filter_info[0]
            result  = _make_grayscale_postage(
                fp, STAMP_OUTPUT_DIR, source_id, ra, dec, flt,
                thumb_size_arcsec=thumb_size,
                scalebar_size_arcsec=DEFAULT_SCALEBAR_SIZE_ARCSEC,
                pmin=pmin, pmax=pmax, show_dual=True,
                z_source=z_source, z_deflector=z_deflector,
            )

        elif n_filters == 2:
            blue_flt, blue_fp = filter_info[0]   # shorter wavelength = blue
            red_flt,  red_fp  = filter_info[1]
            result = _make_2color_postage(
                red_fp, blue_fp, STAMP_OUTPUT_DIR, source_id, ra, dec,
                red_flt, blue_flt,
                thumb_size_arcsec=thumb_size,
                scalebar_size_arcsec=DEFAULT_SCALEBAR_SIZE_ARCSEC,
                pmin=pmin, pmax=pmax,
                **rgb_limits,
                z_source=z_source, z_deflector=z_deflector,
            )

        else:   # 3+ filters: RGB
            blue_flt,  blue_fp  = filter_info[0]
            green_flt, green_fp = filter_info[n_filters // 2]
            red_flt,   red_fp   = filter_info[-1]
            result = _make_rgb_postage(
                red_fp, green_fp, blue_fp, STAMP_OUTPUT_DIR, source_id, ra, dec,
                red_flt, green_flt, blue_flt,
                thumb_size_arcsec=thumb_size,
                scalebar_size_arcsec=DEFAULT_SCALEBAR_SIZE_ARCSEC,
                pmin=pmin, pmax=pmax,
                **rgb_limits,
                z_source=z_source, z_deflector=z_deflector,
            )

        if result:
            created_files.append(result)
            print(f"  -> Saved: {result.name}")
        else:
            failed_sources.append((source_id, "image generation failed"))

        # Optional: individual per-filter grayscale stamps
        if CREATE_INDIVIDUAL_GRAYSCALE:
            indiv_dir = STAMP_OUTPUT_DIR / "individual_filters"
            indiv_dir.mkdir(parents=True, exist_ok=True)
            for flt, fp in filter_info:
                r = _make_grayscale_postage(
                    fp, indiv_dir, source_id, ra, dec, flt,
                    thumb_size_arcsec=thumb_size,
                    scalebar_size_arcsec=DEFAULT_SCALEBAR_SIZE_ARCSEC,
                    pmin=pmin, pmax=pmax, show_dual=False,
                    z_source=z_source, z_deflector=z_deflector,
                )
                if r:
                    created_files.append(r)

    print("\n" + "="*70)
    print(f"POSTAGE STAMP SUMMARY  ({time.time() - t_stage:.0f}s total)")
    print("="*70)
    print(f"Created : {len(created_files)} images")
    print(f"Failed  : {len(failed_sources)} sources")
    if failed_sources:
        for src, reason in failed_sources:
            print(f"  {src}: {reason}")
    print(f"Output  : {STAMP_OUTPUT_DIR}")


# ============================================================
# ===                    MAIN                              ===
# ============================================================

_ALL_ARCHIVES = ["jwst", "chandra", "xrism", "xmm", "erosita", "erosmaster", "swift"]


def parse_args():
    p = argparse.ArgumentParser(description="AGEL Archive Cross-Match Pipeline")
    p.add_argument(
        "--archive", dest="archives", nargs="+", metavar="NAME",
        choices=_ALL_ARCHIVES,
        help=(f"One or more archives to process. Choices: {_ALL_ARCHIVES}. "
              f"Defaults to {DEFAULT_ARCHIVES} when omitted."),
    )
    p.add_argument("--skip-crossmatch", action="store_true",
                   help="Skip Stage 1; load existing cross-match CSVs instead")
    p.add_argument("--skip-filter",     action="store_true",
                   help="Skip Stage 2; load existing filtered CSVs")
    p.add_argument("--skip-download",   action="store_true",
                   help="Skip Stage 3 (useful when files are already downloaded)")
    p.add_argument("--skip-stamps",     action="store_true",
                   help="Skip Stage 4 postage stamp generation")
    p.add_argument("--only-stamps",     action="store_true",
                   help="Run Stage 4 only (shorthand for --skip-crossmatch "
                        "--skip-filter --skip-download)")
    p.add_argument("--clean",           action="store_true",
                   help="After Stage 4, delete full-FoV FITS files, keeping only "
                        "the postage stamp PNGs")
    return p.parse_args()


def run_clean():
    """
    Delete full field-of-view FITS files from JWST_DATA_DIR and TEMP_DIR,
    leaving only the postage-stamp PNGs in STAMP_OUTPUT_DIR.
    """
    print("\n" + "="*70)
    print("CLEAN: removing full-FoV FITS files")
    print("="*70)

    deleted = 0
    for directory in (JWST_DATA_DIR, TEMP_DIR):
        fits_files = list(directory.glob("*.fits"))
        for f in fits_files:
            f.unlink()
            print(f"  Deleted {f}")
            deleted += 1

    print(f"Removed {deleted} FITS file(s). Postage stamps retained in {STAMP_OUTPUT_DIR}")


def main():
    global ACTIVE_ARCHIVES
    args = parse_args()

    # Resolve which archives to use for this run
    ACTIVE_ARCHIVES = set(args.archives if args.archives else DEFAULT_ARCHIVES)
    print(f"Active archives: {sorted(ACTIVE_ARCHIVES)}")

    if args.only_stamps:
        args.skip_crossmatch = True
        args.skip_filter     = True
        args.skip_download   = True

    # ---- Stage 1: Cross-match ----
    crossmatch_results = None
    if not args.skip_crossmatch:
        crossmatch_results = run_crossmatch()
    else:
        print("\n[Skipping Stage 1 – loading existing cross-match CSVs]")
        crossmatch_results = load_crossmatch_results()

    # ---- Stage 2: Filter ----
    filtered = None
    if not args.skip_filter and crossmatch_results is not None:
        filtered = run_filter(crossmatch_results)
    elif args.skip_filter:
        print("\n[Skipping Stage 2 – loading existing filtered CSVs]")
        filtered = {}
        if "jwst" in ACTIVE_ARCHIVES:
            p = MAIN_DIR / "jwst" / "jwst_filtered.csv"
            filtered["jwst"] = pd.read_csv(p) if p.exists() else None
        for key in _active_heasarc_archives():
            p = MAIN_DIR / key / f"{key}_filtered.csv"
            filtered[key] = pd.read_csv(p) if p.exists() else None

    # ---- Stage 3: Download ----
    if not args.skip_download and filtered is not None:
        run_download(filtered)
    elif args.skip_download:
        print("\n[Skipping Stage 3 – assuming files already downloaded]")

    # ---- Stage 4: Postage stamps ----
    if not args.skip_stamps:
        run_postage_stamps()
    else:
        print("\n[Skipping Stage 4 – postage stamp generation]")

    # ---- Clean: remove full-FoV FITS ----
    if args.clean:
        run_clean()

    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print("="*70)
    print(f"All output in: {MAIN_DIR}")


if __name__ == "__main__":
    main()
