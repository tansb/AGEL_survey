"""Instrument constants and default masks.

Everything here is a *documented default* that the user can override from the
GUI / config file; nothing is silently hard-wired into the algorithms.
"""
from __future__ import annotations

# ---------------------------------------------------------------------------
# KCWI / KCRM resolving power by grating and slicer (Keck instrument pages).
# R = lambda / FWHM  (so the FWHM in Angstrom is lambda_central / R).
# These are *nominal* values; arc-lamp or sky-line measurements of your own
# reduction are better, and can be typed into the GUI directly.
# ---------------------------------------------------------------------------
GRATING_R_TABLE = {
    # blue arm (KCWI)
    "BL":  {"Large": 900,  "Medium": 1800, "Small": 3600},
    "BM":  {"Large": 2000, "Medium": 4000, "Small": 8000},
    "BH1": {"Large": 4500, "Medium": 9000, "Small": 18000},
    "BH2": {"Large": 4500, "Medium": 9000, "Small": 18000},
    "BH3": {"Large": 4500, "Medium": 9000, "Small": 18000},
    # red arm (KCRM)
    "RL":  {"Large": 500,  "Medium": 1000, "Small": 2000},
    "RM1": {"Large": 1400, "Medium": 2800, "Small": 5600},
    "RM2": {"Large": 1400, "Medium": 2800, "Small": 5600},
    "RH1": {"Large": 3250, "Medium": 6500, "Small": 13000},
    "RH2": {"Large": 3250, "Medium": 6500, "Small": 13000},
    "RH3": {"Large": 3250, "Medium": 6500, "Small": 13000},
    "RH4": {"Large": 3250, "Medium": 6500, "Small": 13000},
}
BLUE_GRATINGS = ("BL", "BM", "BH1", "BH2", "BH3")
RED_GRATINGS = ("RL", "RM1", "RM2", "RH1", "RH2", "RH3", "RH4")
SLICER_NAMES = ("Large", "Medium", "Small")
DEFAULT_GRATING = {"blue": "BL", "red": "RL"}
DEFAULT_SLICER = "Medium"

# FITS header keywords inspected (in this order) to identify the set-up.
GRATING_KEYS = {"blue": ("BGRATNAM", "GRATNAM", "GRATNAME", "GRATING"),
                "red": ("RGRATNAM", "GRATNAM", "GRATNAME", "GRATING")}
SLICER_KEYS = ("IFUNAM", "SLICER", "SLICENAM", "SLICE")

# Keck II (Maunakea) -- used only for the optional barycentric-correction helper.
KECK_LON_DEG, KECK_LAT_DEG, KECK_HEIGHT_M = -155.4747, 19.8260, 4160.0

# ---------------------------------------------------------------------------
# Spectral regions that are NOT astrophysical signal at any redshift.
# Observed-frame, air wavelengths [Angstrom].  Grouped so the GUI can toggle
# them as physical classes rather than as an opaque list.
# ---------------------------------------------------------------------------
MASK_GROUPS = {
    # Natural night-sky emission lines: sky-subtraction residuals, not source signal.
    "airglow": [
        (5570.0, 5585.0),   # [O I] 5577.34
        (5885.0, 5900.0),   # Na I D airglow 5890/5896
        (6295.0, 6310.0),   # [O I] 6300.30
        (6355.0, 6372.0),   # [O I] 6363.78
    ],
    # The red/NIR OH Meinel-band forest: too crowded to catalogue line by line, and it
    # sits at the very red edge of the KCRM good range (WAVGOOD1 ~ 8940 A).
    "oh_forest": [
        (8750.0, 8950.0),
    ],
    # Artificial light-pollution lines.  Site dependent (weak on Maunakea) -- a
    # precaution rather than a certainty, hence a separate switch.
    "light_pollution": [
        (4043.0, 4051.0),   # Hg I 4046.6
        (4353.0, 4363.0),   # Hg I 4358.3
        (5456.0, 5466.0),   # Hg I 5460.7
        (5679.0, 5692.0),   # Na I / Hg I blend
        (5765.0, 5796.0),   # Hg I 5769.6 / 5790.7
    ],
    # Telluric ABSORPTION bands (no telluric correction is applied by this pipeline).
    "telluric": [
        (6860.0, 6890.0),   # O2 B band
        (7150.0, 7350.0),   # H2O
        (7580.0, 7700.0),   # O2 A band
    ],
}

# ---------------------------------------------------------------------------
# Template-library metadata
# ---------------------------------------------------------------------------
# pPXF's authors treat E-MILES as air-native and the other SPS libraries as
# vacuum-native.  Comparing a galaxy in the wrong frame injects ~ -85..-95 km/s.
SPS_NATIVE_FRAME = {"emiles": "air", "fsps": "vacuum", "galaxev": "vacuum", "xsl": "vacuum"}
SPS_LIBRARIES = ("emiles", "fsps", "galaxev", "xsl")
SPS_DOWNLOAD_BASE = "https://raw.githubusercontent.com/micappe/ppxf_data/main/"

KC96_TEMPLATES = {
    "kc96_starb1": "starb1_template.fits", "kc96_starb2": "starb2_template.fits",
    "kc96_starb3": "starb3_template.fits", "kc96_starb4": "starb4_template.fits",
    "kc96_starb5": "starb5_template.fits", "kc96_starb6": "starb6_template.fits",
    "kc96_elliptical": "elliptical_template.fits", "kc96_bulge": "bulge_template.fits",
    "kc96_s0": "s0_template.fits", "kc96_sa": "sa_template.fits",
    "kc96_sb": "sb_template.fits", "kc96_sc": "sc_template.fits",
}
KC96_LIBRARIES = {
    "kc96_starburst": [f"kc96_starb{i}" for i in range(1, 7)],
    "kc96_galaxy": ["kc96_elliptical", "kc96_bulge", "kc96_s0", "kc96_sa", "kc96_sb", "kc96_sc"],
    "kc96_all": list(KC96_TEMPLATES),
}
# The atlas README states no line-spread function; its native sampling is 5 A/pix,
# so claiming anything finer than 5 A FWHM would be physically implausible.
KC96_ASSUMED_FWHM_A = 5.0
LIBRARY_CHOICES = list(SPS_LIBRARIES) + list(KC96_LIBRARIES) + ["custom"]
