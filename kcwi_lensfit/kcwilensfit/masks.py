"""Wavelength masks (observed frame, AIR) assembled from physical classes of contaminants."""
from __future__ import annotations

from typing import Dict, Iterable, List, Sequence, Tuple

import numpy as np

from . import constants as K
from .lines import LINES, LINES_BY_NAME, observed_wavelengths_air
from .schema import parse_windows
from .wavelength import C_KMS

Window = Tuple[float, float]

# Forbidden / nebular lines that stellar templates cannot reproduce (vacuum rest wavelengths).
FORBIDDEN_NEBULAR = ["[O II] 3727", "[Ne III] 3869", "[O III] 5007", "[N II] 6585", "[S II] 6718", "[O I] 6302"]
BALMER = ["H8 3890", "H-delta", "H-gamma", "H-beta", "H-alpha"]
# emission lines masked in ONE aperture's fit because another plane's source emits them
STRONG_EMISSION = ["Ly-alpha", "C IV 1549 (em)", "He II 1640", "O III] 1663", "C III] 1908",
                   "Mg II 2799 (em)", "[O II] 3727", "[Ne III] 3869", "H-delta", "H-gamma", "H-beta",
                   "[O III] 5007", "[O I] 6302", "H-alpha", "[N II] 6585", "[S II] 6718"]


def default_windows(cfg: dict) -> List[Window]:
    """Sky / telluric windows selected by the config's mask_* switches."""
    out: List[Window] = []
    for group in ("airglow", "oh_forest", "light_pollution", "telluric"):
        if cfg.get(f"mask_{group}", True):
            out += list(K.MASK_GROUPS[group])
    out += parse_windows(cfg.get("extra_mask_windows", ""))
    return out


def line_windows(names: Iterable[str], z: float, half_kms: float, wmin: float, wmax: float) -> List[Window]:
    """Observed-air windows ``+- half_kms`` around every component of the named lines at redshift z."""
    out: List[Window] = []
    for n in names:
        line = LINES_BY_NAME.get(n)
        if line is None:
            continue
        for lam in observed_wavelengths_air(line, z):
            h = lam * half_kms / C_KMS
            if lam + h >= wmin and lam - h <= wmax:
                out.append((float(lam - h), float(lam + h)))
    return out


def merge_windows(windows: Sequence[Window]) -> List[Window]:
    if not windows:
        return []
    w = sorted((min(a, b), max(a, b)) for a, b in windows)
    out = [list(w[0])]
    for lo, hi in w[1:]:
        if lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return [tuple(x) for x in out]


def in_windows(wave: np.ndarray, windows: Sequence[Window]) -> np.ndarray:
    m = np.zeros(wave.shape, bool)
    for lo, hi in windows:
        m |= (wave >= lo) & (wave <= hi)
    return m
