"""Template-library management for pPXF.

Libraries
---------
``emiles / fsps / galaxev / xsl``
    Stellar-population-synthesis grids shipped with the ``ppxf`` data repository
    (Vazdekis+2016; Conroy+2009; Bruzual & Charlot 2003; Verro+2022).  Looked for in the
    installed ppxf package, then in the local cache, and downloaded on first use.
``kc96_*``
    Empirical *galaxy* spectra (Kinney+1996 / Calzetti+1994): starburst, spiral, elliptical
    composites.  These ship with this package.  Their resolution (~5 A) is coarser than KCWI, so
    they are for **redshifts of blue, star-forming sources**, not for velocity dispersions.
``custom``
    Any pPXF ``sps_lib``-format ``.npz`` file.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from importlib import resources
from typing import List, Optional, Tuple
from urllib import request

import numpy as np
from astropy.io import fits

from . import constants as K


def cache_dir() -> str:
    d = os.environ.get("KCWILENSFIT_CACHE", os.path.join(os.path.expanduser("~"), ".cache", "kcwilensfit"))
    os.makedirs(d, exist_ok=True)
    return d


@dataclass
class Library:
    name: str
    path: str
    native_frame: str       # 'air' | 'vacuum'
    is_empirical: bool = False   # True for KC96 (not an SPS age/metallicity grid)

    def metadata(self) -> dict:
        d = np.load(self.path)
        return dict(lam_min=float(d["lam"].min()), lam_max=float(d["lam"].max()),
                    fwhm_median=float(np.median(d["fwhm"])), n_age=len(d["ages"]), n_met=len(d["metals"]))


def sps_file(name: str) -> str:
    """Locate (or download) ``spectra_<name>_9.0.npz``."""
    base = f"spectra_{name}_9.0.npz"
    try:
        import ppxf
        cand = os.path.join(os.path.dirname(ppxf.__file__), "sps_models", base)
        if os.path.isfile(cand):
            return cand
    except ImportError:
        pass
    local = os.path.join(cache_dir(), base)
    if os.path.isfile(local):
        return local
    url = K.SPS_DOWNLOAD_BASE + base
    print(f"[templates] '{base}' not found locally -- downloading from {url} ...", flush=True)
    try:
        request.urlretrieve(url, local + ".part")
        os.replace(local + ".part", local)
    except Exception as e:
        if os.path.exists(local + ".part"):
            os.remove(local + ".part")
        raise RuntimeError(f"Could not download {base} ({e}). Download it manually from {url} and put it "
                           f"in {cache_dir()} (or point KCWILENSFIT_CACHE at its folder).") from e
    return local


def _kc96_fits_path(name: str) -> str:
    fn = K.KC96_TEMPLATES[name]
    return str(resources.files("kcwilensfit") / "data" / "kc96" / fn)


def _read_kc96(name: str) -> Tuple[np.ndarray, np.ndarray]:
    with fits.open(_kc96_fits_path(name)) as h:
        d = h[1].data
        lam = np.asarray(d["WAVELENGTH"], float)
        fl = np.asarray(d["FLUX"], float)
    o = np.argsort(lam)
    return lam[o], fl[o]


def kc96_library_file(lib_name: str) -> str:
    """Build (once) an ``sps_lib``-compatible npz from one KC96 template or a set of them.

    A set becomes one linear-combination library: pPXF fits a non-negative mix, exactly as with a
    real SPS grid.  The fake 'age' axis only labels the templates -- it has no physical meaning.
    """
    names = K.KC96_LIBRARIES.get(lib_name, [lib_name])
    path = os.path.join(cache_dir(), f"{lib_name}_sps.npz")
    if os.path.isfile(path):
        return path
    waves, fluxes = zip(*[_read_kc96(n) for n in names])
    lo, hi = max(w.min() for w in waves), min(w.max() for w in waves)
    lam = waves[0][(waves[0] >= lo) & (waves[0] <= hi)]
    tem = np.zeros((lam.size, len(names), 1))
    for i, (w, f) in enumerate(zip(waves, fluxes)):
        f = np.interp(lam, w, f)
        tem[:, i, 0] = f / np.median(f[f > 0])        # comparable normalisation across templates
    np.savez(path, templates=tem, masses=np.ones((len(names), 1)), ages=np.arange(1, len(names) + 1, dtype=float),
             metals=np.array([0.0]), lam=lam, fwhm=np.full_like(lam, K.KC96_ASSUMED_FWHM_A))
    return path


def resolve_library(name: str, custom_path: str = "", custom_frame: str = "air") -> Library:
    if name == "custom":
        if not custom_path or not os.path.isfile(os.path.expanduser(custom_path)):
            raise FileNotFoundError("sps_library='custom' needs a valid custom_template_path (.npz)")
        return Library("custom", os.path.expanduser(custom_path), custom_frame)
    if name in K.SPS_LIBRARIES:
        return Library(name, sps_file(name), K.SPS_NATIVE_FRAME[name])
    if name in K.KC96_LIBRARIES or name in K.KC96_TEMPLATES:
        # ground-based optical part (the dominant range) of the atlas is air; the IUE-UV splice
        # differs by < 1 A (~0.2 pixel at 5 A sampling) -- sub-resolution.
        return Library(name, kc96_library_file(name), "air", is_empirical=True)
    raise ValueError(f"unknown template library {name!r}")


# UV rest-frame emission lines used as a fall-back gas list when pPXF's own (optical-only) list
# has nothing in the fitted rest range.
UV_GAS_LINES = {"Lya": 1215.67, "NV": 1240.14, "CIV": 1549.48, "HeII": 1640.4, "OIII]": 1666.15,
                "NIII]": 1750.0, "SiIII]": 1892.03, "CIII]": 1908.73, "MgII": 2799.12}
