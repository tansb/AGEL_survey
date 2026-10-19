"""Spectral-line catalogue used for redshift identification.

Design notes (the astronomy behind the table)
---------------------------------------------
* All rest wavelengths are stored in **vacuum**.  Optical lines are conventionally
  tabulated in air, so those are entered in air and converted at import time with
  the same dispersion formula used for the data (:mod:`kcwilensfit.wavelength`).
* Lines are grouped into *multiplets* (``Line``): components that share ONE amplitude with
  nominal relative strengths.  Atomic doublets with an exact branching ratio ([O III] 4960/5008,
  [N II] 6549/6585, [O I]) have ``fixed_ratio=True`` -- genuinely more constraining than
  independent lines.  Doublets whose ratio depends on density / optical depth ([O II], [S II],
  Mg II, C IV, Na D ...) carry a *nominal* ratio (``fixed_ratio=False``): the template is right
  about the wavelength SEPARATION, which is what distinguishes a doublet from a noise spike,
  while a strength mismatch only costs a little sensitivity.
* Every line carries a *kind* because different kinds trace the systemic redshift
  with different reliability:

  ===========  ==================================================  ==================
  kind         physical origin                                     systemic offset
  ===========  ==================================================  ==================
  neb          nebular (HII-region / AGN narrow-line) emission     ~0 (best)
  stel_abs     stellar photospheric absorption (Ca H&K, G, Mg b)   ~0 (best)
  ism_abs      interstellar / outflow absorption (Fe II, Mg II..)  blueshifted, -50..-300 km/s
  res_em       resonant emission (Lya, Mg II em., C IV em.)        redshifted/uncertain, +100..+500 km/s
  ===========  ==================================================  ==================

  The pipeline reports which kinds a redshift is built from, and attaches a
  systematic uncertainty when it relies on outflow-affected lines only.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Sequence, Tuple

import numpy as np

from .wavelength import air_to_vac

#: Typical offset of each line kind relative to the systemic redshift [km/s]
#: (value, 1-sigma scatter).  Used only to annotate results with a systematic
#: error; the fitted redshift itself is the *observed* line redshift.
KIND_OFFSET_KMS = {
    "neb": (0.0, 30.0),
    "stel_abs": (0.0, 30.0),
    "ism_abs": (-150.0, 150.0),
    "res_em": (200.0, 250.0),
}

#: kinds that are (and are not) emission features
EMISSION_KINDS = ("neb", "res_em")
ABSORPTION_KINDS = ("stel_abs", "ism_abs")


@dataclass(frozen=True)
class Line:
    name: str
    comps: Tuple[Tuple[float, float], ...]  # ((vacuum rest wavelength, relative strength), ...)
    kind: str
    fixed_ratio: bool = True   # True: one shared amplitude with the given relative strengths
    note: str = ""

    @property
    def is_emission(self) -> bool:
        return self.kind in EMISSION_KINDS

    @property
    def lam_ref(self) -> float:
        """Strength-weighted mean vacuum wavelength."""
        w = np.array([c[1] for c in self.comps])
        lam = np.array([c[0] for c in self.comps])
        return float(np.sum(w * lam) / np.sum(w))


def _air(lam):  # optical lines: entered in air, stored in vacuum
    return float(air_to_vac(lam))


def _build() -> List[Line]:
    L: List[Line] = []
    add = lambda *a, **k: L.append(Line(*a, **k))

    # ---------------- emission: nebular ----------------
    add("[O II] 3727", ((3727.092, 1.0), (3729.875, 1.25)), "neb", fixed_ratio=False,
        note="3727+3729 doublet; ratio density-dependent (nominal 1:1.25)")
    add("[Ne III] 3869", ((3869.86, 1.0), (3968.59, 0.30)), "neb", note="3869 + 3968 (fixed 1:0.30)")
    add("H8 3890", ((3890.17, 1.0),), "neb")
    add("H-delta", ((4102.89, 1.0),), "neb")
    add("H-gamma", ((4341.69, 1.0),), "neb")
    add("H-beta", ((4862.68, 1.0),), "neb")
    add("[O III] 5007", ((4960.30, 1.0 / 2.98), (5008.24, 1.0)), "neb", note="4960+5008 (fixed 1:2.98)")
    add("H-alpha", ((6564.61, 1.0),), "neb")
    add("[N II] 6585", ((6549.86, 1.0 / 3.05), (6585.27, 1.0)), "neb", note="6549+6585 (fixed 1:3.05)")
    add("[S II] 6718", ((6718.29, 1.0), (6732.67, 0.85)), "neb", fixed_ratio=False,
        note="6718+6733 doublet; ratio density-dependent")
    add("[O I] 6302", ((6302.05, 1.0), (6365.54, 1.0 / 3.0)), "neb", note="6300+6364 (fixed 3:1)")
    # rest-UV nebular / wind lines of high-z star-forming galaxies
    add("C III] 1908", ((1906.68, 1.0), (1908.73, 0.7)), "neb", fixed_ratio=False,
        note="1907/1909 blend")
    add("O III] 1663", ((1660.81, 0.4), (1666.15, 1.0)), "neb")
    add("He II 1640", ((1640.42, 1.0),), "neb")
    add("[Ne IV]/C II] 2326", ((2326.0, 1.0),), "neb", note="C II] 2326 blend")
    # ---------------- emission: resonant ----------------
    add("Ly-alpha", ((1215.67, 1.0),), "res_em", note="asymmetric, redshifted peak")
    add("C IV 1549 (em)", ((1548.19, 1.0), (1550.77, 0.5)), "res_em", note="AGN / young stars")
    add("N V 1240 (em)", ((1238.82, 1.0), (1242.80, 0.5)), "res_em", note="AGN")
    add("Mg II 2799 (em)", ((2796.35, 1.0), (2803.53, 0.7)), "res_em", fixed_ratio=False,
        note="often partly absorbed")

    # ---------------- absorption: ISM / outflow ----------------
    for nm, lam in [("Si II 1260", 1260.42), ("O I 1302", 1302.17), ("Si II 1304", 1304.37),
                    ("C II 1334", 1334.53), ("Si II 1527", 1526.71), ("Fe II 1608", 1608.45),
                    ("Al II 1670", 1670.79), ("Fe II 2344", 2344.21), ("Fe II 2374", 2374.46),
                    ("Fe II 2382", 2382.77), ("Fe II 2586", 2586.65), ("Fe II 2600", 2600.17),
                    ("Mg I 2852", 2852.96)]:
        add(nm, ((lam, 1.0),), "ism_abs")
    add("Si IV 1394", ((1393.76, 1.0),), "ism_abs")
    add("Si IV 1403", ((1402.77, 1.0),), "ism_abs")
    add("C IV 1549 (abs)", ((1548.19, 1.0), (1550.77, 0.5)), "ism_abs", fixed_ratio=False)
    add("Al III 1855", ((1854.72, 1.0), (1862.79, 0.5)), "ism_abs", fixed_ratio=False)
    add("Mg II 2799 (abs)", ((2796.35, 1.0), (2803.53, 0.8)), "ism_abs", fixed_ratio=False,
        note="blueshifted in outflows")

    # ---------------- absorption: stellar photospheric (optical, entered in air) ----------------
    add("Ca II K", ((_air(3933.66), 1.0),), "stel_abs")
    add("Ca II H", ((_air(3968.47), 1.0),), "stel_abs", note="blends with H-epsilon")
    add("G-band 4304", ((_air(4304.4), 1.0),), "stel_abs")
    add("Mg b 5175", ((_air(5167.3), 1.0), (_air(5172.7), 1.0), (_air(5183.6), 1.0)), "stel_abs",
        fixed_ratio=True, note="unresolved triplet")
    add("Fe 5270", ((_air(5269.5), 1.0),), "stel_abs")
    add("Na D 5893", ((_air(5889.95), 1.0), (_air(5895.92), 0.8)), "stel_abs", fixed_ratio=False,
        note="may be ISM / airglow-contaminated")
    add("Ca II 8498", ((_air(8498.02), 1.0),), "stel_abs")
    add("Ca II 8542", ((_air(8542.09), 1.0),), "stel_abs")
    add("Ca II 8662", ((_air(8662.14), 1.0),), "stel_abs")
    return L


LINES: List[Line] = _build()
LINES_BY_NAME: Dict[str, Line] = {l.name: l for l in LINES}


def select_lines(kinds: Sequence[str] | None = None, emission: bool | None = None,
                 names: Sequence[str] | None = None) -> List[Line]:
    """Filter the catalogue by kind / emission-vs-absorption / explicit names."""
    out = []
    for l in LINES:
        if kinds is not None and l.kind not in kinds:
            continue
        if emission is not None and l.is_emission != emission:
            continue
        if names is not None and l.name not in names:
            continue
        out.append(l)
    return out


def default_line_names(mode: str = "both") -> List[str]:
    """Names used by the blind line search for ``mode`` in {'emission','absorption','both'}."""
    mode = mode.lower()
    out = []
    for l in LINES:
        if mode == "emission" and not l.is_emission:
            continue
        if mode == "absorption" and l.is_emission:
            continue
        out.append(l.name)
    return out


def observed_wavelengths_air(line: Line, z: float) -> np.ndarray:
    """Observed AIR wavelengths of a line's components at redshift z.

    (1+z) scales *vacuum* wavelengths exactly; the result is then converted to air
    because the pipeline works on an air wavelength axis.
    """
    from .wavelength import vac_to_air
    lam_vac = np.array([c[0] for c in line.comps]) * (1.0 + z)
    return vac_to_air(lam_vac)


def expected_lines(z: float, wmin: float, wmax: float, names: Sequence[str] | None = None
                   ) -> List[Tuple[Line, float]]:
    """Lines whose strongest component falls inside [wmin, wmax] at redshift z."""
    out = []
    for l in (LINES if names is None else [LINES_BY_NAME[n] for n in names if n in LINES_BY_NAME]):
        strongest = int(np.argmax([c[1] for c in l.comps]))
        lam = observed_wavelengths_air(l, z)[strongest]
        if wmin <= lam <= wmax:
            out.append((l, float(lam)))
    return out
