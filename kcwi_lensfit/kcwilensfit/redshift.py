"""Per-aperture redshift / kinematics determination and compound-lens consistency checks.

The :class:`RedshiftEngine` owns the fitting *policy* -- which method for which kind of object, how
other planes' emission is masked, which cross-checks are run -- and uses the numerical building
blocks in :mod:`ppxf_fit` (stellar continuum) and :mod:`linefit` (emission/absorption lines).

Policy, in astronomer's terms
-----------------------------
* **Deflector** (massive early-type galaxy): pPXF stellar-continuum fit.  Scan chi2(z) with sigma
  fixed -> pick a minimum -> free (V, sigma) fit -> bootstrap error -> degree/library sweep for the
  systematic -> blue-only vs red-only consistency -> independent stellar-absorption-line check.
  Emission lines from *other* planes that land in its window are masked.
* **Source** (lensed star-forming galaxy): matched-filter line search (emission, absorption or
  both) -> joint line fit with per-class velocity offsets; the continuum template fit is available
  as a second, independent estimate.  A single isolated line is never called secure.
* **Lensing geometry** is enforced as a *consistency test*, never imposed: a source must lie behind
  the primary deflector; images grouped as one source must agree; a compound lens must show a
  sensible redshift ladder.  Violations are reported loudly, not silently "fixed".
"""
from __future__ import annotations

import os
import warnings
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import ppxf_fit as pf
from . import linefit as lf
from .apertures import Aperture
from .extraction import ExtractedSpectrum
from .lines import default_line_names, select_lines
from .masks import (STRONG_EMISSION, default_windows, line_windows, merge_windows)
from .schema import parse_float_list, parse_name_list, parse_range
from .templates import resolve_library
from .wavelength import C_KMS, dz_to_kms, vac_to_air

ProgressCB = Optional[Callable[[float, str], None]]
Window = Tuple[float, float]


class Cancelled(Exception):
    """Raised by a progress callback to abort a long computation."""


@dataclass
class ApertureResult:
    id: str
    role: str
    group: str
    mode: str                               # 'continuum' | 'lines' | 'both'
    z: float = float("nan")
    z_err: float = float("nan")             # adopted statistical error
    z_sys_err: float = 0.0                  # systematic (sweep / line-class offsets), in z
    status: str = "none"                    # secure | tentative | unconstrained | failed | none
    method: str = ""                        # which estimate is adopted: 'lines' | 'continuum'
    sigma: float = float("nan")
    sigma_err: float = float("nan")
    sigma_sys: float = float("nan")
    kinematics_ok: bool = False
    snr_per_A: float = float("nan")
    chi2: float = float("nan")
    flags: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)
    detail: Dict[str, object] = field(default_factory=dict)       # JSON-able extras
    artifacts: Dict[str, object] = field(default_factory=dict)    # arrays / fit objects for plots

    def row(self) -> dict:
        return dict(id=self.id, role=self.role, group=self.group, method=self.method, status=self.status,
                    z=self.z, z_err=self.z_err, z_sys_err=self.z_sys_err, sigma_kms=self.sigma,
                    z_sys_estimate=self.detail.get("z_sys_estimate", float("nan")),
                    sigma_err_kms=self.sigma_err, sigma_sys_kms=self.sigma_sys,
                    kinematics_reliable=self.kinematics_ok, snr_per_A=self.snr_per_A, chi2_red=self.chi2,
                    flags="; ".join(self.flags))


@dataclass
class Candidate:
    z: float
    score: float                 # line score S(z), or contrast for continuum scans
    kind: str                    # 'lines' | 'continuum'
    info: str = ""
    below_deflector: bool = False
    significant: Optional[bool] = None   # beats the noise-only maximum of the whole scan


@dataclass
class ScanResult:
    aid: str
    kind: str
    z: np.ndarray
    value: np.ndarray            # chi2 (continuum) or score (lines)
    candidates: List[Candidate]
    extra: dict = field(default_factory=dict)


class RedshiftEngine:
    def __init__(self, cfg: dict, spectra: Dict[str, ExtractedSpectrum], apertures: Dict[str, Aperture],
                 fwhm: Dict[str, float], arms_available: Sequence[str], log: Callable[[str], None] = print):
        self.cfg = cfg
        self.spectra = spectra
        self.apertures = apertures
        self.fwhm = fwhm
        self.arms_available = list(arms_available)
        self.log = log
        self.results: Dict[str, ApertureResult] = {}
        self.last_mask_info: dict = dict(masked=[], skipped_same_plane=[])

    # ------------------------------------------------------------------ helpers
    def _arms(self) -> List[str]:
        sel = self.cfg.get("fit_arms", "both")
        arms = [a for a in (["blue", "red"] if sel == "both" else [sel]) if a in self.arms_available]
        if not arms:
            raise ValueError(f"fit_arms='{sel}' but that arm has no data")
        return arms

    def _role(self, aid: str) -> str:
        sp = self.spectra[aid]
        return sp.role

    def effective_mode(self, aid: str) -> str:
        ap = self.apertures.get(aid)
        mode = ap.mode if ap is not None else "auto"
        if mode == "skip":
            return "skip"
        if mode != "auto":
            return mode
        role = self._role(aid)
        return {"deflector": "continuum", "source": "lines", "other": "both"}.get(role, "both")

    def z_range(self, aid: str) -> Tuple[float, float]:
        ap = self.apertures.get(aid)
        role = self._role(aid)
        if role == "deflector":
            lo, hi = parse_range(self.cfg["deflector_z_range"])
        elif role == "source":
            lo, hi = parse_range(self.cfg["source_z_range"])
        else:
            lo = parse_range(self.cfg["deflector_z_range"])[0]
            hi = parse_range(self.cfg["source_z_range"])[1]
        if ap is not None:
            if ap.z_min is not None:
                lo = ap.z_min
            if ap.z_max is not None:
                hi = ap.z_max
        return float(lo), float(hi)

    def deflector_z(self) -> Optional[float]:
        zs = [r.z for r in self.results.values() if r.role == "deflector" and np.isfinite(r.z)
              and r.status in ("secure", "tentative")]
        return min(zs) if zs else None

    def _members(self, aid: str) -> set:
        """Spectra belonging to the same physical plane as ``aid`` (same group / the combined entry)."""
        sp = self.spectra[aid]
        grp = sp.group or ""
        mem = {aid}
        if grp:
            mem |= {k for k, s in self.spectra.items() if (s.group or "") == grp}
            mem.add(f"group:{grp}")
        if aid.startswith("group:"):
            g = aid[6:]
            mem |= {k for k, s in self.spectra.items() if (s.group or "") == g}
        return mem

    def static_windows(self) -> List[Window]:
        return default_windows(self.cfg)

    def other_plane_windows(self, aid: str, wmin: float, wmax: float, z_self: Optional[float] = None) -> List[Window]:
        """Emission-line masks for the redshifts of every OTHER plane with a trusted result.

        A plane is NOT masked when it is the same plane as this aperture:
          * members of the same group (declared images of one source), or
          * any other region whose redshift lies within ``same_plane_tol_kms`` of ``z_self``.
        The second rule matters when two regions are not grouped but may be the same source (that is
        what you are testing): masking region A's lines while fitting region B would erase the very
        lines that show B sits at A's redshift.  ``z_self`` is only known once a redshift guess exists,
        so scans over unknown z pass ``None`` and callers skip other-plane masking there altogether.
        The decision is recorded in ``self.last_mask_info`` for the fit summary.
        """
        info = dict(masked=[], skipped_same_plane=[])
        self.last_mask_info = info
        if not self.cfg.get("mask_other_planes_lines", True):
            return []
        half = float(self.cfg.get("other_line_mask_kms", 600.0))
        tol = float(self.cfg.get("same_plane_tol_kms", 1500.0))
        mem = self._members(aid)
        out: List[Window] = []
        seen = []
        for k, r in self.results.items():
            if k in mem or r.status not in ("secure", "tentative") or not np.isfinite(r.z):
                continue
            if z_self is not None and abs(C_KMS * (r.z - z_self) / (1.0 + z_self)) < tol:
                info["skipped_same_plane"].append((k, float(r.z)))
                continue
            if any(abs(r.z - s_) < 2e-3 for s_ in seen):
                continue
            seen.append(r.z)
            info["masked"].append((k, float(r.z)))
            names = list(STRONG_EMISSION)
            lf_det = r.detail.get("lines_detected")
            if lf_det:
                names = sorted(set(names) | set(lf_det))
            out += line_windows(names, r.z, half, wmin, wmax)
        return merge_windows(out)

    def _settings(self, aid: str, scan: bool = False, **over) -> pf.FitSettings:
        cfg = self.cfg
        role = self._role(aid)
        libname = cfg["sps_library"] if role == "deflector" else cfg["source_sps_library"]
        lib = resolve_library(libname, cfg.get("custom_template_path", ""), cfg.get("custom_template_frame", "air"))
        age = parse_range(cfg["age_range"]) if not lib.is_empirical else None
        met = parse_range(cfg["metal_range"]) if not lib.is_empirical else None
        st = pf.FitSettings(library=lib, moments=int(cfg["moments"]), degree=int(cfg["degree"]),
                            mdegree=int(cfg["mdegree"]), age_range=age, metal_range=met,
                            sigma_init=float(cfg["sigma_init"]), vel_bound=float(cfg["vel_bound_kms"]),
                            mask_nebular=cfg["mask_nebular"] if role == "deflector" else "none",
                            sigma_fixed=None)
        rr = str(cfg.get("fit_rest_range", "") or "").strip()
        if rr and role == "deflector":
            st.rest_range = parse_range(rr)
        for k, v in over.items():
            setattr(st, k, v)
        return st

    def _line_names(self) -> List[str]:
        return default_line_names(self.cfg["line_sets"])

    def _progress_scaled(self, progress: ProgressCB, lo: float, hi: float) -> ProgressCB:
        if progress is None:
            return None
        return lambda f, m: progress(lo + (hi - lo) * f, m)

    # ------------------------------------------------------------------ scans
    def scan(self, aid: str, kind: Optional[str] = None, progress: ProgressCB = None) -> List[ScanResult]:
        """Run the blind redshift scan(s) for an aperture; returns one ScanResult per method."""
        mode = self.effective_mode(aid)
        kinds = {"continuum": ["continuum"], "lines": ["lines"], "both": ["lines", "continuum"]}[mode]
        if kind:
            kinds = [kind]
        out = []
        for i, k in enumerate(kinds):
            p = self._progress_scaled(progress, i / len(kinds), (i + 1) / len(kinds))
            out.append(self._scan_continuum(aid, p) if k == "continuum" else self._scan_lines(aid, p))
        return out

    def _windows_for(self, aid: str, arms, z_self: Optional[float] = None, use_other: bool = True) -> List[Window]:
        """Static sky/telluric windows + (optionally) other planes' emission lines.

        ``use_other=False`` for any search over an UNKNOWN redshift (it must be free to find a line at
        another plane's redshift); with a known guess pass ``z_self`` so the same plane is not masked."""
        sp = self.spectra[aid]
        wmin = min(sp.arms[a].wave.min() for a in arms)
        wmax = max(sp.arms[a].wave.max() for a in arms)
        self.last_mask_info = dict(masked=[], skipped_same_plane=[])
        other = self.other_plane_windows(aid, wmin, wmax, z_self) if use_other else []
        return merge_windows(self.static_windows() + other)

    def _scan_lines(self, aid: str, progress: ProgressCB = None) -> ScanResult:
        cfg = self.cfg
        arms = self._arms()
        sp = self.spectra[aid]
        win = self._windows_for(aid, arms, use_other=False)      # unknown z: never mask other planes' lines
        ls = lf.prepare_line_spectrum(sp, arms, self.fwhm, win, float(cfg["line_continuum_window_A"]))
        zlo, zhi = self.z_range(aid)
        step = float(cfg["z_step_lines"])
        # uniform in ln(1+z): constant velocity step
        z_grid = np.arange(zlo, zhi + step / 2, step)
        names = self._line_names()
        sig = float(cfg["line_sigma_kms"])
        scan = lf.scan_lines(ls, names, z_grid, sig, float(cfg["detect_snr"]), None if progress is None
                             else self._progress_scaled(progress, 0.0, 0.5))
        cands = lf.find_line_candidates(scan)
        thr = None
        nsim = int(cfg.get("n_null_sims", 0))
        if nsim > 0:
            thr = lf.null_score_threshold(ls, names, z_grid, sig, nsim, progress=None if progress is None
                                          else self._progress_scaled(progress, 0.5, 1.0))
        zdef = self.deflector_z()
        out_c = []
        for c in cands:
            info = ", ".join(f"{n} ({s:.1f}σ)" for n, s in c["lines"][:5])
            below = (self._role(aid) == "source" and zdef is not None and c["z"] <= zdef)
            sig_ok = None if thr is None else bool(c["score"] > thr)
            tag = "" if sig_ok is None else (" — SIGNIFICANT" if sig_ok else " — consistent with noise")
            out_c.append(Candidate(z=c["z"], score=c["score"], kind="lines",
                                   info=f"{c['n_det']} line(s) ≥{cfg['detect_snr']:.0f}σ: {info}{tag}",
                                   below_deflector=below, significant=sig_ok))
        res = ScanResult(aid=aid, kind="lines", z=z_grid, value=scan["score"], candidates=out_c,
                         extra=dict(n_det=scan["n_det"], ls=ls, scan=scan, windows=win, null_threshold=thr))
        return res

    def _scan_continuum(self, aid: str, progress: ProgressCB = None) -> ScanResult:
        cfg = self.cfg
        arms = self._arms()
        sp = self.spectra[aid]
        st = self._settings(aid)
        frame = st.library.native_frame
        # unknown z: other planes' lines are masked only for deflectors (protecting sigma); a source/other
        # region may itself lie at another plane's redshift, so its scan must not be blinded
        win = self._windows_for(aid, arms, use_other=(self._role(aid) == "deflector"))
        zlo, zhi = self.z_range(aid)
        step = float(cfg["z_step_continuum"])
        z_grid = np.arange(zlo, zhi + step / 2, step)
        sigmas = self.scan_sigmas()
        fs = pf.build_fit_spectrum(sp, arms, self.fwhm, frame, win)
        found = []                       # one dict per (sigma, candidate)
        curves = []
        nS = len(sigmas)
        for si, sg in enumerate(sigmas):
            st_s = pf.FitSettings(**{**st.__dict__, "sigma_fixed": sg})
            p1 = self._progress_scaled(progress, si / nS, (si + 0.7) / nS)
            scan = pf.scan_continuum_z(sp, arms, self.fwhm, frame, win, z_grid, st_s, progress=p1)
            curves.append((sg, scan["chi2"]))
            coarse = pf.find_chi2_candidates(scan["z"], scan["chi2"], scan["z_fit"], n_max=6)
            # refine the best few minima with the full pPXF machinery (sigma fixed, V free)
            ref_st = pf.FitSettings(**{**st.__dict__, "has_gas": False, "moments": 2, "vel_bound": 700.0,
                                       "sigma_fixed": sg})
            for i, c in enumerate(coarse):
                try:
                    r = pf.fit_at_z(fs, c["z_grid"], ref_st)
                    found.append(dict(z=r["z"], chi2=r["chi2"], contrast=c["contrast"], sigma=sg))
                except Exception:
                    found.append(dict(z=c["z"], chi2=np.nan, contrast=c["contrast"], sigma=sg))
                if progress:
                    progress((si + 0.7 + 0.3 * (i + 1) / max(len(coarse), 1)) / nS, f"refining candidate {i + 1}/{len(coarse)} (sigma {sg:.0f})")
        # merge candidates that are the same redshift (within 1500 km/s) across the sigma scans
        merged: List[dict] = []
        for c in sorted(found, key=lambda d: (np.inf if not np.isfinite(d["chi2"]) else d["chi2"])):
            m = next((m for m in merged if abs(C_KMS * (c["z"] - m["z"]) / (1 + m["z"])) < 1500.0), None)
            if m is None:
                merged.append(dict(c, sigmas=[c["sigma"]]))
            elif c["sigma"] not in m["sigmas"]:
                m["sigmas"].append(c["sigma"])
        # candidates found for EVERY sigma are robust to that assumption: list them first
        merged.sort(key=lambda d: (-len(d["sigmas"]), np.inf if not np.isfinite(d["chi2"]) else d["chi2"]))
        zdef = self.deflector_z()
        out_c = []
        for c in merged:
            robust = (len(c["sigmas"]) == nS) if nS > 1 else None
            tag = "" if robust is None else (f"; present at all {nS} sigma values (stable)" if robust else
                                             f"; only at sigma={','.join(f'{x:.0f}' for x in c['sigmas'])} -> depends on the assumed sigma, treat with caution")
            out_c.append(Candidate(z=c["z"], score=c["contrast"], kind="continuum",
                                   info=f"pPXF chi2/dof={c['chi2']:.3f}; scan chi2 is {c['contrast']:.1f}x lower than a typical trial z{tag}",
                                   below_deflector=(self._role(aid) == "source" and zdef is not None and c["z"] <= zdef),
                                   significant=robust))
        return ScanResult(aid=aid, kind="continuum", z=z_grid, value=curves[0][1], candidates=out_c,
                          extra=dict(scan=dict(z=z_grid, chi2=curves[0][1]), curves=curves, windows=win,
                                     library=st.library.name, sigmas=sigmas))

    def scan_sigmas(self) -> List[float]:
        """Velocity dispersions [km/s] held fixed in the continuum scan (one scan per value)."""
        vals = [v for v in parse_float_list(self.cfg.get("scan_sigma_kms", "")) if v > 0]
        return vals or [220.0]

    @staticmethod
    def best_candidate(scan: ScanResult, allow_below_deflector: bool = False) -> Optional[Candidate]:
        for c in scan.candidates:
            if c.below_deflector and not allow_below_deflector:
                continue
            return c
        return None

    # ------------------------------------------------------------------ fits
    def fit(self, aid: str, z_guess: Optional[float] = None, mode: Optional[str] = None,
            scans: Optional[Sequence[ScanResult]] = None, progress: ProgressCB = None) -> ApertureResult:
        """Full measurement for one aperture; ``z_guess`` from the user or a scan candidate."""
        mode = mode or self.effective_mode(aid)
        sp = self.spectra[aid]
        res = ApertureResult(id=aid, role=sp.role, group=sp.group, mode=mode)
        if mode == "skip":
            res.notes.append("skipped by user")
            self.results[aid] = res
            return res
        ap = self.apertures.get(aid)
        if z_guess is None and ap is not None and ap.z_known is not None:
            z_guess = ap.z_known
        scans = list(scans or [])
        try:
            if mode in ("lines", "both"):
                zg = z_guess
                if zg is None:
                    sc = next((s for s in scans if s.kind == "lines"), None) or \
                        self._scan_lines(aid, self._progress_scaled(progress, 0, 0.4))
                    scans.append(sc)
                    best = self.best_candidate(sc)
                    if best is None:
                        res.flags.append("line scan found no candidate that is physically allowed for this role")
                    else:
                        zg = best.z
                if zg is not None:
                    self._fit_lines(res, aid, zg, scans)
            if mode in ("continuum", "both"):
                zg = z_guess
                if zg is None and res.status in ("secure", "tentative") and np.isfinite(res.z):
                    zg = res.z
                if zg is None:
                    sc = next((s for s in scans if s.kind == "continuum"), None) or \
                        self._scan_continuum(aid, self._progress_scaled(progress, 0.4, 0.7))
                    scans.append(sc)
                    best = self.best_candidate(sc)
                    if best is not None:
                        zg = best.z
                if zg is not None:
                    self._fit_continuum(res, aid, zg, progress=self._progress_scaled(progress, 0.7, 1.0))
                else:
                    res.flags.append("continuum scan found no usable candidate")
        except Cancelled:
            raise
        except Exception as e:
            import traceback
            tb = traceback.extract_tb(e.__traceback__)
            where = f" (in {os.path.basename(tb[-1].filename)}:{tb[-1].lineno}, {tb[-1].name})" if tb else ""
            res.status = "failed"
            res.flags.append(f"fit failed: {type(e).__name__}: {e}{where}")
            res.detail["error_traceback"] = "".join(traceback.format_exception(type(e), e, e.__traceback__))
            self.log(f"[{aid}] FAILED: {type(e).__name__}: {e}{where}")
            if os.environ.get("KCWILENSFIT_DEBUG"):
                self.log(res.detail["error_traceback"])
        res.artifacts["scans"] = scans
        self._finalise_status(res)
        try:
            res.detail["provenance_text"] = self.provenance(res, aid, scans)
            for line in res.detail["provenance_text"]:
                self.log("    " + line)
        except Exception as e:                                   # never let reporting break a fit
            res.notes.append(f"could not assemble the fit summary: {e}")
        self.results[aid] = res
        return res

    # ---- lines
    def _fit_lines(self, res: ApertureResult, aid: str, z_guess: float, scans) -> None:
        cfg = self.cfg
        arms = self._arms()
        sp = self.spectra[aid]
        sc = next((s for s in scans if s.kind == "lines"), None)
        if sc is not None:
            ls = sc.extra["ls"]
        else:
            win = self._windows_for(aid, arms, use_other=False)
            ls = lf.prepare_line_spectrum(sp, arms, self.fwhm, win, float(cfg["line_continuum_window_A"]))
        names = self._line_names()
        sig0 = float(cfg["line_sigma_kms"])
        z0 = lf.refine_peak(ls, names, z_guess, sigma_v=sig0)
        fit = lf.refine_lines(ls, z0, names, sigma_v0=sig0, select_snr=2.0)
        res.artifacts["ls"] = ls
        if fit is None:
            res.flags.append(f"no line significant (>2σ) at z≈{z0:.4f}")
            res.detail["line_fit"] = None
            if res.status in ("none",):
                res.status = "failed"
            return
        res.artifacts["line_fit"] = fit
        det = [r for r in fit["lines"] if r["snr"] >= float(cfg["detect_snr"])]
        res.detail["lines"] = [{k: (v if not isinstance(v, np.floating) else float(v)) for k, v in r.items()}
                               for r in fit["lines"]]
        res.detail["lines_detected"] = [r["name"] for r in det]
        res.detail["line_sigma_kms"] = fit["sigma_v"]
        res.detail["line_class_offsets_kms"] = fit["class_offsets"]
        res.detail["line_z_consistency_chi2"] = fit["z_consistency_chi2"]
        if "z_sys_estimate" in fit:
            res.detail["z_sys_estimate"] = fit["z_sys_estimate"]
        res.detail["z_lines"] = fit["z"]
        res.detail["z_lines_err"] = fit["z_err"]
        n_ind = len(det)
        thr = sc.extra.get("null_threshold") if sc is not None else None
        score_here = float(lf.scan_lines(ls, names, np.array([fit["z"]]), sig0)["score"][0])
        res.detail["line_score"] = score_here
        res.detail["line_score_null_p95_max"] = thr
        res.z, res.z_err, res.method = fit["z"], fit["z_err"], "lines"
        res.z_sys_err = fit.get("z_sys_err", 0.0)
        res.detail["line_kinds"] = fit["kinds_used"]
        if not np.isfinite(fit["z_err"]) or fit["z_err"] > 0.02:
            res.status = "unconstrained"
            res.flags.append("the joint line fit is degenerate (redshift error is not finite or > 0.02) -- "
                             "no reliable line redshift")
        elif thr is not None and score_here <= thr:
            res.status = "tentative"
            res.flags.append(f"line score {score_here:.1f} does not exceed the noise-only maximum of the scan "
                             f"({thr:.1f}): this 'detection' is consistent with chance given the number of trials")
        elif n_ind >= int(cfg["min_lines_secure"]):
            res.status = "secure"
        elif n_ind >= 1:
            res.status = "tentative"
            lam_f = det[0]["lam_obs"]
            idt = lf.identification_table(ls, lam_f, "emission" if det[0]["emission"] else "absorption", sig0)
            res.detail["identification_table"] = idt[:8]
            res.flags.append(f"only {n_ind} line detected -> competing identifications are listed in the report; "
                             "a single line is not proof of a redshift")
        else:
            res.status = "unconstrained"
            res.flags.append("no line reaches the detection threshold in the joint fit")
        if fit.get("z_consistency_chi2") is not None and np.isfinite(fit["z_consistency_chi2"]) and fit["n_consistency"] >= 2:
            ndof = fit["n_consistency"] - 1
            if fit["z_consistency_chi2"] > 3.0 * ndof + 6.0:
                res.flags.append(f"individual line redshifts are mutually inconsistent "
                                 f"(chi2={fit['z_consistency_chi2']:.1f} for {ndof} dof) -- check the identifications")
        if not fit["have_ref"] and det:
            if "z_sys_estimate" in fit:
                res.notes.append(
                    f"no nebular/stellar (systemic) line detected: z is defined by the strongest class "
                    f"('{fit.get('reference_kind')}'), which typically sits {fit['typical_offset_kms']:+.0f}±{fit['sys_err_kms']:.0f} km/s "
                    f"from the systemic redshift (systemic estimate z≈{fit['z_sys_estimate']:.5f}); this systematic is attached")
        zdef = self.deflector_z()
        if res.role == "source" and zdef is not None and res.z <= zdef:
            res.flags.append(f"z={res.z:.5f} is NOT above the deflector's z={zdef:.5f}: it cannot be a lensed image "
                             "of a source behind that deflector (interloper, or wrong line identification)")

    # ---- continuum
    def _fit_continuum(self, res: ApertureResult, aid: str, z_guess: float, progress: ProgressCB = None) -> None:
        cfg = self.cfg
        arms = self._arms()
        sp = self.spectra[aid]
        st = self._settings(aid)
        frame = st.library.native_frame
        win = self._windows_for(aid, arms, z_self=z_guess)
        res.detail["mask_info"] = dict(self.last_mask_info)
        fit_scale = bool(cfg["fit_red_scale"]) and len(arms) == 2
        extra_bad: List[Window] = []
        fit = pf.fit_joint(sp, arms, self.fwhm, frame, win, z_guess, st, fit_scale=fit_scale, extra_bad=extra_bad)
        # re-centre once so |V| is small and the template mapping is accurate
        if abs(fit["v"]) > 150.0:
            fit = pf.fit_joint(sp, arms, self.fwhm, frame, win, fit["z"], st, fit_scale=fit_scale,
                               scale_init=fit["scale_r"])
        if fit["at_bound"]:
            res.flags.append("velocity landed on the search bound -- not converged; try another redshift guess")
        res.artifacts["ppxf"] = fit
        fs = fit["fs"]
        gp_mask = np.zeros(fs.lam.size, bool)
        gp_mask[fit["goodpixels"]] = True
        snr = pf.snr_per_angstrom(fs, gp_mask)
        res.snr_per_A = snr
        res.chi2 = fit["chi2"]
        z_err = fit["z_err"]
        sig_err = fit["sigma_err"]
        boot = None
        nb = int(cfg["n_bootstrap"])
        if nb > 0:
            boot = pf.bootstrap_fit(fit, n=nb, method=cfg["bootstrap_method"], seed=1,
                                    progress=self._progress_scaled(progress, 0.0, 0.5))
            res.artifacts["bootstrap"] = boot
            if boot["ok"]:
                if np.isfinite(boot["z_err"]) and boot["z_err"] > 3 * fit["z_err"]:
                    res.notes.append("bootstrap z error >> formal error: the formal error is optimistic here "
                                     "(continuum/template degeneracy or imperfect noise model)")
                z_err = max(z_err, boot["z_err"]) if np.isfinite(boot["z_err"]) else z_err
                sig_err = boot["sigma_err"]
            else:
                res.flags.append(f"only {boot['n_success']}/{boot['n_requested']} bootstrap fits converged; "
                                 "falling back to formal errors")
        # robustness sweep -> systematic
        degs = [int(d) for d in parse_float_list(cfg["sweep_degrees"])]
        libs = parse_name_list(cfg["sweep_libraries"])
        sweep = None
        if degs or libs:
            custom = (cfg.get("custom_template_path", ""), cfg.get("custom_template_frame", "air"))
            sweep = pf.systematic_sweep(sp, arms, self.fwhm, win, fit["z"], st, lambda lib: lib.native_frame,
                                        degrees=degs, libraries=libs, fit_scale=fit_scale, scale_init=fit["scale_r"],
                                        custom=custom, progress=self._progress_scaled(progress, 0.5, 0.8))
            res.detail["sweep"] = sweep["rows"]
        sigma_sys = sweep.get("sigma_sys", float("nan")) if sweep else float("nan")
        z_sys = sweep.get("z_sys", 0.0) if sweep else 0.0
        # blue-only vs red-only consistency
        if cfg.get("crosscheck_arms", True) and len(arms) == 2:
            xc = {}
            for arm in arms:
                try:
                    r1 = pf.fit_joint(sp, [arm], self.fwhm, frame, win, fit["z"], st, fit_scale=False)
                    xc[arm] = dict(z=r1["z"], z_err=r1["z_err"], sigma=r1["sigma"], sigma_err=r1["sigma_err"],
                                   snr=pf.snr_per_angstrom(r1["fs"]))
                except Exception as e:
                    xc[arm] = dict(error=str(e))
            res.detail["arm_check"] = xc
            if all("error" not in v for v in xc.values()):
                a, b = xc[arms[0]], xc[arms[1]]
                d = abs(a["z"] - b["z"])
                e = np.hypot(a["z_err"], b["z_err"])
                nsig = d / e if e > 0 else np.inf
                res.detail["arm_check_nsigma"] = float(nsig)
                if nsig > 3.0:
                    res.flags.append(f"blue-only and red-only redshifts differ by {nsig:.1f}σ "
                                     f"(Δz={a['z'] - b['z']:+.5f}) -- calibration, template mismatch or contamination?")
        # kinematics reliability
        ok = True
        if not np.isfinite(snr) or snr < float(cfg["snr_min_kinematics"]):
            ok = False
            res.flags.append(f"S/N per Å = {snr:.1f} < {cfg['snr_min_kinematics']:.0f}: velocity dispersion is unreliable "
                             "(redshift may still be fine)")
        if fit.get("sigma_at_bound"):
            ok = False
            res.flags.append(f"sigma={fit['sigma']:.0f} km/s sits on the fit's allowed range (20-600 km/s): not converged / "
                             "unconstrained -- usually too little S/N or a contaminated aperture")
        if fit["sigma"] < 0.5 * fs.velscale:
            ok = False
            res.flags.append(f"sigma={fit['sigma']:.0f} km/s is below half a pixel ({0.5 * fs.velscale:.0f} km/s): unresolved")
        if st.library.is_empirical:
            ok = False
            res.notes.append("empirical templates (~5 Å) cannot constrain sigma at KCWI resolution; sigma not reported as reliable")
        if fit.get("rest_range_ignored"):
            res.notes.append("the data do not cover the requested rest-frame fitting window; the whole spectrum was used "
                             "(sigma then depends on the template resolution in the rest-UV)")
        if fit["degraded_frac"] > 0.01:
            res.notes.append(f"galaxy spectrum was degraded to the template resolution over {fit['degraded_frac']:.0%} of the "
                             "fitted range (noise kept at its nominal level, i.e. conservative); sigma is relative to that matched resolution")
        res.detail.update(ppxf_z=fit["z"], ppxf_z_err=z_err, ppxf_z_formal=fit["z_err"], ppxf_v=fit["v"],
                          ppxf_sigma=fit["sigma"], ppxf_sigma_err=sig_err, scale_r=fit["scale_r"],
                          n_good_pix=fit["n_good"], velscale=fs.velscale, library=st.library.name,
                          mdegree=st.mdegree, degree=st.degree, frame=frame)
        # redshift: independent line check for deflectors (stellar absorption), if requested
        if res.role == "deflector":
            try:
                lines_ok = self._deflector_line_check(res, aid, fit["z"], arms)
            except Exception as e:
                res.notes.append(f"stellar-line cross-check not possible: {e}")
        adopt_continuum = (res.method != "lines") or res.status not in ("secure",)
        if adopt_continuum:
            res.z, res.z_err, res.method = fit["z"], z_err, "continuum"
            res.z_sys_err = float(z_sys)
            res.sigma, res.sigma_err, res.sigma_sys = fit["sigma"], sig_err, sigma_sys
            res.kinematics_ok = ok
            res.status = "secure" if (not fit["at_bound"]) else "tentative"
            if fit["chi2"] > 3.0:
                res.flags.append(f"chi2/dof = {fit['chi2']:.2f}: poor fit (template mismatch, contamination or "
                                 "under-estimated noise); errors have been rescaled by sqrt(chi2)")
        else:
            res.detail["continuum_agrees_with_lines_nsigma"] = self._agree(res.z, res.z_err, fit["z"], z_err)
            res.sigma, res.sigma_err, res.sigma_sys = fit["sigma"], sig_err, sigma_sys
            res.kinematics_ok = ok
            nsig = res.detail["continuum_agrees_with_lines_nsigma"]
            if nsig > 3.0:
                res.flags.append(f"continuum fit z={fit['z']:.5f} disagrees with the line redshift by {nsig:.1f}σ")
        zdef = self.deflector_z()
        if res.role == "source" and zdef is not None and res.z <= zdef:
            res.flags.append(f"z={res.z:.5f} is NOT above the deflector's z={zdef:.5f}")

    @staticmethod
    def _agree(z1, e1, z2, e2) -> float:
        e = np.hypot(e1, e2)
        return float(abs(z1 - z2) / e) if e > 0 else float("inf")

    def _deflector_line_check(self, res: ApertureResult, aid: str, z: float, arms) -> bool:
        """Independent z from stellar absorption lines (Ca H&K, G, Mg b, Na D, Ca triplet)."""
        names = [l.name for l in select_lines(kinds=("stel_abs",))]
        sp = self.spectra[aid]
        win = self._windows_for(aid, arms, use_other=False)
        ls = lf.prepare_line_spectrum(sp, arms, self.fwhm, win, float(self.cfg["line_continuum_window_A"]))
        fit = lf.refine_lines(ls, z, names, sigma_v0=float(np.mean(self.scan_sigmas())), select_snr=2.0,
                              sigma_range=(100.0, 450.0))
        if fit is None:
            res.notes.append("no stellar absorption line significant for the independent line check")
            return False
        det = [r for r in fit["lines"] if r["snr"] >= 3.0]
        res.detail["stellar_line_check"] = dict(z=fit["z"], z_err=fit["z_err"], n_lines=len(det),
                                                lines=[r["name"] for r in det])
        res.artifacts["stellar_line_fit"] = fit
        nsig = self._agree(fit["z"], fit["z_err"], z, res.detail.get("ppxf_z_formal", 1e-4))
        res.detail["stellar_line_check"]["nsigma_vs_ppxf"] = nsig
        if nsig > 3.0 and len(det) >= 2:
            res.flags.append(f"stellar-line redshift ({fit['z']:.5f}) differs from the pPXF value by {nsig:.1f}σ")
        return True

    # ------------------------------------------------------------------ provenance
    def provenance(self, res: ApertureResult, aid: str, scans) -> List[str]:
        """Plain-language record of WHAT was fitted and HOW, for the console, the GUI and the report."""
        cfg = self.cfg
        sp = self.spectra[aid]
        arms = self._arms()
        L: List[str] = []
        L.append(f"FIT SUMMARY for '{aid}' [{sp.role}], mode '{res.mode}':")
        # ---- data used
        cov = []
        for a in arms:
            w = sp.arms[a].wave[sp.arms[a].good]
            cov.append(f"{a} {w.min():.0f}-{w.max():.0f} Å" if w.size else f"{a} (no valid pixels)")
        which = {"both": "BOTH arms (joint fit)", "blue": "BLUE arm only", "red": "RED arm only"}[cfg.get("fit_arms", "both")]
        L.append(f"  arms fitted: {which}; usable coverage: " + "; ".join(cov)
                 + (f"; cut-over between arms at {res.artifacts['ppxf']['fs'].cutover:.0f} Å" if "ppxf" in res.artifacts and len(arms) == 2 else ""))
        wins = self.static_windows()
        mi = res.detail.get("mask_info") or {}
        om = ", ".join(f"{k} (z={z_:.4f})" for k, z_ in mi.get("masked", []))
        osk = ", ".join(f"{k} (z={z_:.4f})" for k, z_ in mi.get("skipped_same_plane", []))
        L.append(f"  masked: {len(wins)} sky/telluric windows; other planes' emission lines masked in the continuum fit: "
                 + (om or ("none" if cfg.get("mask_other_planes_lines", True) else "off (setting disabled)"))
                 + (f"; NOT masked because same plane as this aperture (within {cfg.get('same_plane_tol_kms', 1500.0):g} km/s): {osk}" if osk else "")
                 + ("; line searches never mask other planes" if res.mode in ("lines", "both") else ""))
        # ---- methods
        kinds = {sc.kind for sc in (scans or [])}
        if "lines" in kinds or "line_fit" in res.artifacts:
            ls_names = self._line_names()
            sc = next((s_ for s_ in scans if s_.kind == "lines"), None)
            zlo, zhi = self.z_range(aid)
            nsim = int(cfg.get("n_null_sims", 0))
            L.append(f"  [lines] matched-filter line search, {cfg['line_sets']} lines ({len(ls_names)} catalogue multiplets), "
                     f"z {zlo:g}-{zhi:g} step {cfg['z_step_lines']:g}, intrinsic width {cfg['line_sigma_kms']:g} km/s"
                     + (f", significance from {nsim} noise simulations" if nsim else ", no noise-simulation significance test"))
            lfit = res.artifacts.get("line_fit")
            if lfit is not None:
                det = [r["name"] for r in lfit["lines"] if r["snr"] >= float(cfg["detect_snr"])]
                L.append(f"          joint line fit: z={lfit['z']:.5f}; detected (>= {cfg['detect_snr']:g}σ): {', '.join(det) or 'none'}; "
                         f"classes used: {', '.join(lfit['kinds_used']) or '—'}")
        if "continuum" in kinds or "ppxf" in res.artifacts:
            sc = next((s_ for s_ in scans if s_.kind == "continuum"), None)
            fit = res.artifacts.get("ppxf")
            lib = fit["settings"].library if fit is not None else self._settings(aid).library
            if sc is not None:
                zlo, zhi = self.z_range(aid)
                L.append(f"  [continuum scan] high-pass non-negative template fit, library '{lib.name}', z {zlo:g}-{zhi:g} step "
                         f"{cfg['z_step_continuum']:g}, sigma held fixed at {', '.join(f'{x:g}' for x in sc.extra.get('sigmas', []))} km/s; "
                         "top minima refined with full pPXF")
            if fit is not None:
                st = fit["settings"]
                fs = fit["fs"]
                gp = fit["goodpixels"]
                lam = fs.lam[gp]
                air = (lambda x: float(vac_to_air(x))) if fs.frame == "vacuum" else (lambda x: float(x))
                blue_px = gp[~fs.is_red[gp]]
                red_px = gp[fs.is_red[gp]]
                parts = []
                if blue_px.size:
                    parts.append(f"blue {air(fs.lam[blue_px].min()):.0f}-{air(fs.lam[blue_px].max()):.0f} Å")
                if red_px.size:
                    parts.append(f"red {air(fs.lam[red_px].min()):.0f}-{air(fs.lam[red_px].max()):.0f} Å")
                rest = fit["lam_rest"][gp]
                missing_arms = [a for a, px in (("blue", blue_px), ("red", red_px)) if a in arms and px.size == 0]
                if missing_arms:
                    msg = (f"no {'/'.join(missing_arms)}-arm pixels survive the rest-frame window / masks, so the fit "
                           f"used the {'red' if 'blue' in missing_arms else 'blue'} arm only although both arms are selected")
                    L.append("  NOTE: " + msg)
                    res.notes.append(msg)
                ntemp = fit["templates"].shape[1]
                L.append(f"  [pPXF fit] library '{st.library.name}' ({st.library.native_frame}-native, {ntemp} templates"
                         + ("" if st.library.is_empirical else f"; ages {st.age_range[0]:g}-{st.age_range[1]:g} Gyr, [M/H] {st.metal_range[0]:g}-{st.metal_range[1]:g}")
                         + f"), template file {os.path.basename(st.library.path)}")
                L.append(f"          pixels fitted: {gp.size} (observed {', '.join(parts)}; rest {rest.min():.0f}-{rest.max():.0f} Å)"
                         + ("; requested rest window not covered -> whole spectrum used" if fit.get("rest_range_ignored") else
                            (f"; rest window {st.rest_range[0]:g}-{st.rest_range[1]:g} Å" if st.rest_range else "")))
                L.append(f"          moments={st.moments}, multiplicative poly degree={st.mdegree}, additive degree={st.degree}, "
                         f"red/blue scale {'fitted (x%.3f)' % fit['scale_r'] if len(arms) == 2 and cfg.get('fit_red_scale', True) else 'not fitted'}, "
                         f"nebular masking '{st.mask_nebular}'")
        if "bootstrap" in res.artifacts:
            b = res.artifacts["bootstrap"]
            L.append(f"  [errors] {b['method']} bootstrap, {b['n_success']}/{b['n_requested']} fits converged")
        elif res.mode in ("continuum", "both") and "ppxf" in res.artifacts:
            L.append("  [errors] formal pPXF errors only (bootstrap off)")
        if res.detail.get("sweep"):
            L.append("  [systematics] sweep: " + ", ".join(r["label"] for r in res.detail["sweep"]))
        if res.detail.get("arm_check"):
            L.append("  [cross-check] blue-only vs red-only fits: " + ", ".join(res.detail["arm_check"]))
        if res.detail.get("stellar_line_check"):
            c = res.detail["stellar_line_check"]
            L.append(f"  [cross-check] stellar absorption-line fit: z={c['z']:.5f}, {c['n_lines']} lines ({', '.join(c['lines'])})")
        for r in sp.meta.get("subtracted", []):
            L.append(f"  [extraction] manual subtraction: - {r['scale']:.3f} x '{r['contaminator']}'" + (" (auto)" if r["auto"] else ""))
        L.append(f"  adopted: z = {res.z:.5f} ± {res.z_err:.5f} from the {res.method or 'n/a'} method ({res.status})")
        return L

    # ------------------------------------------------------------------ status
    def _finalise_status(self, res: ApertureResult) -> None:
        if res.status == "none" and np.isfinite(res.z):
            res.status = "tentative"

    # ------------------------------------------------------------------ batch
    def run_all(self, progress: ProgressCB = None) -> Dict[str, ApertureResult]:
        """Non-interactive pass: sources first (their lines are masked in the deflector fit), then deflectors,
        then a second pass over everything so each fit sees the final masks."""
        ids = [k for k in self.spectra if self.effective_mode(k) != "skip"]
        order = sorted(ids, key=lambda k: {"source": 0, "other": 1, "deflector": 2}.get(self._role(k), 3))
        n = len(order)
        for passno in (1, 2):
            if passno == 2:
                # second pass only matters if some result changed the masks of another plane
                if n < 2:
                    break
            for i, k in enumerate(order):
                self.log(f"[pass {passno}] {k} ({self._role(k)}) ...")
                if progress:
                    progress(((passno - 1) * n + i) / (2 * n), f"pass {passno}: {k}")
                prev = self.results.get(k)
                z_prev = prev.z if prev is not None and prev.status in ("secure", "tentative") else None
                ap = self.apertures.get(k)
                zg = None
                if passno == 2 and z_prev is not None:
                    zg = z_prev
                if self.cfg.get("auto_accept", True):
                    self.fit(k, z_guess=zg)
                r = self.results.get(k)
                if r is not None:
                    self.log(f"    -> status={r.status} z={r.z:.5f}±{r.z_err:.5f} via {r.method}"
                             + (f"; sigma={r.sigma:.0f} km/s" if np.isfinite(r.sigma) else ""))
        return self.results

    # ------------------------------------------------------------------ lens-system consistency
    def system_checks(self) -> Tuple[List[dict], List[str]]:
        """Redshift ladder + lensing-geometry consistency messages."""
        rows = []
        for k, r in self.results.items():
            if np.isfinite(r.z) and r.status in ("secure", "tentative"):
                rows.append(dict(id=k, role=r.role, group=r.group, z=r.z, z_err=r.z_err, status=r.status, method=r.method))
        rows.sort(key=lambda d: d["z"])
        msgs: List[str] = []
        defl = [r for r in rows if r["role"] == "deflector"]
        srcs = [r for r in rows if r["role"] == "source"]
        if defl and srcs:
            zd = min(r["z"] for r in defl)
            for s in srcs:
                if s["z"] <= zd:
                    msgs.append(f"GEOMETRY: source '{s['id']}' (z={s['z']:.4f}) is not behind the primary deflector "
                                f"(z={zd:.4f}).")
                elif dz_to_kms(s["z"] - zd, zd) < 2000.0:
                    msgs.append(f"NOTE: source '{s['id']}' is within 2000 km/s of the deflector -- a physical companion "
                                "rather than a lensed background galaxy?")
        # compound-lens ladder
        planes = []
        for r in rows:
            if not planes or abs(r["z"] - planes[-1][0]) > 0.003:
                planes.append([r["z"], [r]])
            else:
                planes[-1][1].append(r)
        if len({round(p[0], 2) for p in planes}) >= 3 and defl and srcs:
            msgs.append("COMPOUND LENS: %d distinct redshift planes found: %s." %
                        (len(planes), ", ".join(f"z={p[0]:.4f} [{'/'.join(sorted({x['role'] for x in p[1]}))}]" for p in planes)))
        # images of one source must agree
        groups: Dict[str, list] = {}
        for k, r in self.results.items():
            if r.group and not k.startswith("group:") and np.isfinite(r.z) and r.status in ("secure", "tentative"):
                groups.setdefault(r.group, []).append(r)
        for g, members in groups.items():
            if len(members) < 2:
                continue
            disagree = False
            for i in range(len(members)):
                for j in range(i + 1, len(members)):
                    a, b = members[i], members[j]
                    e = float(np.sqrt(a.z_err ** 2 + b.z_err ** 2 + a.z_sys_err ** 2 + b.z_sys_err ** 2))
                    nsig = abs(a.z - b.z) / e if e > 0 else np.inf
                    if nsig > 3.0:
                        disagree = True
                        msgs.append(f"IMAGES: '{a.id}' (z={a.z:.5f}) and '{b.id}' (z={b.z:.5f}) share group '{g}' but "
                                    f"differ by {nsig:.1f}σ -- not the same source, or a contaminated spectrum.")
            comb = self.results.get(f"group:{g}")
            if disagree and comb is not None and comb.status == "secure":
                comb.status = "tentative"
                flag = (f"status capped at 'tentative': member images of group '{g}' do not agree on a redshift, so the "
                        "combined spectrum mixes different sources/noise")
                if flag not in comb.flags:
                    comb.flags.append(flag)
                msgs.append(f"GROUP: combined spectrum 'group:{g}' downgraded to tentative because its member images disagree.")
        # regions at the SAME redshift but not grouped: possibly one source (or a compound's images)
        cand = [r for r in rows if r["role"] in ("source", "other") and not r["id"].startswith("group:")]
        for i in range(len(cand)):
            for j in range(i + 1, len(cand)):
                a, b = cand[i], cand[j]
                if a["group"] and a["group"] == b["group"]:
                    continue
                ra, rb = self.results[a["id"]], self.results[b["id"]]
                e = float(np.sqrt(ra.z_err ** 2 + rb.z_err ** 2 + ra.z_sys_err ** 2 + rb.z_sys_err ** 2))
                dv = abs(C_KMS * (ra.z - rb.z) / (1.0 + ra.z))
                if (e > 0 and abs(ra.z - rb.z) / e < 3.0) or dv < 300.0:
                    msgs.append(f"SAME PLANE?: '{a['id']}' (z={ra.z:.5f}) and '{b['id']}' (z={rb.z:.5f}) are not grouped but agree "
                                f"within {dv:.0f} km/s -- possibly images of one source; if you confirm that, give them the same "
                                "group to combine them.")
        return rows, msgs
