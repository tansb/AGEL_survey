"""The pipeline object: holds the state and exposes one method per stage.

Stages (each can be run from the GUI, from Python, or all at once with :meth:`run_all`)::

    load()          cubes -> air-wavelength arrays, footprints, white-light images, registration
    estimate_sky()  source-free box -> background, noise, noise correlation
    set_apertures() / apertures_from_regions() / suggest_apertures()
    extract()       aperture sums or spectral deblending -> per-aperture spectra (+ combined images)
    fit_all() / scan() / fit_aperture()   redshifts and kinematics
    summarize()     redshift ladder + lensing-geometry consistency
    save()          spectra, tables, figures, report, the exact settings used
"""
from __future__ import annotations

import copy
import json
import os
import time
import traceback
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import numpy as np

from . import cubeio, noise as nz
from .apertures import (Aperture, aperture_mask, detect_peaks, overlap_fraction, parse_ds9_regions)
from .cubeio import Arm, load_arm, white_light
from .extraction import (ArmSpectrum, ExtractedSpectrum, deblend, extract_sum, inverse_variance_combine,
                         parse_subtraction_rules, reject_spikes, subtract_spectrum)
from .profiles import Component, ProfileFit, fit_profiles
from .redshift import ApertureResult, Cancelled, RedshiftEngine, ScanResult
from .registration import Registration, build_registration
from .schema import (FIELD_BY_KEY, defaults, normalise_aperture, parse_range, validate)

ProgressCB = Optional[Callable[[float, str], None]]
FWHM2SIG = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


class PipelineError(RuntimeError):
    pass


class Pipeline:
    def __init__(self, cfg: Optional[dict] = None, log: Optional[Callable[[str], None]] = None):
        self.cfg: dict = defaults()
        if cfg:
            self.update_config(cfg)
        self.log = log or (lambda m: print(m, flush=True))
        self.reset()

    # ------------------------------------------------------------------ state
    def reset(self):
        self.arms: Dict[str, Arm] = {}
        self.whitelight: Dict[str, np.ndarray] = {}
        self.registration: Registration = Registration()
        self.sky: Dict[str, nz.SkyModel] = {}
        self.apertures: List[Aperture] = []
        self.masks: Dict[str, Dict[str, np.ndarray]] = {}       # arm -> id -> fractional mask
        self.spectra: Dict[str, ExtractedSpectrum] = {}
        self.extract_info: dict = {}
        self.profile_fit: Optional[ProfileFit] = None
        self.engine: Optional[RedshiftEngine] = None
        self.results: Dict[str, ApertureResult] = {}
        self.system_rows: List[dict] = []
        self.system_messages: List[str] = []
        self.timings: Dict[str, float] = {}

    def update_config(self, new: dict) -> None:
        for k, v in new.items():
            if k == "apertures":
                self.cfg["apertures"] = [normalise_aperture(a) for a in v]
            else:
                self.cfg[k] = v

    # ------------------------------------------------------------------ helpers
    @property
    def ref_arm(self) -> str:
        """Arm whose pixel grid defines aperture coordinates (blue if present)."""
        return "blue" if "blue" in self.arms else "red"

    @property
    def fwhm(self) -> Dict[str, float]:
        return {k: a.fwhm_A for k, a in self.arms.items()}

    def _need(self, what: str, ok: bool):
        if not ok:
            raise PipelineError(f"Stage '{what}' needs the earlier stages to be run first.")

    def arm_aperture(self, ap: Aperture, arm: str) -> Aperture:
        """Aperture expressed in ``arm``'s pixel frame."""
        a = self.arms.get(arm)
        if ap.has_sky and self.cfg.get("sky_regions_use_wcs", True) and a is not None and a.wcs is not None:
            from .apertures import aperture_in_arm
            ref = self.arms[self.ref_arm]
            return aperture_in_arm(ap, a.wcs, a.pixscale[1], ref.wcs)
        if arm == self.ref_arm:
            return ap
        return self.registration.to_red(ap) if arm == "red" else ap

    def _frame_shift(self, frame_arm: str, arm: str) -> Tuple[float, float, float]:
        """(dx, dy, scale) mapping pixel coordinates of ``frame_arm`` to ``arm``."""
        if frame_arm == arm:
            return 0.0, 0.0, 1.0
        r = self.registration
        if frame_arm == "blue" and arm == "red":
            return r.dx, r.dy, r.scale
        if frame_arm == "red" and arm == "blue":
            return -r.dx / r.scale, -r.dy / r.scale, 1.0 / r.scale
        return 0.0, 0.0, 1.0

    # ------------------------------------------------------------------ stage 1
    def load(self, progress: ProgressCB = None) -> None:
        t0 = time.time()
        problems = [p for p in validate(self.cfg) if "Aperture" not in p]
        if problems:
            raise PipelineError("Configuration problems:\n  - " + "\n  - ".join(problems))
        cfg = self.cfg
        self.reset()
        sib = None
        if cfg.get("red_cube") and cfg.get("blue_cube"):
            from astropy.io import fits as _f
            try:
                sib_red = _f.getheader(os.path.expanduser(cfg["red_cube"]))
                sib_blue = _f.getheader(os.path.expanduser(cfg["blue_cube"]))
            except Exception:
                sib_red = sib_blue = None
        else:
            sib_red = sib_blue = None
        if cfg.get("blue_cube"):
            self.log("Loading blue cube ...")
            self.arms["blue"] = load_arm("blue", cfg["blue_cube"], cfg, cfg.get("blue_var_cube", ""), sibling_header=sib_red)
        if progress:
            progress(0.4, "blue loaded")
        if cfg.get("red_cube"):
            self.log("Loading red cube ...")
            self.arms["red"] = load_arm("red", cfg["red_cube"], cfg, cfg.get("red_var_cube", ""), sibling_header=sib_blue)
        for name, arm in self.arms.items():
            self.log(f"--- {name} arm: cube {arm.flux.shape} (nw,ny,nx), {arm.pixscale[1]:.3f}\"/pix")
            for line in arm.info:
                self.log("    " + line)
            self.whitelight[name] = white_light(arm)
        if len(self.arms) == 2:
            a, b = self.arms["blue"], self.arms["red"]
            if a.flux.shape[1:] != b.flux.shape[1:]:
                self.log("[note] the two cubes have different spatial shapes; registration mode 'wcs' is recommended.")
        self.registration = build_registration(
            self.arms.get("blue"), self.arms.get("red"), self.whitelight.get("blue"), self.whitelight.get("red"),
            cfg["arm_registration"], (cfg["arm_dx"], cfg["arm_dy"]))
        r = self.registration
        if len(self.arms) == 2:
            d = r.diagnostics
            self.log(f"Registration: {r.mode} (dx={r.dx:+.2f}, dy={r.dy:+.2f} pix) -- {r.note}")
            if "wcs_dx" in d:
                self.log(f"    header-WCS offset would be ({d['wcs_dx']:+.2f}, {d['wcs_dy']:+.2f}) pix; "
                         f"cross-correlation ({d.get('xcorr_dx', float('nan')):+.2f}, {d.get('xcorr_dy', float('nan')):+.2f}) "
                         f"pix, sharpness {d.get('xcorr_sharpness', float('nan')):.2f}")
        self.timings["load"] = time.time() - t0
        if progress:
            progress(1.0, "cubes loaded")

    # ------------------------------------------------------------------ stage 2
    def suggest_sky_box(self) -> Dict[str, tuple]:
        self._need("sky", bool(self.arms))
        out = {}
        size = int(self.cfg["sky_box_size"])
        # keep the box away from anything bright in EITHER arm: use the ref arm image masked by bright pixels
        for name, img in self.whitelight.items():
            avoid = None
            out[name] = nz.suggest_sky_box(img, self.arms[name].valid, size, avoid)
        return out

    def estimate_sky(self) -> None:
        self._need("sky", bool(self.arms))
        cfg = self.cfg
        t0 = time.time()
        boxes: Dict[str, tuple] = {}
        if cfg["sky_box_mode"] == "auto":
            boxes = self.suggest_sky_box()
            # the same sky patch must be source-free in both arms: use the ref-arm box for both when registered 1:1
            if len(self.arms) == 2 and self.registration.mode == "identity":
                boxes["red"] = boxes[self.ref_arm]
        else:
            txt_b, txt_r = cfg.get("sky_box_blue", ""), cfg.get("sky_box_red", "")
            for name in self.arms:
                txt = txt_b if name == "blue" else (txt_r or txt_b)
                if not txt:
                    raise PipelineError(f"sky_box_mode='manual' but no box was given for the {name} arm")
                boxes[name] = nz.parse_box(txt)
        for name, arm in self.arms.items():
            self.sky[name] = nz.estimate_sky(arm, boxes[name], cfg["noise_source"], "median",
                                             bool(cfg["correct_noise_correlation"]))
            s = self.sky[name]
            self.log(f"{name} sky box (y0,y1,x0,x1) = {s.box}: {s.n_pix} spaxels; noise from {s.source}; "
                     f"correlation exponent alpha = {s.alpha:.2f} (0.5 = independent), "
                     f"calibration x{s.ratio1:.2f}  ->  a 20-spaxel sum is x{s.correlation_factor(20):.2f} noisier than sqrt(N)")
            med_bg = float(np.nanmedian(s.bg))
            med_nz = float(np.nanmedian(s.noise))
            self.log(f"    background median {med_bg:.3g}, per-spaxel noise median {med_nz:.3g}")
            # sanity: box must be source-free
            wl = self.whitelight[name]
            y0, y1, x0, x1 = s.box
            sub = wl[y0:y1, x0:x1]
            if np.isfinite(sub).any():
                thresh = np.nanmedian(wl) + 3 * 1.4826 * np.nanmedian(np.abs(wl[np.isfinite(wl)] - np.nanmedian(wl)))
                if np.nanmax(sub) > thresh * 3 and np.nanmax(sub) > 5 * np.nanmedian(np.abs(sub)):
                    self.log(f"    [warning] the {name} sky box contains a bright pixel -- check it is source-free.")
        self.timings["sky"] = time.time() - t0

    # ------------------------------------------------------------------ stage 3
    def set_apertures(self, aps: Sequence) -> None:
        out = []
        for a in aps:
            out.append(a if isinstance(a, Aperture) else Aperture.from_dict(a))
        ids = [a.id for a in out]
        if len(ids) != len(set(ids)):
            raise PipelineError("aperture ids must be unique")
        self.apertures = out
        self.cfg["apertures"] = [a.to_dict() for a in out]
        self.masks = {}

    def apertures_from_regions(self, path: str) -> List[str]:
        self._need("apertures", bool(self.arms))
        arm = self.arms[self.ref_arm]
        aps, warns = parse_ds9_regions(os.path.expanduser(path), arm.wcs, arm.pixscale[1])
        return aps, warns

    def suggest_apertures(self, nsigma: float = 5.0, max_peaks: int = 6) -> List[Aperture]:
        self._need("apertures", bool(self.arms))
        arm = self.arms[self.ref_arm]
        img = self.whitelight[self.ref_arm]
        noise = None
        if self.ref_arm in self.sky:
            y0, y1, x0, x1 = self.sky[self.ref_arm].box
            noise = float(1.4826 * np.nanmedian(np.abs(img[y0:y1, x0:x1] - np.nanmedian(img[y0:y1, x0:x1]))))
        return detect_peaks(img, arm.valid, noise, nsigma, max_peaks=max_peaks)

    def build_masks(self) -> None:
        self._need("apertures", bool(self.arms) and bool(self.apertures))
        self.masks = {}
        for name, arm in self.arms.items():
            self.masks[name] = {}
            for ap in self.apertures:
                m = aperture_mask(arm.flux.shape[1:], self.arm_aperture(ap, name))
                self.masks[name][ap.id] = m
        # sanity checks
        for ap in self.apertures:
            for name in self.arms:
                m = self.masks[name][ap.id]
                if m.sum() < 1.0:
                    self.log(f"[warning] aperture '{ap.id}' covers < 1 spaxel in the {name} arm")
                v = self.arms[name].valid
                frac_valid = float((m * v).sum() / max(m.sum(), 1e-9))
                if frac_valid < 0.8:
                    self.log(f"[warning] {100*(1-frac_valid):.0f}% of aperture '{ap.id}' lies outside the {name} IFU footprint")
        if "blue" in self.sky:
            for name in self.arms:
                y0, y1, x0, x1 = self.sky[name].box
                for ap in self.apertures:
                    mm = self.masks[name][ap.id][y0:y1, x0:x1]
                    if mm.sum() > 0:
                        self.log(f"[warning] aperture '{ap.id}' overlaps the {name} sky box -- background subtraction would "
                                 "remove part of its own light. Re-pick the sky box.")
        ids = [a.id for a in self.apertures]
        for i in range(len(ids)):
            for j in range(i + 1, len(ids)):
                f = overlap_fraction(self.masks[self.ref_arm][ids[i]], self.masks[self.ref_arm][ids[j]])
                if f > 0.05:
                    self.log(f"    apertures '{ids[i]}' and '{ids[j]}' overlap ({f:.0%} of the smaller one): "
                             + ("the deblend will separate them." if self.cfg["extraction_method"] == "deblend"
                                else "'sum' extraction does NOT correct this -- consider extraction_method='deblend'."))

    # ------------------------------------------------------------------ stage 4
    def _profile_band(self) -> Tuple[str, Tuple[float, float]]:
        txt = str(self.cfg.get("profile_band_A", "") or "").strip()
        if txt:
            lo, hi = parse_range(txt)
        else:
            arm = "red" if "red" in self.arms else "blue"
            w = self.arms[arm].wave
            lo, hi = float(np.percentile(w, 50)), float(w.max())
        # choose the arm with the largest overlap with the band
        best, bo = None, -1.0
        for name, arm in self.arms.items():
            ov = min(hi, arm.wave.max()) - max(lo, arm.wave.min())
            if ov > bo:
                best, bo = name, ov
        if bo <= 0:
            raise PipelineError(f"profile band {lo:.0f}-{hi:.0f} A lies outside both cubes")
        return best, (lo, hi)

    def extract(self, progress: ProgressCB = None) -> None:
        self._need("extraction", bool(self.sky) and bool(self.apertures))
        cfg = self.cfg
        t0 = time.time()
        self.build_masks()
        method = cfg["extraction_method"]
        self.spectra = {}
        self.extract_info = dict(method=method, flags=[], per_aperture={})
        sub_bg = bool(cfg["subtract_background"])
        corr = bool(cfg["correct_noise_correlation"])
        if method == "sum":
            for k, ap in enumerate(self.apertures):
                arms = {}
                for name, arm in self.arms.items():
                    wimg = self.whitelight[name] if cfg["weighting"] == "intensity" else None
                    if wimg is not None:
                        wimg = wimg - np.nanmedian(self.whitelight[name][self.arms[name].valid])
                    arms[name] = extract_sum(arm, self.sky[name], self.masks[name][ap.id], wimg, sub_bg, corr)
                self.spectra[ap.id] = ExtractedSpectrum(ap.id, ap.role, ap.group, arms, dict(method="sum"))
                if progress:
                    progress((k + 1) / len(self.apertures), f"extracted {ap.id}")
        elif method == "deblend":
            self._extract_deblend(progress)
        else:
            raise PipelineError(f"unknown extraction method {method!r}")
        self._record_counts()
        self._apply_manual_subtractions()
        # post-processing: spike rejection, quality numbers
        for sid, sp in self.spectra.items():
            if cfg.get("reject_spikes", True):
                n = sum(reject_spikes(a) for a in sp.arms.values())
                sp.meta["spikes_rejected"] = n
            sp.meta["snr"] = {a: s.median_snr() for a, s in sp.arms.items()}
            self.extract_info["per_aperture"][sid] = dict(snr=sp.meta["snr"],
                                                         npix={a: s.npix for a, s in sp.arms.items()},
                                                         counts=sp.meta.get("counts"))
        if cfg["combine_groups"]:
            self._combine_groups()
        for sid, sp in self.spectra.items():
            msg = ", ".join(f"{a}: S/N/pix={s:.2f}" for a, s in sp.meta["snr"].items())
            cnt = sp.meta.get("counts") or {a: dict(area=s.npix, touched=int(round(s.npix))) for a, s in sp.arms.items()}
            spax = ", ".join(f"{a}: {c['area']:.1f} spaxels ({c['touched']} touched)" for a, c in cnt.items())
            sub = ""
            for r in sp.meta.get("subtracted", []):
                sub += f"  [minus {r['scale']:.3f} x {r['contaminator']}{' (auto)' if r['auto'] else ''}]"
            self.log(f"  {sid:12s} [{sp.role}]  {spax}  |  {msg}"
                     + (f"  (spikes masked: {sp.meta.get('spikes_rejected', 0)})" if sp.meta.get("spikes_rejected") else "") + sub)
        self.timings["extract"] = time.time() - t0
        self.engine = None

    def _record_counts(self) -> None:
        """Number of spaxels behind every extracted spectrum (area = sum of fractional coverage;
        'touched' = spaxels with any coverage inside the IFU footprint)."""
        for ap in self.apertures:
            sp = self.spectra.get(ap.id)
            if sp is None:
                continue
            cnt = {}
            for name, arm in self.arms.items():
                m = self.masks[name][ap.id]
                cnt[name] = dict(area=float(m.sum()), touched=int(((m > 0.01) & arm.valid).sum()),
                                 inside=float((m * arm.valid).sum()))
            sp.meta["counts"] = cnt

    def _estimate_leakage(self, con_id: str, tgt_id: str) -> Tuple[float, str]:
        """Fraction of the contaminator's aperture flux that falls inside the target aperture, from the
        contaminator's own spatial profile (PSF-convolved Sersic) fitted to a white-light image with
        every OTHER aperture's pixels excluded."""
        aps = {a.id: a for a in self.apertures}
        frame_arm, band = self._profile_band()
        arm = self.arms[frame_arm]
        img = white_light(arm, band[0], band[1])
        y0, y1, x0, x1 = self.sky[frame_arm].box
        sub = img[y0:y1, x0:x1]
        sigma_img = float(1.4826 * np.nanmedian(np.abs(sub - np.nanmedian(sub)))) or 1e-6
        img = img - np.nanmedian(sub)
        valid = arm.valid.copy()
        for k, m in self.masks[frame_arm].items():
            if k != con_id:
                valid &= ~(m > 0.01)              # other regions' light must not drive the profile
        con_fr = self.arm_aperture(aps[con_id], frame_arm)
        con_fr = Aperture.from_dict({**con_fr.to_dict(), "profile": "sersic"})
        psf_sigma = float(self.cfg["seeing_fwhm_arcsec"]) / arm.pixscale[1] * FWHM2SIG
        pfit = fit_profiles(img, valid, [con_fr], psf_sigma, sigma_img, frame_arm, int(self.cfg["deblend_pad_pix"]))
        P = pfit.components[0].render(arm.flux.shape[1:])
        own = float((P * self.masks[frame_arm][con_id]).sum())
        leak = float((P * self.masks[frame_arm][tgt_id]).sum())
        if own <= 0:
            raise PipelineError(f"cannot estimate the leakage of '{con_id}' into '{tgt_id}': its model has no flux in its own aperture")
        p = pfit.components[0].params
        return leak / own, f"Sersic r_e={p['re']:.1f} pix, n={p['n']:.1f}, chi2_red={pfit.chi2red:.1f}"

    def _apply_manual_subtractions(self) -> None:
        rules = parse_subtraction_rules(self.cfg.get("manual_subtractions", ""))
        if not rules:
            return
        ids = [a.id for a in self.apertures]
        for tgt, con, _ in rules:
            for k in (tgt, con):
                if k not in ids:
                    raise PipelineError(f"manual subtraction refers to unknown aperture '{k}' (available: {ids})")
        raw = {k: sp.copy() for k, sp in self.spectra.items()}       # contaminators are used as originally extracted
        if self.cfg["extraction_method"] == "deblend":
            self.log("[note] manual subtraction on top of 'deblend': the overlap is already modelled, so this removes "
                     "only what is left over -- check the result.")
        for tgt, con, scale in rules:
            auto = scale is None
            how = ""
            if auto:
                scale, how = self._estimate_leakage(con, tgt)
            for arm in self.spectra[tgt].arms:
                subtract_spectrum(self.spectra[tgt].arms[arm], raw[con].arms[arm], float(scale))
            rec = dict(contaminator=con, scale=float(scale), auto=auto, how=how)
            self.spectra[tgt].meta.setdefault("subtracted", []).append(rec)
            self.log(f"  manual subtraction: '{tgt}' -= {scale:.3f} x '{con}'" + (f"  (auto: {how})" if auto else ""))
            if scale < 0 or scale > 1.5:
                self.log(f"[warning] subtraction scale {scale:.3f} is outside the physically expected 0-1 range")

    def _extract_deblend(self, progress: ProgressCB = None) -> None:
        cfg = self.cfg
        frame_arm, band = self._profile_band()
        arm = self.arms[frame_arm]
        img = white_light(arm, band[0], band[1])
        valid = arm.valid
        # noise of the band image: scatter of the sky box in that image
        y0, y1, x0, x1 = self.sky[frame_arm].box
        sub = img[y0:y1, x0:x1]
        sigma_img = float(1.4826 * np.nanmedian(np.abs(sub - np.nanmedian(sub))))
        if not np.isfinite(sigma_img) or sigma_img <= 0:
            sigma_img = float(np.nanstd(sub)) or 1e-6
        img_sub = img - np.nanmedian(sub)
        psf_sigma = float(cfg["seeing_fwhm_arcsec"]) / arm.pixscale[1] * FWHM2SIG
        aps_frame = [self.arm_aperture(ap, frame_arm) for ap in self.apertures]
        self.log(f"Fitting spatial profiles on the {band[0]:.0f}-{band[1]:.0f} A white-light image of the {frame_arm} arm "
                 f"(PSF sigma {psf_sigma:.2f} pix) ...")
        pfit = fit_profiles(img_sub, valid, aps_frame, psf_sigma, sigma_img, frame_arm, int(cfg["deblend_pad_pix"]),
                            bool(cfg["deblend_fit_background"]))
        self.profile_fit = pfit
        self.extract_info["profile_chi2red"] = pfit.chi2red
        self.extract_info["flags"] += pfit.flags
        for c in pfit.components:
            if c.kind == "sersic":
                p = c.params
                self.log(f"    {c.ap.id:10s} sersic: r_e={p['re']:.2f} pix ({p['re'] * arm.pixscale[1]:.2f}\"), n={p['n']:.2f}, "
                         f"q={p['q']:.2f}, PA={p['pa']:.0f}")
            else:
                self.log(f"    {c.ap.id:10s} {c.kind}")
        self.log(f"    profile fit chi2_red = {pfit.chi2red:.2f}" + ("" if pfit.success else "  [optimiser did not converge]"))
        for f in pfit.flags:
            self.log("    [warning] " + f)
        out: Dict[str, Dict[str, ArmSpectrum]] = {ap.id: {} for ap in self.apertures}
        diags = {}
        for name, a in self.arms.items():
            shape = a.flux.shape[1:]
            templates, shifts = [], []
            for ap, c in zip(self.apertures, pfit.components):
                a_arm, a_fr = self.arm_aperture(ap, name), self.arm_aperture(ap, frame_arm)
                dx, dy = a_arm.x - a_fr.x, a_arm.y - a_fr.y
                sc = arm.pixscale[1] / a.pixscale[1]
                shifts.append((dx, dy))
                templates.append(c.render(shape, dx, dy, sc))
            dx, dy = float(np.mean([t[0] for t in shifts])), float(np.mean([t[1] for t in shifts]))
            masks = [self.masks[name][ap.id] for ap in self.apertures]
            y_lo, y_hi, x_lo, x_hi = pfit.window
            region = np.zeros(shape, bool)
            ys0, ys1 = max(0, int(round(y_lo + dy)) - 0), min(shape[0], int(round(y_hi + dy)))
            xs0, xs1 = max(0, int(round(x_lo + dx)) - 0), min(shape[1], int(round(x_hi + dx)))
            region[ys0:ys1, xs0:xs1] = True
            specs, diag = deblend(a, self.sky[name], templates, masks, region, bool(cfg["subtract_background"]),
                                  bool(cfg["deblend_fit_background"]), bool(cfg["correct_noise_correlation"]))
            for ap, s in zip(self.apertures, specs):
                out[ap.id][name] = s
            diags[name] = diag
            self.log(f"    {name}: deblend of {len(templates)} components over {diag['n_pix']} spaxels, "
                     f"template condition number {diag['cond']:.1f}")
            if diag["cond"] > 1e4:
                self.extract_info["flags"].append(f"deblend matrix is ill-conditioned in the {name} arm ({diag['cond']:.0f})")
                self.log("    [warning] the spatial templates are nearly degenerate -- spectra may be noisy/anticorrelated")
            if progress:
                progress(0.5 + 0.5 * (list(self.arms).index(name) + 1) / len(self.arms), f"deblended {name}")
        self.extract_info["deblend"] = {k: dict(cond=v["cond"], n_pix=v["n_pix"]) for k, v in diags.items()}
        self.extract_info["resid_img"] = {k: v["resid_img"] for k, v in diags.items()}
        for ap in self.apertures:
            self.spectra[ap.id] = ExtractedSpectrum(ap.id, ap.role, ap.group, out[ap.id],
                                                    dict(method="deblend", profile=pfit.components[
                                                        [a.id for a in self.apertures].index(ap.id)].kind))

    def _combine_groups(self) -> None:
        groups: Dict[str, List[str]] = {}
        for sid, sp in list(self.spectra.items()):
            if sp.group and sp.role == "source":
                groups.setdefault(sp.group, []).append(sid)
        for g, members in groups.items():
            if len(members) < 2:
                continue
            arms = {}
            for name in self.arms:
                arms[name] = inverse_variance_combine([self.spectra[m].arms[name] for m in members])
            cnt = {a: dict(area=sum(self.spectra[m].meta["counts"][a]["area"] for m in members),
                           touched=sum(self.spectra[m].meta["counts"][a]["touched"] for m in members),
                           inside=sum(self.spectra[m].meta["counts"][a]["inside"] for m in members)) for a in self.arms
                   if all("counts" in self.spectra[m].meta for m in members)}
            sp = ExtractedSpectrum(f"group:{g}", "source", g, arms,
                                   dict(method="ivw-combination", members=members, counts=cnt or None,
                                        snr={a: s.median_snr() for a, s in arms.items()}))
            self.spectra[sp.id] = sp
            ind = ", ".join(f"{m}: " + "/".join(f"{self.spectra[m].meta['snr'][a]:.2f}" for a in self.arms) for m in members)
            self.log(f"  combined images {members} of source group '{g}' (S/N per pixel {ind} -> "
                     + "/".join(f"{sp.meta['snr'][a]:.2f}" for a in self.arms) + ")")

    # ------------------------------------------------------------------ stage 5/6
    def make_engine(self) -> RedshiftEngine:
        self._need("fit", bool(self.spectra))
        aps = {a.id: a for a in self.apertures}
        # group-combined spectra inherit the settings of their first member
        for sid, sp in self.spectra.items():
            if sid not in aps and sp.meta.get("members"):
                first = aps[sp.meta["members"][0]]
                aps[sid] = Aperture(id=sid, role="source", group=sp.group, mode=first.mode,
                                    z_known=first.z_known, z_min=first.z_min, z_max=first.z_max)
        self.engine = RedshiftEngine(self.cfg, self.spectra, aps, self.fwhm, list(self.arms), self.log)
        self.engine.results = self.results
        return self.engine

    def fit_all(self, progress: ProgressCB = None) -> Dict[str, ApertureResult]:
        t0 = time.time()
        eng = self.engine or self.make_engine()
        eng.results = self.results
        # members of a combined group are fitted too (independent confirmation); the group entry gets the best S/N
        eng.run_all(progress)
        self.results = eng.results
        self.summarize()
        self.timings["fit"] = time.time() - t0
        return self.results

    def scan(self, aid: str, progress: ProgressCB = None) -> List[ScanResult]:
        eng = self.engine or self.make_engine()
        return eng.scan(aid, progress=progress)

    def fit_aperture(self, aid: str, z_guess: Optional[float] = None, mode: Optional[str] = None, scans=None,
                     progress: ProgressCB = None) -> ApertureResult:
        eng = self.engine or self.make_engine()
        eng.results = self.results
        r = eng.fit(aid, z_guess, mode, scans, progress)
        self.results = eng.results
        self.summarize()
        return r

    def clear_fit(self, aid: str) -> bool:
        """Discard the recorded fit of one spectrum (in place, so the engine sees it too)."""
        found = self.results.pop(aid, None) is not None
        self.summarize()
        return found

    def summarize(self) -> None:
        if self.engine is None:
            return
        self.system_rows, self.system_messages = self.engine.system_checks()

    # ------------------------------------------------------------------ everything
    def run_all(self, progress: ProgressCB = None, save: bool = True) -> Dict[str, ApertureResult]:
        def sub(lo, hi):
            return (lambda f, m: progress(lo + (hi - lo) * f, m)) if progress else None
        self.load(sub(0.0, 0.1))
        self.estimate_sky()
        if not self.apertures:
            aps = self.cfg.get("apertures") or []
            if not aps and self.cfg.get("region_file"):
                new, warns = self.apertures_from_regions(self.cfg["region_file"])
                for w in warns:
                    self.log("[region file] " + w)
                aps = [a.to_dict() for a in new]
                self.log("[note] roles were not specified in the region file: every aperture defaults to 'source' unless "
                         "its label contains 'defl'/'lens'. Edit the config to set roles before trusting the results.")
            if not aps:
                raise PipelineError("No apertures defined. Add them to the config ('apertures') or give a region_file.")
            self.set_apertures(aps)
        self.extract(sub(0.1, 0.2))
        self.fit_all(sub(0.2, 0.95))
        if save:
            from .report import save_outputs
            save_outputs(self, progress=sub(0.95, 1.0))
        if progress:
            progress(1.0, "done")
        return self.results
