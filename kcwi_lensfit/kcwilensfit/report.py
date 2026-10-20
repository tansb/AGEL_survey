"""Writing results: spectra (FITS + ASCII), tables (CSV/JSON), figures and a Markdown report."""
from __future__ import annotations

import csv
import json
import os
import re
from typing import Callable, Optional

import numpy as np
from astropy.io import fits
from astropy.table import Table

from . import __version__
from . import plotting as P
from .schema import save_config


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def _jsonable(o):
    if isinstance(o, (np.floating,)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, float):
        return None if not np.isfinite(o) else o
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return [_jsonable(x) for x in o.tolist()]
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(x) for x in o]
    return o


def save_spectra_fits(pipe, path: str) -> None:
    hdus = [fits.PrimaryHDU()]
    hdus[0].header["CREATOR"] = f"kcwilensfit {__version__}"
    hdus[0].header["TARGET"] = pipe.cfg.get("target_name", "")
    hdus[0].header["COMMENT"] = "WAVE = air wavelength [Angstrom]; FLUX in cube units (flux per spaxel summed); VAR = variance"
    for sid, sp in pipe.spectra.items():
        for arm, s in sp.arms.items():
            cols = [fits.Column(name="WAVE", format="D", array=s.wave),
                    fits.Column(name="FLUX", format="D", array=np.nan_to_num(s.flux, nan=0.0)),
                    fits.Column(name="VAR", format="D", array=np.where(np.isfinite(s.var), s.var, -1.0)),
                    fits.Column(name="GOOD", format="L", array=s.good)]
            h = fits.BinTableHDU.from_columns(cols, name=f"{_safe(sid)}_{arm}".upper()[:60])
            h.header["APID"] = sid
            h.header["ARM"] = arm
            h.header["ROLE"] = sp.role
            h.header["GROUP"] = sp.group
            h.header["NPIX"] = float(s.npix)
            h.header["METHOD"] = str(sp.meta.get("method", ""))
            r = pipe.results.get(sid)
            if r is not None and np.isfinite(r.z):
                h.header["Z"] = float(r.z)
                h.header["Z_ERR"] = float(r.z_err)
            h.header["COMMENT"] = "VAR = -1 marks masked pixels"
            hdus.append(h)
    fits.HDUList(hdus).writeto(path, overwrite=True)


def save_ascii(pipe, folder: str) -> None:
    os.makedirs(folder, exist_ok=True)
    for sid, sp in pipe.spectra.items():
        rows = []
        for arm, s in sp.arms.items():
            t = Table(dict(arm=np.full(s.wave.size, arm), wave=s.wave, flux=s.flux, var=s.var, good=s.good))
            rows.append(t)
        from astropy.table import vstack
        t = vstack(rows)
        t.meta.update(dict(id=sid, role=sp.role, group=sp.group))
        t.write(os.path.join(folder, f"{_safe(sid)}.ecsv"), format="ascii.ecsv", overwrite=True)


def results_rows(pipe):
    return [r.row() for r in pipe.results.values()]


def write_csv(rows, path):
    if not rows:
        return
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if (isinstance(v, float) and not np.isfinite(v)) else v) for k, v in r.items()})


def markdown_report(pipe) -> str:
    cfg = pipe.cfg
    L = []
    L.append(f"# KCWI lens redshift report — {cfg.get('target_name', '')}\n")
    L.append(f"*kcwilensfit {__version__}*\n")
    L.append("## Data and set-up\n")
    for name, a in pipe.arms.items():
        L.append(f"- **{name}**: `{os.path.basename(a.path)}`, {a.wave.min():.0f}–{a.wave.max():.0f} Å (air), "
                 f"native axis {a.native_frame}, instrumental FWHM {a.fwhm_A:.2f} Å, {a.pixscale[1]:.3f}\"/spaxel")
        for line in a.info:
            L.append(f"    - {line}")
    r = pipe.registration
    L.append(f"- registration: **{r.mode}** (dx={r.dx:+.2f}, dy={r.dy:+.2f} pix) — {r.note}")
    for name, s in pipe.sky.items():
        L.append(f"- sky box ({name}) {s.box}: noise α={s.alpha:.2f}")
    L.append(f"- extraction: **{pipe.extract_info.get('method')}**; "
             f"deflector profile χ²ᵣ={pipe.extract_info.get('profile_chi2red', float('nan')):.2f}")
    for f in pipe.extract_info.get("flags", []):
        L.append(f"    - ⚠ {f}")
    L.append("\n## Apertures\n")
    L.append("| id | role | group | shape | x | y | a | b | angle |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for ap in pipe.apertures:
        L.append(f"| {ap.id} | {ap.role} | {ap.group} | {ap.shape} | {ap.x:.2f} | {ap.y:.2f} | {ap.a:.2f} | {ap.b:.2f} | {ap.angle:.0f} |")
    L.append("\n## Extraction\n")
    L.append("| spectrum | role | spaxels (area; touched) per arm | median S/N per pixel | manual subtractions |")
    L.append("|---|---|---|---|---|")
    for sid, sp in pipe.spectra.items():
        cnt = sp.meta.get("counts") or {}
        sp_txt = "; ".join(f"{a}: {c['area']:.1f} ({c['touched']})" for a, c in cnt.items()) or "—"
        sn = "; ".join(f"{a}: {v:.2f}" for a, v in sp.meta.get("snr", {}).items())
        sub = "; ".join(f"− {r['scale']:.3f}×{r['contaminator']}" + (" (auto)" if r["auto"] else "")
                        for r in sp.meta.get("subtracted", [])) or "—"
        L.append(f"| {sid} | {sp.role} | {sp_txt} | {sn} | {sub} |")
    L.append("\n## Redshifts\n")
    L.append("| aperture | role | status | method | z | σ_z (stat) | σ_z (sys) | σ_* [km/s] |")
    L.append("|---|---|---|---|---|---|---|---|")
    for k, rr in pipe.results.items():
        sig = (f"{rr.sigma:.0f} ± {rr.sigma_err:.0f} (stat)" + (f" ± {rr.sigma_sys:.0f} (sys)" if np.isfinite(rr.sigma_sys) else "")
               + ("" if rr.kinematics_ok else " — *unreliable*")) if np.isfinite(rr.sigma) else "—"
        L.append(f"| {k} | {rr.role} | {rr.status} | {rr.method} | {rr.z:.5f} | {rr.z_err:.5f} | {rr.z_sys_err:.5f} | {sig} |")
    L.append("\nAll redshifts are in the frame of the cube wavelengths (heliocentric if the DRP applied its VCORR).\n")
    if pipe.system_messages:
        L.append("## Lens-system consistency\n")
        for m in pipe.system_messages:
            L.append(f"- {m}")
    L.append("\n## Per-aperture notes\n")
    for k, rr in pipe.results.items():
        L.append(f"### {k}\n")
        for f in rr.flags:
            L.append(f"- ⚠ {f}")
        for n in rr.notes:
            L.append(f"- {n}")
        d = rr.detail
        if d.get("provenance_text"):
            L.append("```")
            L += d["provenance_text"]
            L.append("```")
        if d.get("lines"):
            L.append("\n| line | kind | S/N | z_line |\n|---|---|---|---|")
            for ln in d["lines"]:
                if ln["snr"] >= 2.0:
                    L.append(f"| {ln['name']} | {ln['kind']} | {ln['snr']:.1f} | {ln['z_line']:.5f} |")
        if d.get("identification_table"):
            L.append("\nCompeting identifications of the single detected line:\n")
            L.append("| if the line is | z | other lines covered | detected |\n|---|---|---|---|")
            for it in d["identification_table"]:
                L.append(f"| {it['line']} | {it['z']:.4f} | {it['n_other_covered']} | {it['n_other_detected']} |")
        if d.get("sweep"):
            L.append("\nRobustness sweep:\n")
            for row in d["sweep"]:
                if "error" in row:
                    L.append(f"- {row['label']}: failed ({row['error']})")
                else:
                    L.append(f"- {row['label']}: z={row['z']:.5f}, σ={row['sigma']:.1f} km/s, χ²ᵣ={row['chi2']:.2f}")
        if d.get("arm_check"):
            L.append("\nBlue-only vs red-only: " + "; ".join(
                f"{a}: z={v.get('z', float('nan')):.5f}±{v.get('z_err', float('nan')):.5f}" for a, v in d["arm_check"].items()))
        L.append("")
    L.append("## Method summary\n")
    L.append("See `docs/science_notes.md` in the repository for the physical rationale of every default.")
    return "\n".join(L)


def quick_summary(pipe, out_dir: str = "", saved: bool = False) -> str:
    """Short plain-text overview of the current state/results, for review BEFORE saving."""
    L = []
    L.append(f"Target: {pipe.cfg.get('target_name', '')}")
    if out_dir:
        exists = os.path.isdir(out_dir)
        L.append(f"Output directory: {out_dir}" + ("   [saved here]" if saved else
                 ("   (exists; files with the same names will be overwritten)" if exists else "   (will be created on save)")))
    arms = ", ".join(pipe.arms) or "none loaded"
    L.append(f"Arms: {arms}   Apertures: {len(pipe.apertures)}   Spectra: {len(pipe.spectra)}   "
             f"Fitted: {len(pipe.results)}")
    if not pipe.spectra:
        L.append("\nNothing extracted yet -- load cubes, define apertures and run the extraction first.")
        return "\n".join(L)
    unfitted = [k for k in pipe.spectra if k not in pipe.results]
    L.append("")
    if pipe.results:
        L.append(f"{'region':<14}{'role':<10}{'status':<11}{'z':>10}{'+/-':>10}  {'sigma [km/s]':<18}flags")
        for k, r in pipe.results.items():
            sig = f"{r.sigma:.0f} +/- {r.sigma_err:.0f}" + ("" if r.kinematics_ok else " (?)") if np.isfinite(r.sigma) else "-"
            z = f"{r.z:.5f}" if np.isfinite(r.z) else "-"
            ze = f"{r.z_err:.5f}" if np.isfinite(r.z_err) else "-"
            L.append(f"{k:<14}{r.role:<10}{r.status:<11}{z:>10}{ze:>10}  {sig:<18}{'; '.join(r.flags)}")
    else:
        L.append("No fits recorded yet.")
    if unfitted:
        L.append("\nNo recorded fit (left out of the redshift table and the saved spectra plot): " + ", ".join(unfitted))
    if pipe.system_messages:
        L.append("\nLens-system consistency:")
        L += [f"  * {m}" for m in pipe.system_messages]
    return "\n".join(L)


def writable_output_dir(requested: str, target: str = "", log: Callable[[str], None] = print) -> str:
    """Return a directory we can really write to.

    The requested folder is tested by creating (and deleting) a small file. Read-only data folders
    (common when the output folder was accidentally set to the raw-data directory) and
    macOS-protected locations fail this test; the outputs then go to ``~/kcwilensfit_output/<target>``
    and the log says so, instead of the run dying after the (expensive) fits.
    """
    req = os.path.abspath(os.path.expanduser(requested or "kcwilensfit_output"))

    def ok(d):
        try:
            os.makedirs(d, exist_ok=True)
            probe = os.path.join(d, ".kcwilensfit_write_test")
            with open(probe, "w") as fh:
                fh.write("x")
            os.remove(probe)
            return True
        except OSError:
            return False

    if ok(req):
        return req
    fb = os.path.join(os.path.expanduser("~"), "kcwilensfit_output", _safe(target) if target else "run")
    if not ok(fb):
        raise PermissionError(f"Cannot write to '{req}' nor to the fallback '{fb}'. Choose a writable output directory "
                              "(stage 7 or stage 1).")
    log(f"[warning] output directory '{req}' is not writable (read-only data folder or protected location); "
        f"writing to '{fb}' instead. Change 'Output directory' to avoid this message.")
    return fb


def save_outputs(pipe, progress: Optional[Callable[[float, str], None]] = None) -> str:
    out = writable_output_dir(pipe.cfg.get("output_dir"), pipe.cfg.get("target_name", ""), pipe.log)
    pipe.cfg["output_dir"] = out
    cfg_out = dict(pipe.cfg)
    cfg_out["apertures"] = [a.to_dict() for a in pipe.apertures]
    save_config(cfg_out, os.path.join(out, "config_used.json"))
    if pipe.cfg.get("save_fits_bundle", True) and pipe.spectra:
        save_spectra_fits(pipe, os.path.join(out, "spectra.fits"))
    if pipe.cfg.get("save_ascii_spectra", True) and pipe.spectra:
        save_ascii(pipe, os.path.join(out, "spectra"))
    rows = results_rows(pipe)
    write_csv(rows, os.path.join(out, "redshifts.csv"))
    with open(os.path.join(out, "results.json"), "w") as fh:
        json.dump(_jsonable(dict(
            version=__version__, results={k: dict(**r.row(), notes=r.notes, detail=r.detail) for k, r in pipe.results.items()},
            system=dict(ladder=pipe.system_rows, messages=pipe.system_messages))), fh, indent=2, default=str)
    if pipe.cfg.get("save_plots", True):
        fdir = os.path.join(out, "figures")
        os.makedirs(fdir, exist_ok=True)
        save_figures(pipe, fdir, progress)
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write(markdown_report(pipe))
    return out


def save_figures(pipe, fdir: str, progress=None) -> None:
    def save(fig, name):
        if fig is not None:
            fig.savefig(os.path.join(fdir, name), dpi=130)

    if pipe.arms and pipe.whitelight:
        save(P.fig_apertures(pipe) if pipe.apertures else None, "01_apertures.png")
    if pipe.sky:
        save(P.fig_sky(pipe), "02_sky_noise.png")
    save(P.fig_deblend(pipe), "03_deblend.png")
    if pipe.spectra:
        sm = float(pipe.cfg.get("display_smooth_A", 0.0))
        win = None
        if pipe.engine is not None:
            win = pipe.engine.static_windows()
        # only spectra with a recorded fit, excluding regions merged into a combined spectrum
        hide = P.combined_member_ids(pipe)
        fitted = [a for a in pipe.spectra if a in pipe.results and a not in hide]
        if fitted:
            # raw data: sky/telluric windows stay shaded for reference but nothing is blanked or despiked
            save(P.fig_spectra(pipe, fitted, smooth_A=sm, windows=win, despike=False, blank_windows=False),
                 "04_spectra.png")
    k = 0
    for aid, r in pipe.results.items():
        safe = _safe(aid)
        for sc in r.artifacts.get("scans", []):
            save(P.fig_scan([sc]), f"05_scan_{safe}_{sc.kind}.png")
        save(P.fig_continuum_fit(r), f"06_fit_{safe}_continuum.png")
        save(P.fig_lines(r), f"07_lines_{safe}.png")
        save(P.fig_bootstrap(r), f"08_bootstrap_{safe}.png")
        k += 1
        if progress:
            progress(k / max(len(pipe.results), 1), f"figures: {aid}")
    save(P.fig_ladder(pipe.system_rows), "09_redshift_ladder.png")
