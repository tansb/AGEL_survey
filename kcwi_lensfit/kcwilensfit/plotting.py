"""Figures.  Every function returns a ``matplotlib.figure.Figure`` built WITHOUT pyplot, so the
same objects can be saved by the CLI or embedded in the Tk GUI without global-state clashes."""
from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Circle, Ellipse, Rectangle

from .lines import expected_lines
from .masks import in_windows

ROLE_COLOR = {"deflector": "#00bcd4", "source": "#7cfc00", "other": "#ffb300"}


def _scale(img, lo=5, hi=99.5):
    f = img[np.isfinite(img)]
    if f.size == 0:
        return 0, 1
    return np.percentile(f, lo), np.percentile(f, hi)


def draw_aperture(ax, ap, color=None, label=True, lw=1.8):
    """Draw an Aperture (already in this axes' pixel frame)."""
    color = color or ROLE_COLOR.get(ap.role, "w")
    if ap.shape == "circle":
        ax.add_patch(Circle((ap.x, ap.y), ap.a, fill=False, ec=color, lw=lw))
        top = ap.a
    elif ap.shape == "ellipse":
        ax.add_patch(Ellipse((ap.x, ap.y), 2 * ap.a, 2 * ap.b, angle=ap.angle, fill=False, ec=color, lw=lw))
        top = max(ap.a, ap.b)
    else:
        ax.add_patch(Rectangle((ap.x - ap.a / 2, ap.y - ap.b / 2), ap.a, ap.b, angle=ap.angle,
                               rotation_point="center", fill=False, ec=color, lw=lw))
        top = max(ap.a, ap.b) / 2
    if label:
        ax.text(ap.x, ap.y + top + 0.8, ap.id, color=color, fontsize=8, ha="center", va="bottom")


def fig_apertures(pipe, title="Apertures as registered on each cube") -> Figure:
    """White-light image of each arm with the apertures EXACTLY as they will be extracted from it
    (sky-defined regions mapped through each cube's own WCS; others via the blue/red registration),
    the sky box, and the red image as contours on the blue panel for a registration check."""
    arms = list(pipe.arms)
    n = len(arms) + (1 if len(arms) == 2 else 0)
    fig = Figure(figsize=(5.2 * n, 5.2), constrained_layout=True)
    axes = fig.subplots(1, n, squeeze=False)[0]
    for ax, name in zip(axes, arms):
        img = pipe.whitelight[name]
        lo, hi = _scale(img)
        ax.imshow(img, origin="lower", cmap="inferno", vmin=lo, vmax=hi)
        for ap in pipe.apertures:
            draw_aperture(ax, pipe.arm_aperture(ap, name))
        if name in pipe.sky:
            y0, y1, x0, x1 = pipe.sky[name].box
            ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0, fill=False, ec="white", ls="--", lw=1.2))
            ax.text(x0, y1 + 0.3, "sky", color="white", fontsize=8)
        ys, xs = np.nonzero(pipe.arms[name].valid)
        ax.set_xlim(xs.min() - 3, xs.max() + 3)
        ax.set_ylim(ys.min() - 3, ys.max() + 3)
        ax.set_title(f"{name} white light (x,y = this cube's pixels)", fontsize=9)
        ax.set_xlabel("x [pix]"); ax.set_ylabel("y [pix]")
    if len(arms) == 2:
        ax = axes[2]
        b, r = pipe.whitelight["blue"], pipe.whitelight["red"]
        lo, hi = _scale(b)
        ax.imshow(b, origin="lower", cmap="gray", vmin=lo, vmax=hi)
        reg = pipe.registration
        rr = np.nan_to_num(r)
        from scipy.ndimage import shift as ndshift
        rs = ndshift(rr, (-reg.dy, -reg.dx), order=1)       # red resampled into the blue frame
        lv = np.nanpercentile(r[np.isfinite(r)], [70, 85, 95, 99])
        ax.contour(rs, levels=lv, colors="red", linewidths=0.8)
        for ap in pipe.apertures:
            draw_aperture(ax, pipe.arm_aperture(ap, "blue"), color="cyan", label=False, lw=1.0)
            a_red = pipe.arm_aperture(ap, "red")
            # red-arm aperture expressed in the blue frame (undo the registration shift)
            from dataclasses import replace
            draw_aperture(ax, replace(a_red, x=a_red.x - reg.dx, y=a_red.y - reg.dy), color="yellow", label=False, lw=1.0)
        ys, xs = np.nonzero(pipe.arms["blue"].valid)
        ax.set_xlim(xs.min() - 3, xs.max() + 3); ax.set_ylim(ys.min() - 3, ys.max() + 3)
        ax.set_title("registration check: blue (gray), red contours,\nblue apertures (cyan) vs red apertures (yellow)", fontsize=9)
    fig.suptitle(title, fontsize=10)
    return fig


def fig_sky_box(pipe, boxes: Dict[str, tuple], title: str = "Sky box on the field of view") -> Figure:
    """White-light image of each arm with the (proposed or estimated) sky box drawn on the FoV.

    The box must contain no source light: the title reports its brightest pixel relative to the
    field so an obviously contaminated box is visible at a glance.
    """
    arms = list(pipe.arms)
    fig = Figure(figsize=(5.4 * len(arms), 5.4), constrained_layout=True)
    axes = fig.subplots(1, len(arms), squeeze=False)[0]
    for ax, name in zip(axes, arms):
        img = pipe.whitelight[name]
        lo, hi = _scale(img, 5, 99.5)
        ax.imshow(img, origin="lower", cmap="inferno", vmin=lo, vmax=hi)
        box = boxes.get(name)
        sub = ""
        if box is not None:
            y0, y1, x0, x1 = box
            ax.add_patch(Rectangle((x0 - 0.5, y0 - 0.5), x1 - x0, y1 - y0, fill=False, ec="cyan", lw=2))
            ax.text(x0, y1 + 0.5, f"sky {x1 - x0}x{y1 - y0}", color="cyan", fontsize=9)
            cut = img[y0:y1, x0:x1]
            if np.isfinite(cut).any():
                bright = np.nanmax(cut)
                typical = np.nanmedian(img[np.isfinite(img)])
                sub = f"  (box max = {100 * (bright - typical) / max(hi - typical, 1e-30):.0f}% of the display range above the field median)"
        for ap in pipe.apertures:
            draw_aperture(ax, pipe.arm_aperture(ap, name), label=False, lw=1.0)
        ys, xs = np.nonzero(pipe.arms[name].valid)
        ax.set_xlim(xs.min() - 3, xs.max() + 3)
        ax.set_ylim(ys.min() - 3, ys.max() + 3)
        ax.set_title(f"{name} white light{sub}", fontsize=9)
        ax.set_xlabel("x [pix]"); ax.set_ylabel("y [pix]")
    fig.suptitle(title, fontsize=10)
    return fig


def fig_sky(pipe) -> Figure:
    arms = list(pipe.sky)
    fig = Figure(figsize=(11, 2.6 * len(arms) + 1.5), constrained_layout=True)
    gs = fig.add_gridspec(len(arms) + 1, 2, height_ratios=[1] * len(arms) + [0.9])
    for i, name in enumerate(arms):
        s, arm = pipe.sky[name], pipe.arms[name]
        ax = fig.add_subplot(gs[i, 0])
        ax.plot(arm.wave, s.bg, lw=0.6, color="steelblue", label="background (median)")
        ax.plot(arm.wave, s.noise, lw=0.6, color="darkorange", label="noise (per spaxel)")
        ax.axhline(0, color="k", lw=0.4)
        ax.set_title(f"{name}: sky box {s.box}  ({s.n_pix} spaxels)", fontsize=9)
        ax.legend(fontsize=7)
        ax2 = fig.add_subplot(gs[i, 1])
        d = s.scaling or {}
        if "npix" in d:
            ax2.loglog(d["npix"], np.sqrt(d["var_ratio"] * d["npix"]) , "o-", label="measured σ_N/σ_1·√N… (norm.)")
            nn = np.array(d["npix"], float)
            ax2.loglog(nn, np.sqrt(nn), "k--", label="independent: √N")
            ax2.set_xlabel("N spaxels in block"); ax2.set_ylabel("σ_sum / σ_pix")
            ax2.set_title(f"noise correlation: α={s.alpha:.2f}", fontsize=9)
            ax2.legend(fontsize=7)
        else:
            ax2.axis("off")
    axl = fig.add_subplot(gs[len(arms), :])
    axl.axis("off")
    axl.text(0.0, 0.5, "Sky box must contain no source light. Noise correlation α > 0.5 means neighbouring spaxels are "
             "correlated (resampling): errors of summed apertures are inflated by N^(α-0.5).", fontsize=8, wrap=True)
    return fig


def fig_deblend(pipe) -> Optional[Figure]:
    pf = pipe.profile_fit
    if pf is None:
        return None
    arm = pipe.arms[pf.components[0].frame_arm]
    from .cubeio import white_light
    band = pipe._profile_band()[1]
    img = white_light(arm, band[0], band[1])
    img = img - np.nanmedian(img)
    fig = Figure(figsize=(14, 4.4), constrained_layout=True)
    ax = fig.subplots(1, 4)
    y0, y1, x0, x1 = pf.window
    lo, hi = _scale(img[y0:y1, x0:x1], 2, 99.8)
    ext = (x0 - 0.5, x1 - 0.5, y0 - 0.5, y1 - 0.5)
    ax[0].imshow(img[y0:y1, x0:x1], origin="lower", extent=ext, cmap="inferno", vmin=lo, vmax=hi); ax[0].set_title("data")
    ax[1].imshow(pf.model[y0:y1, x0:x1], origin="lower", extent=ext, cmap="inferno", vmin=lo, vmax=hi); ax[1].set_title("spatial model")
    res = (img - pf.model)[y0:y1, x0:x1]
    m = np.nanpercentile(np.abs(res), 99)
    ax[2].imshow(res, origin="lower", extent=ext, cmap="RdBu_r", vmin=-m, vmax=m); ax[2].set_title("data - model")
    for a in ax[:3]:
        for ap in pipe.apertures:
            draw_aperture(a, pipe.arm_aperture(ap, pf.components[0].frame_arm), label=False, lw=1.0)
    ri = pipe.extract_info.get("resid_img", {}).get(pf.components[0].frame_arm)
    if ri is not None:
        r2 = ri[y0:y1, x0:x1]
        m2 = np.nanpercentile(np.abs(r2), 99) if np.isfinite(r2).any() else 1
        ax[3].imshow(r2, origin="lower", extent=ext, cmap="RdBu_r", vmin=-m2, vmax=m2)
        ax[3].set_title("median spectral residual\n(after per-λ deblend)", fontsize=9)
    fig.suptitle(f"Spatial deblending model ({band[0]:.0f}-{band[1]:.0f} Å white light, χ²ᵣ={pf.chi2red:.2f})", fontsize=10)
    return fig


def _line_overlay(ax, z, wmin, wmax, emission=True, absorption=True, fs=6):
    from .lines import LINES
    names_done = set()
    for ln, lam in expected_lines(z, wmin, wmax):
        if ln.is_emission and not emission:
            continue
        if (not ln.is_emission) and not absorption:
            continue
        c = "seagreen" if ln.is_emission else "firebrick"
        ax.axvline(lam, color=c, ls=":" if ln.is_emission else "--", lw=0.6, alpha=0.7)
        ax.text(lam, 0.98 if not ln.is_emission else 0.02, ln.name, color=c, fontsize=fs, rotation=90,
                ha="right", va="top" if not ln.is_emission else "bottom", transform=ax.get_xaxis_transform())


def _clean_for_display(arm, windows=None, despike=True, width=11, nsig=5.0):
    """Plot-only copy of an arm spectrum with known-bad pixels set to NaN.

    Blanked: pixels inside ``windows`` (sky lines / telluric / light pollution / user masks) and,
    if ``despike``, narrow outliers (flux or error spikes a few pixels wide, e.g. the instrumental
    spike near the red end of the blue arm).  The median filter is wider than any resolved line,
    so real emission/absorption features are kept.  The data used for fitting are not modified.
    Returns (cleaned ArmSpectrum, bool array of blanked pixels).
    """
    from scipy.ndimage import binary_dilation, median_filter
    d = arm.copy()
    bad = np.zeros(d.wave.shape, bool)
    if windows:
        bad |= in_windows(d.wave, windows)
    if despike and d.wave.size > 3 * width:
        f = np.where(np.isfinite(d.flux), d.flux, np.nan)
        fz = np.where(np.isfinite(f), f, np.nanmedian(f) if np.isfinite(f).any() else 0.0)
        res = fz - median_filter(fz, size=width, mode="nearest")
        sig = 1.4826 * np.nanmedian(np.abs(res - np.nanmedian(res)))
        spike = np.abs(res) > nsig * max(sig, 1e-30)
        err = np.sqrt(np.where(np.isfinite(d.var) & (d.var > 0), d.var, np.nan))
        if np.isfinite(err).any():
            ez = np.where(np.isfinite(err), err, np.nanmedian(err))
            spike |= ez > 4.0 * median_filter(ez, size=5 * width, mode="nearest")
        bad |= binary_dilation(spike, iterations=1)
    d.flux = np.where(bad, np.nan, d.flux)
    d.var = np.where(bad, np.nan, d.var)
    return d, bad


def combined_member_ids(pipe) -> set:
    """IDs of spectra that were merged into a combined (group) spectrum."""
    return {m for sp in pipe.spectra.values() for m in (sp.meta.get("members") or []) if m in pipe.spectra}


def fig_spectra(pipe, aids: Optional[Sequence[str]] = None, smooth_A: float = 0.0, z_marks: Optional[Dict[str, float]] = None,
                windows=None, despike: bool = True, blank_windows: bool = True) -> Figure:
    from .extraction import smooth_for_display
    if aids is None:
        hide = combined_member_ids(pipe)
        aids = [a for a in pipe.spectra if a not in hide]
    aids = list(aids)
    fig = Figure(figsize=(13, 2.4 * len(aids) + 0.8), constrained_layout=True)
    axes = fig.subplots(len(aids), 1, squeeze=False, sharex=True)[:, 0]
    xlims = []
    for ax, aid in zip(axes, aids):
        sp = pipe.spectra[aid]
        # x-range from the arms actually present, so a single-arm spectrum is not squished by
        # windows (e.g. sky lines) that lie in the missing arm
        wall = np.concatenate([s.wave[s.good] if s.good.any() else s.wave for s in sp.arms.values()])
        xlo, xhi = float(wall.min()), float(wall.max())
        xlims.append((xlo, xhi))
        for lo, hi in (windows or []):
            if hi < xlo or lo > xhi:
                continue
            ax.axvspan(max(lo, xlo), min(hi, xhi), color="0.88", zorder=0)
        clean = {}
        for arm, s in sp.arms.items():
            s, _ = _clean_for_display(s, windows if blank_windows else None, despike)
            clean[arm] = s
            d = smooth_for_display(s, smooth_A)
            ax.plot(d.wave, d.flux, lw=0.6, color="tab:blue" if arm == "blue" else "tab:red", label=arm)
            g = d.good
            if g.any():
                ax.plot(d.wave[g], np.sqrt(d.var[g]), lw=0.4, color="0.5")
        ax.axhline(0, color="k", lw=0.4)
        r = pipe.results.get(aid)
        z = None
        if z_marks and aid in z_marks:
            z = z_marks[aid]
        elif r is not None and np.isfinite(r.z):
            z = r.z
        if z is not None:
            w = np.concatenate([s.wave for s in sp.arms.values()])
            _line_overlay(ax, z, w.min(), w.max())
        ttl = f"{aid} [{sp.role}]"
        if r is not None and np.isfinite(r.z):
            ttl += f"   z = {r.z:.5f} ± {r.z_err:.5f}  ({r.status}, {r.method})"
        ax.set_title(ttl, fontsize=9, loc="left")
        # robust y-limits
        allf = np.concatenate([s.flux[s.good] for s in clean.values() if s.good.any()] or [np.zeros(0)])
        if allf.size:
            lo, hi = np.percentile(allf, [0.5, 99.5])
            pad = 0.3 * (hi - lo)
            ax.set_ylim(lo - pad, hi + pad)
        ax.set_ylabel("flux")
    pad = 0.01 * (max(h for _, h in xlims) - min(l for l, _ in xlims))
    axes[0].set_xlim(min(l for l, _ in xlims) - pad, max(h for _, h in xlims) + pad)
    axes[-1].set_xlabel("observed wavelength (air) [Å]")
    axes[0].legend(fontsize=7, loc="upper right")
    return fig


def fig_scan(scans) -> Figure:
    fig = Figure(figsize=(11, 3.0 * len(scans)), constrained_layout=True)
    axes = fig.subplots(len(scans), 1, squeeze=False)[:, 0]
    for ax, sc in zip(axes, scans):
        curves = sc.extra.get("curves") if sc.kind == "continuum" else None
        if curves and len(curves) > 1:
            for k, (sg, ch) in enumerate(curves):
                ax.plot(sc.z, ch, lw=0.8, label=f"sigma fixed = {sg:.0f} km/s")
            ax.legend(fontsize=7)
        else:
            ax.plot(sc.z, sc.value, lw=0.8, color="k")
        for i, c in enumerate(sc.candidates[:6]):
            ax.axvline(c.z, color="C3" if not c.below_deflector else "0.6", lw=0.8, ls="--")
            ax.text(c.z, 0.97 - 0.07 * (i % 5), f" {i + 1}", color="C3", transform=ax.get_xaxis_transform(), fontsize=8)
        if sc.kind == "lines":
            ax.set_ylabel("line score S(z)")
            ax.set_title(f"{sc.aid}: blind line search (null-subtracted matched-filter score)", fontsize=9)
        else:
            ax.set_ylabel("χ²/N (σ fixed, high-pass)")
            ax.set_title(f"{sc.aid}: continuum template scan (lower is better)", fontsize=9)
        ax.set_xlabel("trial z")
    return fig


def fig_continuum_fit(res) -> Optional[Figure]:
    fit = res.artifacts.get("ppxf")
    if fit is None:
        return None
    pp, fs = fit["pp"], fit["fs"]
    lam = fs.lam
    gp = np.zeros(lam.size, bool)
    gp[fit["goodpixels"]] = True
    fig = Figure(figsize=(13, 5.5), constrained_layout=True)
    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1])
    ax = fig.add_subplot(gs[0]); ax2 = fig.add_subplot(gs[1], sharex=ax)
    ax.plot(lam, fit["flux"], color="k", lw=0.6, label="data")
    ax.plot(lam, pp.bestfit, color="crimson", lw=0.9, label="pPXF best fit")
    if not gp.all():
        ax.plot(lam[~gp], fit["flux"][~gp], ".", ms=1.5, color="0.7", label="masked")
    ax.legend(fontsize=7, ncol=3)
    resid = fit["flux"] - pp.bestfit
    ax2.plot(lam, np.where(gp, resid, np.nan), color="0.3", lw=0.5)
    ax2.axhline(0, color="k", lw=0.5)
    ax2.set_xlabel(f"observed wavelength ({fs.frame}) [Å]"); ax2.set_ylabel("resid")
    ttl = f"{res.id}: z={res.z:.5f}"
    if np.isfinite(res.sigma):
        ttl += f", σ={res.sigma:.0f}±{res.sigma_err:.0f} km/s"
    ttl += f", χ²/dof={fit['chi2']:.2f}, S/N/Å={res.snr_per_A:.1f}"
    ax.set_title(ttl, fontsize=10)
    return fig


def fig_lines(res) -> Optional[Figure]:
    lfit = res.artifacts.get("line_fit")
    ls = res.artifacts.get("ls")
    if lfit is None or ls is None:
        return None
    rows = [r for r in lfit["lines"]]
    rows.sort(key=lambda r: -r["snr"])
    rows = rows[:12]
    n = len(rows)
    ncol = min(4, max(n, 1))
    nrow = int(np.ceil(n / ncol))
    fig = Figure(figsize=(3.2 * ncol, 2.4 * nrow + 0.6), constrained_layout=True)
    axes = np.atleast_1d(fig.subplots(nrow, ncol, squeeze=False)).ravel()
    for ax in axes[n:]:
        ax.axis("off")
    for ax, r in zip(axes, rows):
        lam0 = r["lam_obs"]
        sel = np.abs(ls.wave - lam0) < 40
        ax.step(ls.wave[sel], (ls.flux - ls.cont)[sel], where="mid", color="k", lw=0.6)
        ax.plot(ls.wave[sel], lfit["model"][sel], color="crimson", lw=1)
        ax.axvline(lam0, color="C0", ls=":", lw=0.7)
        ax.set_title(f"{r['name']}  ({r['snr']:.1f}σ)", fontsize=8)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"{res.id}: line fits at z = {lfit['z']:.5f} ± {lfit['z_err']:.5f}", fontsize=10)
    return fig


def fig_bootstrap(res) -> Optional[Figure]:
    b = res.artifacts.get("bootstrap")
    if b is None or not b.get("ok"):
        return None
    fig = Figure(figsize=(9, 3), constrained_layout=True)
    ax = fig.subplots(1, 2)
    ax[0].hist(b["z_samples"], bins=20, color="steelblue"); ax[0].axvline(res.z, color="k", ls="--")
    ax[0].set_xlabel("z"); ax[0].set_title(f"σ_z = {b['z_err']:.2e}", fontsize=9)
    ax[1].hist(b["sigma_samples"], bins=20, color="darkorange"); ax[1].axvline(res.sigma, color="k", ls="--")
    ax[1].set_xlabel("σ [km/s]"); ax[1].set_title(f"σ_σ = {b['sigma_err']:.1f} km/s", fontsize=9)
    return fig


def fig_ladder(rows: List[dict]) -> Optional[Figure]:
    if not rows:
        return None
    fig = Figure(figsize=(8, 0.5 * len(rows) + 1.6), constrained_layout=True)
    ax = fig.subplots()
    cols = ROLE_COLOR
    for i, r in enumerate(rows):
        ax.errorbar(r["z"], i, xerr=max(r["z_err"], 1e-5), fmt="o", color=cols.get(r["role"], "k"),
                    mfc="white" if r["status"] != "secure" else None)
        ax.text(r["z"], i + 0.25, f" {r['id']}", fontsize=8)
    ax.set_yticks([]); ax.set_xlabel("redshift")
    ax.set_title("redshift ladder (open symbol = tentative)", fontsize=9)
    return fig
