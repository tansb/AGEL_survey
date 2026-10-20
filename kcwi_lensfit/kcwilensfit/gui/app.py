"""The graphical interface: one tab per pipeline stage.

Every input is pre-filled with its documented default and carries a visible description.  Heavy
work runs in a background thread so the window stays responsive; progress and messages stream
into the log pane.  Nothing here contains science logic -- it only drives
:class:`kcwilensfit.pipeline.Pipeline`.
"""
from __future__ import annotations

import json
import os
import queue
import sys
import threading
import traceback
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from typing import Callable, Dict, List, Optional

import numpy as np

from .. import __version__
from .. import plotting as P
from ..apertures import Aperture
from ..pipeline import Pipeline, PipelineError
from ..redshift import Cancelled
from ..schema import (APERTURE_COLUMNS, APERTURE_DEFAULTS, APERTURE_MODES, APERTURE_PROFILES, APERTURE_ROLES,
                      APERTURE_SHAPES, FIELDS, STAGES, defaults, load_config, normalise_aperture, save_config,
                      validate)
from .widgets import FigureFrame, Form, GREY, LogPane, remember_dir, start_dir

AP_TREE_COLS = ["id", "role", "group", "shape", "x", "y", "a", "b", "angle", "profile", "mode", "z_known"]


class App(tk.Tk):
    def __init__(self, config_path: Optional[str] = None):
        super().__init__()
        self._base_title = f"kcwilensfit {__version__} — KCWI/KCRM lens redshifts"
        self.title(self._base_title)
        self.geometry("1500x940")
        self.minsize(1100, 700)
        self.msgq: "queue.Queue" = queue.Queue()
        self.pipe = Pipeline(log=lambda m: self.msgq.put(("log", m)))
        self.cancel_flag = threading.Event()
        self.busy = False
        self.config_path: Optional[str] = None
        self.ap_rows: List[dict] = []
        self.scan_cache: Dict[str, list] = {}
        self._build()
        self.after(100, self._poll)
        if config_path:
            self.open_config(config_path)
        else:
            self._populate(defaults())

    # ------------------------------------------------------------------ layout
    def _build(self):
        style = ttk.Style(self)
        try:
            style.theme_use("clam" if sys.platform.startswith("linux") else style.theme_use())
        except tk.TclError:
            pass
        menu = tk.Menu(self)
        fm = tk.Menu(menu, tearoff=0)
        fm.add_command(label="New target (clear all)…", command=self.new_target)
        fm.add_command(label="Open config…", command=lambda: self.open_config(None))
        fm.add_command(label="Save config…", command=self.save_config_dialog)
        fm.add_separator()
        fm.add_command(label="Quit", command=self.destroy)
        menu.add_cascade(label="File", menu=fm)
        rm = tk.Menu(menu, tearoff=0)
        rm.add_command(label="Run ALL stages (non-interactive)", command=self.run_everything)
        rm.add_command(label="Cancel running task", command=self.cancel)
        menu.add_cascade(label="Run", menu=rm)
        hm = tk.Menu(menu, tearoff=0)
        hm.add_command(label="About", command=lambda: messagebox.showinfo(
            "About", f"kcwilensfit {__version__}\nRedshifts and kinematics of KCWI/KCRM strong-lens systems.\n"
                     "Every default is documented in the tab it appears in and in docs/science_notes.md."))
        menu.add_cascade(label="Help", menu=hm)
        self.config(menu=menu)

        paned = ttk.PanedWindow(self, orient="vertical")
        paned.pack(fill="both", expand=True)
        top = ttk.Frame(paned)
        paned.add(top, weight=5)
        self.nb = ttk.Notebook(top)
        self.nb.pack(fill="both", expand=True)
        bottom = ttk.Frame(paned)
        paned.add(bottom, weight=1)
        bar = ttk.Frame(bottom)
        bar.pack(fill="x")
        self.progress = ttk.Progressbar(bar, mode="determinate", maximum=1.0, length=300)
        self.progress.pack(side="left", padx=4, pady=2)
        self.status = ttk.Label(bar, text="Ready.", foreground=GREY)
        self.status.pack(side="left", padx=6)
        ttk.Button(bar, text="Cancel", command=self.cancel).pack(side="right", padx=4)
        ttk.Button(bar, text="Clear all (new target)", command=self.new_target).pack(side="right", padx=4)
        self.log = LogPane(bottom, height=8)
        self.log.pack(fill="both", expand=True)

        self.forms: Dict[str, Form] = {}
        self.figs: Dict[str, FigureFrame] = {}
        stage_fields = {s: [f for f in FIELDS if f.stage == s] for s, _ in STAGES}
        self.tabs = {}
        for key, title in STAGES:
            tab = ttk.Frame(self.nb)
            self.tabs[key] = tab
            self.nb.add(tab, text=title)
        self._tab_data(stage_fields["inputs"])
        self._tab_noise(stage_fields["noise"])
        self._tab_apertures(stage_fields["apertures"])
        self._tab_extraction(stage_fields["extraction"])
        self._tab_fitting(stage_fields["fitting"])
        self._tab_output(stage_fields["output"])
        tn = self.forms["inputs"].rows["target_name"].var
        tn.trace_add("write", lambda *_: self._update_title())
        self._update_title()
        self.nb.bind("<<NotebookTabChanged>>", self._on_tab_changed)

    def _update_title(self):
        name = self.forms["inputs"].rows["target_name"].var.get().strip()
        self.title(f"{self._base_title} — {name}" if name else self._base_title)

    def _on_tab_changed(self, _=None):
        if self.nb.select() == str(self.tabs["output"]):
            self.refresh_output_preview()

    def _split(self, tab, specs, left_w=560):
        left = ttk.Frame(tab, width=left_w)
        left.pack(side="left", fill="y")
        left.pack_propagate(False)
        right = ttk.Frame(tab)
        right.pack(side="left", fill="both", expand=True)
        return left, right

    @staticmethod
    def _btn(parent, text, cmd, **kw):
        b = ttk.Button(parent, text=text, command=cmd, **kw)
        b.pack(side="left", padx=3, pady=3)
        return b

    # ------------------------------------------------------------------ tabs
    def _tab_data(self, specs):
        left, right = self._split(self.tabs["inputs"], specs)
        intro = ttk.Label(left, foreground=GREY, wraplength=520, justify="left",
                          text="Welcome. Choose your reduced KCWI (blue) and/or KCRM (red) cubes, check the instrument "
                               "settings below, then press “Load cubes”. Every field is pre-filled with a documented "
                               "default and described underneath.")
        intro.pack(anchor="w", padx=8, pady=6)
        f = Form(left, specs)
        f.pack(fill="both", expand=True)
        self.forms["inputs"] = f
        bb = ttk.Frame(left); bb.pack(fill="x")
        self._btn(bb, "Load cubes (applies these settings)", self.do_load)
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)
        for key, title in (("data", "White light"), ("instrument", "Registration check")):
            fr = FigureFrame(nb); self.figs[key] = fr; nb.add(fr, text=title)

    def _tab_noise(self, specs):
        left, right = self._split(self.tabs["noise"], specs)
        f = Form(left, specs); f.pack(fill="both", expand=True)
        self.forms["noise"] = f
        bb = ttk.Frame(left); bb.pack(fill="x")
        self._btn(bb, "Preview sky box on FoV", self.preview_sky_box)
        self._btn(bb, "Suggest box → manual fields", self.do_suggest_sky)
        self._btn(bb, "Estimate sky & noise", self.do_sky)
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)
        for key, title in (("noise_box", "Sky box on FoV"), ("noise", "Background & noise")):
            fr = FigureFrame(nb); self.figs[key] = fr; nb.add(fr, text=title)
        self.noise_nb = nb

    def _current_sky_boxes(self):
        """Boxes implied by the form: 'auto' -> the suggested box, 'manual' -> the typed limits."""
        from ..noise import parse_box
        if not self.pipe.arms:
            return {}
        cfg = {k: row.get() for k, row in self.forms["noise"].rows.items()}
        if cfg["sky_box_mode"] == "auto":
            self.pipe.update_config(dict(sky_box_size=cfg["sky_box_size"]))
            boxes = self.pipe.suggest_sky_box()
            if len(self.pipe.arms) == 2 and self.pipe.registration.mode == "identity":
                boxes["red"] = boxes[self.pipe.ref_arm]
            return boxes
        out = {}
        for name in self.pipe.arms:
            txt = cfg["sky_box_blue"] if name == "blue" else (cfg["sky_box_red"] or cfg["sky_box_blue"])
            if txt:
                out[name] = parse_box(txt)
        return out

    def preview_sky_box(self, quiet=False):
        if not self.pipe.arms:
            if not quiet:
                messagebox.showinfo("Load first", "Load the cubes (Data tab) first.")
            return
        try:
            boxes = self._current_sky_boxes()
            self.pipe.apertures = [Aperture.from_dict(a) for a in self.ap_rows]
            fig = P.fig_sky_box(self.pipe, boxes)
        except Exception as e:
            if not quiet:
                messagebox.showerror("Cannot preview the sky box", str(e))
            return
        self.figs["noise_box"].show(fig)
        self.noise_nb.select(self.figs["noise_box"])

    def _tab_apertures(self, specs):
        tab = self.tabs["apertures"]
        left = ttk.Frame(tab, width=640); left.pack(side="left", fill="y"); left.pack_propagate(False)
        right = ttk.Frame(tab); right.pack(side="left", fill="both", expand=True)
        f = Form(left, specs)
        f.pack(fill="x")
        self.forms["apertures"] = f
        ttk.Label(left, foreground=GREY, wraplength=610, justify="left", text=
                  "Define where to extract spectra. 'deflector' = foreground lens galaxy; 'source' = lensed background image "
                  "(give all images of one source the same group); 'other' = companion/interloper. Region files in WCS "
                  "(fk5/icrs) are converted to pixels separately for each cube with that cube's own WCS, and the shapes "
                  "actually used are drawn on the right.").pack(anchor="w", padx=8, pady=4)
        tv = ttk.Treeview(left, columns=AP_TREE_COLS, show="headings", height=8, selectmode="browse")
        for c in AP_TREE_COLS:
            tv.heading(c, text=c)
            tv.column(c, width=48 if c not in ("id", "group", "role", "shape") else 70, anchor="w")
        tv.pack(fill="x", padx=6)
        tv.bind("<<TreeviewSelect>>", self._ap_selected)
        self.ap_tree = tv
        bb = ttk.Frame(left); bb.pack(fill="x")
        self._btn(bb, "Add", self.ap_add)
        self._btn(bb, "Duplicate", self.ap_dup)
        self._btn(bb, "Delete", self.ap_del)
        self._btn(bb, "Load DS9 region file…", self.ap_load_regions)
        bb2 = ttk.Frame(left); bb2.pack(fill="x")
        self._btn(bb2, "Auto-detect sources", self.ap_autodetect)
        self._btn(bb2, "Refresh preview", self.refresh_ap_preview)
        self.click_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(bb2, text="click on the reference (left) image to move the selected aperture",
                        variable=self.click_var).pack(side="left", padx=6)
        # edit form
        ed = ttk.LabelFrame(left, text="Selected aperture (edit, then ‘Apply’)")
        ed.pack(fill="x", padx=6, pady=4)
        self.ap_vars: Dict[str, tk.StringVar] = {}
        self.ap_help = ttk.Label(ed, foreground=GREY, wraplength=600, justify="left", text="(hover/focus a field for its description)")
        r = 0
        choices = {"role": APERTURE_ROLES, "shape": APERTURE_SHAPES, "profile": APERTURE_PROFILES, "mode": APERTURE_MODES}
        cols = [c for c in APERTURE_COLUMNS if c[0] in AP_TREE_COLS]
        for i, (key, label, dflt, helptxt) in enumerate(cols):
            rr, cc = divmod(i, 4)
            fr = ttk.Frame(ed); fr.grid(row=rr, column=cc, padx=4, pady=2, sticky="w")
            ttk.Label(fr, text=label, font=("TkDefaultFont", 9, "bold")).pack(anchor="w")
            v = tk.StringVar(value=str(dflt)); self.ap_vars[key] = v
            if key in choices:
                w = ttk.Combobox(fr, textvariable=v, values=choices[key], state="readonly", width=11)
            else:
                w = ttk.Entry(fr, textvariable=v, width=11)
            w.pack(anchor="w")
            w.bind("<FocusIn>", lambda e, t=helptxt: self.ap_help.configure(text=t))
            w.bind("<Enter>", lambda e, t=helptxt: self.ap_help.configure(text=t))
        self.ap_help.grid(row=10, column=0, columnspan=4, sticky="w", padx=4, pady=4)
        ttk.Button(ed, text="Apply", command=self.ap_apply).grid(row=11, column=0, pady=4, sticky="w", padx=4)
        self.figs["apertures"] = FigureFrame(right)
        self.figs["apertures"].pack(fill="both", expand=True)
        self.figs["apertures"].on_click(self._ap_click)

    def _tab_extraction(self, specs):
        left, right = self._split(self.tabs["extraction"], specs)
        f = Form(left, specs); f.pack(fill="both", expand=True)
        self.forms["extraction"] = f
        bb = ttk.Frame(left); bb.pack(fill="x")
        self._btn(bb, "Run extraction", self.do_extract)
        box = ttk.LabelFrame(right, text="Extracted spectra: how many spaxels, what S/N")
        box.pack(fill="x", padx=4, pady=2)
        cols = ("spectrum", "role", "spaxels_blue", "spaxels_red", "snr_blue", "snr_red", "subtracted")
        heads = ("spectrum", "role", "spaxels (blue)", "spaxels (red)", "S/N/pix (blue)", "S/N/pix (red)", "manual subtraction")
        self.ex_tree = ttk.Treeview(box, columns=cols, show="headings", height=6, selectmode="none")
        for c, h, w in zip(cols, heads, (110, 70, 130, 130, 100, 100, 320)):
            self.ex_tree.heading(c, text=h); self.ex_tree.column(c, width=w, anchor="w")
        self.ex_tree.pack(fill="x", padx=4, pady=3)
        ttk.Label(box, foreground=GREY, wraplength=900, justify="left", text=
                  "spaxels = aperture area as the sum of fractional pixel coverage (the number of spaxels' worth of light summed); "
                  "the count of spaxels touched at all is in brackets. Combined-image spectra ('group:…') add their members.").pack(
            anchor="w", padx=6, pady=(0, 3))
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)
        for key, title in (("ex_spec", "Spectra"), ("ex_deb", "Deblend model"), ("ex_ap", "Apertures")):
            fr = FigureFrame(nb); self.figs[key] = fr; nb.add(fr, text=title)

    def _refresh_ex_table(self):
        self.ex_tree.delete(*self.ex_tree.get_children())
        for sid, sp in self.pipe.spectra.items():
            cnt = sp.meta.get("counts") or {}
            snr = sp.meta.get("snr", {})
            cell = lambda a: (f"{cnt[a]['area']:.1f}  ({cnt[a]['touched']})" if a in cnt else "—")
            sn = lambda a: (f"{snr[a]:.2f}" if a in snr else "—")
            sub = "; ".join(f"− {r['scale']:.3f} × {r['contaminator']}" + (" (auto)" if r["auto"] else "")
                            for r in sp.meta.get("subtracted", [])) or ""
            self.ex_tree.insert("", "end", values=(sid, sp.role, cell("blue"), cell("red"), sn("blue"), sn("red"), sub))

    def _tab_fitting(self, specs):
        tab = self.tabs["fitting"]
        left = ttk.Frame(tab, width=520); left.pack(side="left", fill="y"); left.pack_propagate(False)
        right = ttk.Frame(tab); right.pack(side="left", fill="both", expand=True)
        f = Form(left, specs); f.pack(fill="both", expand=True)
        self.forms["fitting"] = f
        # interactive control panel
        ctl = ttk.LabelFrame(right, text="Fit aperture by aperture (or all at once)")
        ctl.pack(fill="x", padx=4, pady=2)
        row1 = ttk.Frame(ctl); row1.pack(fill="x")
        ttk.Label(row1, text="Spectrum:").pack(side="left", padx=4)
        self.fit_target = tk.StringVar()
        self.fit_combo = ttk.Combobox(row1, textvariable=self.fit_target, state="readonly", width=16)
        self.fit_combo.pack(side="left")
        self.fit_combo.bind("<<ComboboxSelected>>", lambda e: self._refresh_fit_view())
        ttk.Label(row1, text="mode:").pack(side="left", padx=(10, 2))
        self.fit_mode = tk.StringVar(value="auto")
        ttk.Combobox(row1, textvariable=self.fit_mode, values=APERTURE_MODES[:-1], state="readonly", width=10).pack(side="left")
        self._btn(row1, "1. Scan redshift", self.do_scan)
        self._btn(row1, "2. Fit selected candidate / typed z", self.do_fit_selected)
        self._btn(row1, "Auto-fit this", self.do_fit_auto)
        self._btn(row1, "Auto-fit ALL", self.do_fit_all)
        self._btn(row1, "Clear fit (this region)", self.do_clear_fit)
        row2 = ttk.Frame(ctl); row2.pack(fill="x")
        ttk.Label(row2, text="…or type a redshift to refine around:").pack(side="left", padx=4)
        self.z_manual = tk.StringVar()
        ttk.Entry(row2, textvariable=self.z_manual, width=10).pack(side="left")
        ttk.Label(row2, text="Candidates (double-click to choose; red = below deflector, impossible for a source):",
                  foreground=GREY).pack(side="left", padx=10)
        cols = ("n", "method", "z", "score", "info")
        self.cand_tree = ttk.Treeview(ctl, columns=cols, show="headings", height=5, selectmode="browse")
        for c, w in zip(cols, (30, 80, 80, 70, 700)):
            self.cand_tree.heading(c, text=c); self.cand_tree.column(c, width=w, anchor="w")
        self.cand_tree.pack(fill="x", padx=4)
        self.cand_tree.tag_configure("below", foreground="#b00020")
        self.cand_tree.bind("<Double-1>", lambda e: self.do_fit_selected())
        rcols = ("id", "role", "status", "method", "z", "z_err", "sigma", "flags")
        self.res_tree = ttk.Treeview(ctl, columns=rcols, show="headings", height=5)
        for c, w in zip(rcols, (90, 70, 80, 80, 80, 70, 100, 600)):
            self.res_tree.heading(c, text=c); self.res_tree.column(c, width=w, anchor="w")
        self.res_tree.pack(fill="x", padx=4, pady=3)
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)
        sumfr = ttk.Frame(nb); nb.add(sumfr, text="Fit summary")
        self.fit_summary = tk.Text(sumfr, wrap="word", font=("Menlo", 10), height=12)
        self.fit_summary.pack(fill="both", expand=True)
        self.fit_summary.configure(state="disabled")
        for key, title in (("fit_scan", "Scan"), ("fit_cont", "Continuum fit"), ("fit_lines", "Line fits"),
                           ("fit_spec", "Spectra + z"), ("fit_boot", "Bootstrap"), ("fit_ladder", "Redshift ladder")):
            fr = FigureFrame(nb); self.figs[key] = fr; nb.add(fr, text=title)

    def _tab_output(self, specs):
        left, right = self._split(self.tabs["output"], specs)
        dirbox = ttk.LabelFrame(left, text="Output directory")
        dirbox.pack(fill="x", padx=6, pady=4)
        self.out_dir_label = ttk.Label(dirbox, text="", wraplength=520, justify="left")
        self.out_dir_label.pack(anchor="w", padx=6, pady=2)
        ttk.Button(dirbox, text="Change folder…", command=self.change_out_dir).pack(anchor="w", padx=6, pady=(0, 4))
        f = Form(left, specs); f.pack(fill="both", expand=True)
        self.forms["output"] = f
        bb = ttk.Frame(left); bb.pack(fill="x")
        self._btn(bb, "Refresh preview", self.refresh_output_preview)
        self._btn(bb, "Save all outputs", self.do_save)
        self._btn(bb, "Open output folder", self.open_out)
        nb = ttk.Notebook(right); nb.pack(fill="both", expand=True)
        sumfr = ttk.Frame(nb); nb.add(sumfr, text="Summary")
        self.summary_text = tk.Text(sumfr, wrap="none", font=("Menlo", 10), height=12)
        self.summary_text.pack(fill="both", expand=True)
        for key, title in (("out_ap", "Apertures"), ("out_spec", "Spectra + z (fitted regions)")):
            fr = FigureFrame(nb); self.figs[key] = fr; nb.add(fr, text=title)
        repfr = ttk.Frame(nb); nb.add(repfr, text="Full report")
        self.report_text = tk.Text(repfr, wrap="word", font=("Menlo", 10))
        self.report_text.pack(fill="both", expand=True)

    def _planned_out_dir(self) -> str:
        try:
            d = self.forms["inputs"].rows["output_dir"].get()
        except Exception:
            d = self.pipe.cfg.get("output_dir", "")
        return os.path.abspath(os.path.expanduser(d or "kcwilensfit_output"))

    def change_out_dir(self):
        p = filedialog.askdirectory(title="Output directory", initialdir=start_dir(self._planned_out_dir()), mustexist=False)
        if p:
            remember_dir(p)
            self.forms["inputs"].rows["output_dir"].set(p)
            self.refresh_output_preview()

    def refresh_output_preview(self, saved: bool = False):
        """Summary text + the aperture and fitted-spectra plots, shown BEFORE anything is written."""
        from ..report import markdown_report, quick_summary
        out = self._planned_out_dir()
        self.out_dir_label.configure(text=out)
        txt = quick_summary(self.pipe, out, saved=saved)
        self.summary_text.configure(state="normal")
        self.summary_text.delete("1.0", "end"); self.summary_text.insert("1.0", txt)
        self.summary_text.configure(state="disabled")
        if self.pipe.apertures and self.pipe.arms:
            try:
                self.figs["out_ap"].show(P.fig_apertures(self.pipe))
            except Exception as e:
                self.figs["out_ap"].clear(f"Cannot draw apertures: {e}")
        else:
            self.figs["out_ap"].clear("Load cubes and define apertures first.")
        hide = P.combined_member_ids(self.pipe)
        fitted = [a for a in self.pipe.spectra if a in self.pipe.results and a not in hide]
        if fitted:
            try:
                self.figs["out_spec"].show(P.fig_spectra(
                    self.pipe, fitted, float(self.pipe.cfg["display_smooth_A"]),
                    windows=self.pipe.engine.static_windows() if self.pipe.engine else None))
            except Exception as e:
                self.figs["out_spec"].clear(f"Cannot draw spectra: {e}")
        else:
            self.figs["out_spec"].clear("No fitted regions yet.")
        self.report_text.delete("1.0", "end")
        if self.pipe.spectra:
            try:
                self.report_text.insert("1.0", markdown_report(self.pipe))
            except Exception as e:
                self.report_text.insert("1.0", f"(report preview unavailable: {e})")

    # ------------------------------------------------------------------ config <-> widgets
    def _populate(self, cfg: dict):
        for f in self.forms.values():
            f.populate(cfg)
        self.ap_rows = [normalise_aperture(a) for a in cfg.get("apertures", [])]
        self._refresh_ap_tree()

    def _collect(self) -> dict:
        cfg = {}
        for f in self.forms.values():
            cfg.update(f.collect())
        cfg["apertures"] = [dict(a) for a in self.ap_rows]
        cfg["region_file"] = self.pipe.cfg.get("region_file", "")
        return cfg

    def _apply_cfg(self) -> bool:
        try:
            cfg = self._collect()
        except Exception as e:
            messagebox.showerror("Invalid input", str(e))
            return False
        self.pipe.update_config(cfg)
        return True

    def open_config(self, path: Optional[str]):
        path = path or filedialog.askopenfilename(title="Open config", initialdir=start_dir())
        remember_dir(path)
        if not path:
            return
        try:
            cfg = load_config(path)
        except Exception as e:
            messagebox.showerror("Cannot read config", str(e))
            return
        self.config_path = path
        self._populate(cfg)
        self.pipe.update_config(cfg)
        self.write(f"Loaded settings from {path}")

    def new_target(self):
        """Discard all loaded data, results, plots and settings so a new target can be set up."""
        if self.busy:
            messagebox.showinfo("Busy", "Wait for (or cancel) the running task first.")
            return
        if not messagebox.askyesno("Clear all", "Clear all data, results, plots and settings (reset to defaults)?\n"
                                                "Unsaved results will be lost."):
            return
        self.pipe.reset()
        self.pipe.cfg = defaults()
        self.scan_cache.clear()
        self.config_path = None
        self._populate(defaults())
        self.pipe.update_config(defaults())
        for fr in self.figs.values():
            fr.clear("(no figure yet)")
        for tv in (self.ex_tree, self.cand_tree, self.res_tree):
            tv.delete(*tv.get_children())
        self.fit_target.set("")
        self.fit_combo["values"] = []
        for t in (self.fit_summary, self.report_text, self.summary_text):
            was = str(t.cget("state"))
            t.configure(state="normal"); t.delete("1.0", "end"); t.configure(state=was)
        self.progress["value"] = 0
        self.status.configure(text="Ready.")
        self.refresh_output_preview()
        self.write("Cleared all data and settings -- ready for a new target.")

    def save_config_dialog(self):
        if not self._apply_cfg():
            return
        path = filedialog.asksaveasfilename(title="Save config", defaultextension=".json",
                                            initialfile="kcwilensfit_config.json", initialdir=start_dir())
        if path:
            save_config(self.pipe.cfg, path)
            self.write(f"Saved settings to {path}")

    # ------------------------------------------------------------------ background tasks
    def write(self, msg: str):
        self.log.write(msg)

    def _poll(self):
        try:
            while True:
                kind, payload = self.msgq.get_nowait()
                if kind == "log":
                    self.write(payload)
                elif kind == "progress":
                    f, m = payload
                    self.progress["value"] = max(0.0, min(1.0, f))
                    self.status.configure(text=m)
                elif kind == "done":
                    fn, ok, res = payload
                    self.busy = False
                    self.progress["value"] = 0
                    self.status.configure(text="Ready." if ok else "Failed — see log.")
                    if fn:
                        try:
                            fn(ok, res)
                        except Exception:
                            self.write("GUI error:\n" + traceback.format_exc())
        except queue.Empty:
            pass
        self.after(100, self._poll)

    def cancel(self):
        if self.busy:
            self.cancel_flag.set()
            self.write("Cancelling … (finishes the current step)")

    def run_task(self, label: str, work: Callable, done: Optional[Callable] = None):
        if self.busy:
            messagebox.showinfo("Busy", "A task is already running.")
            return
        if not self._apply_cfg():
            return
        self.busy = True
        self.cancel_flag.clear()
        self.status.configure(text=label + " …")

        def progress(f, m):
            if self.cancel_flag.is_set():
                raise Cancelled()
            self.msgq.put(("progress", (f, f"{label}: {m}")))

        def target():
            try:
                res = work(progress)
                self.msgq.put(("done", (done, True, res)))
            except Cancelled:
                self.msgq.put(("log", "Cancelled."))
                self.msgq.put(("done", (None, False, None)))
            except Exception as e:
                tb = traceback.format_exc()
                self.msgq.put(("log", f"FAILED: {type(e).__name__}: {e}"))
                if os.environ.get("KCWILENSFIT_DEBUG"):
                    self.msgq.put(("log", tb))
                self.msgq.put(("done", (None, False, None)))
        threading.Thread(target=target, daemon=True).start()

    # ------------------------------------------------------------------ stage actions
    def do_load(self):
        def done(ok, _):
            self.refresh_ap_preview(quiet=True)
            self._show_load_figs()
            self.preview_sky_box(quiet=True)
            self.write("Cubes loaded. Next: the Noise & sky and Apertures tabs.")
        self.run_task("Loading cubes", lambda pr: self.pipe.load(pr), done)

    def _show_load_figs(self):
        if not self.pipe.arms:
            return
        fig = P.fig_apertures(self.pipe, "White light, both arms (apertures appear once defined)")
        self.figs["data"].show(fig)
        self.figs["instrument"].show(P.fig_apertures(self.pipe, "Registration check"))

    def do_suggest_sky(self):
        if not self.pipe.arms:
            messagebox.showinfo("Load first", "Load the cubes (Data tab) first.")
            return
        if not self._apply_cfg():
            return
        try:
            boxes = self.pipe.suggest_sky_box()
        except Exception as e:
            messagebox.showerror("No sky box", str(e)); return
        f = self.forms["noise"]
        b = boxes.get("blue") or boxes.get("red")
        r = boxes.get("red") or b
        if len(boxes) == 2 and self.pipe.registration.mode == "identity":
            r = b          # same spaxels in both arms: use ONE source-free patch for both
        f.rows["sky_box_blue"].set(",".join(map(str, b)))
        f.rows["sky_box_red"].set(",".join(map(str, r)))
        f.rows["sky_box_mode"].set("manual")
        self.write(f"Suggested sky box (y0,y1,x0,x1): blue {b}, red {r}. Mode switched to 'manual' so you can adjust it.")
        self.preview_sky_box(quiet=True)

    def do_sky(self):
        if not self.pipe.arms:
            messagebox.showinfo("Load first", "Load the cubes (Data tab) first.")
            return
        def work(pr):
            self.pipe.estimate_sky()
        def done(ok, _):
            self.figs["noise"].show(P.fig_sky(self.pipe))
            self.figs["noise_box"].show(P.fig_sky_box(self.pipe, {k: v.box for k, v in self.pipe.sky.items()},
                                                      "Sky box used for the background/noise estimate"))
            self.refresh_ap_preview(quiet=True)
        self.run_task("Estimating sky & noise", work, done)

    # ---- apertures
    def _refresh_ap_tree(self):
        self.ap_tree.delete(*self.ap_tree.get_children())
        for i, a in enumerate(self.ap_rows):
            vals = []
            for c in AP_TREE_COLS:
                v = a.get(c, "")
                vals.append(f"{v:.2f}" if isinstance(v, float) and c in ("x", "y", "a", "b", "angle") else v)
            self.ap_tree.insert("", "end", iid=str(i), values=vals)

    def _sel_index(self) -> Optional[int]:
        s = self.ap_tree.selection()
        return int(s[0]) if s else None

    def _ap_selected(self, _=None):
        i = self._sel_index()
        if i is None:
            return
        a = self.ap_rows[i]
        for k, v in self.ap_vars.items():
            x = a.get(k, "")
            v.set(f"{x:g}" if isinstance(x, float) else str(x))

    def ap_apply(self):
        i = self._sel_index()
        if i is None:
            messagebox.showinfo("Select a row", "Select an aperture in the table first (or use Add).")
            return
        a = dict(self.ap_rows[i])
        for k, v in self.ap_vars.items():
            a[k] = v.get()
        try:
            a = normalise_aperture(a)
        except ValueError as e:
            messagebox.showerror("Invalid value", str(e)); return
        if any(abs(a[k] - self.ap_rows[i][k]) > 1e-9 for k in ("x", "y", "a", "b", "angle")):
            for k in ("ra", "dec", "a_as", "b_as"):
                a[k] = ""             # edited in pixels: the sky definition no longer applies
        if any(j != i and r["id"] == a["id"] for j, r in enumerate(self.ap_rows)):
            messagebox.showerror("Duplicate id", "Aperture ids must be unique."); return
        self.ap_rows[i] = a
        self._refresh_ap_tree()
        self.ap_tree.selection_set(str(i))
        self.refresh_ap_preview(quiet=True)

    def ap_add(self):
        n = len(self.ap_rows) + 1
        a = normalise_aperture(dict(APERTURE_DEFAULTS, id=f"ap{n}", role="deflector" if n == 1 else "source"))
        if self.pipe.arms:
            ref = self.pipe.arms[self.pipe.ref_arm]
            ys, xs = np.nonzero(ref.valid)
            a["x"], a["y"] = float(np.mean(xs)), float(np.mean(ys))
        self.ap_rows.append(a)
        self._refresh_ap_tree()
        self.ap_tree.selection_set(str(len(self.ap_rows) - 1))
        self.refresh_ap_preview(quiet=True)

    def ap_dup(self):
        i = self._sel_index()
        if i is None:
            return
        a = dict(self.ap_rows[i]); a["id"] += "_copy"; a["x"] += 2
        self.ap_rows.append(a); self._refresh_ap_tree(); self.refresh_ap_preview(quiet=True)

    def ap_del(self):
        i = self._sel_index()
        if i is None:
            return
        del self.ap_rows[i]; self._refresh_ap_tree(); self.refresh_ap_preview(quiet=True)

    def ap_load_regions(self):
        if not self.pipe.arms:
            messagebox.showinfo("Load first", "Load the cubes first: WCS regions need the cube WCS.")
            return
        path = filedialog.askopenfilename(title="DS9 region file", initialdir=start_dir())
        if not path:
            return
        remember_dir(path)
        try:
            aps, warns = self.pipe.apertures_from_regions(path)
        except Exception as e:
            messagebox.showerror("Cannot read regions", str(e)); return
        for w in warns:
            self.write("[region file] " + w)
        if not aps:
            messagebox.showinfo("No regions", "No circle/ellipse/box regions were found."); return
        rows = []
        for a in aps:
            d = a.to_dict()
            rows.append(normalise_aperture(d))
        if self.ap_rows and messagebox.askyesno("Replace?", f"Replace the {len(self.ap_rows)} existing apertures with the "
                                                           f"{len(rows)} regions from the file? (No = append)"):
            self.ap_rows = rows
        else:
            used = {a['id'] for a in self.ap_rows}
            for r in rows:
                while r["id"] in used:
                    r["id"] += "_"
                used.add(r["id"])
            self.ap_rows += rows
        sky = sum(1 for r in rows if r.get("ra") not in ("", None))
        self.write(f"Loaded {len(rows)} regions from {os.path.basename(path)} ({sky} defined in sky coordinates -> mapped to each "
                   "cube with its own WCS). Set each aperture's ROLE and GROUP, then check the preview.")
        self.pipe.cfg["region_file"] = path
        self._refresh_ap_tree()
        self.refresh_ap_preview(quiet=True)

    def ap_autodetect(self):
        if not self.pipe.arms:
            messagebox.showinfo("Load first", "Load the cubes first.")
            return
        aps = self.pipe.suggest_apertures()
        if not aps:
            messagebox.showinfo("Nothing found", "No significant peaks above the threshold.")
            return
        self.ap_rows = [normalise_aperture(a.to_dict()) for a in aps]
        self.write(f"Auto-detected {len(aps)} peaks; the brightest is labelled 'deflector' — CHECK roles and radii yourself.")
        self._refresh_ap_tree(); self.refresh_ap_preview(quiet=True)

    def refresh_ap_preview(self, quiet=False):
        if not self.pipe.arms:
            if not quiet:
                messagebox.showinfo("Load first", "Load the cubes first.")
            return
        try:
            self.pipe.apertures = [Aperture.from_dict(a) for a in self.ap_rows]
            self.pipe.update_config(dict(sky_regions_use_wcs=self.forms["apertures"].rows["sky_regions_use_wcs"].get()))
            fig = P.fig_apertures(self.pipe)
        except Exception as e:
            self.write(f"[warning] preview failed: {e}")
            return
        self.figs["apertures"].show(fig)
        self.figs["ex_ap"].show(P.fig_apertures(self.pipe))

    def _ap_click(self, ev):
        if not self.click_var.get() or ev.inaxes is None or ev.xdata is None:
            return
        i = self._sel_index()
        if i is None:
            return
        fig = self.figs["apertures"].canvas.figure
        if ev.inaxes is not fig.axes[0]:
            self.write("Click on the LEFT image (reference arm) to move an aperture.")
            return
        a = dict(self.ap_rows[i])
        a["x"], a["y"] = float(ev.xdata), float(ev.ydata)
        for k in ("ra", "dec", "a_as", "b_as"):
            a[k] = ""
        self.ap_rows[i] = a
        self._refresh_ap_tree(); self.ap_tree.selection_set(str(i)); self._ap_selected()
        self.refresh_ap_preview(quiet=True)

    # ---- extraction
    def do_extract(self):
        if not self.pipe.sky:
            messagebox.showinfo("Noise first", "Run Estimate sky & noise (Noise & sky tab) first."); return
        if not self.ap_rows:
            messagebox.showinfo("Apertures", "Define at least one aperture (Apertures tab) first."); return
        def work(pr):
            self.pipe.set_apertures(self.ap_rows)
            self.pipe.extract(pr)
            self.pipe.results = {}            # earlier fits belong to the old extraction
            self.pipe.system_rows, self.pipe.system_messages = [], []
            self.scan_cache.clear()
        def done(ok, _):
            self.figs["ex_spec"].show(P.fig_spectra(self.pipe, smooth_A=float(self.pipe.cfg["display_smooth_A"])))
            fd = P.fig_deblend(self.pipe)
            if fd is not None:
                self.figs["ex_deb"].show(fd)
            self.figs["ex_ap"].show(P.fig_apertures(self.pipe))
            self._refresh_fit_targets()
            self._refresh_fit_view()
            self._refresh_ex_table()
            self.write("Extraction finished. Next: the Redshift & kinematics tab.")
        self.run_task("Extracting", work, done)

    # ---- fitting
    def _refresh_fit_targets(self):
        ids = list(self.pipe.spectra)
        self.fit_combo["values"] = ids
        if ids and self.fit_target.get() not in ids:
            self.fit_target.set(ids[0])

    def _need_spectra(self) -> bool:
        if not self.pipe.spectra:
            messagebox.showinfo("Extract first", "Run the extraction (Extraction tab) first.")
            return False
        return True

    def do_scan(self):
        if not self._need_spectra():
            return
        aid = self.fit_target.get()
        mode = self.fit_mode.get()          # Tk variables must be read on the main thread
        def work(pr):
            eng = self.pipe.make_engine() if self.pipe.engine is None else self.pipe.engine
            eng.cfg = self.pipe.cfg
            eng.results = self.pipe.results
            if mode in ("auto", "both"):
                scans = eng.scan(aid, progress=pr)
            else:
                scans = eng.scan(aid, kind=("lines" if mode == "lines" else "continuum"), progress=pr)
            return scans
        def done(ok, scans):
            self.scan_cache[aid] = scans
            if self.fit_target.get() == aid:
                self._refresh_fit_view()
            self.write(f"{aid}: scan complete. Pick a candidate (or type a z) and press ‘Fit selected’.")
        self.run_task(f"Scanning {aid}", work, done)

    def do_fit_selected(self):
        if not self._need_spectra():
            return
        aid = self.fit_target.get()
        z = None
        txt = self.z_manual.get().strip()
        if txt:
            try:
                z = float(txt)
            except ValueError:
                messagebox.showerror("Bad redshift", f"'{txt}' is not a number."); return
        else:
            sel = self.cand_tree.selection()
            if sel:
                z = float(sel[0].split(":")[1])
        if z is None:
            messagebox.showinfo("Choose a redshift", "Run ‘Scan redshift’ and select a candidate, or type a redshift."); return
        self._fit(aid, z)

    def do_fit_auto(self):
        if self._need_spectra():
            self._fit(self.fit_target.get(), None)

    def _fit(self, aid, z):
        mode = self.fit_mode.get()          # read Tk variables on the main thread
        def work(pr):
            eng = self.pipe.engine or self.pipe.make_engine()
            eng.cfg = self.pipe.cfg
            return self.pipe.fit_aperture(aid, z, None if mode == "auto" else mode, self.scan_cache.get(aid), pr)
        def done(ok, r):
            self._show_results(aid)
        self.run_task(f"Fitting {aid}", work, done)

    def do_clear_fit(self):
        aid = self.fit_target.get()
        if self.busy:
            messagebox.showinfo("Busy", "Wait for (or cancel) the running task first.")
            return
        if not aid or aid not in self.pipe.results:
            messagebox.showinfo("No fit", f"'{aid}' has no recorded fit to clear.")
            return
        if not messagebox.askyesno("Clear fit", f"Discard the recorded fit for '{aid}'?\n"
                                                "It will be left out of the saved outputs and the redshift ladder."):
            return
        self.pipe.clear_fit(aid)
        self._show_results(aid)
        self.write(f"Cleared the fit for '{aid}'.")

    def do_fit_all(self):
        if not self._need_spectra():
            return
        def work(pr):
            self.pipe.engine = None
            self.pipe.make_engine()
            return self.pipe.fit_all(pr)
        def done(ok, _):
            self._show_results(self.fit_target.get())
            for m in self.pipe.system_messages:
                self.write("  * " + m)
        self.run_task("Fitting all apertures", work, done)

    def _refresh_fit_view(self, _=None):
        """Show ONLY what belongs to the spectrum selected in the 'Spectrum' box: its scan candidates,
        its result row and its figures. Nothing from a previously selected spectrum stays visible."""
        aid = self.fit_target.get()
        self.z_manual.set("")
        # --- candidates + scan figure
        self.cand_tree.delete(*self.cand_tree.get_children())
        scans = self.scan_cache.get(aid)
        r = self.pipe.results.get(aid)
        if not scans and r is not None:
            scans = r.artifacts.get("scans")
        if scans:
            n = 0
            for sc in scans:
                for c in sc.candidates[:8]:
                    n += 1
                    self.cand_tree.insert("", "end", iid=f"{sc.kind}:{c.z:.6f}",
                                          values=(n, sc.kind, f"{c.z:.5f}", f"{c.score:.1f}",
                                                  c.info + ("   [BELOW DEFLECTOR z: impossible for a lensed source]" if c.below_deflector else "")),
                                          tags=("below",) if c.below_deflector else ())
            self.figs["fit_scan"].show(P.fig_scan(scans))
        else:
            self.figs["fit_scan"].clear(f"No scan run yet for '{aid}'.")
        # --- result row (this spectrum only)
        self.res_tree.delete(*self.res_tree.get_children())
        if r is not None:
            sig = (f"{r.sigma:.0f}±{r.sigma_err:.0f}" + ("" if r.kinematics_ok else " (?)")) if np.isfinite(r.sigma) else "—"
            self.res_tree.insert("", "end", iid=aid, values=(aid, r.role, r.status, r.method, f"{r.z:.5f}",
                                                             f"{r.z_err:.5f}", sig, "; ".join(r.flags)))
        txt = "\n".join(r.detail.get("provenance_text", [])) if r is not None else f"'{aid}' has not been fitted yet."
        self.fit_summary.configure(state="normal")
        self.fit_summary.delete("1.0", "end")
        self.fit_summary.insert("1.0", txt or "(no summary available)")
        self.fit_summary.configure(state="disabled")
        none_msg = f"'{aid}' has not been fitted yet."
        for key, fn in (("fit_cont", P.fig_continuum_fit), ("fit_lines", P.fig_lines), ("fit_boot", P.fig_bootstrap)):
            fig = fn(r) if r is not None else None
            if fig is not None:
                self.figs[key].show(fig)
            else:
                self.figs[key].clear(none_msg if r is None else "Not produced for this spectrum/method.")
        if aid in self.pipe.spectra:
            self.figs["fit_spec"].show(P.fig_spectra(self.pipe, [aid], float(self.pipe.cfg["display_smooth_A"]),
                                                     windows=self.pipe.engine.static_windows() if self.pipe.engine else None))
        else:
            self.figs["fit_spec"].clear("Run the extraction first.")

    def _show_results(self, aid):
        """After a fit: refresh the view of the selected spectrum, update the (system-wide) ladder, log notes."""
        if aid in self.pipe.spectra:
            self.fit_target.set(aid)
        self._refresh_fit_view()
        for m in self.pipe.system_messages:
            self.write("  * " + m)
        fig = P.fig_ladder(self.pipe.system_rows)
        if fig is not None:
            self.figs["fit_ladder"].show(fig)
        else:
            self.figs["fit_ladder"].clear("No fitted redshifts yet.")
        r = self.pipe.results.get(aid)
        if r is not None:
            for f in r.flags:
                self.write(f"[{aid}] ⚠ {f}")
            for n in r.notes:
                self.write(f"[{aid}] {n}")

    # ---- output
    def do_save(self):
        if not self.pipe.spectra:
            messagebox.showinfo("Nothing to save", "Run the extraction (and fits) first."); return
        from ..report import save_outputs
        def work(pr):
            return save_outputs(self.pipe, pr)
        def done(ok, out):
            self.write(f"Outputs written to {out}")
            for st in ("inputs", "output"):
                if "output_dir" in self.forms[st].rows:
                    self.forms[st].rows["output_dir"].set(out)
            self.refresh_output_preview(saved=True)
        self.run_task("Saving", work, done)

    def open_out(self):
        out = os.path.abspath(os.path.expanduser(self.pipe.cfg.get("output_dir", ".")))
        if not os.path.isdir(out):
            messagebox.showinfo("Not yet", "The output folder does not exist yet."); return
        import subprocess
        subprocess.Popen(["open" if sys.platform == "darwin" else ("explorer" if os.name == "nt" else "xdg-open"), out])

    def run_everything(self):
        if not self._apply_cfg():
            return
        def work(pr):
            self.pipe.reset()
            self.pipe.update_config(dict(apertures=self.ap_rows))
            self.pipe.set_apertures(self.ap_rows)
            ap = list(self.pipe.apertures)
            self.pipe.run_all(pr)
        def done(ok, _):
            self._show_load_figs()
            self._refresh_fit_targets()
            self.figs["ex_spec"].show(P.fig_spectra(self.pipe))
            self._refresh_ex_table()
            self._show_results(self.fit_target.get())
        self.run_task("Running all stages", work, done)


def launch(config_path: Optional[str] = None):
    app = App(config_path)
    app.mainloop()


if __name__ == "__main__":
    launch(sys.argv[1] if len(sys.argv) > 1 else None)
