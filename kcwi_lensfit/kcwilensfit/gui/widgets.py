"""Reusable Tk widgets: described input fields, scrollable forms, embedded matplotlib figures."""
from __future__ import annotations

import os
import tkinter as tk
from tkinter import filedialog, ttk
from typing import Any, Callable, Dict, List, Optional

from ..schema import FieldSpec, coerce

GREY = "#5a5a5a"


_LAST_DIR = {"path": ""}


def remember_dir(path: str) -> None:
    """Remember the folder of a chosen file/folder so the NEXT dialog opens there."""
    if not path:
        return
    p = os.path.expanduser(str(path))
    d = p if os.path.isdir(p) else os.path.dirname(p)
    if d and os.path.isdir(d):
        _LAST_DIR["path"] = d


def start_dir(current: str = "") -> str:
    """Where a dialog should open: the folder already in the field, else the last-used folder, else home."""
    p = os.path.expanduser(str(current or ""))
    if p:
        d = p if os.path.isdir(p) else os.path.dirname(p)
        if d and os.path.isdir(d):
            return d
    return _LAST_DIR["path"] or os.path.expanduser("~")


def safe_filetypes(patterns):
    """File-dialog filters that cannot crash Tk.

    macOS Tk 8.6 aborts the whole process (NSInvalidArgumentException) for filters containing
    compound extensions such as '*.fits.gz' or space-joined pattern strings, so on macOS no filter
    is used at all; elsewhere only simple single-dot patterns are passed, as a tuple.
    """
    import sys
    if sys.platform == "darwin" or not patterns:
        return None
    simple = tuple(p for p in patterns if p.count(".") == 1 and " " not in p)
    return [("Files", simple), ("All files", "*")] if simple else None


class ScrollFrame(ttk.Frame):
    """A frame whose content scrolls vertically (mouse wheel / trackpad aware)."""

    def __init__(self, master, **kw):
        super().__init__(master, **kw)
        self.canvas = tk.Canvas(self, highlightthickness=0, borderwidth=0)
        self.vbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.inner = ttk.Frame(self.canvas)
        self.inner.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))
        self._win = self.canvas.create_window((0, 0), window=self.inner, anchor="nw")
        self.canvas.bind("<Configure>", lambda e: self.canvas.itemconfigure(self._win, width=e.width))
        self.canvas.configure(yscrollcommand=self.vbar.set)
        self.canvas.pack(side="left", fill="both", expand=True)
        self.vbar.pack(side="right", fill="y")
        for w in (self.canvas, self.inner):
            w.bind("<Enter>", self._bind_wheel)
            w.bind("<Leave>", self._unbind_wheel)

    def _bind_wheel(self, _):
        self.canvas.bind_all("<MouseWheel>", self._on_wheel)
        self.canvas.bind_all("<Button-4>", lambda e: self.canvas.yview_scroll(-2, "units"))
        self.canvas.bind_all("<Button-5>", lambda e: self.canvas.yview_scroll(2, "units"))

    def _unbind_wheel(self, _):
        self.canvas.unbind_all("<MouseWheel>")
        self.canvas.unbind_all("<Button-4>")
        self.canvas.unbind_all("<Button-5>")

    def _on_wheel(self, e):
        d = e.delta
        step = -1 * (d // 120 if abs(d) >= 120 else (1 if d > 0 else -1))
        self.canvas.yview_scroll(int(step) * 2, "units")


class FieldRow(ttk.Frame):
    """One setting: label, input widget pre-filled with its default, unit and a visible description."""

    def __init__(self, master, spec: FieldSpec, on_change: Optional[Callable[[], None]] = None):
        super().__init__(master)
        self.spec = spec
        self._on_change = on_change
        self.columnconfigure(1, weight=1)
        ttk.Label(self, text=spec.label, font=("TkDefaultFont", 10, "bold")).grid(row=0, column=0, sticky="w", padx=(0, 8))
        self.var: tk.Variable
        k = spec.kind
        if k == "bool":
            self.var = tk.BooleanVar(value=bool(spec.default))
            w = ttk.Checkbutton(self, variable=self.var, text="enabled", command=self._changed)
            w.grid(row=0, column=1, sticky="w")
        elif k == "choice":
            self.var = tk.StringVar(value=str(spec.default))
            w = ttk.Combobox(self, textvariable=self.var, values=list(spec.choices or []), state="readonly", width=24)
            w.grid(row=0, column=1, sticky="w")
            w.bind("<<ComboboxSelected>>", lambda e: self._changed())
        elif k == "text":
            self.var = None
            w = tk.Text(self, height=3, width=40, wrap="word")
            w.insert("1.0", str(spec.default))
            w.grid(row=0, column=1, sticky="ew")
            w.bind("<FocusOut>", lambda e: self._changed())
            self._text = w
        else:
            self.var = tk.StringVar(value=self._fmt(spec.default))
            if k in ("path", "dir"):
                fr = ttk.Frame(self)
                fr.grid(row=0, column=1, sticky="ew")
                fr.columnconfigure(0, weight=1)
                e = ttk.Entry(fr, textvariable=self.var)
                e.grid(row=0, column=0, sticky="ew")
                ttk.Button(fr, text="Browse…", width=8, command=self._browse).grid(row=0, column=1, padx=(4, 0))
            else:
                e = ttk.Entry(self, textvariable=self.var, width=26)
                e.grid(row=0, column=1, sticky="w")
            e.bind("<FocusOut>", lambda ev: self._changed())
            e.bind("<Return>", lambda ev: self._changed())
        if spec.unit:
            ttk.Label(self, text=spec.unit, foreground=GREY).grid(row=0, column=2, sticky="w", padx=(6, 0))
        self.desc = ttk.Label(self, text=spec.help, foreground=GREY, wraplength=520, justify="left",
                              font=("TkDefaultFont", 9))
        self.desc.grid(row=1, column=0, columnspan=3, sticky="w", pady=(1, 8))
        self.bind("<Configure>", self._rewrap)

    def _rewrap(self, e):
        self.desc.configure(wraplength=max(240, e.width - 20))

    @staticmethod
    def _fmt(v):
        if isinstance(v, float):
            return f"{v:g}"
        return "" if v is None else str(v)

    def _browse(self):
        init = start_dir(self.var.get())
        if self.spec.kind == "dir":
            p = filedialog.askdirectory(title=self.spec.label, initialdir=init, mustexist=False)
        else:
            ft = safe_filetypes(self.spec.filetypes)
            kw = {"filetypes": ft} if ft else {}
            p = filedialog.askopenfilename(title=self.spec.label, initialdir=init, **kw)
        if p:
            remember_dir(p)
            self.var.set(p)
            self._changed()

    def _changed(self):
        if self._on_change:
            self._on_change()

    def get(self) -> Any:
        if self.spec.kind == "text":
            raw = self._text.get("1.0", "end").strip()
        else:
            raw = self.var.get()
        if self.spec.kind in ("int", "float") and str(raw).strip() == "":
            raw = self.spec.default
        return coerce(self.spec, raw)

    def set(self, value) -> None:
        if self.spec.kind == "text":
            self._text.delete("1.0", "end")
            self._text.insert("1.0", str(value or ""))
        elif self.spec.kind == "bool":
            self.var.set(bool(value))
        else:
            self.var.set(self._fmt(value))
            if self.spec.kind == "path" and value:       # input files only: an output folder must not hijack the start dir
                remember_dir(str(value))


class Form(ttk.Frame):
    """All fields of one stage, grouped under sub-headings, with a collapsible 'advanced' block."""

    def __init__(self, master, specs: List[FieldSpec], on_change=None):
        super().__init__(master)
        self.rows: Dict[str, FieldRow] = {}
        self.sf = ScrollFrame(self)
        self.sf.pack(fill="both", expand=True)
        body = self.sf.inner
        body.columnconfigure(0, weight=1)
        group = None
        adv_frame = None
        adv_specs = [s for s in specs if s.advanced]
        r = 0
        for s in specs:
            if s.advanced:
                continue
            if s.group != group:
                group = s.group
                if group:
                    ttk.Separator(body).grid(row=r, column=0, sticky="ew", pady=(8, 2)); r += 1
                    ttk.Label(body, text=group, font=("TkDefaultFont", 11, "bold"), foreground="#1f4e79").grid(
                        row=r, column=0, sticky="w", pady=(0, 4)); r += 1
            fr = FieldRow(body, s, on_change)
            fr.grid(row=r, column=0, sticky="ew", padx=6)
            self.rows[s.key] = fr
            r += 1
        if adv_specs:
            self._adv_open = tk.BooleanVar(value=False)
            ttk.Separator(body).grid(row=r, column=0, sticky="ew", pady=(10, 2)); r += 1
            self._adv_btn = ttk.Checkbutton(body, text=f"Show {len(adv_specs)} advanced settings (the defaults suit most data)",
                                            variable=self._adv_open,
                                            command=self._toggle_adv)
            self._adv_btn.grid(row=r, column=0, sticky="w"); r += 1
            self._adv_frame = ttk.Frame(body)
            self._adv_frame.columnconfigure(0, weight=1)
            self._adv_row = r
            rr = 0
            g2 = None
            for s in adv_specs:
                if s.group != g2:
                    g2 = s.group
                    if g2:
                        ttk.Label(self._adv_frame, text=g2 + " (advanced)", font=("TkDefaultFont", 10, "bold"),
                                  foreground="#1f4e79").grid(row=rr, column=0, sticky="w", pady=(6, 2)); rr += 1
                fr = FieldRow(self._adv_frame, s, on_change)
                fr.grid(row=rr, column=0, sticky="ew", padx=6)
                self.rows[s.key] = fr
                rr += 1

    def _toggle_adv(self):
        if self._adv_open.get():
            self._adv_frame.grid(row=self._adv_row, column=0, sticky="ew")
        else:
            self._adv_frame.grid_forget()

    def collect(self) -> Dict[str, Any]:
        out = {}
        for k, row in self.rows.items():
            out[k] = row.get()
        return out

    def populate(self, cfg: Dict[str, Any]) -> None:
        for k, row in self.rows.items():
            if k in cfg:
                row.set(cfg[k])


class FigureFrame(ttk.Frame):
    """Embeds a matplotlib Figure (created without pyplot) with the standard navigation toolbar."""

    def __init__(self, master):
        super().__init__(master)
        self.canvas = None
        self.toolbar = None
        self.placeholder = ttk.Label(self, text="(no figure yet)", foreground=GREY)
        self.placeholder.pack(expand=True)
        self._click_cb = None

    def show(self, fig) -> None:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
        for w in self.winfo_children():
            w.destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self)
        self.toolbar = NavigationToolbar2Tk(self.canvas, self, pack_toolbar=False)
        self.toolbar.update()
        self.toolbar.pack(side="bottom", fill="x")
        self.canvas.get_tk_widget().pack(side="top", fill="both", expand=True)
        if self._click_cb:
            self.canvas.mpl_connect("button_press_event", self._click_cb)
        self.canvas.draw()

    def clear(self, message: str = "(nothing to show for this spectrum yet)") -> None:
        """Remove the current figure and show a placeholder."""
        for w in self.winfo_children():
            w.destroy()
        self.canvas = None
        self.toolbar = None
        ttk.Label(self, text=message, foreground=GREY).pack(expand=True)

    def on_click(self, cb) -> None:
        self._click_cb = cb
        if self.canvas is not None:
            self.canvas.mpl_connect("button_press_event", cb)


class LogPane(ttk.Frame):
    def __init__(self, master, height=9):
        super().__init__(master)
        self.text = tk.Text(self, height=height, wrap="word", state="disabled", font=("Menlo", 9) if os.name != "nt" else ("Consolas", 9))
        sb = ttk.Scrollbar(self, command=self.text.yview)
        self.text.configure(yscrollcommand=sb.set)
        self.text.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        self.text.tag_configure("warn", foreground="#b36b00")
        self.text.tag_configure("err", foreground="#b00020")

    def write(self, msg: str) -> None:
        self.text.configure(state="normal")
        tag = "err" if ("FAILED" in msg or "Error" in msg) else ("warn" if ("[warning]" in msg or "WARNING" in msg or "⚠" in msg) else None)
        self.text.insert("end", msg + "\n", tag)
        self.text.see("end")
        self.text.configure(state="disabled")
