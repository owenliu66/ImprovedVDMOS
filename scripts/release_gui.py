"""Release GUI (session 30, Phase E): five columns over the package API.

  1. FILES   : pick CSVs, preset, ref temp, Parse (ingest + report)
  2. DATA    : parsed tables + rows (displayNum-formatted, crash fixed)
  3. PARAMS  : free-checkbox parameter tables (channel + diode), fit
               buttons -- channel@ref, diode, TC; double-click a value
               cell to edit any parameter by hand; Plot diode / Plot TC
               -- fits run on a worker thread so the window stays live
  4. CAPS    : capacitance points + ramp-sum resolution knob (n_sub)
  5. EXPORT  : plots (IV measured-vs-model, caps), optional recovery
               extraction (LTspice grid -- slow), Save subckt: folder
               dialog, writes ONLY <part>_fast.subckt (+ _full.subckt
               when recovery knobs exist) -- no caps/cds artifacts

Thin shell: every number comes from vdmos_fit (ingest/presets/stages/
release/export_subckt); nothing is re-implemented here.

    python scripts/release_gui.py
"""

from __future__ import annotations

import math
import queue
import sys
import threading
import tkinter as tk
import traceback
from pathlib import Path
from tkinter import filedialog, messagebox, ttk
from tkinter import font as tkfont

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np  # noqa: E402
from matplotlib.figure import Figure  # noqa: E402

from vdmos_fit.diode import vdmos_diode  # noqa: E402
from vdmos_fit.display import displayNum, parseEng  # noqa: E402
from vdmos_fit.export_subckt import RecovParams  # noqa: E402
from vdmos_fit.ingest import (  # noqa: E402
    CapTable,
    IVTable,
    QgTable,
    RdsSpec,
    RecTable,
    VthSpec,
    load_tables,
)
from vdmos_fit.model import vdmos_id  # noqa: E402
from vdmos_fit.params import BOUNDS, Params  # noqa: E402
from vdmos_fit.presets import PRESETS  # noqa: E402
from vdmos_fit.release import (  # noqa: E402
    caps_model_from_tables,
    export_part,
    merge_ivs,
    recov_from_table,
)
from vdmos_fit.stages import (  # noqa: E402
    channel_rows,
    diode_rows,
    fit_channel_ref,
    fit_channel_tc,
    fit_diode_stage,
    model_ron,
    model_vth,
    resolve_test_current,
)

ROW_CAP = 400
# Params fields that belong to the diode tree (everything else in the
# channel tree); the GUI shows ALL model parameters under both presets,
# marking preset-unused ones as locked rather than hiding them (user,
# session 30).
DIODE_FIELDS = frozenset({"i_s", "k_n", "i_sr", "k_nr", "i_kf", "a_ikf",
                          "r_diode", "a_rdiode", "b_rdiode", "t_nom",
                          "k_xti", "v_eg"})
COLS = {
    "IV": ("V_gs", "V_ds", "I_ds", "T_j"),
    "CAP": ("V", "C"),
    "QG": ("V_gs", "Q_g", "condition"),
    "RDS": ("V_gs", "R_ds(on)", "I_test", "T_j", "condition"),
    "VTH": ("T_j", "V_th", "I_test"),
    "REC": ("Q_rr", "t_rr", "I_RM", "I_sd", "di/dt"),
}


def table_kind(t):
    for kind, cls in (("IV", IVTable), ("CAP", CapTable), ("QG", QgTable),
                      ("RDS", RdsSpec), ("VTH", VthSpec), ("REC", RecTable)):
        if isinstance(t, cls):
            return kind
    return "?"


def table_rows(t):
    if isinstance(t, IVTable):
        return list(zip(t.vgs, t.vds, t.ids, t.tj, strict=True))[:ROW_CAP]
    if isinstance(t, CapTable):
        return list(zip(t.v, t.c, strict=True))[:ROW_CAP]
    if isinstance(t, QgTable):
        return list(zip(t.vgs, t.qg, t.condition, strict=True))[:ROW_CAP]
    if isinstance(t, RdsSpec):
        return list(zip(t.vgs, t.rds, t.i_test, t.tj, t.condition,
                        strict=True))[:ROW_CAP]
    if isinstance(t, VthSpec):
        return list(zip(t.tj, t.vth, t.i_test, strict=True))[:ROW_CAP]
    if isinstance(t, RecTable):
        return list(zip(t.qrr, t.trr, t.irm, t.i_f, t.didt,
                        strict=True))[:ROW_CAP]
    return []


def table_nrows(t):
    return len(next(iter(t.__dict__.values()), ()))


class App:
    def __init__(self, root: tk.Tk):
        self.root = root
        root.title("ImprovedVDMOS release pipeline")
        # Row height from the ACTUAL font metrics: the default Treeview
        # rowheight clips text under Windows display scaling (session 30
        # field report).
        style = ttk.Style(root)
        fam = style.lookup("Treeview", "font") or "TkDefaultFont"
        linespace = tkfont.nametofont(fam).metrics("linespace")
        rh = max(int(style.lookup("Treeview", "rowheight") or 0),
                 int(linespace * 1.55))
        style.configure("Treeview", rowheight=rh)
        style.configure("Treeview.Heading", rowpadding=4)
        self.q: queue.Queue = queue.Queue()
        self.paths: list[str] = []
        self.tables: list = []
        self.report: list[str] = []
        # Live parameter state: the parameter TREES are the source of
        # truth (double-click a value cell to edit by hand, session 30);
        # fits write back into these dicts.
        self.ch_vals: dict[str, float] = {}
        self.dio_vals: dict[str, float] = {}
        # per-parameter optimisation limits, shown and editable in the
        # trees (defaults from params.BOUNDS; session 30)
        self.ch_min: dict[str, float] = {}
        self.ch_max: dict[str, float] = {}
        self.dio_min: dict[str, float] = {}
        self.dio_max: dict[str, float] = {}
        self.ch_fitted = False
        self.dio_fitted = False
        self.tc: dict = {}
        self.recov: RecovParams | None = None
        self.caps_model = None
        self.rmse: dict = {}
        self.free_ch: set[str] = set()
        self.free_dio: set[str] = set()
        self.busy = False

        self.preset_name = tk.StringVar(value="trench")
        self.ref_temp = tk.StringVar(value="25")
        self.part_name = tk.StringVar(value="MYFET")
        self.n_sub = tk.IntVar(value=12)
        self.status = tk.StringVar(value="pick CSVs and Parse")
        # plot axis controls (session 31): blank range = autoscale, scale
        # "auto" = keep the plot's own default (e.g. caps are log-y).
        self.ax_xmin = tk.StringVar(value="")
        self.ax_xmax = tk.StringVar(value="")
        self.ax_ymin = tk.StringVar(value="")
        self.ax_ymax = tk.StringVar(value="")
        self.ax_xscale = tk.StringVar(value="auto")
        self.ax_yscale = tk.StringVar(value="auto")

        panes = ttk.Panedwindow(root, orient=tk.HORIZONTAL)
        panes.pack(fill=tk.BOTH, expand=True)
        self._build_files(panes)
        self._build_data(panes)
        self._build_params(panes)
        self._build_caps(panes)
        self._build_export(panes)
        ttk.Label(root, textvariable=self.status,
                  anchor="w").pack(fill=tk.X, side=tk.BOTTOM)
        self._apply_preset()
        root.after(120, self._poll)

    # ---- column 1: files ----------------------------------------------------
    def _build_files(self, panes):
        f = ttk.Labelframe(panes, text="1 Files")
        panes.add(f, weight=1)
        self.file_list = tk.Listbox(f, height=14, width=34,
                                    exportselection=False)
        self.file_list.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        row = ttk.Frame(f)
        row.pack(fill=tk.X, padx=4)
        ttk.Button(row, text="Add CSVs",
                   command=self._add_files).pack(side=tk.LEFT)
        ttk.Button(row, text="Remove",
                   command=self._del_file).pack(side=tk.LEFT, padx=4)
        row2 = ttk.Frame(f)
        row2.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(row2, text="preset").pack(side=tk.LEFT)
        cb = ttk.Combobox(row2, textvariable=self.preset_name, width=8,
                          values=("trench", "sj", "sic"), state="readonly")
        cb.pack(side=tk.LEFT, padx=4)
        cb.bind("<<ComboboxSelected>>", lambda e: self._apply_preset())
        ttk.Label(row2, text="ref °C").pack(side=tk.LEFT)
        ttk.Entry(row2, textvariable=self.ref_temp, width=6).pack(
            side=tk.LEFT, padx=4)
        ttk.Button(f, text="Parse", command=self._parse).pack(
            fill=tk.X, padx=4, pady=4)
        self.report_txt = tk.Text(f, height=10, width=34, state="disabled",
                                  wrap="word")
        self.report_txt.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

    def _add_files(self):
        paths = filedialog.askopenfilenames(
            initialdir=str(ROOT / "data"), filetypes=[("CSV", "*.csv")])
        self.paths.extend(p for p in paths if p not in self.paths)
        self._refresh_file_list()

    def _del_file(self):
        sel = list(self.file_list.curselection())
        for i in reversed(sel):
            del self.paths[i]
        self._refresh_file_list()

    def _refresh_file_list(self):
        self.file_list.delete(0, tk.END)
        for p in self.paths:
            self.file_list.insert(tk.END, Path(p).name)

    def _parse(self):
        if not self.paths:
            messagebox.showerror("No files", "Add CSVs first.")
            return
        try:
            self.tables, self.report = load_tables(
                self.paths, float(self.ref_temp.get()))
        except Exception as exc:
            messagebox.showerror("Parse failed", str(exc))
            return
        self._log(self.report)
        self.tbl_tree.delete(*self.tbl_tree.get_children())
        for t in self.tables:
            self.tbl_tree.insert(
                "", tk.END, values=(table_kind(t), Path(t.source).name,
                                    table_nrows(t),
                                    "; ".join(getattr(t, "notes", ()))))
        self._build_caps_model()
        self.status.set(f"parsed {len(self.tables)} tables")

    # ---- column 2: parsed data ----------------------------------------------
    def _build_data(self, panes):
        f = ttk.Labelframe(panes, text="2 Parsed data")
        panes.add(f, weight=2)
        self.tbl_tree = ttk.Treeview(f, columns=("kind", "file", "rows",
                                                 "notes"),
                                     show="headings", height=8)
        for c, w in (("kind", 45), ("file", 150), ("rows", 55),
                     ("notes", 170)):
            self.tbl_tree.heading(c, text=c)
            self.tbl_tree.column(c, width=w)
        self.tbl_tree.pack(fill=tk.X, padx=4, pady=4)
        self.tbl_tree.bind("<<TreeviewSelect>>", self._show_rows)
        self.row_tree = ttk.Treeview(f, show="headings", height=18)
        self.row_tree.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)

    def _show_rows(self, _event=None):
        sel = self.tbl_tree.selection()
        if not sel:
            return
        idx = self.tbl_tree.index(sel[0])
        t = self.tables[idx]
        cols = COLS[table_kind(t)]
        self.row_tree.delete(*self.row_tree.get_children())
        self.row_tree["columns"] = cols
        for c in cols:
            self.row_tree.heading(c, text=c)
            self.row_tree.column(c, width=80)
        for r in table_rows(t):
            self.row_tree.insert("", tk.END, values=tuple(
                displayNum(v) if isinstance(v, float) else str(v)
                for v in r))

    # ---- column 3: parameters + fits -----------------------------------------
    def _build_params(self, panes):
        f = ttk.Labelframe(panes, text="3 Parameters & fits")
        panes.add(f, weight=2)
        self.ch_tree = ttk.Treeview(
            f, columns=("param", "min", "value", "max", "free"),
            show="headings", height=14)
        for c, w in (("param", 80), ("min", 65), ("value", 90),
                     ("max", 65), ("free", 40)):
            self.ch_tree.heading(c, text=c)
            self.ch_tree.column(c, width=w)
        self.ch_tree.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.ch_tree.bind("<Button-1>",
                          lambda e: self._toggle(e, self.ch_tree,
                                                 self.free_ch, "ch"))
        self.ch_tree.bind("<Double-1>",
                          lambda e: self._edit_cell(e, self.ch_tree, {
                              2: self.ch_min, 3: self.ch_vals,
                              4: self.ch_max}))
        row = ttk.Frame(f)
        row.pack(fill=tk.X, padx=4)
        self.btn_ch = ttk.Button(row, text="Fit channel @ref",
                                 command=self._fit_channel)
        self.btn_ch.pack(side=tk.LEFT)
        self.btn_tc = ttk.Button(row, text="Fit TC", command=self._fit_tc)
        self.btn_tc.pack(side=tk.LEFT, padx=4)
        ttk.Button(row, text="Plot TC",
                   command=self._plot_tc).pack(side=tk.LEFT, padx=4)
        self.dio_tree = ttk.Treeview(
            f, columns=("param", "min", "value", "max", "free"),
            show="headings", height=8)
        for c, w in (("param", 80), ("min", 65), ("value", 90),
                     ("max", 65), ("free", 40)):
            self.dio_tree.heading(c, text=c)
            self.dio_tree.column(c, width=w)
        self.dio_tree.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        self.dio_tree.bind("<Button-1>",
                           lambda e: self._toggle(e, self.dio_tree,
                                                  self.free_dio, "dio"))
        self.dio_tree.bind("<Double-1>",
                           lambda e: self._edit_cell(e, self.dio_tree, {
                               2: self.dio_min, 3: self.dio_vals,
                               4: self.dio_max}))
        row2 = ttk.Frame(f)
        row2.pack(fill=tk.X, padx=4, pady=(0, 4))
        self.btn_dio = ttk.Button(row2, text="Fit diode",
                                  command=self._fit_diode)
        self.btn_dio.pack(side=tk.LEFT)
        ttk.Button(row2, text="Plot diode",
                   command=self._plot_diode).pack(side=tk.LEFT, padx=4)
        self.fit_lbl = tk.StringVar(value="")
        ttk.Label(row2, textvariable=self.fit_lbl).pack(side=tk.LEFT,
                                                        padx=8)

    def _apply_preset(self):
        p = PRESETS[self.preset_name.get()]
        self.free_ch = set(p.channel_free25)
        self.free_dio = set(p.diode_free)
        base = Params().replace(**{**p.channel_fixed, **p.channel_start})
        vals = base.as_dict()
        if "v_ovmax" not in self.free_ch:
            vals["v_ovmax"] = 0.0   # user rule: charge cap OFF unless fitted
        self.ch_vals = {k: float(vals[k]) for k in vals
                        if k not in DIODE_FIELDS}
        self.ch_min, self.ch_max = self._default_bounds(self.ch_vals, BOUNDS)
        self._sync_ktheta_min()
        self.ch_fitted = False
        # show EVERY model parameter (dataclass order); the preset only
        # decides which are free -- unused ones stay visible, locked
        self._fill_tree(self.ch_tree, list(self.ch_vals), self.ch_vals,
                        self.free_ch, self.ch_min, self.ch_max)
        dbase = Params().replace(**{**p.diode_fixed, **p.diode_start})
        dvals = dbase.as_dict()
        self.dio_vals = {k: float(dvals[k]) for k in dvals
                         if k in DIODE_FIELDS}
        # the diode tree shows the EFFECTIVE limits: BOUNDS widened by
        # the preset's diode_bounds (what fit_diode_stage actually uses)
        dtab = {**BOUNDS, **p.diode_bounds}
        self.dio_min, self.dio_max = self._default_bounds(self.dio_vals,
                                                          dtab)
        self.dio_fitted = False
        self._fill_tree(self.dio_tree, list(self.dio_vals), self.dio_vals,
                        self.free_dio, self.dio_min, self.dio_max)

    @staticmethod
    def _default_bounds(keys, table):
        lo = {k: float(table.get(k, (-math.inf, math.inf))[0])
              for k in keys}
        hi = {k: float(table.get(k, (-math.inf, math.inf))[1])
              for k in keys}
        return lo, hi

    def _sync_ktheta_min(self):
        """k_theta's safe lower bound depends on the knee form:
        - mobility form (k_kneeform = 0): the knee-divisor zero -1/v_ovmax
          (a negative k_theta widens the knee; below this 1 + k_theta*v_ov_eff2
          <= 0 at the max overdrive and v_knee blows up). DYNAMIC in v_ovmax.
          With the charge cap off (v_ovmax = 0) the bound is unbounded, so keep
          the static guard.
        - linear knee-scale form (k_kneeform = 1): v_knee = k_theta*v_ov_eff2,
          so k_theta is a non-negative multiplier -> lower bound 0 (negative
          would floor the knee at v_kneemin, i.e. switch the knee off)."""
        if "k_theta" not in self.ch_min:
            return
        if self.ch_vals.get("k_kneeform", 0.0) >= 0.5:
            self.ch_min["k_theta"] = 0.0
        else:
            vom = self.ch_vals.get("v_ovmax", 0.0)
            if vom > 0.0:
                self.ch_min["k_theta"] = -1.0 / vom
        for row in self.ch_tree.get_children():
            if self.ch_tree.item(row, "values")[0] == "k_theta":
                self.ch_tree.set(row, "min", displayNum(self.ch_min["k_theta"]))
                break

    def _fill_tree(self, tree, keys, values, free, lo, hi):
        tree.delete(*tree.get_children())
        for k in keys:
            tree.insert("", tk.END,
                        values=(k, displayNum(lo[k]),
                                displayNum(values.get(k, 0.0)),
                                displayNum(hi[k]),
                                "X" if k in free else ""))

    def _toggle(self, event, tree, free, _which):
        row = tree.identify_row(event.y)
        col = tree.identify_column(event.x)
        if not row or col != "#5":
            return
        key = tree.item(row, "values")[0]
        if key in free:
            free.discard(key)
            tree.set(row, "free", "")
        else:
            free.add(key)
            tree.set(row, "free", "X")

    def _update_tree(self, tree, vals: dict[str, float]):
        for row in tree.get_children():
            k = tree.item(row, "values")[0]
            tree.set(row, "value", displayNum(vals[k]))

    def _edit_cell(self, event, tree, colmaps: dict[int, dict[str, float]]):
        """Manual edit (session 30): double-click a min / value / max
        cell, type a number -- engineering suffixes accepted (2.483m,
        -123.5k, 1e-3, -inf/inf for the limits) -- Enter commits,
        Escape cancels."""
        row = tree.identify_row(event.y)
        col = tree.identify_column(event.x)
        if not row or col not in ("#2", "#3", "#4"):
            return
        vals = colmaps[int(col[1:])]
        key = tree.item(row, "values")[0]
        bbox = tree.bbox(row, col)
        if not bbox:
            return
        x, y, w, h = bbox
        ent = tk.Entry(tree)
        ent.place(x=x, y=y, width=w, height=h)
        ent.insert(0, tree.set(row, col))
        ent.focus_set()
        ent.select_range(0, tk.END)

        def commit(_=None):
            txt = ent.get()
            ent.destroy()
            try:
                val = parseEng(txt)
            except ValueError:
                messagebox.showerror(
                    "Bad number", f"{txt!r} is not a number (try e.g. "
                    "2.483m, -123.5k, 1e-3, or inf/-inf for limits)")
                return
            vals[key] = val
            tree.set(row, col, displayNum(val))
            # editing v_ovmax or the knee-form selector moves k_theta's dynamic
            # lower bound (-1/v_ovmax for the mobility form, 0 for the scale form)
            if key in ("v_ovmax", "k_kneeform") and "k_theta" in colmaps[2]:
                self._sync_ktheta_min()

        ent.bind("<Return>", commit)
        ent.bind("<FocusOut>", commit)
        ent.bind("<Escape>", lambda e: ent.destroy())

    def _bounds(self, keys, lo: dict[str, float], hi: dict[str, float]):
        """Validate + assemble bounds_overrides for fit(); raises
        ValueError on an inverted box."""
        out = {}
        for k in keys:
            a, b = lo[k], hi[k]
            if not a < b:
                raise ValueError(
                    f"{k}: min {displayNum(a)} must be below max "
                    f"{displayNum(b)}")
            out[k] = (a, b)
        return out

    def _channel_params(self) -> Params:
        return Params().replace(**self.ch_vals)

    def _diode_params(self) -> Params:
        return Params().replace(**self.dio_vals)

    def _iv(self):
        ivs = [t for t in self.tables if isinstance(t, IVTable)]
        return merge_ivs(ivs) if ivs else None

    def _specs(self):
        rds = next((t for t in self.tables if isinstance(t, RdsSpec)), None)
        vth = next((t for t in self.tables if isinstance(t, VthSpec)), None)
        return rds, vth

    def _tc_rds(self):
        """The Rds(on) table that spans the MOST temperatures -- the TC data.
        The TC fit and the TC plot must use the SAME table. The first-parsed
        RdsSpec is often the 25 C-only chord file, which cannot constrain the
        drain-resistance temperature coefficients (a_rd, b_rd) at all: at
        T = ref every row has dT = 0, so r_d(T) = r_d25 regardless of a_rd/
        b_rd and they drift to their bounds (session 32: TC RMSE read 2.1 %
        while the high-T Rds(on) curves were garbage)."""
        cands = [t for t in self.tables if isinstance(t, RdsSpec)]
        if not cands:
            return None
        return max(cands, key=lambda t: len(set(t.tj.tolist())))

    def _fit_channel(self):
        iv = self._iv()
        if iv is None:
            messagebox.showerror("No data", "Parse an IV CSV first.")
            return
        rds, _ = self._specs()
        p = PRESETS[self.preset_name.get()]
        if "r_sk" in self.free_ch and "r_s" not in self.free_ch:
            messagebox.showwarning(
                "Kelvin note",
                "r_sk is optimised as the ratio r_sk/r_s (the kelvin tap is "
                "a segment OF the source path). With r_s locked, r_sk can "
                "only scale that locked value -- free r_s too to fit the "
                "kelvin split itself.")
        free = tuple(sorted(self.free_ch))
        try:
            bounds = self._bounds(self.ch_vals, self.ch_min, self.ch_max)
        except ValueError as exc:
            messagebox.showerror("Bad limits", str(exc))
            return

        def work():
            res = fit_channel_ref(
                iv, rds, p, ref_temp=float(self.ref_temp.get()),
                free_override=free, bounds_override=bounds,
                start_override=dict(self.ch_vals),
                progress=lambda m: self.q.put(("progress",
                                               "channel@ref " + m)))
            return ("ch", res)
        self._run(work)

    def _fit_diode(self):
        ivs = [t for t in self.tables if isinstance(t, IVTable)]
        if not ivs:
            messagebox.showerror("No data", "Parse an IV CSV first.")
            return
        p = PRESETS[self.preset_name.get()]
        free = tuple(sorted(self.free_dio)) or None
        try:
            bounds = self._bounds(self.dio_vals, self.dio_min, self.dio_max)
        except ValueError as exc:
            messagebox.showerror("Bad limits", str(exc))
            return

        def work():
            res = fit_diode_stage(
                ivs, p, free_override=free, bounds_override=bounds,
                progress=lambda m: self.q.put(("progress", "diode " + m)))
            return ("dio", res)
        self._run(work)

    def _fit_tc(self):
        rds = self._tc_rds()
        vth = next((t for t in self.tables if isinstance(t, VthSpec)), None)
        if rds is None and vth is None:
            messagebox.showerror("No specs", "Parse an RDS or VTH CSV.")
            return
        if rds is not None and len(set(rds.tj.tolist())) < 2:
            messagebox.showwarning(
                "Rds(on) is single-temperature",
                "The Rds(on) data spans only one temperature, so the drain-"
                "resistance temperature coefficients (a_rd, b_rd) cannot be "
                "identified from it -- they will sit at their bounds and the "
                "high-temperature Rds(on) curves will be unreliable. Parse a "
                "multi-temperature Rds(on) CSV (Rds(on) vs V_gs at several "
                "T_j) for a meaningful TC fit.")
        p = PRESETS[self.preset_name.get()]

        def work():
            res = fit_channel_tc(
                self._channel_params(), rds, vth, p,
                progress=lambda m: self.q.put(("progress", "TC " + m)))
            return ("tc", res)
        self._run(work)

    # ---- column 4: caps -------------------------------------------------------
    def _build_caps(self, panes):
        f = ttk.Labelframe(panes, text="4 Caps")
        panes.add(f, weight=2)
        self.cap_tree = ttk.Treeview(f, columns=("branch", "v_ds", "cap"),
                                     show="headings", height=20)
        for c, w in (("branch", 60), ("v_ds", 80), ("cap", 90)):
            self.cap_tree.heading(c, text=c)
            self.cap_tree.column(c, width=w)
        self.cap_tree.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
        row = ttk.Frame(f)
        row.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(row, text="export resolution n_sub").pack(side=tk.LEFT)
        ttk.Spinbox(row, from_=5, to=40, textvariable=self.n_sub,
                    width=5).pack(side=tk.LEFT, padx=6)
        # recovery extraction lives with the caps data (user, session 30)
        self.btn_recov = ttk.Button(f, text="Extract recovery (slow)",
                                    command=self._extract_recov)
        self.btn_recov.pack(fill=tk.X, padx=4, pady=(0, 4))

    def _fill_caps_tree(self):
        self.cap_tree.delete(*self.cap_tree.get_children())
        for t in self.tables:
            if isinstance(t, CapTable):
                for v, c in zip(t.v, t.c, strict=True):
                    self.cap_tree.insert("", tk.END,
                                         values=(t.kind, displayNum(v),
                                                 displayNum(c)))

    def _build_caps_model(self):
        caps = [t for t in self.tables if isinstance(t, CapTable)]
        if not caps:
            self.caps_model = None
            return
        qgs = [t for t in self.tables if isinstance(t, QgTable)]
        self.caps_model, _ = caps_model_from_tables(caps, qgs, r_ds_on=0.0)
        self._fill_caps_tree()

    # ---- column 5: plots + export ---------------------------------------------
    def _build_export(self, panes):
        f = ttk.Labelframe(panes, text="5 Plots & export")
        panes.add(f, weight=3)
        row = ttk.Frame(f)
        row.pack(fill=tk.X, padx=4, pady=4)
        ttk.Label(row, text="part").pack(side=tk.LEFT)
        ttk.Entry(row, textvariable=self.part_name, width=16).pack(
            side=tk.LEFT, padx=4)
        row2 = ttk.Frame(f)
        row2.pack(fill=tk.X, padx=4)
        ttk.Button(row2, text="Plot IV",
                   command=self._plot_iv).pack(side=tk.LEFT)
        ttk.Button(row2, text="Plot caps",
                   command=self._plot_caps).pack(side=tk.LEFT, padx=4)
        self.btn_save = ttk.Button(row2, text="Save subckt",
                                   command=self._save)
        self.btn_save.pack(side=tk.LEFT, padx=4)
        # axis range + scale controls (session 31): blank = autoscale,
        # scale "auto" = keep the plot's own default. Engineering suffixes
        # accepted in the range boxes. Applies to the current plot.
        row3 = ttk.Frame(f)
        row3.pack(fill=tk.X, padx=4, pady=(4, 0))
        for lbl, var in (("x0", self.ax_xmin), ("x1", self.ax_xmax)):
            ttk.Label(row3, text=lbl).pack(side=tk.LEFT)
            ttk.Entry(row3, textvariable=var, width=7).pack(side=tk.LEFT,
                                                            padx=(1, 4))
        ttk.Combobox(row3, textvariable=self.ax_xscale, width=6,
                     values=("auto", "linear", "log"),
                     state="readonly").pack(side=tk.LEFT, padx=4)
        for lbl, var in (("y0", self.ax_ymin), ("y1", self.ax_ymax)):
            ttk.Label(row3, text=lbl).pack(side=tk.LEFT)
            ttk.Entry(row3, textvariable=var, width=7).pack(side=tk.LEFT,
                                                            padx=(1, 4))
        ttk.Combobox(row3, textvariable=self.ax_yscale, width=6,
                     values=("auto", "linear", "log"),
                     state="readonly").pack(side=tk.LEFT, padx=4)
        ttk.Button(row3, text="Apply axes",
                   command=self._reapply_axes).pack(side=tk.LEFT, padx=4)
        self.fig = Figure(figsize=(5, 4), dpi=90)
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
        self.canvas = FigureCanvasTkAgg(self.fig, master=f)
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True,
                                         padx=4, pady=4)

    def _apply_axes(self):
        """Apply the axis range/scale controls to every axis on the figure.
        An empty range entry autoscales that side; a filled entry pins it.
        Scale "auto" = keep the plot's own default (caps default to log-y).
        Engineering suffixes accepted in the boxes."""
        def lim(var):
            t = var.get().strip()
            if not t:
                return None
            try:
                return parseEng(t)
            except Exception:
                return None
        xs, ys = self.ax_xscale.get(), self.ax_yscale.get()

        def auto(a, axis):
            # re-enable autoscale (a prior manual set_xlim disabled it),
            # recompute the data limits, then apply -- so an empty box
            # restores autoscale instead of leaving the old range stuck.
            a.autoscale(enable=True, axis=axis)
            a.relim()
            a.autoscale_view(scalex=(axis == "x"), scaley=(axis == "y"))

        for ax in self.fig.axes:
            if xs == "log":
                ax.set_xscale("log")
            elif xs == "linear":
                ax.set_xscale("linear")
            if ys == "log":
                ax.set_yscale("log")
            elif ys == "linear":
                ax.set_yscale("linear")
            xl, xh = lim(self.ax_xmin), lim(self.ax_xmax)
            yl, yh = lim(self.ax_ymin), lim(self.ax_ymax)
            if xs == "log":
                xl = None if (xl is not None and xl <= 0) else xl
                xh = None if (xh is not None and xh <= 0) else xh
            if ys == "log":
                yl = None if (yl is not None and yl <= 0) else yl
                yh = None if (yh is not None and yh <= 0) else yh
            # autoscale each axis, then pin only the sides the user filled in
            auto(ax, "x")
            xlo, xhi = ax.get_xlim()
            ax.set_xlim(xl if xl is not None else xlo,
                        xh if xh is not None else xhi)
            auto(ax, "y")
            ylo, yhi = ax.get_ylim()
            ax.set_ylim(yl if yl is not None else ylo,
                        yh if yh is not None else yhi)

    def _reapply_axes(self):
        if self.fig.axes:
            self._apply_axes()
            self.canvas.draw()

    def _plot_iv(self):
        iv = self._iv()
        if iv is None:
            return
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        p25 = self._channel_params()
        D = channel_rows(iv, float(self.ref_temp.get()))
        for vg in sorted(set(D.v_gs.tolist())):
            m = D.v_gs == vg
            ax.plot(D.v_ds[m], D.i_d[m], ".", label=f"meas {vg:g} V")
            v = sorted(set(D.v_ds[m].tolist()))
            ax.plot(v, [float(vdmos_id(x, vg, p25)) for x in v], "-")
        ax.set_xlabel("v_DS [V]")
        ax.set_ylabel("i_D [A]")
        ax.legend(fontsize=7)
        self._apply_axes()
        self.canvas.draw()

    def _plot_caps(self):
        if self.caps_model is None:
            return
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        for t in self.tables:
            if isinstance(t, CapTable):
                ax.plot(t.v, t.c * 1e12, ".", label=f"meas {t.kind}")
        # The EXPORTED model recombined into the TERMINAL capacitance the
        # subckt actually realizes (session 30 crash-course: branch plots
        # were misleading -- the gd branch is indexed by v_gd ~= -V_ds,
        # the gs/ds tables by V_ds):
        #   Ciss = c_gs(v) + c_gd(-v)   Coss = c_ds(v) + c_gd(-v)
        #   Crss = c_gd(-v)             (branch-measured kinds directly)
        m = self.caps_model

        def exported(kind, v):
            if kind == "iss":
                return m.c_gs(v) + m.c_gd(-v)
            if kind == "oss":
                return m.c_ds(v) + m.c_gd(-v)
            if kind == "rss":
                return m.c_gd(-v)
            if kind == "gs":
                return m.c_gs(v)
            if kind == "ds":
                return m.c_ds(v)
            if kind == "gd":
                return m.c_gd(v)
            return None

        for t in self.tables:
            if isinstance(t, CapTable):
                v = np.linspace(float(min(t.v)), float(max(t.v)), 200)
                y = exported(t.kind, v)
                if y is not None:
                    ax.plot(v, y * 1e12, "-",
                            label=f"model {t.kind} (exported)")
        ax.set_yscale("log")
        ax.set_xlabel("V_DS [V]")
        ax.set_ylabel("C [pF]")
        ax.legend(fontsize=7)
        self._apply_axes()
        self.canvas.draw()

    def _plot_diode(self):
        """Measured diode rows (gate-off reverse IV) vs the diode card
        currently in the tree (fitted or hand-edited)."""
        self.fig.clear()
        ax = self.fig.add_subplot(111)
        ds = [diode_rows(t) for t in self.tables if isinstance(t, IVTable)]
        ds = [d for d in ds if d is not None]
        if not ds:
            messagebox.showinfo(
                "No diode rows", "The parsed IV tables carry no gate-off "
                "reverse rows (I_ds < 0 with V_gs <= 0).")
            self.canvas.draw()
            return
        temps: set[float] = set()
        for d in ds:
            temps |= set(d.temp_c.tolist())
            for T in sorted(set(d.temp_c.tolist())):
                m = d.temp_c == T
                ax.plot(d.v_sd[m], d.i_sd[m], ".", label=f"meas {T:g} C")
        p = self._diode_params()
        vs = np.linspace(0.2, float(max(d.v_sd.max() for d in ds)), 120)
        for T in sorted(temps):
            ax.plot(vs, [float(vdmos_diode(v, p, t_c=T)) for v in vs],
                    "-", label=f"model {T:g} C")
        ax.set_xlabel("v_SD forward [V]")
        ax.set_ylabel("i_SD forward [A]")
        ax.legend(fontsize=7)
        self._apply_axes()
        self.canvas.draw()

    def _plot_tc(self):
        """Measured Rds(on)/Vth against the datasheet table's own x column,
        with the TC law applied to the channel parameters currently in the
        tree. The Rds(on) panel uses the Rds(on) table that SPANS
        temperatures (the TC data), not the 25 C chord file, and plots
        against that CSV's first data column (V_gs for an Rds(on)-vs-V_gs
        table), one series per temperature. The Vth panel separates the two
        test currents instead of mixing them."""
        rds = self._tc_rds()
        vth = next((t for t in self.tables if isinstance(t, VthSpec)), None)
        if rds is None and vth is None:
            messagebox.showerror("No specs", "Parse an RDS or VTH CSV.")
            return
        p25 = self._channel_params()
        tcd = {n: float(self.tc.get(n, 0.0))
               for n in ("s_vth", "beta_p", "a_rd", "b_rd")}
        self.fig.clear()
        nup = 2 if (rds is not None and vth is not None) else 1
        if rds is not None:
            try:
                itest = resolve_test_current(rds, None)
            except ValueError as exc:
                messagebox.showerror("Test current", str(exc))
                return
            ax1 = self.fig.add_subplot(nup, 1, 1)
            xkind = rds.x_kind
            if xkind == "vgs":
                xarr, xlab = rds.vgs, "V_gs [V]"
                sarr, sfmt = rds.tj, "T_j {v:g} C"
            elif xkind == "ids":
                xarr, xlab = rds.i_test, "I_test [A]"
                sarr, sfmt = rds.vgs, "V_gs {v:g} V"
            else:  # tj, or unrecognised -> legacy T_j axis
                xarr, xlab = rds.tj, "T_j [C]"
                sarr, sfmt = rds.vgs, "V_gs {v:g} V"
            y = rds.rds * 1e3

            def model1(s, xv):
                if xkind == "vgs":            # s = temperature, vary v_gs
                    it = float(itest[sarr == s][0])
                    return model_ron(p25, float(s), float(xv), it, tcd) * 1e3
                if xkind == "ids":            # s = v_gs, vary test current
                    T = float(rds.tj[sarr == s][0])
                    return model_ron(p25, T, float(s), float(xv), tcd) * 1e3
                it = float(itest[sarr == s][0])   # s = v_gs, vary T_j
                return model_ron(p25, float(xv), float(s), it, tcd) * 1e3

            for s in sorted(set(sarr.tolist())):
                m = sarr == s
                ax1.plot(xarr[m], y[m], ".",
                         label=f"meas {sfmt.format(v=s)}")
                xg = np.linspace(float(xarr[m].min()),
                                 float(xarr[m].max()), 40)
                ys = []
                for xv in xg:
                    try:
                        ys.append(model1(s, xv))
                    except ValueError:
                        ys.append(np.nan)
                ax1.plot(xg, ys, "-", label=f"model {sfmt.format(v=s)}")
            ax1.set_xlabel(xlab)
            ax1.set_ylabel("Rds(on) [mOhm]")
            ax1.legend(fontsize=6)
        if vth is not None:
            ax2 = self.fig.add_subplot(nup, 1, nup)
            xkind = vth.x_kind
            if xkind == "ids":            # Vth CSV leads with I_d
                xarr, xlab = vth.i_test, "I_test [A]"
                sarr, sfmt = vth.tj, "T_j {v:g} C"
            else:                         # T_j (the datasheet layout)
                xarr, xlab = vth.tj, "T_j [C]"
                sarr, sfmt = vth.i_test, "I_test {v:g} A"
            for s in sorted(set(sarr.tolist())):
                m = sarr == s
                ax2.plot(xarr[m], vth.vth[m], ".",
                         label=f"meas {sfmt.format(v=s)}")
                xg = np.linspace(float(xarr[m].min()),
                                 float(xarr[m].max()), 40)
                ys = []
                for xv in xg:
                    try:
                        if xkind == "ids":
                            ys.append(model_vth(p25, float(s), float(xv), tcd))
                        else:
                            ys.append(model_vth(p25, float(xv), float(s), tcd))
                    except ValueError:
                        ys.append(np.nan)
                ax2.plot(xg, ys, "-", label=f"model {sfmt.format(v=s)}")
            ax2.set_xlabel(xlab)
            ax2.set_ylabel("Vth [V]")
            ax2.legend(fontsize=7)
        self.fig.tight_layout()
        self._apply_axes()
        self.canvas.draw()

    def _extract_recov(self):
        rectbl = next((t for t in self.tables
                       if isinstance(t, RecTable)), None)
        if rectbl is None:
            messagebox.showerror(
                "No REC table", "Recovery needs a CSV with Q_rr, t_rr, "
                "I_RM, I_sd, di/dt columns.")
            return
        dio = self._diode_params()

        def work():
            rec, rms, rtext = recov_from_table(
                rectbl, self.part_name.get(), dio,
                progress=lambda m: self.q.put(("progress", "recovery "
                                               + m)))
            return ("recov", (rec, rms, rtext))
        self._run(work)

    def _save(self):
        if self.caps_model is None:
            messagebox.showerror("No caps", "Parse a CAP CSV first.")
            return
        p25 = self._channel_params()
        if "v_ovmax" not in self.free_ch:
            p25 = p25.replace(v_ovmax=0.0)  # user rule: cap off unless fitted
        dio = self._diode_params()
        # re-splice c_gdmax with the fitted on-resistance (release order)
        caps = [t for t in self.tables if isinstance(t, CapTable)]
        qgs = [t for t in self.tables if isinstance(t, QgTable)]
        model, _ = caps_model_from_tables(
            caps, qgs, r_ds_on=p25.r_d + p25.r_s)
        warnings = list(self.report)
        if not self.ch_fitted:
            warnings.append("channel parameters are preset/manual values, "
                            "not fitted")
        if not self.dio_fitted:
            warnings.append("diode card is preset/manual values, not fitted")
        # user rule 2026-11-09: pick the destination folder via a dialog;
        # save ONLY the device subckts (FAST always, FULL when recovery
        # knobs exist) -- the caps subckt and cds expression are embedded
        # inline in them, no standalone artifacts.
        part = self.part_name.get()
        folder = filedialog.askdirectory(
            parent=self.root,
            title=f"Choose folder for the {part} subckt files")
        if not folder:
            self.status.set("save cancelled")
            return
        folder = Path(folder)
        # overwrite guard (session 30 field report: a GUI Save clobbered
        # the official ISC040 artifacts when the part name collided)
        targets = [folder / f"{part}_fast.subckt"]
        if self.recov is not None:
            targets.append(folder / f"{part}_full.subckt")
        clash = [t for t in targets if t.exists()]
        if clash and not messagebox.askyesno(
                "Overwrite?",
                "These existing files will be OVERWRITTEN:\n"
                + "\n".join(t.name for t in clash)
                + "\n\nContinue? (If this is an official generated part, "
                "regenerate it with scripts/gen_full_subckt.py instead.)"):
            return
        res, files = export_part(
            folder, part, p25, dio, self.recov, model,
            n_sub=int(self.n_sub.get()), rmse=self.rmse,
            warnings=warnings, preset=self.preset_name.get(),
            ref_temp=float(self.ref_temp.get()), tc=dict(self.tc),
            caps_files=False, both_variants=True)
        self._log(res.warnings)
        self.status.set("wrote " + ", ".join(Path(f).name for f in files))
        messagebox.showinfo("Saved", "\n".join(files))

    # ---- worker plumbing -------------------------------------------------------
    def _run(self, work):
        if self.busy:
            messagebox.showinfo("Busy", "A fit is already running.")
            return
        self.busy = True
        for b in (self.btn_ch, self.btn_tc, self.btn_dio, self.btn_recov):
            b.state(["disabled"])
        self.status.set("running (worker thread)...")

        def target():
            try:
                self.q.put(("ok", work()))
            except Exception:
                self.q.put(("err", traceback.format_exc()))
        threading.Thread(target=target, daemon=True).start()

    def _poll(self):
        try:
            while True:
                kind, payload = self.q.get_nowait()
                if kind == "progress":
                    self.fit_lbl.set(payload)
                    continue
                self.busy = False
                for b in (self.btn_ch, self.btn_tc, self.btn_dio,
                          self.btn_recov):
                    b.state(["!disabled"])
                if kind == "err":
                    self.status.set("fit failed")
                    messagebox.showerror("Failed", payload)
                    continue
                what, res = payload
                if what == "ch":
                    self.ch_fitted = True
                    for k in self.ch_vals:
                        self.ch_vals[k] = float(getattr(res.params, k))
                    self._update_tree(self.ch_tree, self.ch_vals)
                    self.rmse["channel_ref"] = res.rmse_rel
                    self.fit_lbl.set(f"channel RMSE "
                                     f"{res.rmse_rel * 100:.2f}%")
                    self.status.set("channel ref fit done")
                elif what == "dio":
                    if res is None:
                        self.fit_lbl.set("no diode rows for this preset")
                    else:
                        self.dio_fitted = True
                        for k in self.dio_vals:
                            self.dio_vals[k] = float(
                                getattr(res.params, k))
                        self._update_tree(self.dio_tree, self.dio_vals)
                        self.rmse["diode"] = res.rmse_rel
                        self.fit_lbl.set(f"diode RMSE "
                                         f"{res.rmse_rel * 100:.2f}%")
                elif what == "tc":
                    p = PRESETS[self.preset_name.get()]
                    self.tc = {k: res[k] for k in p.tc_free}
                    self.rmse["tc"] = res["rmse_overall"]
                    self.fit_lbl.set(f"TC RMSE "
                                     f"{res['rmse_overall'] * 100:.2f}%")
                elif what == "recov":
                    self.recov, rms, _ = res
                    # extract() returns rms_rel ALREADY in percent
                    # (session-30 display bug: it was multiplied by 100
                    # again -> "103.55 %" for a 1.04 % fit)
                    if rms is not None:
                        self.fit_lbl.set(f"recovery RMS {rms:.3f} %")
                        self._log([
                            f"recovery fit: tau={displayNum(self.recov.tau)},"
                            f" t_minor={displayNum(self.recov.t_minor)},"
                            f" t_pole={displayNum(self.recov.t_pole)},"
                            f" RMS {rms:.3f} %"])
                    else:
                        self.fit_lbl.set("recovery extracted")
                        self._log([
                            "recovery fit: tau="
                            f"{displayNum(self.recov.tau)}, t_minor="
                            f"{displayNum(self.recov.t_minor)}, t_pole="
                            f"{displayNum(self.recov.t_pole)}"])
        except queue.Empty:
            pass
        self.root.after(120, self._poll)

    def _log(self, lines):
        self.report_txt.configure(state="normal")
        for ln in lines:
            self.report_txt.insert(tk.END, ln + "\n")
        self.report_txt.configure(state="disabled")


def main() -> None:
    root = tk.Tk()
    App(root)
    root.mainloop()


if __name__ == "__main__":
    main()
