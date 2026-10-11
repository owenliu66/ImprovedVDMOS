"""Pop-up digitizer window: datasheet crop + Cartesian overlay + marker table."""

from __future__ import annotations

import tkinter as tk
from tkinter import filedialog, messagebox, ttk

from matplotlib.backend_bases import MouseEvent
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

from .calibration import CalibrationPicker
from .corners import CornerFinder
from .export import default_data_dir, export_csv
from .grid import COLORS, LINE_STYLES, LOG_MINOR_PRESETS, GridState, apply_axes
from .markers import MarkerLayer
from .snapping import CurveSnapper
from .table import MarkerTable
from .transform import SCALES, AxisTransform, format_value


class DigitizerWindow(tk.Toplevel):
    """One pop-up per selected datasheet region.

    Internal coordinates are crop pixels (u right, v down); datasheet units come
    from :class:`AxisTransform` and only appear in tick labels, the table and
    the CSV export.
    """

    def __init__(self, master, image, page_index: int, region_id: int) -> None:
        super().__init__(master)
        self.image = image
        self.width, self.height = image.size
        self.page_index = page_index
        self.region_id = region_id
        self.title(f"Curve digitizer - page {page_index + 1}, region {region_id}")

        self.snapper = CurveSnapper(image, threshold=0.55, radius=12.0)
        self.corners = CornerFinder(image, threshold=0.55)
        self.grid_state = GridState()
        # Pixel anchors of the two typed values.  The default is the crop edge
        # pair (the original edge-based calibration); two-point calibration
        # replaces them with the clicked tick positions.  _apply_axes() always
        # rebuilds the transforms from these, so "Apply grid" can never
        # silently discard a calibration.
        self.x_anchors: tuple[float, float] = (0.0, float(self.width))
        self.y_anchors: tuple[float, float] = (float(self.height), 0.0)
        self.xform_x = AxisTransform(*self.x_anchors, 0, 10)
        self.xform_y = AxisTransform(*self.y_anchors, 0, 10)

        # Pack order fixes space allocation: chrome first, plot last (expand).
        self._build_toolbar()
        self._build_axis_panel()
        self._build_grid_panel()
        self._build_hint()
        self._build_table_parent()
        self._build_plot()

        self.markers = MarkerLayer(self.ax, self.canvas, self.snapper, host=self)
        self.calib = CalibrationPicker(self.ax, self.canvas, self.corners, host=self)
        self.table = MarkerTable(self.table_parent, host=self)
        self._apply_axes()
        self._fit_to_screen()
        self.protocol("WM_DELETE_WINDOW", self.destroy)

    def _fit_to_screen(self) -> None:
        """Open at a usable size instead of the (possibly huge) requested size."""
        self.update_idletasks()
        screen_w = max(self.winfo_screenwidth(), 800)
        screen_h = max(self.winfo_screenheight(), 600)
        width = min(self.winfo_reqwidth(), int(0.94 * screen_w))
        height = min(self.winfo_reqheight(), int(0.94 * screen_h))
        self.geometry(f"{width}x{height}")
        self.minsize(640, 420)

    # ----------------------------------------------------------------- UI
    def _build_toolbar(self) -> None:
        bar = ttk.Frame(self, padding=(6, 4))
        bar.pack(side="top", fill="x")

        ttk.Button(bar, text="Export CSV", command=self.export).pack(side="left")

        self.snap_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(bar, text="Snap to curve", variable=self.snap_var).pack(
            side="left", padx=(14, 4)
        )
        ttk.Label(bar, text="radius px").pack(side="left")
        self.radius_var = tk.DoubleVar(value=12.0)
        ttk.Spinbox(bar, from_=2, to=60, increment=1, textvariable=self.radius_var, width=5).pack(
            side="left", padx=(2, 10)
        )

        ttk.Label(bar, text="dark threshold").pack(side="left")
        self.threshold_var = tk.DoubleVar(value=0.55)
        ttk.Scale(
            bar,
            from_=0.05,
            to=0.95,
            variable=self.threshold_var,
            command=lambda _v: self._schedule_rebuild(),
        ).pack(side="left", padx=(2, 10))

        self.mask_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(
            bar, text="show curve mask", variable=self.mask_var, command=self._toggle_mask
        ).pack(side="left", padx=(0, 10))

        self.suppress_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            bar,
            text="ignore printed grid/axes",
            variable=self.suppress_var,
            command=self._schedule_rebuild,
        ).pack(side="left", padx=(0, 10))

        ttk.Button(bar, text="Clear markers", command=self.clear_markers).pack(side="left")
        self.cursor_label = ttk.Label(bar, text="", width=34)
        self.cursor_label.pack(side="right")

    def _build_axis_panel(self) -> None:
        panel = ttk.LabelFrame(self, text="Axes (two known points)", padding=6)
        panel.pack(side="top", fill="x", padx=6, pady=2)

        row = ttk.Frame(panel)
        row.pack(side="top", fill="x")
        self.x_min = self._entry(row, "X at pt 1", "0")
        self.x_max = self._entry(row, "X at pt 2", "10")
        self.x_tick = self._entry(row, "X tick", "1")
        self.x_scale = self._scale(row, "X scale", "linear")

        self.y_min = self._entry(row, "Y at pt 1", "0")
        self.y_max = self._entry(row, "Y at pt 2", "10")
        self.y_tick = self._entry(row, "Y tick", "1")
        self.y_scale = self._scale(row, "Y scale", "linear")

        row2 = ttk.Frame(panel)
        row2.pack(side="top", fill="x", pady=(4, 0))
        self.calib_x_btn = ttk.Button(row2, text="Calibrate X", command=lambda: self._arm("x"))
        self.calib_x_btn.pack(side="left", padx=(0, 4))
        self.calib_y_btn = ttk.Button(row2, text="Calibrate Y", command=lambda: self._arm("y"))
        self.calib_y_btn.pack(side="left", padx=(0, 12))

        self.calib_snap_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(
            row2, text="snap clicks to tick corners", variable=self.calib_snap_var
        ).pack(side="left")

        ttk.Button(row2, text="Apply axes", command=self._apply_axes_clicked).pack(
            side="right", padx=(12, 0)
        )

    def _build_grid_panel(self) -> None:
        panel = ttk.LabelFrame(self, text="Grid overlay", padding=6)
        panel.pack(side="top", fill="x", padx=6, pady=2)

        self.grid_major = tk.BooleanVar(value=True)
        self.grid_minor = tk.BooleanVar(value=False)
        ttk.Checkbutton(panel, text="major", variable=self.grid_major).pack(side="left")
        ttk.Checkbutton(panel, text="minor", variable=self.grid_minor).pack(side="left", padx=4)

        self.grid_style = tk.StringVar(value="dashed")
        self._combo(panel, "style", self.grid_style, list(LINE_STYLES), 8)
        self.grid_color = tk.StringVar(value="#3b7dd8")
        self._combo(panel, "color", self.grid_color, list(COLORS), 9)

        self.grid_alpha = tk.DoubleVar(value=0.55)
        ttk.Label(panel, text="alpha").pack(side="left", padx=(8, 2))
        ttk.Scale(panel, from_=0.05, to=1.0, variable=self.grid_alpha, length=80).pack(side="left")

        self.grid_sub_x = tk.IntVar(value=5)
        ttk.Label(panel, text="X minor div").pack(side="left", padx=(8, 2))
        ttk.Spinbox(panel, from_=2, to=20, textvariable=self.grid_sub_x, width=4).pack(side="left")

        self.grid_sub_y = tk.IntVar(value=5)
        ttk.Label(panel, text="Y minor div").pack(side="left", padx=(8, 2))
        ttk.Spinbox(panel, from_=2, to=20, textvariable=self.grid_sub_y, width=4).pack(side="left")

        self.grid_log_minor_x = tk.StringVar(value="2..9")
        self._combo(panel, "X log minor", self.grid_log_minor_x, list(LOG_MINOR_PRESETS), 6)

        self.grid_log_minor_y = tk.StringVar(value="2..9")
        self._combo(panel, "Y log minor", self.grid_log_minor_y, list(LOG_MINOR_PRESETS), 6)

        ttk.Button(panel, text="Apply grid", command=self._apply_axes).pack(
            side="left", padx=(12, 4)
        )

    def _build_hint(self) -> None:
        self.hint_default = (
            "click plot: add marker    drag: move (snaps to curve)    hold X: horizontal only    "
            "hold Y: vertical only    select + Del: remove    double-click cell: edit"
        )
        self.hint_label = ttk.Label(self, text=self.hint_default, padding=(6, 2))
        self.hint_label.pack(side="bottom", fill="x")

    def _build_table_parent(self) -> None:
        self.table_parent = ttk.LabelFrame(self, text="Markers", padding=4)
        self.table_parent.pack(side="bottom", fill="x", padx=6, pady=2)

    def _build_plot(self) -> None:
        # The canvas requests figsize*dpi, so an uncapped figure would make the
        # window taller than the screen and push the table off-screen.
        max_px_w, max_px_h = 900, 470
        scale = min(max_px_w / self.width, max_px_h / self.height, 1.0)
        fig_w = max(self.width * scale, 320.0) / 100.0
        fig_h = max(self.height * scale, 240.0) / 100.0

        self.figure = Figure(figsize=(fig_w, fig_h), dpi=100)
        self.ax = self.figure.add_axes([0.09, 0.10, 0.90, 0.88])
        self.ax.imshow(
            self.image,
            extent=(0, self.width, self.height, 0),
            interpolation="auto",
            zorder=1,
        )
        self.ax.set_xlim(0, self.width)
        self.ax.set_ylim(self.height, 0)
        self.mask_artist = self.ax.imshow(
            self.snapper.mask,
            cmap="gray",
            extent=(0, self.width, self.height, 0),
            alpha=0.45,
            zorder=2,
        )
        self.mask_artist.set_visible(False)

        self.canvas = FigureCanvasTkAgg(self.figure, master=self)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.mpl_connect("motion_notify_event", self._on_cursor)

    def _entry(self, parent, label, initial):
        ttk.Label(parent, text=label).pack(side="left", padx=(6, 2))
        var = tk.StringVar(value=initial)
        entry = ttk.Entry(parent, textvariable=var, width=8)
        entry.pack(side="left")
        entry.bind("<Return>", lambda _e: self._apply_axes())
        return var

    def _scale(self, parent, label, initial):
        return self._combo(parent, label, tk.StringVar(value=initial), list(SCALES), 7)

    def _combo(self, parent, label, var, values, width):
        ttk.Label(parent, text=label).pack(side="left", padx=(8, 2))
        box = ttk.Combobox(parent, textvariable=var, values=values, width=width, state="readonly")
        box.pack(side="left")
        return var

    # ------------------------------------------------------- host interface
    def snapping_enabled(self) -> bool:
        return bool(self.snap_var.get())

    def snap_radius(self) -> float:
        try:
            return float(self.radius_var.get())
        except (TypeError, tk.TclError, ValueError):
            return 12.0

    def corner_snap_enabled(self) -> bool:
        return bool(self.calib_snap_var.get())

    def on_marker_added(self, _marker) -> None:
        self.table.refresh()

    def on_marker_moved(self, marker) -> None:
        self.table.update_row(marker.index)

    def on_marker_removed(self) -> None:
        self.table.refresh()

    def on_selection_changed(self) -> None:
        self.table.select_row(self.markers.selected)

    def marker_count(self) -> int:
        return len(self.markers.markers)

    def row_text(self, index: int) -> tuple[str, str, str, str]:
        marker = self.markers.markers[index]
        return (
            str(marker.index + 1),
            format_value(self.xform_x.to_data(marker.u)),
            format_value(self.xform_y.to_data(marker.v)),
            marker.condition,
        )

    def on_cell_edited(self, index: int, column: str, text: str) -> None:
        if not 0 <= index < len(self.markers.markers):
            return
        marker = self.markers.markers[index]
        text = text.strip()
        if column == "condition":
            marker.condition = text
            self.table.update_row(index)
            return
        try:
            value = float(text.replace(",", "."))
        except ValueError:
            messagebox.showerror("Invalid value", f"'{text}' is not a number", parent=self)
            self.table.update_row(index)
            return
        transform = self.xform_x if column == "x" else self.xform_y
        pixel = transform.to_pix(value)
        if pixel != pixel:  # NaN: e.g. a non-positive value on a log axis
            messagebox.showerror(
                "Out of range",
                f"{value} cannot be placed on a log {column.upper()} axis",
                parent=self,
            )
            self.table.update_row(index)
            return
        # Typed coordinates are authoritative: no snapping on table edits.
        u, v = (pixel, marker.v) if column == "x" else (marker.u, pixel)
        self.markers.set_pixel(index, u, v)

    def on_row_selected(self, index: int | None) -> None:
        self.markers.select(index)

    def on_delete_requested(self, index: int) -> None:
        self.markers.remove(index)

    # --------------------------------------------------- calibration interface
    def _arm(self, axis: str) -> None:
        self.calib.arm(axis)

    def on_calibration_armed(self, axis: str) -> None:
        self.markers.set_input_enabled(False)
        self.canvas.get_tk_widget().focus_set()  # mpl key events need canvas focus
        self.hint_label.config(
            text=f"calibrating {axis.upper()}: click the tick that carries "
            f"{'X at pt 1' if axis == 'x' else 'Y at pt 1'}, then the second one "
            "(Esc cancels; only the "
            f"{axis.upper() if axis == 'x' else 'Y'} coordinate of each click is used)"
        )

    def on_calibration_cancelled(self) -> None:
        self.markers.set_input_enabled(True)
        self.hint_label.config(text=self.hint_default)
        self.cursor_label.config(text="calibration cancelled")

    def on_calibration_status(self, text: str) -> None:
        self.cursor_label.config(text=text)

    def on_calibration_complete(self, axis: str, point1, point2) -> None:
        if axis == "x":
            self.x_anchors = (float(point1[0]), float(point2[0]))
        else:
            self.y_anchors = (float(point1[1]), float(point2[1]))
        self.markers.set_input_enabled(True)
        self.hint_label.config(text=self.hint_default)
        self._apply_axes()
        self.cursor_label.config(
            text=f"{axis.upper()} calibrated at pixels "
            f"({self.x_anchors[0] if axis == 'x' else self.y_anchors[0]:.1f}) -> "
            f"({self.x_anchors[1] if axis == 'x' else self.y_anchors[1]:.1f})"
        )

    # ------------------------------------------------------------- actions
    def _apply_axes_clicked(self) -> None:
        """Apply-axes button: the green calibration marks have done their job."""
        self.calib.clear_anchors()
        self._apply_axes()

    def _apply_axes(self) -> None:
        try:
            x_min, x_max = float(self.x_min.get()), float(self.x_max.get())
            y_min, y_max = float(self.y_min.get()), float(self.y_max.get())
        except ValueError:
            messagebox.showerror("Invalid axis range", "Min/max must be numbers", parent=self)
            return
        if x_min == x_max or y_min == y_max:
            messagebox.showerror("Invalid axis range", "Min and max must differ", parent=self)
            return
        try:
            self.xform_x = AxisTransform.from_points(
                *self.x_anchors, x_min, x_max, self.x_scale.get()
            )
            self.xform_y = AxisTransform.from_points(
                *self.y_anchors, y_min, y_max, self.y_scale.get()
            )
        except ValueError as exc:
            messagebox.showerror("Invalid axis", str(exc), parent=self)
            return

        self.grid_state.x_interval = self._interval(self.x_tick)
        self.grid_state.y_interval = self._interval(self.y_tick)
        self.grid_state.major = bool(self.grid_major.get())
        self.grid_state.minor = bool(self.grid_minor.get())
        self.grid_state.style = self.grid_style.get()
        self.grid_state.color = self.grid_color.get()
        self.grid_state.alpha = float(self.grid_alpha.get())
        self.grid_state.x_minor_divisions = int(self.grid_sub_x.get())
        self.grid_state.y_minor_divisions = int(self.grid_sub_y.get())
        self.grid_state.x_log_minor = self.grid_log_minor_x.get()
        self.grid_state.y_log_minor = self.grid_log_minor_y.get()

        apply_axes(self.ax, self.xform_x, self.xform_y, self.grid_state)
        self.canvas.draw_idle()
        self.table.refresh()

    @staticmethod
    def _interval(var) -> float | None:
        try:
            value = float(var.get())
        except (ValueError, tk.TclError):
            return None
        return value if value > 0 else None

    def _schedule_rebuild(self) -> None:
        if getattr(self, "_threshold_job", None):
            self.after_cancel(self._threshold_job)
        self._threshold_job = self.after(250, self._rebuild_snapper)

    def _rebuild_snapper(self) -> None:
        self._threshold_job = None
        threshold = float(self.threshold_var.get())
        self.snapper.rebuild(
            threshold=threshold,
            suppress_printed_grid=bool(self.suppress_var.get()),
        )
        # Corner detection needs its own mask (frames and ticks are kept), but
        # it shares the dark-threshold slider so both agree on "what is ink".
        self.corners.rebuild(threshold=threshold)
        self.mask_artist.set_data(self.snapper.mask)
        self._toggle_mask()

    def _toggle_mask(self) -> None:
        self.mask_artist.set_visible(bool(self.mask_var.get()))
        self.canvas.draw_idle()

    def _on_cursor(self, event: MouseEvent) -> None:
        if event.inaxes is not self.ax or event.xdata is None or event.ydata is None:
            self.cursor_label.config(text="")
            return
        if self.calib.armed():
            # While calibrating, data values would be misleading: the axis is
            # exactly what is not yet defined.  Show raw crop pixels instead.
            self.cursor_label.config(
                text=f"calibrating {self.calib.axis.upper()}:   "
                f"u = {event.xdata:.1f} px   v = {event.ydata:.1f} px"
            )
            return
        self.cursor_label.config(
            text=f"x = {format_value(self.xform_x.to_data(event.xdata))}   "
            f"y = {format_value(self.xform_y.to_data(event.ydata))}"
        )

    def clear_markers(self) -> None:
        if not self.markers.markers:
            return
        if messagebox.askyesno("Clear markers", "Delete all markers?", parent=self):
            self.markers.markers.clear()
            self.markers.selected = None
            self.markers._rebuild_artists()
            self.table.refresh()
            self.canvas.draw_idle()

    def export(self) -> None:
        if not self.markers.markers:
            messagebox.showinfo("Export CSV", "No markers to export.", parent=self)
            return
        path = filedialog.asksaveasfilename(
            parent=self,
            title="Save digitized curve as CSV",
            initialdir=str(default_data_dir()),
            initialfile=f"page{self.page_index + 1}_region{self.region_id}.csv",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
        )
        if not path:
            return
        try:
            written = export_csv(path, self.markers.markers, self.xform_x, self.xform_y)
        except OSError as exc:
            messagebox.showerror("Export failed", str(exc), parent=self)
            return
        self.cursor_label.config(text=f"saved {written.name}")
