"""Marker model and mouse/keyboard interaction on the digitizer plot.

Plot data coordinates *are* crop-pixel coordinates, so a marker stores (u, v)
in pixels and the host converts to datasheet units for display and export.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from matplotlib.axes import Axes
from matplotlib.backend_bases import KeyEvent, MouseEvent
from matplotlib.lines import Line2D

PALETTE = (
    "#e8332a",
    "#1f77b4",
    "#2ca02c",
    "#ff7f0e",
    "#9467bd",
    "#17becf",
    "#bcbd22",
    "#8c564b",
)

GRAB_TOLERANCE_PX = 10.0


@dataclass
class Marker:
    u: float
    v: float
    condition: str = ""
    index: int = field(default=0)  # 0-based, displayed as index + 1


class MarkerLayer:
    """Owns the marker artists plus add / drag / axis-lock / delete handling.

    The ``host`` object must expose ``snapping_enabled()``, ``snap_radius()``,
    ``on_marker_added()``, ``on_marker_moved()``, ``on_marker_removed()`` and
    ``on_selection_changed()``.
    """

    def __init__(self, ax: Axes, canvas, snapper, host) -> None:
        self.ax = ax
        self.canvas = canvas
        self.snapper = snapper
        self.host = host
        self.markers: list[Marker] = []
        self.selected: int | None = None
        self._lock: str | None = None
        self._drag: int | None = None
        self._lines: list[Line2D] = []
        self._labels = []

        canvas.mpl_connect("button_press_event", self._on_press)
        canvas.mpl_connect("motion_notify_event", self._on_motion)
        canvas.mpl_connect("button_release_event", self._on_release)
        canvas.mpl_connect("key_press_event", self._on_key_press)
        canvas.mpl_connect("key_release_event", self._on_key_release)

    # ------------------------------------------------------------- geometry
    def limits(self) -> tuple[float, float, float, float]:
        x0, x1 = self.ax.get_xlim()
        y0, y1 = self.ax.get_ylim()  # y is inverted: y0 = H (bottom), y1 = 0 (top)
        return min(x0, x1), max(x0, x1), min(y0, y1), max(y0, y1)

    def _clamp(self, u: float, v: float) -> tuple[float, float]:
        x_lo, x_hi, y_lo, y_hi = self.limits()
        return min(max(u, x_lo), x_hi), min(max(v, y_lo), y_hi)

    def _display(self, marker: Marker) -> tuple[float, float]:
        return self.ax.transData.transform((marker.u, marker.v))

    # ------------------------------------------------------------- mutation
    def add(self, u: float, v: float, snap: bool | None = None) -> Marker:
        u, v = self._clamp(u, v)
        if snap is None:
            snap = self.host.snapping_enabled()
        if snap:
            point = self.snapper.snap(u, v, None, self.host.snap_radius())
            if point is not None:
                u, v = point
        marker = Marker(u=u, v=v, index=len(self.markers))
        self.markers.append(marker)
        self.selected = marker.index
        self._rebuild_artists()
        self.host.on_marker_added(marker)
        self.host.on_selection_changed()
        self.canvas.draw_idle()
        return marker

    def remove_selected(self) -> None:
        if self.selected is None:
            return
        del self.markers[self.selected]
        for new_index, marker in enumerate(self.markers):
            marker.index = new_index
        self.selected = None if not self.markers else min(self.selected, len(self.markers) - 1)
        self._rebuild_artists()
        self.host.on_marker_removed()
        self.host.on_selection_changed()
        self.canvas.draw_idle()

    def remove(self, index: int) -> None:
        if 0 <= index < len(self.markers):
            self.selected = index
            self.remove_selected()

    def select(self, index: int | None) -> None:
        if index is None or 0 <= index < len(self.markers):
            self.selected = index
            self._restyle()
            self.canvas.draw_idle()

    def set_pixel(self, index: int, u: float, v: float, redraw: bool = True) -> None:
        if not 0 <= index < len(self.markers):
            return
        marker = self.markers[index]
        marker.u, marker.v = self._clamp(u, v)
        self._update_artist(index)
        self.host.on_marker_moved(marker)
        if redraw:
            self.canvas.draw_idle()

    # ------------------------------------------------------------ rendering
    def _rebuild_artists(self) -> None:
        for line in self._lines:
            line.remove()
        for text in self._labels:
            text.remove()
        self._lines, self._labels = [], []
        for index in range(len(self.markers)):
            line = self.ax.plot([], [], ls="", marker="o", ms=7, zorder=6, picker=True)[0]
            text = self.ax.annotate(
                "",
                xy=(0, 0),
                xytext=(7, -7),
                textcoords="offset points",
                fontsize=7,
                zorder=7,
            )
            self._lines.append(line)
            self._labels.append(text)
            self._update_artist(index)
        self._restyle()

    def _update_artist(self, index: int) -> None:
        marker = self.markers[index]
        line = self._lines[index]
        line.set_data([marker.u], [marker.v])
        text = self._labels[index]
        text.set_text(str(marker.index + 1))
        text.set_position((marker.u, marker.v))

    def _restyle(self) -> None:
        for index, _marker in enumerate(self.markers):
            color = PALETTE[index % len(PALETTE)]
            line = self._lines[index]
            is_selected = index == self.selected
            line.set_color(color)
            line.set_markerfacecolor("#ffffff" if is_selected else color)
            line.set_markeredgecolor(color if is_selected else "#202020")
            line.set_markeredgewidth(2.2 if is_selected else 1.0)
            line.set_markersize(10 if is_selected else 7)
            self._labels[index].set_color(color)

    # ------------------------------------------------------------- events
    def _hit_marker(self, event: MouseEvent) -> int | None:
        best, best_dist = None, GRAB_TOLERANCE_PX
        for marker in self.markers:
            dx, dy = self._display(marker)
            dist = ((dx - event.x) ** 2 + (dy - event.y) ** 2) ** 0.5
            if dist <= best_dist:
                best, best_dist = marker.index, dist
        return best

    def _on_press(self, event: MouseEvent) -> None:
        if event.inaxes is not self.ax or event.button != 1:
            return
        hit = self._hit_marker(event)
        if hit is None:
            marker = self.add(event.xdata, event.ydata)
            self._drag = marker.index
        else:
            self.selected = hit
            self._restyle()
            self.host.on_selection_changed()
            self._drag = hit
            self.canvas.draw_idle()

    def _on_motion(self, event: MouseEvent) -> None:
        if self._drag is None or event.inaxes is not self.ax:
            return
        index = self._drag
        marker = self.markers[index]
        u, v = event.xdata, event.ydata
        if self._lock == "x":
            v = marker.v
        elif self._lock == "y":
            u = marker.u
        if self.host.snapping_enabled():
            point = self.snapper.snap(u, v, self._lock, self.host.snap_radius())
            if point is not None:
                u, v = point
        self.set_pixel(index, u, v)

    def _on_release(self, event: MouseEvent) -> None:
        self._drag = None

    def _on_key_press(self, event: KeyEvent) -> None:
        if event.key in ("x", "X"):
            self._lock = "x"
        elif event.key in ("y", "Y"):
            self._lock = "y"
        elif event.key in ("delete", "backspace"):
            self.remove_selected()

    def _on_key_release(self, event: KeyEvent) -> None:
        if event.key in ("x", "X") and self._lock == "x":
            self._lock = None
        elif event.key in ("y", "Y") and self._lock == "y":
            self._lock = None
