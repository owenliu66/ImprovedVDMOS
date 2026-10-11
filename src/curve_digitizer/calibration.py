"""Click-two-points axis calibration on the digitizer plot.

The datasheet crop is calibrated by clicking the two printed ticks whose data
values are already typed in the axis panel.  While armed, this picker owns the
mouse: it draws crosshair and connector artists, optionally snaps each click
onto the tick/axis corner found by :class:`~curve_digitizer.corners.CornerFinder`,
validates the pair, and hands the two pixel anchors to the host window.

The two anchor pixels are directional: point 1 maps to the first typed value
and point 2 to the second, so clicking right-to-left or bottom-to-top simply
inverts the axis the same way the transform already does.
"""

from __future__ import annotations

import math

from matplotlib.axes import Axes
from matplotlib.backend_bases import MouseEvent

MIN_SEPARATION_PX = 5.0  # two clicks closer than this cannot define an axis
COLOR = "#2ca02c"


class CalibrationPicker:
    """State machine: idle -> point 1 -> point 2 -> idle.

    The ``host`` object must expose ``corner_snap_enabled()``,
    ``on_calibration_armed(axis)``, ``on_calibration_cancelled()``,
    ``on_calibration_status(text)`` and
    ``on_calibration_complete(axis, point1, point2)``.
    """

    def __init__(self, ax: Axes, canvas, corners, host) -> None:
        self.ax = ax
        self.canvas = canvas
        self.corners = corners
        self.host = host
        self.axis: str | None = None
        self.pending: list[tuple[float, float]] = []
        self._pending_artists: list = []
        self._anchor_artists: dict[str, list] = {"x": [], "y": []}

        canvas.mpl_connect("button_press_event", self._on_press)
        canvas.mpl_connect("key_press_event", self._on_key_press)

    # ------------------------------------------------------------------ state
    def armed(self) -> bool:
        return self.axis is not None

    def arm(self, axis: str) -> None:
        """Start calibrating ``axis``; pressing the same button again cancels."""
        if self.axis == axis:
            self.cancel()
            return
        self._clear(self._pending_artists)
        self.axis = axis
        self.pending = []
        self.host.on_calibration_armed(axis)

    def cancel(self) -> None:
        if self.axis is None:
            return
        self.axis = None
        self.pending = []
        self._clear(self._pending_artists)
        self.host.on_calibration_cancelled()
        self.canvas.draw_idle()

    def destroy(self) -> None:
        for artists in self._anchor_artists.values():
            self._clear(artists)
        self._clear(self._pending_artists)

    def clear_anchors(self) -> None:
        """Remove the persistent green anchor markers from both axes."""
        for artists in self._anchor_artists.values():
            self._clear(artists)
        self.canvas.draw_idle()

    # ---------------------------------------------------------------- artists
    def _clear(self, artists: list) -> None:
        for artist in artists:
            artist.remove()
        artists.clear()

    def _draw_pending(self) -> None:
        self._clear(self._pending_artists)
        if len(self.pending) == 2:
            (u1, v1), (u2, v2) = self.pending
            self._pending_artists.append(
                self.ax.plot([u1, u2], [v1, v2], ls="--", lw=1.0, color=COLOR, zorder=7)[0]
            )
        for k, (u, v) in enumerate(self.pending):
            self._pending_artists.append(
                self.ax.plot([u], [v], ls="", marker="+", ms=16, mew=2.0, color=COLOR, zorder=8)[0]
            )
            self._pending_artists.append(
                self.ax.annotate(
                    f"{k + 1}",
                    xy=(u, v),
                    xytext=(9, 9),
                    textcoords="offset points",
                    fontsize=8,
                    color=COLOR,
                    zorder=8,
                )
            )

    def _promote(self, axis: str) -> None:
        """Keep the accepted pair visible until the next calibration."""
        self._clear(self._anchor_artists[axis])
        self._anchor_artists[axis], self._pending_artists = self._pending_artists, []

    # ----------------------------------------------------------------- events
    def _on_press(self, event: MouseEvent) -> None:
        if self.axis is None or event.inaxes is not self.ax or event.button != 1:
            return
        if event.xdata is None or event.ydata is None:
            return
        u, v = float(event.xdata), float(event.ydata)
        note = ""
        if self.host.corner_snap_enabled():
            hit = self.corners.find(u, v, self.axis)
            if hit is not None:
                u, v = hit.u, hit.v
                note = f"   snapped to {hit.kind} ({hit.offset:.1f} px)"
            else:
                note = "   no tick corner found - raw click used"

        if self.pending:
            (u1, v1) = self.pending[0]
            if math.hypot(u - u1, v - v1) < MIN_SEPARATION_PX:
                self.host.on_calibration_status(
                    f"points too close ({math.hypot(u - u1, v - v1):.1f} px) - click point 2 again"
                )
                self.pending = [(u1, v1)]
                self._draw_pending()
                self.canvas.draw_idle()
                return
            self.pending.append((u, v))
            self._draw_pending()
            self._promote(self.axis)
            axis, points = self.axis, tuple(self.pending)
            self.axis = None
            self.pending = []
            self.host.on_calibration_complete(axis, points[0], points[1])
            self.canvas.draw_idle()
            return

        self.pending = [(u, v)]
        self._draw_pending()
        self.host.on_calibration_status(
            f"point 1 at u = {u:.1f}  v = {v:.1f}{note}   -   click point 2 (Esc cancels)"
        )
        self.canvas.draw_idle()

    def _on_key_press(self, event) -> None:
        if self.axis is not None and event.key in ("escape", "Escape"):
            self.cancel()
