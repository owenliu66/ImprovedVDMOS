"""Two-point calibration: corner snapping, anchored transforms, per-axis minors.

The corner tests are pure NumPy/PIL on synthetic "datasheets" (a printed frame
with ticks), so they run without a display.  The window test drives Tk without
``mainloop`` and skips itself if no display is available.
"""

from __future__ import annotations

import numpy as np
import pytest
from matplotlib.figure import Figure

from curve_digitizer.corners import CornerFinder
from curve_digitizer.grid import GridState, apply_axes
from curve_digitizer.transform import LINEAR, LOG, AxisTransform

Image = pytest.importorskip("PIL.Image", reason="Pillow ships with matplotlib")


# --------------------------------------------------------------- synthetic sheet
def _line(arr, x0, y0, x1, y1, thick=2):
    steps = int(3 * max(abs(x1 - x0), abs(y1 - y0))) + 1
    for k in range(steps + 1):
        t = k / steps
        i = int(round(y0 + t * (y1 - y0)))
        j = int(round(x0 + t * (x1 - x0)))
        arr[i : i + thick, j : j + thick] = 0


def _sheet(
    w=400,
    h=300,
    frame_row=None,
    frame_col=None,
    thickness=2,
    ticks=(),
    tick_len=9,
    tick_thick=2,
    decoys=(),
):
    """White sheet with a printed axis frame, ticks crossing it, and decoys."""
    arr = np.full((h, w), 255, dtype=np.uint8)
    if frame_row is not None:
        arr[frame_row : frame_row + thickness, :] = 0
    if frame_col is not None:
        arr[:, frame_col : frame_col + thickness] = 0
    if frame_row is not None:
        for x in ticks:
            arr[frame_row - tick_len : frame_row + thickness + tick_len, x : x + tick_thick] = 0
    if frame_col is not None:
        for y in ticks:
            arr[y : y + tick_thick, frame_col - tick_len : frame_col + thickness + tick_len] = 0
    for decoy in decoys:
        _line(arr, *decoy)
    return Image.fromarray(arr, mode="L")


# ---------------------------------------------------------------- corner finding
@pytest.mark.parametrize(
    ("thickness", "tick_thick", "expect"),
    [(1, 1, (120.5, 250.5)), (2, 2, (121.0, 251.0)), (3, 3, (121.5, 251.5))],
)
def test_corner_vertex_is_line_centre(thickness, tick_thick, expect):
    """Symmetric printed bands give an exact centreline vertex, any thickness.

    Display convention: pixel ``j`` spans ``[j, j+1]``, so a band covering
    columns 120-121 centres on 121.0 and a 1-pixel band on 120.5.
    """
    img = _sheet(frame_row=250, thickness=thickness, ticks=(120,), tick_thick=tick_thick)
    hit = CornerFinder(img).find(118.0, 247.0, "x")
    assert hit is not None and hit.kind == "corner"
    assert hit.u == pytest.approx(expect[0], abs=0.05)
    assert hit.v == pytest.approx(expect[1], abs=0.05)


def test_corner_snap_moves_click_toward_the_tick():
    img = _sheet(frame_row=250, ticks=(120,))
    hit = CornerFinder(img).find(116.0, 244.0, "x")
    assert hit is not None
    assert (hit.u, hit.v) == pytest.approx((121.0, 251.0), abs=0.05)
    assert hit.offset == pytest.approx(np.hypot(121.0 - 116.0, 251.0 - 244.0), abs=0.05)


def test_y_axis_corner_uses_the_vertical_frame():
    img = _sheet(frame_col=60, ticks=(150,))
    hit = CornerFinder(img).find(63.0, 147.0, "y")
    assert hit is not None and hit.kind == "corner"
    assert (hit.u, hit.v) == pytest.approx((61.0, 151.0), abs=0.05)


def test_steep_decoy_is_rejected_by_the_right_angle_test():
    """A 70-degree stroke crosses the frame but is not perpendicular to it."""
    img = _sheet(
        frame_row=250,
        ticks=(),
        decoys=[(196, 232, 204, 254, 2)],  # ~70 degrees from horizontal
    )
    assert CornerFinder(img).find(200.0, 246.0, "x") is None


def test_shallow_curve_crossing_is_rejected():
    img = _sheet(
        frame_row=250,
        ticks=(),
        decoys=[(180, 244, 220, 258, 2)],  # ~19 degrees, a curve leaving saturation
    )
    assert CornerFinder(img).find(200.0, 250.0, "x") is None


def test_lone_tick_fallback_when_the_frame_is_cropped_out():
    """No axis line in the window: snap u to the tick stroke, keep the raw v."""
    img = _sheet()
    arr = np.array(img)
    arr[100:116, 120:122] = 0
    hit = CornerFinder(Image.fromarray(arr, mode="L")).find(123.0, 105.0, "x")
    assert hit is not None and hit.kind == "tick"
    assert hit.u == pytest.approx(121.0, abs=0.05)
    assert hit.v == pytest.approx(105.0, abs=1e-9)


def test_no_ink_near_the_click_falls_back_to_the_raw_point():
    img = _sheet(frame_row=250, ticks=(120,))
    assert CornerFinder(img).find(30.0, 30.0, "x") is None


def test_corner_mask_keeps_the_frame_that_the_snapper_deletes():
    """The two masks must stay separate: frames/ticks are noise for snapping."""
    img = _sheet(frame_row=250, ticks=(120,))
    finder = CornerFinder(img)
    assert finder.mask[250:252, 200].any()  # frame survives
    assert finder.mask[241:250, 120:122].any()  # tick survives


# ------------------------------------------------------- anchored transforms
def test_from_points_round_trip_linear():
    xf = AxisTransform.from_points(121.0, 381.0, 0.0, 40.0)
    assert xf.to_data(121.0) == pytest.approx(0.0)
    assert xf.to_data(381.0) == pytest.approx(40.0)
    assert xf.to_data(251.0) == pytest.approx(20.0)
    for u in (121.0, 140.3, 251.0, 300.7, 381.0):
        assert xf.to_pix(xf.to_data(u)) == pytest.approx(u)


def test_from_points_log_axis_fixes_decade_scale():
    xf = AxisTransform.from_points(50.0, 300.0, 1e-3, 1e1, LOG)
    assert xf.to_data(50.0) == pytest.approx(1e-3)
    assert xf.to_data(300.0) == pytest.approx(1e1)
    assert xf.to_data(175.0) == pytest.approx(0.1)  # midpoint of 4 decades
    for u in (50.0, 112.5, 175.0, 237.5, 300.0):
        assert xf.to_pix(xf.to_data(u)) == pytest.approx(u)


def test_from_points_direction_follows_click_order():
    """Clicking right-to-left or bottom-to-top reverses the axis, as printed."""
    xf = AxisTransform.from_points(380.0, 10.0, 0.0, 40.0)
    assert xf.to_data(380.0) == pytest.approx(0.0)
    assert xf.to_data(10.0) == pytest.approx(40.0)
    yf = AxisTransform.from_points(280.0, 10.0, 0.0, 10.0)
    assert yf.to_data(280.0) == pytest.approx(0.0)
    assert yf.to_data(10.0) == pytest.approx(10.0)


def test_from_points_rejects_two_identical_pixels():
    with pytest.raises(ValueError):
        AxisTransform.from_points(120.0, 120.0, 0.0, 10.0)


def test_edge_calibration_is_the_special_case_of_two_points():
    """The original edge-based workflow is unchanged: anchors at 0 and W."""
    edge = AxisTransform(0, 400, 0, 10, LINEAR)
    two = AxisTransform.from_points(0.0, 400.0, 0.0, 10.0)
    for u in (0.0, 37.5, 200.0, 399.0):
        assert two.to_data(u) == pytest.approx(edge.to_data(u))


# --------------------------------------------------------- per-axis minor ticks
def _figure_axes():
    fig = Figure(figsize=(4, 3))
    return fig, fig.add_axes([0.1, 0.1, 0.8, 0.8])


def test_minor_divisions_are_independent_per_axis():
    xf = AxisTransform(0, 100, 0, 10, LINEAR)
    yf = AxisTransform(100, 0, 0, 10, LINEAR)
    _fig, ax = _figure_axes()
    grid = GridState(minor=True, x_minor_divisions=2, y_minor_divisions=5)
    apply_axes(ax, xf, yf, grid)
    assert len(ax.get_xticks(minor=True)) == 5  # 1 major subdivision per interval
    assert len(ax.get_yticks(minor=True)) == 20  # 4 per interval


def test_log_minor_mantissas_are_independent_per_axis():
    xf = AxisTransform(0, 100, 1e-3, 1e1, LOG)
    yf = AxisTransform(100, 0, 1e-3, 1e1, LOG)
    _fig, ax = _figure_axes()
    grid = GridState(minor=True, x_log_minor="2,5", y_log_minor="2..9")
    apply_axes(ax, xf, yf, grid)
    assert len(ax.get_xticks(minor=True)) == 8  # 4 decades x {2, 5}
    assert len(ax.get_yticks(minor=True)) == 32  # 4 decades x {2..9}


def test_log_minors_survive_a_subdivision_of_one():
    """The guard binds as ``not majors or (subdivisions <= 1 and LINEAR)``."""
    xf = AxisTransform(0, 100, 1e-3, 1e1, LOG)
    majors = xf.major_ticks(1)
    assert len(xf.minor_ticks(majors, 1, (2, 5))) == 8


# ----------------------------------------------------------------- window glue
class _Click:
    """Stand-in for a matplotlib MouseEvent aimed at the plot."""

    def __init__(self, ax, x, y, button=1):
        self.inaxes = ax
        self.button = button
        self.xdata = x
        self.ydata = y


class _Key:
    def __init__(self, key):
        self.key = key


def _make_window():
    import tkinter as tk

    from curve_digitizer.digitizer_window import DigitizerWindow

    try:
        root = tk.Tk()
    except tk.TclError:  # pragma: no cover - headless CI
        pytest.skip("no Tk display available")
    root.withdraw()
    # frame at row 250 with ticks at x = 120 and x = 380, i.e. 0 and 40 A
    img = _sheet(frame_row=250, ticks=(120, 380))
    return root, DigitizerWindow(root, img, 0, 1)


@pytest.fixture()
def window():
    root, win = _make_window()
    try:
        yield win
    finally:
        win.update_idletasks()  # flush mpl's queued Tk idle draw before teardown
        win.destroy()
        root.destroy()


def test_two_point_calibration_replaces_the_edge_anchors(window):
    win = window
    win.x_min.set("0")
    win.x_max.set("40")
    win._arm("x")
    assert win.calib.armed() and win.markers._input_enabled is False

    # A stray click while armed must not create a marker.
    win.markers._on_press(_Click(win.ax, 50.0, 100.0))
    assert win.marker_count() == 0

    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win.calib._on_press(_Click(win.ax, 378.0, 247.0))

    assert win.x_anchors == pytest.approx((121.0, 381.0))
    assert not win.calib.armed() and win.markers._input_enabled is True
    assert win.xform_x.to_data(121.0) == pytest.approx(0.0)
    assert win.xform_x.to_data(381.0) == pytest.approx(40.0)
    assert win.xform_x.to_data(251.0) == pytest.approx(20.0)
    # Y is untouched: still the crop-edge calibration.
    assert win.y_anchors == pytest.approx((300.0, 0.0))


def test_apply_axes_keeps_the_calibration_anchors(window):
    """Re-applying the grid must not silently revert to edge calibration."""
    win = window
    win._arm("x")
    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win.calib._on_press(_Click(win.ax, 378.0, 247.0))
    win._apply_axes()  # e.g. the user pressed "Apply grid"
    assert win.x_anchors == pytest.approx((121.0, 381.0))
    assert win.xform_x.to_data(251.0) == pytest.approx(5.0)  # values are 0..10


def test_clicking_the_same_calibrate_button_cancels(window):
    win = window
    win._arm("x")
    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win._arm("x")
    assert not win.calib.armed()
    assert win.markers._input_enabled is True
    assert win.x_anchors == pytest.approx((0.0, 400.0))  # unchanged


def test_escape_cancels_and_restores_marker_input(window):
    win = window
    win._arm("y")
    win.calib._on_key_press(_Key("escape"))
    assert not win.calib.armed()
    assert win.markers._input_enabled is True


def test_two_clicks_too_close_are_not_accepted(window):
    win = window
    win._arm("x")
    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win.calib._on_press(_Click(win.ax, 120.0, 249.0))
    assert win.calib.armed()  # still waiting for a usable second point
    assert len(win.calib.pending) == 1
    assert win.x_anchors == pytest.approx((0.0, 400.0))


def test_corner_snap_can_be_turned_off(window):
    win = window
    win.calib_snap_var.set(False)
    win._arm("x")
    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win.calib._on_press(_Click(win.ax, 378.0, 247.0))
    assert win.x_anchors == pytest.approx((118.0, 378.0))  # raw clicks


def test_apply_axes_button_clears_the_green_calibration_marks(window):
    win = window
    win._arm("x")
    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win.calib._on_press(_Click(win.ax, 378.0, 247.0))
    # The accepted pair stays visible right after calibration as feedback.
    assert win.calib._anchor_artists["x"]

    win._apply_axes_clicked()  # the "Apply axes" button
    assert not win.calib._anchor_artists["x"]
    assert not win.calib._anchor_artists["y"]
    # Anchors themselves survive: only the artists are removed.
    assert win.x_anchors == pytest.approx((121.0, 381.0))


def test_markers_added_after_calibration_report_calibrated_values(window):
    win = window
    win.snap_var.set(False)
    win.x_min.set("0")
    win.x_max.set("40")
    win._arm("x")
    win.calib._on_press(_Click(win.ax, 118.0, 247.0))
    win.calib._on_press(_Click(win.ax, 378.0, 247.0))

    win.markers._on_press(_Click(win.ax, 251.0, 150.0))
    assert win.marker_count() == 1
    assert win.row_text(0)[1] == "20"  # x = 20 A at the pixel midpoint

