"""Adjustable Cartesian grid overlay drawn on top of the datasheet crop.

The grid is generated from the same :class:`~curve_digitizer.transform.AxisTransform`
that produces the tick labels, so editing the axis range, the tick interval or
the linear/log choice re-derives every grid line in data coordinates.  It never
touches the underlying image; it is purely an alignment aid for checking the
calibration against the datasheet's own printed grid.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from .transform import AxisTransform

LINE_STYLES = {
    "solid": "-",
    "dashed": "--",
    "dotted": ":",
    "dash-dot": "-.",
}

COLORS = ("#3b7dd8", "#888888", "#d83b3b", "#2ca02c", "#f0a020", "#ffffff")

# Mantissa sets offered for log minor grid lines.
LOG_MINOR_PRESETS = {
    "none": (),
    "2,5": (2, 5),
    "2..9": (2, 3, 4, 5, 6, 7, 8, 9),
}


@dataclass
class GridState:
    """User-visible grid configuration."""

    major: bool = True
    minor: bool = False
    style: str = "dashed"
    color: str = "#3b7dd8"
    alpha: float = 0.55
    width: float = 0.7
    minor_subdivisions: int = 5
    log_minor: str = "2..9"
    x_interval: float | None = field(default=None)
    y_interval: float | None = field(default=None)


def apply_axes(
    ax,
    xform_x: AxisTransform,
    xform_y: AxisTransform,
    grid: GridState,
) -> tuple[list[float], list[float]]:
    """Install data-valued ticks, labels and the grid overlay on ``ax``.

    Returns the major tick pixel positions (x, y) so callers can reuse them.
    """
    majors_x = xform_x.major_ticks(grid.x_interval)
    majors_y = xform_y.major_ticks(grid.y_interval)
    log_mantissas = LOG_MINOR_PRESETS.get(grid.log_minor, ())
    minors_x = xform_x.minor_ticks(majors_x, grid.minor_subdivisions, log_mantissas)
    minors_y = xform_y.minor_ticks(majors_y, grid.minor_subdivisions, log_mantissas)

    pix_x = [xform_x.to_pix(v) for v in majors_x]
    pix_y = [xform_y.to_pix(v) for v in majors_y]
    ax.set_xticks(pix_x)
    ax.set_xticklabels([_label(xform_x, v) for v in majors_x])
    ax.set_yticks(pix_y)
    ax.set_yticklabels([_label(xform_y, v) for v in majors_y])

    if grid.minor:
        ax.set_xticks([xform_x.to_pix(v) for v in minors_x], minor=True)
        ax.set_yticks([xform_y.to_pix(v) for v in minors_y], minor=True)
    else:
        ax.set_xticks([], minor=True)
        ax.set_yticks([], minor=True)

    style = LINE_STYLES.get(grid.style, "--")
    # matplotlib warns if line properties are passed for a hidden grid, so the
    # properties are only supplied when that grid layer is actually visible.
    if grid.major:
        ax.grid(
            which="major",
            visible=True,
            color=grid.color,
            linestyle=style,
            linewidth=grid.width,
            alpha=grid.alpha,
        )
    else:
        ax.grid(which="major", visible=False)
    if grid.minor:
        ax.grid(
            which="minor",
            visible=True,
            color=grid.color,
            linestyle=style,
            linewidth=max(grid.width * 0.6, 0.2),
            alpha=grid.alpha * 0.5,
        )
    else:
        ax.grid(which="minor", visible=False)
    ax.tick_params(axis="both", which="major", labelsize=8)
    ax.tick_params(axis="both", which="minor", length=2, labelsize=7)
    return pix_x, pix_y


def _label(xform: AxisTransform, value: float) -> str:
    from .transform import format_value

    return format_value(value)
