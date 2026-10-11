"""Sub-pixel 90-degree corner detection for two-point axis calibration.

A printed tick meets the axis line at a right angle, so the most accurate
anchor for a calibration click is the *vertex* where the tick stroke crosses
the axis line, not the raw cursor position.  ``CornerFinder`` finds that vertex
inside a small window around the click:

1. locate the axis line: the window row with the longest horizontal dark run;
2. label the dark components *outside* that band and keep the ones touching it
   -- those are the ticks (and any gridline crossing the frame);
3. reject components whose principal direction is not perpendicular to the
   axis line (within ``ANGLE_TOL_DEG`` of 90 degrees).  This is what discards a
   curve leaving saturation at a shallow angle, arrowheads and text;
4. return the (tick centre, axis centre) vertex.

Sub-pixel positions come from *grayscale* threshold crossings, not from the
binary mask: the band edges are located by linearly interpolating the intensity
profile across the threshold, so an anti-aliased printed line is placed between
pixel centres.  Positions are reported in the display convention used by
``imshow(extent=(0, W, H, 0))``, where pixel column ``j`` spans ``[j, j+1]`` and
therefore a band covering columns ``c0..c1`` has its centre at ``(c0+c1+1)/2``.

Calibrating the Y axis is the same search on a transposed mask and gray array,
so one implementation covers both.

This class deliberately builds its *own* dark-pixel mask.  ``CurveSnapper``
removes long straight rows/columns (frames, gridlines) and small components
(isolated ticks) because they are noise for trace snapping, but they are
exactly what corner detection needs, so the two masks are kept separate.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy import ndimage

WINDOW_RADIUS = 14  # half-size of the search window around the click, in px
TICK_MIN_PX = 3  # fewest dark pixels that can be a tick stroke
ANGLE_TOL_DEG = 12.0  # how far from 90 degrees a stroke may tilt and still count
MIN_LINE_FRAC = 0.5  # axis line must span this fraction of the window width
MAX_LINE_THICKNESS = 6  # a wider dark band is a filled block, not a printed line
EDGE_MARGIN = 3  # rows/columns of gray context probed around a band

AXIS_X, AXIS_Y = "x", "y"
_STRUCTURE = np.ones((3, 3), dtype=int)


@dataclass(frozen=True)
class CornerHit:
    """Snapped vertex plus how far it moved the raw click."""

    u: float
    v: float
    kind: str  # "corner" = tick x axis line, "tick" = lone stroke, no line found
    offset: float  # distance moved from the raw click, in px


def _longest_run(flags: np.ndarray) -> tuple[int, int, int]:
    """Start, end and length of the longest True run in a 1-D boolean array."""
    idx = np.nonzero(flags)[0]
    if idx.size == 0:
        return -1, -1, 0
    breaks = np.nonzero(np.diff(idx) > 1)[0]
    starts = np.concatenate(([idx[0]], idx[breaks + 1]))
    ends = np.concatenate((idx[breaks], [idx[-1]]))
    lengths = ends - starts + 1
    k = int(np.argmax(lengths))
    return int(starts[k]), int(ends[k]), int(lengths[k])


class CornerFinder:
    """Finds the tick/axis corner nearest a click, in crop-pixel coordinates."""

    def __init__(self, image, threshold: float = 0.55) -> None:
        self.gray = np.asarray(image.convert("L"), dtype=np.float32) / 255.0
        self.height, self.width = self.gray.shape
        self.threshold = float(threshold)
        self.rebuild()

    def rebuild(self, threshold: float | None = None) -> None:
        """Rebuild the raw dark mask (no line suppression, no size filter)."""
        if threshold is not None:
            self.threshold = float(threshold)
        self.mask = self.gray < self.threshold

    # ------------------------------------------------------------------ find
    def find(
        self,
        u: float,
        v: float,
        axis: str,
        radius: int = WINDOW_RADIUS,
    ) -> CornerHit | None:
        """Corner nearest the click, or ``None`` (caller keeps the raw click).

        ``axis`` is the axis being calibrated: for "x" the axis line is the
        horizontal frame the ticks hang from; for "y" the search runs on a
        transposed mask so the same code finds a vertical line with horizontal
        ticks.
        """
        if axis == AXIS_X:
            return self._find(self.mask, self.gray, float(u), float(v), int(radius))
        hit = self._find(self.mask.T, self.gray.T, float(v), float(u), int(radius))
        if hit is None:
            return None
        return CornerHit(u=hit.v, v=hit.u, kind=hit.kind, offset=hit.offset)

    def _find(
        self,
        mask: np.ndarray,
        gray: np.ndarray,
        u: float,
        v: float,
        radius: int,
    ) -> CornerHit | None:
        h, w = mask.shape
        j0 = max(int(round(u)) - radius, 0)
        j1 = min(int(round(u)) + radius, w - 1)
        i0 = max(int(round(v)) - radius, 0)
        i1 = min(int(round(v)) + radius, h - 1)
        if j1 <= j0 or i1 <= i0:
            return None
        win = mask[i0 : i1 + 1, j0 : j1 + 1]
        gray_win = gray[i0 : i1 + 1, j0 : j1 + 1]
        cu = u - j0

        band = self._axis_band(win)
        strokes = self._strokes(win, band)
        threshold = self.threshold
        if not strokes:
            return None

        if band is not None:
            lo, hi = band
            v_axis = self._axis_centre(win, gray_win, band, cu, threshold)
            best: tuple[float, float] | None = None
            for stroke in strokes:
                u_tick = self._tick_centre(gray_win, stroke, threshold)
                if u_tick is None:
                    continue
                dist = abs(u_tick - cu)
                if best is None or dist < best[0]:
                    best = (dist, u_tick)
            if best is None:
                return None
            uu, vv = j0 + best[1], i0 + v_axis
            return CornerHit(uu, vv, "corner", math.hypot(uu - u, vv - v))

        # Fallback: a lone perpendicular stroke (the crop cut the frame away).
        best = None
        for stroke in strokes:
            u_tick = self._tick_centre(gray_win, stroke, threshold)
            if u_tick is None:
                continue
            dist = abs(u_tick - cu)
            if best is None or dist < best[0]:
                best = (dist, u_tick)
        if best is None:
            return None
        uu = j0 + best[1]
        return CornerHit(uu, v, "tick", abs(uu - u))

    # ------------------------------------------------------------- geometry
    @staticmethod
    def _axis_band(win: np.ndarray) -> tuple[int, int] | None:
        """Row range of the longest horizontal dark line, or None."""
        need = MIN_LINE_FRAC * win.shape[1]
        best_row, best_len = -1, 0
        for i in range(win.shape[0]):
            _start, _end, length = _longest_run(win[i])
            if length > best_len:
                best_row, best_len = i, length
        if best_row < 0 or best_len < need:
            return None
        lo = hi = best_row
        while lo - 1 >= 0 and win[lo - 1].sum() >= need and hi - lo + 2 <= MAX_LINE_THICKNESS:
            lo -= 1
        while hi + 1 < win.shape[0] and win[hi + 1].sum() >= need:
            if hi + 2 - lo > MAX_LINE_THICKNESS:
                break
            hi += 1
        return lo, hi

    @staticmethod
    def _strokes(win: np.ndarray, band: tuple[int, int] | None) -> list[_Stroke]:
        """Perpendicular dark strokes, as (row, col) pixel arrays.

        With a band found, the band rows are removed first so each tick becomes
        its own component; only components that actually touch the band count.
        Direction is the principal axis of the component, which is why a single
        short run cannot masquerade as a tick: the whole stroke is measured.
        """
        body = win.copy()
        if band is not None:
            body[band[0] : band[1] + 1, :] = False
        labels, count = ndimage.label(body, structure=_STRUCTURE)
        if count == 0:
            return []

        out: list[np.ndarray] = []
        for k in range(1, count + 1):
            ys, xs = np.nonzero(labels == k)
            if ys.size < TICK_MIN_PX:
                continue
            if band is not None:
                lo, hi = band
                touches = np.any(ys == lo - 1) or np.any(ys == hi + 1)
                if not touches:
                    continue
            if not _perpendicular(ys, xs):
                continue
            out.append(_Stroke(ys, xs))
        return out

    @staticmethod
    def _axis_centre(
        win: np.ndarray,
        gray_win: np.ndarray,
        band: tuple[int, int],
        cu: float,
        threshold: float,
    ) -> float:
        """Sub-pixel centre of the axis band, measured at a tick-free column."""
        lo, hi = band
        col = int(round(cu))
        for offset in (0, 6, -6, 12, -12, 20, -20):
            j = col + offset
            if not 0 <= j < win.shape[1]:
                continue
            if not win[lo : hi + 1, j].all():
                continue
            above = lo - 2 >= 0 and not win[lo - 2, j]
            below = hi + 2 < win.shape[0] and not win[hi + 2, j]
            if above and below:
                return _centre_from_profile(gray_win[:, j], lo, hi, threshold)
        return 0.5 * (lo + hi) + 1.0  # quantized fallback: band centre

    @staticmethod
    def _tick_centre(
        gray_win: np.ndarray, stroke: _Stroke, threshold: float
    ) -> float | None:
        """Sub-pixel centre of the tick stroke, probed at its far end."""
        probe = stroke.probe_row(gray_win.shape[0])
        if probe is None:
            return None
        return _centre_from_profile(
            gray_win[probe, :], stroke.xs.min(), stroke.xs.max(), threshold
        )


def _perpendicular(ys: np.ndarray, xs: np.ndarray) -> bool:
    """True when the stroke's principal direction is within ANGLE_TOL_DEG of 90."""
    r = ys.astype(np.float64)
    c = xs.astype(np.float64)
    r -= r.mean()
    c -= c.mean()
    cov = np.array([[float(c @ c), float(c @ r)], [float(r @ c), float(r @ r)]])
    evals, evecs = np.linalg.eigh(cov)
    dx, dy = evecs[:, int(np.argmax(evals))]
    angle = math.degrees(math.atan2(abs(float(dy)), abs(float(dx))))
    return abs(angle - 90.0) <= ANGLE_TOL_DEG


class _Stroke:
    """One candidate tick: its pixels plus where to probe it horizontally."""

    def __init__(self, ys: np.ndarray, xs: np.ndarray) -> None:
        self.ys = ys
        self.xs = xs

    def probe_row(self, height: int) -> int | None:
        """Row at the stroke's far end, where its width is measured."""
        row = int(self.ys.max()) if self.ys.mean() >= height / 2 else int(self.ys.min())
        row = min(max(row, 0), height - 1)
        if not np.any(self.ys == row):
            return None
        return row


def _crossing(samples: np.ndarray, i: int, j: int, threshold: float) -> float:
    """Fractional index where an intensity profile crosses ``threshold``."""
    if i < 0:
        return 0.0
    if j >= samples.size:
        return float(samples.size - 1)
    g_i, g_j = float(samples[i]), float(samples[j])
    if g_i == g_j:
        return float(j)
    t = (threshold - g_i) / (g_j - g_i)
    return float(i) + min(max(t, 0.0), 1.0)


def _centre_from_profile(
    profile: np.ndarray,
    lo: int,
    hi: int,
    threshold: float,
) -> float:
    """Centre of a dark band in *display* coordinates from a gray profile.

    ``lo``/``hi`` are the integer pixel indices of the band found on the binary
    mask; the value is refined by interpolating the grayscale crossing on both
    sides, then shifted by half a pixel because pixel ``j`` spans ``[j, j+1]``
    on screen.  A band at columns 120-121 therefore centres on 121.0.
    """
    lo = max(int(lo), 0)
    hi = min(int(hi), profile.size - 1)
    dark = profile < threshold

    # Walk outward from the known band to the first light pixel on each side.
    left = lo
    while left - 1 >= 0 and left - 1 >= lo - EDGE_MARGIN and dark[left - 1]:
        left -= 1
    right = hi
    while right + 1 < profile.size and right + 1 <= hi + EDGE_MARGIN and dark[right + 1]:
        right += 1

    edge_lo = _crossing(profile, left - 1, left, threshold)
    edge_hi = _crossing(profile, right, right + 1, threshold)
    return 0.5 * (edge_lo + edge_hi) + 0.5
