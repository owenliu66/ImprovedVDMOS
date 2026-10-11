"""Pixel <-> data coordinate transforms for a digitized datasheet plot.

The digitizer keeps the datasheet image geometrically faithful: the matplotlib
axes live in *crop pixel* space (u to the right, v downward) and the image is
drawn with ``imshow(extent=(0, W, H, 0))``.  Data values only ever appear as
tick labels, produced by the :class:`AxisTransform` below.

This matters because a log axis printed on a datasheet already has its decade
spacing compressed into pixels; warping the image onto a log-scaled matplotlib
axis would distort it.  Keeping the image in pixel space and mapping pixels ->
data with a linear or log law is the physically correct model.
"""

from __future__ import annotations

import math

LINEAR = "linear"
LOG = "log"
SCALES = (LINEAR, LOG)


def format_value(value: float) -> str:
    """Compact, human-readable label for an axis tick or table cell."""
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e5 or magnitude < 1e-3:
        return f"{value:.3g}"
    return f"{value:.6g}"


def nice_step(span: float, target: int = 8) -> float:
    """Round tick step for an axis span (1 / 2 / 2.5 / 5 decades)."""
    if span <= 0:
        return 1.0
    raw = span / max(target, 1)
    magnitude = 10.0 ** math.floor(math.log10(raw))
    for mult in (1.0, 2.0, 2.5, 5.0, 10.0):
        if raw <= mult * magnitude:
            return mult * magnitude
    return 10.0 * magnitude


class AxisTransform:
    """Maps one axis between crop-pixel space and datasheet data space.

    ``pix_lo``/``pix_hi`` are fixed by the crop geometry (for X: 0 and W; for Y:
    H and 0, because pixel rows grow downward).  ``data_lo``/``data_hi`` are the
    axis min/max the user enters at those two edges.
    """

    def __init__(
        self,
        pix_lo: float,
        pix_hi: float,
        data_lo: float,
        data_hi: float,
        scale: str = LINEAR,
    ) -> None:
        if pix_hi == pix_lo:
            raise ValueError("pixel extent is degenerate")
        if scale not in SCALES:
            raise ValueError(f"unknown scale {scale!r}")
        if scale == LOG and (data_lo <= 0 or data_hi <= 0):
            raise ValueError("log axis requires strictly positive min and max")
        self.pix_lo = float(pix_lo)
        self.pix_hi = float(pix_hi)
        self.data_lo = float(data_lo)
        self.data_hi = float(data_hi)
        self.scale = scale
        self.lo = min(data_lo, data_hi)
        self.hi = max(data_lo, data_hi)
        self._log_lo = math.log10(self.lo) if self.lo > 0 else None
        self._log_hi = math.log10(self.hi) if self.hi > 0 else None

    @classmethod
    def from_points(
        cls,
        pix1: float,
        pix2: float,
        data1: float,
        data2: float,
        scale: str = LINEAR,
    ) -> AxisTransform:
        """Build the map from two clicked pixels and their known data values.

        Directional: ``pix1`` carries ``data1`` and ``pix2`` carries ``data2``,
        so clicking the ticks in reverse order simply reverses the axis.  For a
        log axis the two points fix the decade scale and offset, which is
        exactly the pixels -> ``log10(data)`` law implemented below.
        """
        if pix1 == pix2:
            raise ValueError("the two calibration points are at the same pixel")
        return cls(pix1, pix2, data1, data2, scale)

    # ------------------------------------------------------------------ maps
    def to_data(self, pix: float) -> float:
        t = (pix - self.pix_lo) / (self.pix_hi - self.pix_lo)
        if self.scale == LOG:
            a, b = math.log10(self.data_lo), math.log10(self.data_hi)
            return 10.0 ** (a + t * (b - a))
        return self.data_lo + t * (self.data_hi - self.data_lo)

    def to_pix(self, data: float) -> float:
        if self.scale == LOG:
            if data <= 0:
                return float("nan")
            a, b = math.log10(self.data_lo), math.log10(self.data_hi)
            t = (math.log10(data) - a) / (b - a)
        else:
            t = (data - self.data_lo) / (self.data_hi - self.data_lo)
        return self.pix_lo + t * (self.pix_hi - self.pix_lo)

    # ----------------------------------------------------------------- ticks
    def major_ticks(self, interval: float | None = None) -> list[float]:
        """Data values of the major ticks inside [lo, hi].

        Linear: ``interval`` is the data-space step.  Log: ``interval`` is the
        number of decades between labelled ticks (1 = every decade).
        """
        span = self.hi - self.lo
        if span <= 0:
            return [self.lo]

        if self.scale == LOG:
            step = 1 if interval is None else max(1, int(round(interval)))
            start = math.floor(self._log_lo)
            stop = math.ceil(self._log_hi)
            values = []
            for k in range(start, stop + 1, step):
                value = 10.0**k
                if self.lo * 0.999999 <= value <= self.hi * 1.000001:
                    values.append(value)
            return values or [self.lo, self.hi]

        step = nice_step(span) if not interval or interval <= 0 else float(interval)
        if span / step > 400:  # guard against a pathological tick step
            step = nice_step(span)
        first = math.ceil(self.lo / step - 1e-9) * step
        values = []
        k = 0
        tol = 1e-9 * max(1.0, abs(self.hi))
        while first + k * step <= self.hi + tol and len(values) < 400:
            values.append(first + k * step)
            k += 1
        return values

    def minor_ticks(
        self,
        majors: list[float],
        subdivisions: int = 5,
        log_mantissas: tuple[int, ...] = (2, 3, 4, 5, 6, 7, 8, 9),
    ) -> list[float]:
        """Data values of the minor ticks (never including the majors)."""
        if not majors or subdivisions <= 1 and self.scale == LINEAR:
            return []

        if self.scale == LOG:
            if not log_mantissas:
                return []
            values = []
            for k in range(math.floor(self._log_lo) - 1, math.ceil(self._log_hi) + 1):
                for m in log_mantissas:
                    value = m * 10.0**k
                    if self.lo < value < self.hi:
                        values.append(value)
            return values

        if len(majors) < 2:
            return []
        step = majors[1] - majors[0]
        if step <= 0:
            return []
        sub = step / subdivisions
        values = []
        for major in majors:
            for j in range(1, subdivisions):
                value = major + j * sub
                if self.lo < value < self.hi:
                    values.append(value)
        return values
