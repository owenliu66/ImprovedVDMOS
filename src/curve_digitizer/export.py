"""CSV export of digitized markers."""

from __future__ import annotations

import csv
from pathlib import Path

HEADER = ("index", "x", "y", "condition")


def default_data_dir() -> Path:
    """Projects/ImprovedVDMOS/data, resolved from this package's location."""
    here = Path(__file__).resolve()
    candidate = here.parents[2] / "data"  # src/curve_digitizer/<file> -> project root
    if candidate.parent.name.lower() == "improvedvdmos":
        candidate.mkdir(parents=True, exist_ok=True)
        return candidate
    fallback = Path.cwd() / "data"
    fallback.mkdir(parents=True, exist_ok=True)
    return fallback


def export_csv(path: str | Path, markers, xform_x, xform_y) -> Path:
    """Write ``index,x,y,condition`` rows in datasheet data units."""
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(HEADER)
        for marker in markers:
            writer.writerow(
                [
                    marker.index + 1,
                    f"{xform_x.to_data(marker.u):.10g}",
                    f"{xform_y.to_data(marker.v):.10g}",
                    marker.condition,
                ]
            )
    return target
