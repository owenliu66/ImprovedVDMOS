#!/usr/bin/env python
"""Thin entry point: fit a measured CSV and write plots + fitted parameters.

All logic lives in the `vdmos_fit` package; this file only wires arguments
to library calls. Run as:  python scripts/fit_datasheet.py data/IPW065N04N.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

# Allow running without `pip install -e .` (development convenience).
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from vdmos_fit import DEFAULT, fit, identifiability, load_csv, sensitivity_report  # noqa: E402
from vdmos_fit.plotting import overlay_data  # noqa: E402


def main(csv_path: str) -> int:
    data = load_csv(csv_path)
    print(f"{len(data)} points from {csv_path}")

    # Two-stage: relative weighting first (captures subthreshold + full range),
    # then refine with ampere weighting so the Ron corner dominates.
    stage1 = fit(data, start=DEFAULT, weight="relative")
    stage2 = fit(data, start=stage1.params, weight="abs")

    print("\n--- stage 1 (relative) ---")
    print(stage1.summary())
    print("\n--- stage 2 (ampere-weighted) ---")
    print(stage2.summary())

    print("\npre-fit sensitivity (what this data can actually see):")
    for k, v in sensitivity_report(DEFAULT, data).items():
        print(f"  {k:<8s} {v:10.3g}")
    print()
    print(identifiability(stage2).summary())

    out = Path("out")
    out.mkdir(exist_ok=True)
    overlay_data(data, stage2.params).figure.savefig(out / "fit_vs_data.png", dpi=150)
    (out / "fitted_params.txt").write_text(stage2.summary(), encoding="utf-8")
    print(f"\nwrote {out / 'fit_vs_data.png'} and {out / 'fitted_params.txt'}")
    return 0 if stage2.success else 1


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("usage: python scripts/fit_datasheet.py <data.csv>")
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
