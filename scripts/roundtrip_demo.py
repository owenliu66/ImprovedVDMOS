#!/usr/bin/env python
"""Synthetic round-trip: generate data from known parameters, refit, compare.

This is the regression that proves the fitter is unbiased before you point it
at real datasheet curves. Run:  python scripts/roundtrip_demo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import numpy as np  # noqa: E402

from vdmos_fit import (  # noqa: E402
    DEFAULT,
    Params,
    fit,
    identifiability,
    make_synthetic,
    sensitivity_report,
)
from vdmos_fit.plotting import overlay_data  # noqa: E402

TRUTH = Params(
    kp=30.0, vth0=3.1, theta=0.02, vsat=1.2, eta=0.03, lam=0.015, gk=0.3, vthk=1.5
)


def main() -> int:
    vgs = np.linspace(2.0, 14.0, 25)
    vds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    data = make_synthetic(TRUTH, vgs, vds, noise=0.01, seed=7)

    print("pre-fit sensitivity (mean |d ln Id / d ln p|):")
    for k, v in sensitivity_report(DEFAULT, data).items():
        print(f"  {k:<8s} {v:10.3g}")

    res = fit(data, start=DEFAULT)
    print()
    print(res.summary())

    print("\nrecovered vs truth:")
    print(f"{'param':<10s} {'truth':>12s} {'fitted':>12s} {'error %':>10s}")
    for name in ("kp", "vth0", "theta", "vsat", "eta", "lam", "gk", "vthk"):
        want, got = getattr(TRUTH, name), getattr(res.params, name)
        err = 100.0 * (got - want) / want
        print(f"{name:<10s} {want:12.4g} {got:12.4g} {err:10.2f}")

    print()
    print(identifiability(res).summary())
    print(
        "\nnote: a large vthk error here is expected, not a fitter bug. The kink"
        "\ntransition (vthk = 1.5 V, kk = 0.1 V) lies below the lowest sampled Vgs,"
        "\nso this grid carries no information about it. Add dense Vgs points through"
        "\nthe kink region, or fix vthk from measurement, before trusting it."
    )

    out = Path("out")
    out.mkdir(exist_ok=True)
    overlay_data(data, res.params, title="Round-trip fit").figure.savefig(
        out / "roundtrip.png", dpi=150
    )
    print(f"\nwrote {out / 'roundtrip.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
