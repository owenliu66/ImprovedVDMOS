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
    k_p=30.0, v_th0=3.1, k_theta=0.02, v_vsat=1.2, k_eta=0.03, k_lambda=0.015, g_k=0.3, v_thk=1.5
)


def main() -> int:
    v_gs = np.linspace(2.0, 14.0, 25)
    v_ds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    data = make_synthetic(TRUTH, v_gs, v_ds, noise=0.01, seed=7)

    print("pre-fit sensitivity (mean |d ln i_d / d ln p|):")
    for k, v in sensitivity_report(DEFAULT, data).items():
        print(f"  {k:<12s} {v:10.3g}")

    res = fit(data, start=DEFAULT)
    print()
    print(res.summary())

    print("\nrecovered vs truth:")
    print(f"{'param':<12s} {'truth':>12s} {'fitted':>12s} {'error %':>10s}")
    for name in ("k_p", "v_th0", "k_theta", "v_vsat", "k_eta", "k_lambda", "g_k", "v_thk"):
        want, got = getattr(TRUTH, name), getattr(res.params, name)
        err = 100.0 * (got - want) / want
        print(f"{name:<12s} {want:12.4g} {got:12.4g} {err:10.2f}")

    print()
    print(identifiability(res).summary())
    print(
        "\nnote: a large v_thk error here is expected, not a fitter bug. The kink"
        "\ntransition (v_thk = 1.5 V, v_kovsoft = 0.1 V) lies below the lowest sampled v_gs,"
        "\nso this grid carries no information about it. Add dense v_gs points through"
        "\nthe kink region, or fix v_thk from measurement, before trusting it."
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
