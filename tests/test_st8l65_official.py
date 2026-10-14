"""Pins the OFFICIAL ST8L65N044M9 record (vdmos_fit.reference).

Promoted 2026-10-13: fit A of scripts/fit_st8l65_nd.py -- 10-param with the
drain-knee softness exponent `nd` free, fitted on the full 41-point measured
set, 3.79 % relative RMSE. These tests do NOT re-run the fit (too slow for a
unit suite); they evaluate the recorded parameter vector against the measured
CSV and against the model's physical invariants, so a model/bounds/weighting
change that silently moves the official fit fails here instead of quietly
invalidating every number quoted in README.md and STATUS.md.

Deliberately NOT pinned: the individual knee parameters (kp, kcap, m, nd, vk,
lam are entangled at |corr| up to 0.999 -- only the SHAPE is identifiable) and
lam/vk sensitivity (sigma = inf is expected: the data stops at Vds = 20 V).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from vdmos_fit import (
    BOUNDS,
    ST8L65_FIT_A_FREE,
    ST8L65_OFFICIAL,
    ST8L65_OFFICIAL_RMSE_REL,
    load_csv,
    vdmos_id,
)

CSV = Path(__file__).resolve().parents[1] / "data" / "ST8L65N044M9_IV.csv"

# Measured anchors from data/README.md / the datasheet.
RON_CHORD_MOHM = 34.5      # measured Ids chord at Vgs = 10 V
PLATEAU_A = 240.2          # Vgs = 10 V, Vds = 19.7 V
VTH_UA = 250.0             # Vgs = Vds = 3.7 V (Vgs(th) spec condition)


@pytest.fixture(scope="module")
def data():
    return load_csv(CSV)


def test_official_record_reproduces_on_measured_data(data):
    """The recorded parameters still reproduce the published 3.79 % RMSE."""
    model = vdmos_id(data.vds, data.vgs, ST8L65_OFFICIAL)
    rel = (model - data.ids) / data.ids
    rmse = float(np.sqrt(np.mean(rel**2)))
    assert rmse == pytest.approx(ST8L65_OFFICIAL_RMSE_REL, abs=2e-3)
    assert rmse < 0.040, "official fit regressed"


def test_official_record_hits_the_measured_anchors(data):
    """Ron chord, saturation plateau and the Vgs(th) point stay in tolerance.

    rd + rs is a LOWER bound on the Rds(on) chord (the rest of the chord is
    the intrinsic channel), so it is checked from below only.
    """
    p = ST8L65_OFFICIAL
    assert (p.rd + p.rs) * 1e3 < RON_CHORD_MOHM
    assert (p.rd + p.rs) * 1e3 == pytest.approx(31.7, abs=1.0)

    plateau = float(vdmos_id(19.7, 10.0, p))
    assert plateau == pytest.approx(PLATEAU_A, rel=0.035)

    vth = float(vdmos_id(3.7, 3.7, p)) * 1e6
    assert vth == pytest.approx(VTH_UA, rel=0.06)

    # The Vgs = 5 curve is unreproducible by any C1 model (data/README.md);
    # every other curve must be tracked to within a few percent.
    rel = (vdmos_id(data.vds, data.vgs, p) - data.ids) / data.ids
    for vg in data.vgs_grid():
        if np.isclose(vg, 5.0):
            continue
        r = float(np.sqrt(np.mean(rel[np.isclose(data.vgs, vg)] ** 2)))
        assert r < 0.045, f"per-curve RMSE regressed at Vgs = {vg:g} V"


def test_official_parameters_are_inside_the_global_bounds():
    """A record that sits on a bound is not a physical result.

    lam is the one exception the dataset forces (it wants zero CLM and pins
    at the 1e-12 floor); everything else, nd included, must be interior.
    """
    for name in ST8L65_FIT_A_FREE:
        if name == "lam":
            continue
        lo, hi = BOUNDS[name]
        val = getattr(ST8L65_OFFICIAL, name)
        assert lo < val < hi, f"{name} = {val:g} rides the ({lo}, {hi}) bound"


def test_official_nd_is_interior_and_not_the_legacy_sentinel():
    """nd = 0.436 is a fitted value, not the nd = 0 tanh sentinel and not a
    bound-pin. The data wants a knee softer than the legacy tanh (2.8) and
    softer than the retired 0.5 guard allowed."""
    nd = ST8L65_OFFICIAL.nd
    assert BOUNDS["nd"][0] < nd < 1.0
    assert nd != pytest.approx(0.0)
    assert nd != pytest.approx(2.8, abs=0.2), "nd regressed to the legacy tanh"


def test_official_model_is_finite_and_monotone(data):
    """No NaN/inf and no negative gds at the official set.

    The only negative gds anywhere is finite-difference noise at the ~1e-14 A
    subthreshold floor (measured min -3.1e-12 S at Vgs = -0.8 / Vds = 95.7 V),
    so the tolerance is far below any physical slope, not a free pass for the
    vsat-divisor dip (which cannot occur here: vsat is fixed at 1e4).
    """
    p = ST8L65_OFFICIAL
    vds = np.linspace(0.0, 200.0, 401)
    for vg in np.linspace(-2.0, 16.0, 10):
        ids = np.asarray(vdmos_id(vds, vg, p))
        assert np.all(np.isfinite(ids))
        assert np.gradient(ids, vds).min() > -1e-6, f"gds went negative at Vgs={vg:g}"

    assert float(vdmos_id(0.0, 0.0, p)) == 0.0
    assert float(vdmos_id(20.0, -2.0, p)) < 1e-12   # off state, no leak element
