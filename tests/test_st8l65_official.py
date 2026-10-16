"""Pins the OFFICIAL ST8L65N044M9 record (vdmos_fit.reference).

Promoted 2026-10-13: fit A of scripts/fit_st8l65_nd.py -- 10-param with the
drain-knee softness exponent `k_kneesharp` free, fitted on the full 41-point measured
set, 3.79 % relative RMSE. These tests do NOT re-run the fit (too slow for a
unit suite); they evaluate the recorded parameter vector against the measured
CSV and against the model's physical invariants, so a model/bounds/weighting
change that silently moves the official fit fails here instead of quietly
invalidating every number quoted in README.md and STATUS.md.

Deliberately NOT pinned: the individual knee parameters (k_p, v_ovmax, k_ovsharp,
k_kneesharp, v_dbsat, k_lambda are entangled at |corr| up to 0.999 -- only the
SHAPE is identifiable) and
k_lambda/v_dbsat sensitivity (sigma = inf is expected: the data stops at v_ds = 20 V).
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
RON_CHORD_MOHM = 34.5      # measured i_d chord at v_gs = 10 V
PLATEAU_A = 240.2          # v_gs = 10 V, v_ds = 19.7 V
VTH_UA = 250.0             # v_gs = v_ds = 3.7 V (v_gs(th) spec condition)


@pytest.fixture(scope="module")
def data():
    return load_csv(CSV)


def test_official_record_reproduces_on_measured_data(data):
    """The recorded parameters still reproduce the published 3.79 % RMSE."""
    model = vdmos_id(data.v_ds, data.v_gs, ST8L65_OFFICIAL)
    rel = (model - data.i_d) / data.i_d
    rmse = float(np.sqrt(np.mean(rel**2)))
    assert rmse == pytest.approx(ST8L65_OFFICIAL_RMSE_REL, abs=2e-3)
    assert rmse < 0.040, "official fit regressed"


def test_official_record_hits_the_measured_anchors(data):
    """Ron chord, saturation plateau and the v_gs(th) point stay in tolerance.

    r_d + r_s is a LOWER bound on the Rds(on) chord (the rest of the chord is
    the intrinsic channel), so it is checked from below only.
    """
    p = ST8L65_OFFICIAL
    assert (p.r_d + p.r_s) * 1e3 < RON_CHORD_MOHM
    assert (p.r_d + p.r_s) * 1e3 == pytest.approx(31.7, abs=1.0)

    plateau = float(vdmos_id(19.7, 10.0, p))
    assert plateau == pytest.approx(PLATEAU_A, rel=0.035)

    v_th = float(vdmos_id(3.7, 3.7, p)) * 1e6
    assert v_th == pytest.approx(VTH_UA, rel=0.06)

    # The v_gs = 5 curve is unreproducible by any C1 model (data/README.md);
    # every other curve must be tracked to within a few percent.
    rel = (vdmos_id(data.v_ds, data.v_gs, p) - data.i_d) / data.i_d
    for vg in data.v_gs_grid():
        if np.isclose(vg, 5.0):
            continue
        r = float(np.sqrt(np.mean(rel[np.isclose(data.v_gs, vg)] ** 2)))
        assert r < 0.045, f"per-curve RMSE regressed at v_gs = {vg:g} V"


def test_official_parameters_are_inside_the_global_bounds():
    """A record that sits on a bound is not a physical result.

    k_lambda is the one exception the dataset forces (it wants zero CLM and pins
    at the 1e-12 floor); everything else, k_kneesharp included, must be interior.
    """
    for name in ST8L65_FIT_A_FREE:
        if name == "k_lambda":
            continue
        lo, hi = BOUNDS[name]
        val = getattr(ST8L65_OFFICIAL, name)
        assert lo < val < hi, f"{name} = {val:g} rides the ({lo}, {hi}) bound"


def test_official_nd_is_interior_and_not_the_legacy_sentinel():
    """k_kneesharp = 0.436 is a fitted value, not the k_kneesharp = 0 tanh sentinel and not a
    bound-pin. The data wants a knee softer than the legacy tanh (2.8) and
    softer than the retired 0.5 guard allowed."""
    k_kneesharp = ST8L65_OFFICIAL.k_kneesharp
    assert BOUNDS["k_kneesharp"][0] < k_kneesharp < 1.0
    assert k_kneesharp != pytest.approx(0.0)
    assert k_kneesharp != pytest.approx(2.8, abs=0.2), "k_kneesharp regressed to the legacy tanh"


def test_official_model_is_finite_and_monotone(data):
    """No NaN/inf and no negative gds at the official set.

    The only negative gds anywhere is finite-difference noise at the ~1e-14 A
    subthreshold floor (measured min -3.1e-12 S at v_gs = -0.8 / v_ds = 95.7 V),
    so the tolerance is far below any physical slope, not a free pass for the
    v_vsat-divisor dip (which cannot occur here: v_vsat is fixed at 1e4).
    """
    p = ST8L65_OFFICIAL
    v_ds = np.linspace(0.0, 200.0, 401)
    for vg in np.linspace(-2.0, 16.0, 10):
        i_d = np.asarray(vdmos_id(v_ds, vg, p))
        assert np.all(np.isfinite(i_d))
        assert np.gradient(i_d, v_ds).min() > -1e-6, f"gds went negative at v_gs={vg:g}"

    assert float(vdmos_id(0.0, 0.0, p)) == 0.0
    assert float(vdmos_id(20.0, -2.0, p)) < 1e-12   # off state, no leak element
