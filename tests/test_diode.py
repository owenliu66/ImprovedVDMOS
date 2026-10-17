"""Body-diode tests: LTspice-faithful equation, series-R solver, sentinel
bit-identity, integration into vdmos_id, extraction round trip, and the
pinned ST8L65 diode record.

Equation reference: LTspice diode DC branch (Is/N/Isr/Nr/Ikf/Rs + Tnom/Xti/Eg),
verified against LTspice 26.0.1 in test_ltspice_cross_validation (skipped when
LTspice is not installed).
"""

from __future__ import annotations

import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import brentq

from vdmos_fit import (
    DEFAULT,
    ST8L65_DIODE,
    ST8L65_DIODE_FREE,
    Params,
    diode_sensitivity_report,
    fit_diode,
    identifiability,
    load_diode_csv,
    make_synthetic_diode,
    vdmos_diode,
    vdmos_id,
)
from vdmos_fit.diode import _diode_current, _scaled_saturation
from vdmos_fit.model import _channel_current, _solve_package_r

ROOT = Path(__file__).resolve().parents[1]
DIODE_CSV = ROOT / "data" / "ST8L65N044M9_diode.csv"

# A representative "on" set (close to the ST8L65 record, r_diode nonzero).
P_ON = DEFAULT.replace(i_s=3.7e-6, k_n=2.0, i_kf=17.5, r_diode=1.6e-3,
                       v_eg=1.10, k_xti=3.0, t_nom=25.0)

Q = 1.602176634e-19
K_B = 1.380649e-23


def naive_diode(v_di, p: Params, t_c: float = 25.0) -> float:
    """Straightforward float64 evaluation of the documented equation.

    Only valid where exp() stays tame (|v_di| well under ~1.2 V) -- used for
    hand-value checks, never for solver brackets.
    """
    t_k = t_c + 273.15
    vt = K_B * t_k / Q
    i_s_t = _scaled_saturation(p.i_s, p.k_n, p, np.array([t_k]))[0]
    x = np.clip(v_di / (p.k_n * vt), -700.0, 700.0)  # bracket-endpoint safety
    j = i_s_t * np.expm1(x)
    if p.i_kf > 0.0 and j > 0.0:
        j = j / np.sqrt(1.0 + j / p.i_kf)
    if p.i_sr > 0.0:
        i_sr_t = _scaled_saturation(p.i_sr, p.k_nr, p, np.array([t_k]))[0]
        j += i_sr_t * np.expm1(np.clip(v_di / (p.k_nr * vt), -700.0, 700.0))
    return j


# ---------------------------------------------------------------------------
# Sentinel / bit-identity (protects the 103-test channel anchor set)
# ---------------------------------------------------------------------------

def test_off_sentinel_bit_identical():
    """Diode-off vdmos_id == the pre-diode model, bit for bit, on v_ds >= 0."""
    vd = np.linspace(0.0, 20.0, 41)
    vg = np.linspace(-2.0, 16.0, 31)
    VD, VG = np.meshgrid(vd, vg)
    assert np.array_equal(_channel_current(VD, VG, DEFAULT), vdmos_id(VD, VG, DEFAULT))
    pkg = DEFAULT.replace(r_d=2e-3, r_s=1e-3)
    ref = _solve_package_r(VD, VG, pkg.r_d, pkg.r_s,
                           lambda d, g: _channel_current(d, g, pkg))
    assert np.array_equal(ref, vdmos_id(VD, VG, pkg))


def test_channel_clamped_zero_for_negative_vds():
    """v_ds <= 0 with the diode off -> exactly zero (no garbage channel tail)."""
    out = vdmos_id(np.linspace(-20.0, -0.01, 25), np.linspace(-2, 16, 25), DEFAULT)
    assert np.all(out == 0.0)


# ---------------------------------------------------------------------------
# The equation itself
# ---------------------------------------------------------------------------

def test_explicit_equation_hand_values():
    """_diode_current matches a naive evaluation, both branches, at 25 degC."""
    p = P_ON.replace(i_sr=1e-8, k_nr=2.0, r_diode=0.0)
    vs = np.array([-1.0, -0.3, 0.0, 0.05, 0.3, 0.55, 0.8, 1.05])
    mine = _diode_current(vs, p, np.full(vs.shape, 298.15))
    ref = np.array([naive_diode(v, p) for v in vs])
    np.testing.assert_allclose(mine, ref, rtol=1e-9, atol=1e-15)
    assert mine[vs == 0.0][0] == 0.0  # exact zero at v_di = 0 (bracket proof)


def test_high_injection_asymptote():
    """J >> i_kf: i_diffuse -> sqrt(i_s(T) * expm1(x) * i_kf) (sqrt compression)."""
    p = DEFAULT.replace(i_s=1e-3, k_n=1.0, i_kf=1e-3)  # J ~ amps >> 1 mA
    t_k = np.array([298.15])
    v = 0.6
    i_model = float(_diode_current(np.array([v]), p, t_k)[0])
    x = v / (K_B * 298.15 / Q)
    i_asym = np.sqrt(p.i_s * np.expm1(x) * p.i_kf)
    # exact relation I = J/sqrt(1+J/i_kf); asymptote within 1 % when J/i_kf ~ 1e3
    assert i_model == pytest.approx(i_asym, rel=0.01)
    # and the model is BELOW the uncompressed J (compression is one-sided)
    assert i_model < p.i_s * np.expm1(x)


def test_temperature_law_hand_values():
    """i_s(T) law vs hand-computed factors at -55/25/150 degC (divisor k_n)."""
    p = DEFAULT.replace(i_s=1.0, k_n=2.0, k_xti=3.0, v_eg=1.10, t_nom=25.0)
    t_nom_k = 298.15
    for t_c in (-55.0, 25.0, 150.0):
        t_k = t_c + 273.15
        manual = (t_k / t_nom_k) ** (3.0 / 2.0) * np.exp(
            1.10 * (Q / K_B) / 2.0 * (1.0 / t_nom_k - 1.0 / t_k))
        got = float(_scaled_saturation(1.0, 2.0, p, np.array([t_k]))[0])
        assert got == pytest.approx(manual, rel=1e-9)
    # negative dV/dT: forward drop at fixed current falls as T rises
    i_test = 20.0

    def v_at(t):
        return brentq(lambda v: float(vdmos_diode(v, P_ON, t_c=t)) - i_test, 0.1, 1.5)

    assert v_at(150.0) < v_at(25.0) < v_at(-55.0)


# ---------------------------------------------------------------------------
# The r_diode series solver
# ---------------------------------------------------------------------------

def test_series_r_solver_vs_brentq():
    """Bracketed Newton vs independent brentq roots, both polarities."""
    p = P_ON

    def root(v_term):
        def f(v_di):
            return naive_diode(v_di, p.replace(r_diode=0.0)) - (v_term - v_di) / p.r_diode

        lo = min(0.0, v_term / p.r_diode)
        hi = max(0.0, v_term / p.r_diode)
        x = brentq(f, lo, hi, xtol=1e-14, rtol=8.9e-16)
        return (v_term - x) / p.r_diode  # current, not internal voltage

    vs = np.array([-1.0, -0.01, 0.0, 0.05, 0.4, 0.73, 1.05])
    sol = vdmos_diode(vs, p)
    ref = np.array([root(v) for v in vs])
    msk = vs != 0.0
    # atol floor: the solver converges on |g| <= 1e-12 A, so relative accuracy
    # at microampere reverse currents is bounded by that absolute tolerance.
    np.testing.assert_allclose(sol[msk], ref[msk], rtol=1e-6, atol=1e-12)
    assert sol[vs == 0.0][0] == 0.0


def test_monotone_and_reverse_floor():
    """Monotone in v_sd, finite on [-2, 2] V, reverse floor = i_s(T) + i_sr(T)."""
    p = P_ON.replace(i_sr=1e-9)
    vs = np.linspace(-2.0, 2.0, 401)
    i = vdmos_diode(vs, p, t_c=25.0)
    assert np.all(np.isfinite(i))
    # non-decreasing: strictly rising in forward, FLAT at reverse saturation
    # (physical), up to the solver's 1e-12 A absolute tolerance
    assert np.all(np.diff(i) >= -1e-12)
    t_k = np.array([298.15])
    floor = float(_scaled_saturation(p.i_s, p.k_n, p, t_k)[0]) \
        + float(_scaled_saturation(p.i_sr, p.k_nr, p, t_k)[0])
    # at -2 V reverse, expm1(-2/(k_n*VT)) = -1 to ~1e-17
    assert float(vdmos_diode(-2.0, p, t_c=25.0)) == pytest.approx(-floor, rel=1e-9)


def test_overflow_safety():
    """No NaN/inf at absurd voltages and extreme bound combinations."""
    p = DEFAULT.replace(i_s=1e-22, k_n=1.0, i_kf=1e-3)
    assert np.isfinite(vdmos_diode(3.0, p))
    assert np.isfinite(vdmos_diode(100.0, p))
    extreme = p.replace(v_eg=2.0, k_xti=10.0)
    assert np.isfinite(vdmos_diode(100.0, extreme, t_c=150.0))
    assert np.isfinite(vdmos_diode(-1e3, extreme, t_c=-55.0))


# ---------------------------------------------------------------------------
# Integration into vdmos_id
# ---------------------------------------------------------------------------

def test_reverse_polarity_carried_by_diode():
    """vdmos_id(v_ds < 0) == -vdmos_diode(-v_ds), independent of v_gs."""
    for v_gs in (-2.0, 0.0, 10.0):
        lhs = vdmos_id(-0.8, v_gs, P_ON)
        rhs = -vdmos_diode(0.8, P_ON)
        assert lhs == pytest.approx(rhs, rel=1e-12)


def test_forward_floor_is_diode_reverse_saturation():
    """At v_ds > 0 the diode adds its (tiny) reverse saturation to i_d."""
    p = P_ON.replace(r_d=0.0, r_s=0.0)
    i_on = float(vdmos_id(20.0, 10.0, p))
    i_off = float(vdmos_id(20.0, 10.0, p.replace(i_s=0.0, i_sr=0.0)))
    # diode at v_sd = -20 V carries -i_s (expm1 -> -1), so i_d gains +i_s
    assert i_on - i_off == pytest.approx(p.i_s, rel=1e-3)


def test_straddling_grid_with_package_r():
    """Mixed-sign v_ds grids with r_d/r_s: finite, no RuntimeError, continuous."""
    p = P_ON.replace(r_d=24e-3, r_s=7e-3)
    vd = np.array([-1.0, -1e-6, 0.0, 1e-6, 0.5, 20.0])
    out = vdmos_id(vd, np.zeros_like(vd), p)
    assert np.all(np.isfinite(out))
    assert out[0] < 0.0 < out[4]
    # continuity across v_ds = 0 (channel 0, diode continuous)
    a = float(vdmos_id(-1e-9, 0.0, p))
    b = float(vdmos_id(1e-9, 0.0, p))
    assert abs(a - b) < 1e-6


# ---------------------------------------------------------------------------
# Data I/O and extraction
# ---------------------------------------------------------------------------

def test_load_diode_csv():
    data = load_diode_csv(DIODE_CSV)
    assert len(data) == 30
    assert data.temp_grid() == [-55.0, 25.0, 150.0]
    assert data.i_sd.min() > 4.0 and data.i_sd.max() < 56.0
    assert data.label == "ST8L65N044M9_diode"


def test_diode_round_trip_recovers_params():
    """fit_diode recovers known diode params from noiseless synthetic data."""
    truth = DEFAULT.replace(i_s=3.7e-6, k_n=2.0, i_kf=15.0, r_diode=1.5e-3,
                            v_eg=1.10, k_xti=3.0, t_nom=25.0)
    data = make_synthetic_diode(truth, np.linspace(5.0, 55.0, 10),
                                [-55.0, 25.0, 150.0])
    assert len(data) == 30
    start = DEFAULT.replace(k_n=2.0, i_s=1e-8, i_kf=3.0, r_diode=0.003,
                            v_eg=0.8, k_xti=3.0, t_nom=25.0)
    res = fit_diode(data, start=start, free=("i_s", "r_diode", "v_eg", "i_kf"),
                    bounds_overrides={"r_diode": (0.0, 0.020)})
    p = res.params
    assert res.rmse_rel < 0.01
    assert p.r_diode == pytest.approx(1.5e-3, rel=0.05)
    assert p.i_kf == pytest.approx(15.0, rel=0.10)
    assert p.v_eg == pytest.approx(1.10, rel=0.05)
    # i_s is the weak prefactor: order-of-magnitude recovery at best
    assert p.i_s == pytest.approx(3.7e-6, rel=0.9)
    ident = identifiability(res, ("i_s", "r_diode", "v_eg", "i_kf"))
    assert ident.rank == 4


def test_diode_round_trip_noisy():
    """1 % multiplicative noise on v_sd still lands within a few percent."""
    truth = DEFAULT.replace(i_s=3.7e-6, k_n=2.0, i_kf=15.0, r_diode=1.5e-3,
                            v_eg=1.10, k_xti=3.0, t_nom=25.0)
    data = make_synthetic_diode(truth, np.linspace(5.0, 55.0, 10),
                                [-55.0, 25.0, 150.0], noise=0.01, seed=3)
    start = DEFAULT.replace(k_n=2.0, i_s=1e-8, i_kf=3.0, r_diode=0.003,
                            v_eg=0.8, k_xti=3.0, t_nom=25.0)
    res = fit_diode(data, start=start, free=("i_s", "r_diode", "v_eg", "i_kf"),
                    bounds_overrides={"r_diode": (0.0, 0.020)})
    # 1 % noise on v_sd is ~7 mV -- through the exponential that is a ~15-20 %
    # current residual. The VOLTAGE data is what's clean; RMSE on current is
    # noise-amplified by design (same reason relative weighting matters).
    assert res.rmse_rel < 0.25
    assert res.params.r_diode == pytest.approx(1.5e-3, rel=0.2)
    assert res.params.v_eg == pytest.approx(1.10, rel=0.15)


def test_diode_sensitivity_report_live_columns():
    """All four fit-A parameters have nonzero pre-fit sensitivity at the record."""
    data = load_diode_csv(DIODE_CSV)
    sens = diode_sensitivity_report(ST8L65_DIODE, data, ST8L65_DIODE_FREE)
    assert all(v > 1e-3 for v in sens.values())


# ---------------------------------------------------------------------------
# Pinned record (mirrors test_st8l65_official for the channel)
# ---------------------------------------------------------------------------

def test_st8l65_diode_record_reproduces():
    """The promoted record still reproduces the measured CSV within tolerance."""
    data = load_diode_csv(DIODE_CSV)
    model = vdmos_diode(data.v_sd, ST8L65_DIODE, t_c=data.temp_c)
    rel = (model - data.i_sd) / data.i_sd
    assert np.sqrt(np.mean(rel**2)) < 0.049, "official diode fit regressed"
    # per-temperature sanity
    for t in (-55.0, 25.0, 150.0):
        msk = np.isclose(data.temp_c, t)
        assert np.sqrt(np.mean(rel[msk] ** 2)) < 0.08
    # interior values, physical reverse floor, sensed v_eg
    assert 0.0 < ST8L65_DIODE.r_diode < 0.020
    assert abs(float(vdmos_diode(-1.0, ST8L65_DIODE))) < 1e-5
    assert 0.9 < ST8L65_DIODE.v_eg < 1.3
    assert ST8L65_DIODE.i_sr == 0.0 and ST8L65_DIODE.k_n == 2.0


# ---------------------------------------------------------------------------
# LTspice cross-validation (skipped when LTspice is not installed)
# ---------------------------------------------------------------------------

def _ltspice_exe() -> str | None:
    for name in ("LTspice", "ltspice"):
        found = shutil.which(name)
        if found:
            return found
    candidates = [
        Path(r"C:\Program Files\ADI\LTspice\LTspice.exe"),
        Path.home() / r"AppData\Local\Programs\ADI\LTspice\LTspice.exe",
    ]
    for c in candidates:
        if c.exists():
            return str(c)
    return None


LTSPICE = _ltspice_exe()


def _run_ltspice_meas(cir_text: str, tmp_path: Path) -> dict[str, float]:
    cir = tmp_path / "xcheck.cir"
    cir.write_text(cir_text)
    subprocess.run([LTSPICE, "-b", "-w", str(cir)], cwd=tmp_path,
                   capture_output=True, timeout=120)
    log = (tmp_path / "xcheck.log").read_text()
    # LTspice .meas log lines look like:  m03: I(V1) =-3.44551587105 at 0.3
    out: dict[str, float] = {}
    for m in re.finditer(r"^(m\d+):\s*I\(V1\)\s*=\s*([-\d.e+]+)", log, re.MULTILINE):
        out[m.group(1)] = float(m.group(2))
    assert out, f"no .meas results parsed from {tmp_path / 'xcheck.log'}"
    return out


@pytest.mark.skipif(LTSPICE is None, reason="LTspice not installed")
def test_ltspice_cross_validation(tmp_path):
    """Diffusion + Ikf + Rs + temperature law match LTspice's own solver.

    Scope note: the Isr branch is EXCLUDED -- LTspice's implementation of its
    own documented recombination term deviates by ~1-17 % below 0.3 V
    (measured 2026-10; see diode.py docstring). Our data lives at 0.56-1.05 V
    where i_sr is unsensed anyway.
    """
    p = DEFAULT.replace(i_s=1e-6, k_n=1.5, i_kf=5.0, r_diode=1e-3,
                        k_xti=3.0, v_eg=1.11, t_nom=25.0)
    card = ("D(Is=1e-6 N=1.5 Ikf=5 Rs=1e-3 Xti=3 Eg=1.11 Tnom=25)")
    biases = (0.3, 0.6, 0.9)
    names = {0.3: "m03", 0.6: "m06", 0.9: "m09"}

    # 25 degC: diffusion + high-injection + Rs
    cir = (
        "* xcheck 25C\n"
        f".model bodydio {card}\n"
        "D1 a 0 bodydio\nV1 a 0 dc 0\n"
        ".dc V1 0 1 0.05\n"
        + "".join(f".meas dc {names[v]} FIND I(V1) AT={v}\n" for v in biases)
        + ".temp 25\n.end\n"
    )
    meas = _run_ltspice_meas(cir, tmp_path)
    for v in biases:
        lt = abs(meas[names[v]])
        ours = float(vdmos_diode(v, p, t_c=25.0))
        assert ours == pytest.approx(lt, rel=2e-3), f"25 C @ {v} V: {ours} vs {lt}"

    # -55 degC: pins the Xti/Eg temperature law (divisor N) against LTspice
    cir_lo = (
        "* xcheck -55C\n"
        f".model bodydio {card}\n"
        "D1 a 0 bodydio\nV1 a 0 dc 0\n"
        ".dc V1 0 1 0.05\n"
        + "".join(f".meas dc {names[v]} FIND I(V1) AT={v}\n" for v in biases)
        + ".temp -55\n.end\n"
    )
    meas_lo = _run_ltspice_meas(cir_lo, tmp_path)
    for v in biases:
        lt = abs(meas_lo[names[v]])
        ours = float(vdmos_diode(v, p, t_c=-55.0))
        assert ours == pytest.approx(lt, rel=2e-2), f"-55 C @ {v} V: {ours} vs {lt}"
