"""LTspice parity bench for the non-collapsed channel expression (the gap the
export docstring flagged). We translate the emitted B-source expression into
NumPy and check it reproduces model.py _channel_current exactly, for the
families the collapsed SJ export never exercised: the DIBL-shift soft cap,
the linear knee-scale form (k_kneeform=1), and the negative-k_theta mobility
widening with its divisor guard. This is what catches the exported-subckt
"saturation pinned at 950 A" class of bug (missing DIBL cap -> v_th runs away
negative -> overdrive saturates at the charge cap for every gate)."""

from __future__ import annotations

import re

import numpy as np
import pytest

from vdmos_fit.export_subckt import (
    _SGN,
    channel_expr,
    channel_expr_auto,
    channel_expr_full,
    device_common,
)
from vdmos_fit.model import _channel_current
from vdmos_fit.params import Params
from vdmos_fit.stages import params_at


def _eval_expr(expr: str, v_ds, v_gs):
    """Evaluate the emitted SPICE B-source expression in NumPy. Valid for
    v_ds >= 0 (the mirror sign factor _SGN is 1 there, so replace it)."""
    s = expr.replace(_SGN, "1")
    s = s.replace("V(dx,sx)", "(_vd)").replace("V(G,sx)", "(_vg)")
    for fn, npfn in (("max", "np.maximum"), ("tanh", "np.tanh"),
                     ("log", "np.log"), ("exp", "np.exp"),
                     ("abs", "np.abs"), ("pow", "np.power")):
        s = re.sub(rf"\b{fn}\(", npfn + "(", s)
    return eval(s, {"np": np, "_vd": v_ds, "_vg": v_gs})  # noqa: S307


def _parity(p: Params):
    v_ds = np.linspace(0.0, 100.0, 51)
    v_gs = np.linspace(0.0, 15.0, 31)
    vd, vg = np.meshgrid(v_ds, v_gs)
    expr = channel_expr_auto(p)
    got = _eval_expr(expr, vd, vg)
    ref = _channel_current(vd, vg, p)
    assert np.all(np.isfinite(got))
    atol = 1e-9 * float(np.max(np.abs(ref)))
    assert np.allclose(got, ref, rtol=1e-6, atol=atol)


def test_parity_dibl_cap_scale_knee():
    """The user's velocity-saturation + DIBL-cap + linear knee-scale fit."""
    _parity(Params(k_p=280.0, v_th0=4.45, v_sub=0.244, v_ovmax=7.9,
                   k_ovsharp=3.0, k_beta_trap=0.0, k_theta=0.95,
                   k_kneeform=1.0, v_vsat=11.3, k_kneesharp=1.63, k_eta=0.11,
                   v_diblmax=7.0, k_diblsharp=0.5, k_lambda=0.0, g_k=0.0,
                   v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.31, v_kneemin=1e-3,
                   r_d=0.0, r_s=0.0, r_sk=0.0))


def test_parity_dibl_cap_negative_ktheta_mobility():
    """Mobility knee with NEGATIVE k_theta (wider knee): the divisor
    1 + k_theta*v_ov_eff2 crosses zero on the grid, so both model and export
    must take the floored branch identically."""
    _parity(Params(k_p=260.0, v_th0=4.45, v_sub=0.244, v_ovmax=8.2,
                   k_ovsharp=3.3, k_beta_trap=0.0, k_theta=-0.2,
                   k_kneeform=0.0, v_vsat=2.3, k_kneesharp=1.83, k_eta=0.11,
                   v_diblmax=7.0, k_diblsharp=0.5, k_lambda=0.0, g_k=0.0,
                   v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.25, v_kneemin=1e-3,
                   r_d=0.0, r_s=0.0, r_sk=0.0))


def test_parity_full_family_clm_kink_trap():
    """CLM, trench kink and trap-filling gain together (mobility knee, cap
    on) -- the full transcription path."""
    _parity(Params(k_p=240.0, v_th0=4.3, v_sub=0.24, v_ovmax=6.0,
                   k_ovsharp=4.0, k_beta_trap=0.05, k_theta=0.02,
                   k_kneeform=0.0, v_vsat=40.0, k_kneesharp=2.0, k_eta=0.05,
                   v_diblmax=9.0, k_diblsharp=3.0, k_lambda=0.004, g_k=0.5,
                   v_thk=1.6, v_kovsoft=0.1, v_dbsat=0.6, v_kneemin=1e-3,
                   r_d=0.0, r_s=0.0, r_sk=0.0))


def test_auto_routes_cap_and_kneeform_to_full():
    """channel_expr_auto must NOT collapse a fit that has the DIBL cap or the
    scale knee on, even with k_theta/k_beta_trap/k_lambda/g_k all zero --
    otherwise the collapsed expression (no cap) silently drops the cap."""
    base = dict(k_theta=0.0, k_beta_trap=0.0, k_lambda=0.0, g_k=0.0)
    cap = Params(**{**base, "v_diblmax": 7.0})
    assert channel_expr_auto(cap) == channel_expr_full(cap)  # routed to full
    knee = Params(**{**base, "k_kneeform": 1.0, "k_theta": 0.9})
    assert channel_expr_auto(knee) == channel_expr_full(knee)
    plain = Params(**base)                                  # all-off collapsed
    assert channel_expr_auto(plain) == channel_expr(plain)


# ---- temperature coefficients (the `.temp` gap: subckt was 25 C-only) ------

_TC = {"s_vth": -0.00748634, "beta_p": 2.2004, "a_rd": 0.0056206,
       "b_rd": 2.00079e-05}
_PBASE = dict(k_p=222.59, v_th0=4.28667, v_sub=0.244088, v_ovmax=5.77459,
              k_ovsharp=2.19649, k_kneesharp=1.78701, k_eta=0.0313065,
              v_dbsat=0.52529, r_d=0.00248295, v_vsat=10.0, v_kneemin=1e-3)


def test_tc_off_is_byte_identical():
    """With no TC coefficients the export must be byte-identical to the
    25 C-only pinned files (no .param block, literal k_p/v_th0/r_d)."""
    p = Params(**_PBASE)
    zero = {n: 0.0 for n in _TC}
    assert channel_expr_auto(p) == channel_expr_auto(p, None)
    assert channel_expr_auto(p) == channel_expr_auto(p, zero)
    assert device_common(p, "caps") == device_common(p, "caps", None)
    assert device_common(p, "caps") == device_common(p, "caps", zero)
    assert ".param kp" not in device_common(p, "caps", zero)


def test_tc_param_formulas_match_model():
    """The emitted .param kp/vth0t/rdt expressions (referencing LTspice TEMP)
    must equal the model's TC law (stages.params_at) at an arbitrary temp."""
    p = Params(**_PBASE)
    T = 100.0
    pt = params_at(p, T, _TC)
    vals = dict(re.findall(r"\.param\s+(\w+)=\{([^}]*)\}",
                           device_common(p, "caps", _TC)))
    assert set(vals) >= {"kp", "vth0t", "rdt"}

    def ev(e):
        e = re.sub(r"\bpow\(", "np.power(", e)
        return eval(e, {"np": np, "TEMP": T})  # noqa: S307

    assert ev(vals["kp"]) == pytest.approx(pt.k_p, rel=1e-9)
    assert ev(vals["vth0t"]) == pytest.approx(pt.v_th0, rel=1e-9)
    assert ev(vals["rdt"]) == pytest.approx(pt.r_d, rel=1e-9)


def test_tc_channel_expr_wiring_parity():
    """channel_expr_auto(p, tc) places kp in the gain slot and vth0t in the
    threshold slot: substituting the TC-applied k_p(T)/v_th0(T) reproduces
    _channel_current at that temperature."""
    p = Params(**{**_PBASE, "k_theta": 0.95, "k_kneeform": 1.0,
                  "v_diblmax": 7.0, "k_diblsharp": 0.5})
    T = 100.0
    pt = params_at(p, T, _TC)
    expr = channel_expr_auto(p, _TC)
    assert "kp" in expr and "vth0t" in expr
    s = expr.replace(_SGN, "1")
    s = s.replace("V(dx,sx)", "(_vd)").replace("V(G,sx)", "(_vg)")
    s = re.sub(r"\bvth0t\b", "(_vth0t)", s)
    s = re.sub(r"\bkp\b", "(_kp)", s)
    for fn, npfn in (("max", "np.maximum"), ("tanh", "np.tanh"),
                     ("log", "np.log"), ("exp", "np.exp"),
                     ("abs", "np.abs"), ("pow", "np.power")):
        s = re.sub(rf"\b{fn}\(", npfn + "(", s)
    v_ds = np.linspace(0.0, 100.0, 41)
    v_gs = np.linspace(0.0, 15.0, 21)
    vd, vg = np.meshgrid(v_ds, v_gs)
    got = eval(s, {"np": np, "_vd": vd, "_vg": vg,  # noqa: S307
                   "_kp": pt.k_p, "_vth0t": pt.v_th0})
    ref = _channel_current(vd, vg, pt)
    atol = 1e-9 * float(np.max(np.abs(ref)))
    assert np.allclose(got, ref, rtol=1e-6, atol=atol)
