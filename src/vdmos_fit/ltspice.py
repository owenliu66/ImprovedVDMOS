"""LTspice export of the capacitance delta (CAP_SPEC validation step 5).

The charge-master model exports through LTspice's general nonlinear
capacitor -- the one primitive that IS a charge master (help, "C.
Capacitor": "instead of specifying the capacitance, one writes an
expression for the charge... LTspice symbolically differentiates it"):

    .subckt <name> G D S
    Cgs G S Q=<q_gs(x)>      ; x = voltage across the device (n1 - n2)
    Cgd G D Q=<q_gd(x)>
    Cds D S Q=<q_ds(x)>
    .ends

Each q_k(x) is a SUM OF RECTIFIED RAMPS over the sampled charge table:

    q(x) = q_0 + sum_k  m_k * max(0, x - x_k),   m_k = dq/dv on segment k

LTspice's symbolic derivative is then exact:

    dq/dx = sum_k m_k * u(x - x_k)   ->  piecewise-CONSTANT C staircase

no numerical d() of a B source, no discontinuous-source chatter. The
staircase error against the PCHIP model C is bounded by the sample
density (``n_sub`` subdivisions per digitised node interval; default 12
keeps it under ~2 % of the local C variation, well below digitiser
noise).

Clamp policy comes for free: outside the sampled span the ramp sum
continues with the END-SEGMENT slope m -- i.e. C stays flat at the last
digitised node and q keeps growing linearly, exactly the CapModel
policy, for all x. (A ``table()`` inside Q= compiles in LTspice 26.0.1
but its symbolic derivative is broken -- measured 2026-10-21, garbage
currents ~1e11 A -- hence the ramp-sum form. A ``table()`` in a B
source is rejected outright: "Function not supported in this context".)

Sign convention (verified against LTspice 26.0.1): Q= is the charge on
the FIRST node's plate; the current INTO the first terminal is dq/dt;
driving the delta with a voltage source Vxx gives I(Vxx) = -dq/dt.

Parity is checked by driving each branch pair with PWL ramps and
comparing I through the drive sources against the Python model
(tests/test_caps.py, skipped when LTspice is not installed).
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np

from .caps import CapModel

__all__ = ["cap_charge_tables", "export_cap_subckt", "write_cap_subckt"]

# Default sample density: subdivisions of each digitised node interval.
# 12 keeps the piecewise-constant-C staircase error under ~2 % of the local
# C variation between nodes (half-variation bound at segment midpoints).
N_SUB = 12


def _refine(nodes: np.ndarray, n_sub: int) -> np.ndarray:
    """Native nodes + n_sub subdivisions of every interval (sorted, unique)."""
    out = [nodes[0]]
    for a, b in zip(nodes[:-1], nodes[1:], strict=True):
        out.extend(np.linspace(a, b, n_sub + 1)[1:])
    return np.unique(np.array(out))


def _grid_gd(model: CapModel, n_sub: int) -> np.ndarray:
    """v_gd sample grid: quadratic head + rss nodes + exp branch."""
    rss_v = model.tables.rss.v
    neg = -_refine(np.concatenate(([0.0], rss_v)), n_sub)[::-1]  # v_gd <= 0 side
    p = model.params
    tau = 1.0 / p.k_agd if np.isfinite(p.k_agd) else rss_v[-1]
    v_hi = min(max(3.0 * tau, 25.0), 0.5 * rss_v[-1])
    pos = np.linspace(0.0, v_hi, 150)[1:]
    return np.unique(np.concatenate((neg, pos)))


def _grid_table(model: CapModel, branch: str, n_sub: int) -> np.ndarray:
    table = model.c_gs_table if branch == "gs" else model.c_ds_table
    return _refine(table.v, n_sub)


def cap_charge_tables(
    model: CapModel, *, n_sub: int = N_SUB
) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """(v, q) sample pairs per branch, sampled from the model's own q().

    q_gd(0) = 0, q_gs(0) = q_ds(0) = 0 -- the same references as CapModel.
    The end segments' slopes carry the clamp-C policy beyond the sampled
    span (the ramp-sum export extrapolates them linearly, matching the
    model), so no extension points are needed.
    """
    vg = _grid_gd(model, n_sub)
    gs = _grid_table(model, "gs", n_sub)
    ds = _grid_table(model, "ds", n_sub)
    qs = {
        "gd": (vg, np.ravel(model.q_gd(vg))),
        "gs": (gs, np.ravel(model.q_gs(gs))),
        "ds": (ds, np.ravel(model.q_ds(ds))),
    }
    for name, (v, q) in qs.items():
        if np.any(np.diff(v) <= 0) or np.any(np.diff(q) <= 0):
            raise ValueError(f"{name}: charge table not strictly increasing")
    return qs


def _rampsum_expr(v: np.ndarray, q: np.ndarray) -> str:
    """Incremental ramp-sum form of the piecewise-linear q(x):

        q(x) = q_0 + m_1*(x - x_0) + sum_k (m_{k+1} - m_k)*max(0, x - x_k)

    with m_k the charge-table segment slopes. LTspice's symbolic dQ/dx is
    then the staircase C(x) = m_1 + sum_k dm_k*u(x-x_k): below x_0 C stays
    at m_1 and above the last knot at m_last -- the clamp-C policy, free.
    (Absolute-slope ramps would sum slopes instead of stepping them --
    caught by LTspice parity, 2026-10-21.)
    """
    slopes = np.diff(q) / np.diff(v)
    x0 = v[0]
    shift0 = f"x+{-x0:.10g}" if x0 < 0 else f"x-{x0:.10g}"
    expr = f"{q[0]:.10g}+{slopes[0]:.10g}*({shift0})"
    for xk, dm in zip(v[1:-1], np.diff(slopes), strict=True):
        if dm == 0.0:
            continue
        shift = f"x+{-xk:.10g}" if xk < 0 else f"x-{xk:.10g}"
        expr += f"+{dm:+.10g}*max(0,{shift})"
    return expr


def _wrap_netlist(prefix: str, expr: str, width: int = 90) -> str:
    """One netlist statement, split into leading-'+' continuation lines.

    Splits only at '+' operators of the ramp sum, so every continuation
    rejoins with an explicit operator regardless of how LTspice joins
    lines.
    """
    tokens = re.findall(r"[^+]*\+?", prefix + expr)
    chunks: list[str] = []
    cur = ""
    for tok in tokens:
        if len(cur) + len(tok) > width and cur:
            chunks.append(cur)
            cur = tok
        else:
            cur += tok
    if cur:
        chunks.append(cur)
    return "\n".join(c if i == 0 else "+ " + c for i, c in enumerate(chunks)) + "\n"


def export_cap_subckt(
    model: CapModel, name: str = "caps", *, n_sub: int = N_SUB
) -> str:
    """The three-capacitor delta as an LTspice .subckt (pins G D S)."""
    qs = cap_charge_tables(model, n_sub=n_sub)
    p = model.params
    header = (
        f"* ImprovedVDMOS capacitance delta -- generated by vdmos_fit.ltspice\n"
        f"* source: {model.tables.label or 'CapModel'}  "
        f"c_gd0={p.c_gd0:.6g} c_gdmax={p.c_gdmax:.6g} k_agd={p.k_agd:.6g}\n"
        f"* charge-master nonlinear caps (Q= ramp-sum of sampled q(v);\n"
        f"* LTspice's symbolic dQ/dx gives a piecewise-CONSTANT C staircase,\n"
        f"* {n_sub} samples/digitised interval; end-segment slope = clamp-C policy)\n"
        f".subckt {name} G D S\n"
    )
    body = (
        _wrap_netlist("Cgs G S Q=", _rampsum_expr(*qs["gs"]))
        + _wrap_netlist("Cgd G D Q=", _rampsum_expr(*qs["gd"]))
        + _wrap_netlist("Cds D S Q=", _rampsum_expr(*qs["ds"]))
    )
    return header + body + f".ends {name}\n"


def write_cap_subckt(
    model: CapModel, path: str | Path, name: str = "caps", *, n_sub: int = N_SUB
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(export_cap_subckt(model, name, n_sub=n_sub), encoding="utf-8")
    return path
