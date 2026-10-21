"""Charge-based capacitance model -- the settled design of CAP_SPEC.md.

Three-capacitor delta (G-S, G-D, D-S) per device, built from digitised
datasheet C-V curves measured at V_GS = 0 (off-state):

    c_gd = rss   (off-state table for v_gd < 0 + exponential-saturation
                  branch for v_gd >= 0, spliced C0/C1 at v_gd = 0)
    c_gs = iss - rss   (plain table)
    c_ds = oss - rss   (plain table)

Charge is master: q(v) = int C dv is the stored quantity, C = dq/dv by
construction (exact by the FTC -- no behavioural C(v)dv/dt terms). Tables are
PCHIP on (v, ln C): positive by construction, monotonicity-preserving through
the super-junction cliff, C1 in v. NOT log-log: the splice sits at v_gd = 0
and log 0 = -inf. Out-of-range policy: clamp C at the nearest node (q grows
linearly beyond).

The c_gd splice: the datasheet Crss(V_DS) at V_GS = 0 IS the v_gd < 0 branch
with v_ds = -v_gd. Near the origin the table branch uses the first-4-point
quadratic anchor C(v) = c_gd0 - s_gd0*v + b2*v^2 (v = |v_gd|), so C_gd(0) =
c_gd0 and dC_gd/dv_gd(0-) = s_gd0 EXACTLY; the exponential branch

    C_gd(v_gd) = c_gdmax - (c_gdmax - c_gd0) * exp(-k_agd * v_gd),
    k_agd = s_gd0 / (c_gdmax - c_gd0)

matches both for ANY c_gdmax > c_gd0 (C1 is automatic -- the tanh form of
LTspice would additionally require c_gdmax < 2*c_gd0, violated by SiC and
trench parts; see CAP_SPEC "Rejected alternatives"). c_gdmax is the single
tunable knob, extracted from the gate-charge curve by extract_cmax().

Naming (names.xlsx): c_ capacitance [F], q_ charge [C]; splice parameters
c_gd0, c_gdmax [F], k_agd [1/V] (unit exception like the k_lambda trio).
"""

from __future__ import annotations

import math
from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.interpolate import PchipInterpolator
from scipy.optimize import brentq

from .data import CapCurveSet, QgCurveSet

__all__ = [
    "CapTables",
    "OriginAnchor",
    "CapParams",
    "CapModel",
    "build_cap_model",
    "CmaxFit",
    "extract_cmax",
]

# 5-point Gauss-Legendre nodes/weights on [-1, 1]; exact for the per-segment
# cubics of ln C to far better than digitiser noise.
_GL_X, _GL_W = leggauss(5)

# Origin-anchor convention: quadratic fit to the FIRST N points (lowest v_ds)
# of the rss branch. The window is part of the convention (CAP_SPEC): the SJ
# slope is window-sensitive because the plateau->cliff inflection is real.
ANCHOR_N = 4


def _dedupe(v: np.ndarray, c: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Strictly increasing v; duplicate v nodes collapse by averaging c."""
    v, c = np.asarray(v, float), np.asarray(c, float)
    order = np.argsort(v)
    v, c = v[order], c[order]
    keep = np.concatenate(([True], np.diff(v) > 0))
    return v[keep], c[keep]


class _BranchTable:
    """One C(v) curve: PCHIP on (v, ln C) + cumulative charge at the nodes.

    q_at(v) integrates with the clamp policy (flat C outside the node range,
    so q grows linearly). All quantities SI (V, F, C).
    """

    def __init__(self, v: np.ndarray, c: np.ndarray) -> None:
        v, c = _dedupe(v, c)
        if len(v) < 2:
            raise ValueError("need at least 2 distinct nodes")
        if np.any(c <= 0):
            raise ValueError("capacitance nodes must be positive")
        self.v, self.c = v, c
        self._p = PchipInterpolator(v, np.log(c))
        seg = np.empty(len(v) - 1)
        for i in range(len(v) - 1):
            a, b = v[i], v[i + 1]
            x = 0.5 * (b - a) * _GL_X + 0.5 * (a + b)
            seg[i] = np.sum(0.5 * (b - a) * _GL_W * np.exp(self._p(x)))
        self._cumq = np.concatenate(([0.0], np.cumsum(seg)))

    def c_at(self, v: np.ndarray | float) -> np.ndarray:
        vq = np.clip(np.asarray(v, float), self.v[0], self.v[-1])
        return np.exp(self._p(vq))

    def q_at(self, v: np.ndarray | float) -> np.ndarray:
        """Integral of C from v[0] to v (clamp policy; negative below v[0])."""
        vq = np.asarray(v, float)
        scalar = vq.ndim == 0
        vq = np.atleast_1d(vq)
        out = np.empty(vq.shape)
        # inside range: node cumsum + partial-segment GL integral
        inside = (vq >= self.v[0]) & (vq <= self.v[-1])
        vi = vq[inside]
        i = np.clip(np.searchsorted(self.v, vi) - 1, 0, len(self.v) - 2)
        a = self.v[i]
        mid = 0.5 * (a + vi)
        half = 0.5 * (vi - a)
        x = mid[:, None] + half[:, None] * _GL_X
        part = np.sum(half[:, None] * _GL_W * np.exp(self._p(x)), axis=1)
        out[inside] = self._cumq[i] + part
        # outside: flat-C clamp -> linear growth
        lo = vq < self.v[0]
        out[lo] = self.c[0] * (vq[lo] - self.v[0])
        hi = vq > self.v[-1]
        out[hi] = self._cumq[-1] + self.c[-1] * (vq[hi] - self.v[-1])
        return out[0] if scalar else out

    def integral(self, a: float, b: float) -> float:
        return float(self.q_at(b) - self.q_at(a))


@dataclass(frozen=True)
class CapTables:
    """The three digitised branches as charge-ready tables (SI units)."""

    iss: _BranchTable
    oss: _BranchTable
    rss: _BranchTable
    label: str = ""

    @classmethod
    def from_curve_set(cls, data: CapCurveSet) -> CapTables:
        missing = {"iss", "oss", "rss"} - set(data.branches())
        if missing:
            raise ValueError(f"{data.label or 'cap data'} missing branch(es) {sorted(missing)}")
        return cls(
            iss=_BranchTable(*data.branch_xy("iss")),
            oss=_BranchTable(*data.branch_xy("oss")),
            rss=_BranchTable(*data.branch_xy("rss")),
            label=data.label,
        )


@dataclass(frozen=True)
class OriginAnchor:
    """Quadratic fit C_rss(v) = c_gd0 - s_gd0*v + b2*v^2 near the origin.

    Convention: first ANCHOR_N points of the rss branch (lowest v_ds).
    s_gd0 > 0 is dC_gd/dv_gd at the splice from the table side (the fit
    slope dC_rss/dv_ds is its negative). b2 is diagnostic only.
    """

    c_gd0: float
    s_gd0: float
    b2: float
    n_points: int = ANCHOR_N


def origin_anchor(rss_v: np.ndarray, rss_c: np.ndarray, n: int = ANCHOR_N) -> OriginAnchor:
    if n < 3:
        raise ValueError("anchor window needs >= 3 points")
    if len(rss_v) < n:
        raise ValueError(f"rss branch has only {len(rss_v)} points < window {n}")
    b2, s, c0 = np.polyfit(rss_v[:n], rss_c[:n], 2)  # C = c0 + s*v + b2*v^2
    if c0 <= 0 or s >= 0:
        raise ValueError(f"implausible origin anchor: c0={c0:g}, s={s:g} (want c0>0, s<0)")
    return OriginAnchor(c_gd0=float(c0), s_gd0=float(-s), b2=float(b2), n_points=n)


@dataclass(frozen=True)
class CapParams:
    """Splice parameters of c_gd(v_gd). c_gd0/s_gd0 come from the table
    anchor; c_gdmax is the single tunable (seed c_gs(0), bracket
    [c_gd0, 2*c_gs(0)], extracted from Qg -- CAP_SPEC "Cmax protocol").
    c_gdmax == c_gd0 is the flat branch (k_agd -> inf limit)."""

    c_gd0: float
    s_gd0: float
    c_gdmax: float

    def __post_init__(self) -> None:
        if self.c_gd0 <= 0 or self.s_gd0 <= 0:
            raise ValueError("c_gd0 and s_gd0 must be positive")
        if self.c_gdmax < self.c_gd0:
            raise ValueError("c_gdmax < c_gd0 (flat branch is c_gdmax == c_gd0)")

    @property
    def k_agd(self) -> float:
        if self.c_gdmax <= self.c_gd0:
            return math.inf
        return self.s_gd0 / (self.c_gdmax - self.c_gd0)


class CapModel:
    """c_gd(v_gd), c_gs(v_gs), c_ds(v_ds) and their charges (SI).

    q_gd(0) = 0 reference; q_gs/q_ds integrate from 0 with the clamp policy.
    The rss table branch of c_gd uses the origin quadratic below the first
    digitised node (exact splice slope) and PCHIP above; the tiny (digitiser-
    scale) value mismatch at the first node is documented in CAP_SPEC.
    """

    def __init__(self, tables: CapTables, anchor: OriginAnchor, params: CapParams) -> None:
        self.tables = tables
        self.anchor = anchor
        self.params = params
        # c_gs / c_ds on the union of branch nodes (difference of the two
        # PCHIP interpolants evaluated there, then PCHIP-lnC of the difference)
        self.c_gs_table = self._diff_table(tables.iss, tables.rss)
        self.c_ds_table = self._diff_table(tables.oss, tables.rss)
        if self.c_gs_table.c.min() <= 0 or self.c_ds_table.c.min() <= 0:
            raise ValueError("decomposition produced non-positive capacitance")

    @staticmethod
    def _diff_table(a: _BranchTable, b: _BranchTable) -> _BranchTable:
        nodes = np.unique(np.concatenate((a.v, b.v)))
        c = a.c_at(nodes) - b.c_at(nodes)
        return _BranchTable(nodes, c)

    # -- rss branch with origin-quadratic head -------------------------------

    def _rss_head(self, v: np.ndarray) -> np.ndarray:
        """C_rss(v) for v >= 0: quadratic below first node, PCHIP above."""
        v = np.atleast_1d(np.asarray(v, float))
        a = self.anchor
        v0 = self.tables.rss.v[0]
        out = np.empty_like(v)
        near = v < v0
        out[near] = a.c_gd0 - a.s_gd0 * v[near] + a.b2 * v[near] ** 2
        out[~near] = self.tables.rss.c_at(v[~near])
        return out

    def _q_rss_from0(self, v: np.ndarray | float) -> np.ndarray:
        """Integral of the head-extended C_rss from 0 to v (v >= 0 side)."""
        v = np.atleast_1d(np.asarray(v, float))
        a = self.anchor
        v0 = self.tables.rss.v[0]

        # analytic quadratic integral 0 -> min(v, v0)
        def quad(x: np.ndarray) -> np.ndarray:
            return a.c_gd0 * x - a.s_gd0 * x**2 / 2 + a.b2 * x**3 / 3

        q_head_at_v0 = float(quad(np.array([v0]))[0])
        t = self.tables.rss
        above = q_head_at_v0 + (t.q_at(np.clip(v, v0, None)) - t.q_at(v0))
        return np.where(v <= v0, quad(v), above)

    # -- c_gd / q_gd ----------------------------------------------------------

    def c_gd(self, v_gd: np.ndarray | float) -> np.ndarray:
        v = np.atleast_1d(np.asarray(v_gd, float))
        out = np.empty_like(v)
        neg = v < 0
        out[neg] = self._rss_head(-v[neg])
        p = self.params
        pos = ~neg
        if math.isinf(p.k_agd):
            out[pos] = p.c_gd0
        else:
            out[pos] = p.c_gdmax - (p.c_gdmax - p.c_gd0) * np.exp(-p.k_agd * v[pos])
        return out

    def q_gd(self, v_gd: np.ndarray | float) -> np.ndarray:
        """Charge on the G-D capacitor, q_gd(0) = 0."""
        v = np.atleast_1d(np.asarray(v_gd, float))
        out = np.empty_like(v)
        neg = v < 0
        out[neg] = -self._q_rss_from0(-v[neg])
        p = self.params
        pos = ~neg
        if math.isinf(p.k_agd):
            out[pos] = p.c_gd0 * v[pos]
        else:
            out[pos] = p.c_gdmax * v[pos] + (p.c_gdmax - p.c_gd0) * (
                np.exp(-p.k_agd * v[pos]) - 1
            ) / p.k_agd
        return out

    # -- c_gs / c_ds ----------------------------------------------------------

    def c_gs(self, v: np.ndarray | float) -> np.ndarray:
        return self.c_gs_table.c_at(v)

    def q_gs(self, v: np.ndarray | float) -> np.ndarray:
        return self.c_gs_table.q_at(v) - self.c_gs_table.q_at(0.0)

    def c_ds(self, v: np.ndarray | float) -> np.ndarray:
        return self.c_ds_table.c_at(v)

    def q_ds(self, v: np.ndarray | float) -> np.ndarray:
        return self.c_ds_table.q_at(v) - self.c_ds_table.q_at(0.0)

    def c_gs0(self) -> float:
        """c_gs at v_gs = 0 (clamp value at/below the first union node)."""
        return float(np.ravel(self.c_gs_table.c_at(0.0))[0])


def build_cap_model(
    data: CapCurveSet,
    *,
    c_gdmax: float | None = None,
    anchor_n: int = ANCHOR_N,
) -> CapModel:
    """Tables + origin anchor + params; c_gdmax defaults to the seed c_gs(0)."""
    tables = CapTables.from_curve_set(data)
    rss_v, rss_c = data.branch_xy("rss")
    rss_v, rss_c = _dedupe(rss_v, rss_c)
    anchor = origin_anchor(rss_v, rss_c, anchor_n)
    # seed needs c_gs(0): build the difference table first via a throwaway
    # model is circular; compute directly:
    gs = CapModel._diff_table(tables.iss, tables.rss)
    seed = float(gs.c_at(0.0).ravel()[0])
    params = CapParams(
        c_gd0=anchor.c_gd0,
        s_gd0=anchor.s_gd0,
        c_gdmax=float(c_gdmax) if c_gdmax is not None else seed,
    )
    return CapModel(tables, anchor, params)


# ---------------------------------------------------------------------------
# Cmax extraction from the gate-charge curve (CAP_SPEC "Cmax protocol")
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class CmaxFit:
    """Per-condition extraction record (all SI)."""

    i_d: float
    v_ds: float
    v_gp_lo: float          # plateau start voltage  (v_gd = V1 - V_DS)
    v_gp_hi: float          # plateau end voltage    (v_gd = 0 here)
    qgs_meas: float
    qgs_model: float
    qgd_meas: float         # plateau charge
    qgd_model: float        # table-branch prediction (NO c_gdmax term)
    plateau_deficit: float  # meas - model; > 0 = channel charge missing
    q_seg3_meas: float
    cgs_seg3_model: float
    remainder: float        # q_seg3_meas - cgs_seg3_model = q_exp target
    c_gdmax_raw: float | None   # None when not bracketed by the solve
    c_gdmax: float              # clamped into [c_gd0, 2*c_gs(0)]
    clamp: str | None           # "LOWER"/"UPPER"/None


def _q_exp_branch(model: CapModel, v: float, c_gdmax: float) -> float:
    """Charge of the exp branch 0..v with a GIVEN c_gdmax (C)."""
    a = model.anchor
    if c_gdmax <= a.c_gd0:
        return a.c_gd0 * v
    k = a.s_gd0 / (c_gdmax - a.c_gd0)
    return c_gdmax * v + (c_gdmax - a.c_gd0) * (math.exp(-k * v) - 1) / k


def extract_cmax(
    model: CapModel,
    qg: QgCurveSet,
    *,
    r_ds_on: float,
    conditions: Sequence[tuple[float, float]] | None = None,
) -> tuple[CmaxFit, ...]:
    """Extract c_gdmax per Qg test condition (CAP_SPEC, corrected split).

    The plateau ENDS at v_gd = 0 (v_ds falls to v_gs), so the plateau charge
    is a pure table-branch prediction -- it VALIDATES the tables but carries
    no c_gdmax information. c_gdmax is extracted from the post-plateau
    segment (v_gd > 0):

        remainder = (Q_end - Q2) - int_{V2}^{V_end} c_gs(v) dv
        remainder = q_exp(V_end - i_d*r_ds_on; c_gdmax)   -> solve

    r_ds_on: datasheet Rds(on) used for the ohmic tail voltage i_d*r_ds_on
    (the v_gd endpoint is V_end - i_d*r_ds_on, not V_end).
    """
    conds = list(conditions) if conditions is not None else qg.conditions()
    a = model.anchor
    c_gs0 = model.c_gs0()
    lo, hi = a.c_gd0, 2.0 * c_gs0
    if hi <= lo:
        raise ValueError("degenerate c_gdmax bracket")
    out: list[CmaxFit] = []
    for i_d, v_ds in conds:
        q, v = qg.condition_xy(i_d, v_ds)
        if q[0] > 0:  # prepend the physical origin corner (Q=0 at V=0)
            q, v = np.concatenate(([0.0], q)), np.concatenate(([0.0], v))
        if len(q) < 4:
            raise ValueError(f"condition (i_d={i_d}, v_ds={v_ds}) has < 3 segments")
        # plateau = max dQ/dV segment
        dq, dv = np.diff(q), np.diff(v)
        plat = int(np.argmax(dq / np.maximum(dv, 1e-12)))
        if plat >= len(dq) - 1:
            raise ValueError("plateau is the last segment -- malformed Qg curve")
        q1, v1 = q[plat], v[plat]          # plateau start
        q2, v2 = q[plat + 1], v[plat + 1]  # plateau end (v_gd = 0)
        q_end, v_end = q[-1], v[-1]

        qgs_model = float(model.q_gs(v1))
        # plateau: table branch over v_ds in [0, V_DS - V1] + c_gs over V1..V2
        qgd_model = float(
            model._q_rss_from0(np.array([v_ds - v1])).ravel()[0]
            + (model.q_gs(v2) - model.q_gs(v1))
        )
        cgs_seg3 = float(model.q_gs(v_end) - model.q_gs(v2))
        remainder = (q_end - q2) - cgs_seg3
        arg = v_end - i_d * r_ds_on
        if arg <= 0:
            raise ValueError("V_end below the ohmic tail voltage -- check r_ds_on")

        f = lambda x, _arg=arg, _rem=remainder: _q_exp_branch(model, _arg, x) - _rem  # noqa: E731
        f_lo, f_hi = f(lo), f(hi)
        if f_lo >= 0:      # even the flat branch overshoots -> clamp lower
            raw, c_gdmax, clamp = None, lo, "LOWER"
        elif f_hi <= 0:    # cannot reach the measured charge -> clamp upper
            raw, c_gdmax, clamp = None, hi, "UPPER"
        else:
            raw = brentq(f, lo, hi, xtol=1e-18, rtol=1e-12)
            c_gdmax, clamp = raw, None
        out.append(
            CmaxFit(
                i_d=i_d, v_ds=v_ds, v_gp_lo=v1, v_gp_hi=v2,
                qgs_meas=float(q1), qgs_model=qgs_model,
                qgd_meas=float(q2 - q1), qgd_model=qgd_model,
                plateau_deficit=float((q2 - q1) - qgd_model),
                q_seg3_meas=float(q_end - q2), cgs_seg3_model=cgs_seg3,
                remainder=remainder, c_gdmax_raw=raw, c_gdmax=c_gdmax,
                clamp=clamp,
            )
        )
    return tuple(out)
