"""Measured / digitised I-V datasets and synthetic test data."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .model import vdmos_id
from .params import Params

__all__ = [
    "CurveSet",
    "DiodeCurveSet",
    "CapCurveSet",
    "QgCurveSet",
    "load_csv",
    "load_diode_csv",
    "load_cap_csv",
    "load_qg_csv",
    "save_csv",
    "save_diode_csv",
    "make_synthetic",
    "make_synthetic_diode",
    "example_dataset",
]

# On-disk column names are an EXTERNAL contract (the curve digitizer exports
# them and the digitised datasheet CSVs already use them), so they keep the
# old SPICE-style spelling even though the in-code attributes follow the
# model naming convention. See README, "Naming convention".
_CSV_COLUMNS: tuple[str, ...] = ("v_gs", "v_ds", "i_d")
_CSV_HEADERS: dict[str, str] = {"v_gs": "vgs", "v_ds": "vds", "i_d": "ids"}


@dataclass(frozen=True)
class CurveSet:
    """Flat arrays of (v_gs, v_ds, i_d) triplets -- one row per measured point.

    Flat form (rather than 2-D grids) is deliberate: datasheet curves digitised
    from graphs are ragged, and different curves can have different lengths.
    """

    v_gs: np.ndarray
    v_ds: np.ndarray
    i_d: np.ndarray
    label: str = ""

    def __post_init__(self) -> None:
        n = len(self.v_gs)
        if not (len(self.v_ds) == len(self.i_d) == n):
            raise ValueError("v_gs, v_ds and i_d must have the same length")

    def __len__(self) -> int:
        return len(self.v_gs)

    def subset(self, mask: np.ndarray) -> CurveSet:
        return CurveSet(self.v_gs[mask], self.v_ds[mask], self.i_d[mask], self.label)

    def v_gs_grid(self, tol: float = 1e-9) -> list[float]:
        return sorted({round(float(v), 9) for v in self.v_gs})

    def v_ds_grid(self, tol: float = 1e-9) -> list[float]:
        return sorted({round(float(v), 9) for v in self.v_ds})


def load_csv(path: str | Path) -> CurveSet:
    """Read a CSV with header `vgs,vds,ids` (column order flexible, case-insensitive).

    The header spellings are the legacy SPICE-style names -- see
    :data:`_CSV_COLUMNS`.
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    cols = {k.strip().lower(): k for k in rows[0] if k}
    for need in _CSV_HEADERS.values():
        if need not in cols:
            raise ValueError(f"{path} missing column '{need}'; header is {list(rows[0])}")
    get = lambda r, c: float(r[cols[c]])  # noqa: E731
    return CurveSet(
        np.array([get(r, "vgs") for r in rows]),
        np.array([get(r, "vds") for r in rows]),
        np.array([get(r, "ids") for r in rows]),
        label=path.stem,
    )


def save_csv(data: CurveSet, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([_CSV_HEADERS[c] for c in _CSV_COLUMNS])
        w.writerows(zip(data.v_gs, data.v_ds, data.i_d, strict=True))


@dataclass(frozen=True)
class DiodeCurveSet:
    """Flat arrays of (i_sd, v_sd, temp_c) triplets -- one row per measured point.

    Body-diode forward data, source-referenced: ``i_sd`` [A] flows source ->
    drain (anode = source), ``v_sd`` [V] is the source-drain drop (positive
    forward), ``temp_c`` [degC] the measurement temperature. The model's
    drain-referenced sign convention is ``v_ds = -v_sd``; the conversion
    lives in the model/fitter, not here -- the CSV keeps the datasheet's
    own spelling.
    """

    i_sd: np.ndarray
    v_sd: np.ndarray
    temp_c: np.ndarray
    label: str = ""

    def __post_init__(self) -> None:
        n = len(self.i_sd)
        if not (len(self.v_sd) == len(self.temp_c) == n):
            raise ValueError("i_sd, v_sd and temp_c must have the same length")

    def __len__(self) -> int:
        return len(self.i_sd)

    def subset(self, mask: np.ndarray) -> DiodeCurveSet:
        return DiodeCurveSet(self.i_sd[mask], self.v_sd[mask], self.temp_c[mask], self.label)

    def temp_grid(self) -> list[float]:
        return sorted({round(float(t), 9) for t in self.temp_c})


def load_diode_csv(path: str | Path) -> DiodeCurveSet:
    """Read a diode CSV with header `I_sd,v_sd,T` (an `index` column is ignored).

    Column match is case-insensitive by header, like :func:`load_csv`.
    `T` is degrees Celsius (datasheet convention).
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    cols = {k.strip().lower(): k for k in rows[0] if k}
    need = {"i_sd": "I_sd", "v_sd": "v_sd", "t": "T"}
    missing = [disp for key, disp in need.items() if key not in cols]
    if missing:
        raise ValueError(f"{path} missing column(s) {missing}; header is {list(rows[0])}")
    get = lambda r, c: float(r[cols[c]])  # noqa: E731
    return DiodeCurveSet(
        np.array([get(r, "i_sd") for r in rows]),
        np.array([get(r, "v_sd") for r in rows]),
        np.array([get(r, "t") for r in rows]),
        label=path.stem,
    )


def save_diode_csv(data: DiodeCurveSet, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["index", "I_sd", "v_sd", "T"])
        w.writerows(
            (i + 1, a, b, t)
            for i, (a, b, t) in enumerate(
                zip(data.i_sd, data.v_sd, data.temp_c, strict=True)
            )
        )


_CAP_BRANCHES: frozenset[str] = frozenset({"iss", "oss", "rss"})


@dataclass(frozen=True)
class CapCurveSet:
    """Flat arrays of (v_ds, c, branch) -- digitised capacitance curves.

    One row per digitised point; ``branch`` is the datasheet curve name
    ("iss", "oss", "rss"), lower-case. Capacitance ``c`` is stored in FARADS;
    the on-disk header spells the column `Cap (pF)` (datasheet units) and
    :func:`load_cap_csv` converts. The delta decomposition contract
    (c_gd = rss, c_gs = iss - rss, c_ds = oss - rss) lives in
    :mod:`vdmos_fit.caps`, not here -- this class just carries the curves.
    See CAP_SPEC.md.
    """

    v_ds: np.ndarray
    c: np.ndarray
    branch: np.ndarray  # string array, values in _CAP_BRANCHES
    label: str = ""

    def __post_init__(self) -> None:
        n = len(self.v_ds)
        if not (len(self.c) == len(self.branch) == n):
            raise ValueError("v_ds, c and branch must have the same length")
        bad = set(self.branch) - _CAP_BRANCHES
        if bad:
            raise ValueError(f"unknown capacitance branch(es) {sorted(bad)}")

    def __len__(self) -> int:
        return len(self.v_ds)

    def branches(self) -> list[str]:
        return sorted({str(b) for b in self.branch})

    def branch_mask(self, branch: str) -> np.ndarray:
        return self.branch == branch

    def branch_xy(self, branch: str) -> tuple[np.ndarray, np.ndarray]:
        """(v_ds, c) of one branch, sorted by v_ds."""
        m = self.branch_mask(branch)
        v, c = self.v_ds[m], self.c[m]
        order = np.argsort(v)
        return v[order], c[order]


def load_cap_csv(path: str | Path) -> CapCurveSet:
    """Read a capacitance CSV `index,V_ds,Cap (pF),Type` (Type in iss/oss/rss).

    Column match is case-insensitive by header, like :func:`load_csv`; the
    `index` column is ignored. `Cap (pF)` is datasheet picofarad and is
    converted to farads.
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    cols = {k.strip().lower(): k for k in rows[0] if k}
    need = {"v_ds": "V_ds", "cap (pf)": "Cap (pF)", "type": "Type"}
    missing = [disp for key, disp in need.items() if key not in cols]
    if missing:
        raise ValueError(f"{path} missing column(s) {missing}; header is {list(rows[0])}")
    get = lambda r, c: float(r[cols[c]])  # noqa: E731
    return CapCurveSet(
        np.array([get(r, "v_ds") for r in rows]),
        np.array([get(r, "cap (pf)") for r in rows]) * 1e-12,
        np.array([str(r[cols["type"]]).strip().lower() for r in rows]),
        label=path.stem,
    )


def save_cap_csv(data: CapCurveSet, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["index", "V_ds", "Cap (pF)", "Type"])
        w.writerows(
            (i + 1, v, c * 1e12, b)
            for i, (v, c, b) in enumerate(
                zip(data.v_ds, data.c, data.branch, strict=True)
            )
        )


@dataclass(frozen=True)
class QgCurveSet:
    """Flat arrays of (q_g, v_gs, i_d, v_ds) -- digitised gate-charge corners.

    Datasheet gate-charge specs are given as V_GS vs Q_g at a fixed test
    condition (I_d, V_DS); digitised curves are polyline corner points, so
    each row is a corner and consecutive corners of one condition define a
    segment (segment 1 = input charging, the max-dQ/dV segment = the Miller
    plateau, the last = post-plateau). ``q_g`` is stored in COULOMBS (the
    on-disk header spells `Q_g(nC)`); ``i_d``/``v_ds`` are parsed per row
    from the `condition` column (e.g. "I_d=29A, V_ds=400V").
    """

    q_g: np.ndarray
    v_gs: np.ndarray
    i_d: np.ndarray
    v_ds: np.ndarray
    label: str = ""

    def __post_init__(self) -> None:
        n = len(self.q_g)
        if not (len(self.v_gs) == len(self.i_d) == len(self.v_ds) == n):
            raise ValueError("q_g, v_gs, i_d and v_ds must have the same length")

    def __len__(self) -> int:
        return len(self.q_g)

    def conditions(self) -> list[tuple[float, float]]:
        """Distinct (i_d, v_ds) test conditions, sorted by v_ds."""
        return sorted({(float(a), float(b)) for a, b in zip(self.i_d, self.v_ds, strict=True)})

    def condition_mask(self, i_d: float, v_ds: float) -> np.ndarray:
        return (np.isclose(self.i_d, i_d)) & (np.isclose(self.v_ds, v_ds))

    def condition_xy(self, i_d: float, v_ds: float) -> tuple[np.ndarray, np.ndarray]:
        """(q_g, v_gs) corners of one test condition, sorted by q_g."""
        m = self.condition_mask(i_d, v_ds)
        q, v = self.q_g[m], self.v_gs[m]
        order = np.argsort(q)
        return q[order], v[order]


_QG_COND_KEYS = ("i_d", "v_ds")


def _parse_qg_condition(text: str) -> tuple[float, float]:
    """Parse 'I_d=29A, V_ds=400V' -> (29.0, 400.0). Tolerant of spacing."""
    out: dict[str, float] = {}
    for part in text.split(","):
        if "=" not in part:
            continue
        key, val = part.split("=", 1)
        key = key.strip().lower()
        val = val.strip().rstrip("aAvV ")
        if key in _QG_COND_KEYS:
            out[key] = float(val)
    if set(out) != set(_QG_COND_KEYS):
        raise ValueError(f"condition {text!r} missing i_d/v_ds fields")
    return out["i_d"], out["v_ds"]


def load_qg_csv(path: str | Path) -> QgCurveSet:
    """Read a gate-charge CSV `index,Q_g(nC),V_gs,condition`.

    `Q_g(nC)` is datasheet nanocoulomb and is converted to coulombs;
    `condition` carries the test bias ("I_d=..A, V_ds=..V") and is parsed
    per row so one file can hold several test curves. Header match ignores
    case and spaces (the digitizer writes `Q_g(nC)`).
    """
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    cols = {k.strip().lower().replace(" ", ""): k for k in rows[0] if k}
    need = {"q_g(nc)": "Q_g(nC)", "v_gs": "V_gs", "condition": "condition"}
    missing = [disp for key, disp in need.items() if key not in cols]
    if missing:
        raise ValueError(f"{path} missing column(s) {missing}; header is {list(rows[0])}")
    q, v, i_d, v_ds = [], [], [], []
    for r in rows:
        q.append(float(r[cols["q_g(nc)"]]) * 1e-9)
        v.append(float(r[cols["v_gs"]]))
        a, b = _parse_qg_condition(r[cols["condition"]])
        i_d.append(a)
        v_ds.append(b)
    return QgCurveSet(np.array(q), np.array(v), np.array(i_d), np.array(v_ds), label=path.stem)


def save_qg_csv(data: QgCurveSet, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["index", "Q_g(nC)", "V_gs", "condition"])
        w.writerows(
            (i + 1, q * 1e9, v, f"I_d={a:g}A, V_ds={b:g}V")
            for i, (q, v, a, b) in enumerate(
                zip(data.q_g, data.v_gs, data.i_d, data.v_ds, strict=True)
            )
        )


def make_synthetic(
    truth: Params,
    v_gs_vals: np.ndarray,
    v_ds_vals: np.ndarray,
    *,
    noise: float = 0.01,
    seed: int = 0,
) -> CurveSet:
    """Full grid of (v_gs, v_ds) points from `truth`, with multiplicative noise.

    Noise is multiplicative and applied *without* clipping: near-zero currents
    are handled in the residual denominator (`i_d_floor` in
    :func:`vdmos_fit.fit.residuals`), not in the data -- rationale in README,
    "Synthetic data". Used to verify the fitter recovers known parameters.
    """
    vg, vd = np.meshgrid(
        np.asarray(v_gs_vals, float), np.asarray(v_ds_vals, float), indexing="ij"
    )
    vg, vd = vg.ravel(), vd.ravel()
    i_d = vdmos_id(vd, vg, truth)
    rng = np.random.default_rng(seed)
    i_d = i_d * (1.0 + noise * rng.standard_normal(i_d.size))
    return CurveSet(vg, vd, i_d, label="synthetic")


def example_dataset(noise: float = 0.01, seed: int = 1) -> CurveSet:
    """Typical power-FET static grid: transfer curves at several v_ds plus output curves."""
    v_gs = np.linspace(2.0, 14.0, 25)
    v_ds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    truth = Params(
        k_p=30.0,
        v_th0=3.1,
        k_theta=0.02,
        v_vsat=1.2,
        k_eta=0.03,
        k_lambda=0.015,
        g_k=0.3,
        v_thk=1.5,
    )
    return make_synthetic(truth, v_gs, v_ds, noise=noise, seed=seed)


def make_synthetic_diode(
    truth: Params,
    i_sd_vals: np.ndarray,
    temps_c: np.ndarray,
    *,
    noise: float = 0.0,
    seed: int = 0,
) -> DiodeCurveSet:
    """Diode forward points: solve v_sd(i_sd) from `truth` at each temperature.

    The datasheet diode spec is current-referenced (V at given I), so the
    generator inverts the model per point (monotone -> unique root) rather
    than sweeping voltage. Used to verify the diode fitter recovers known
    parameters; noise is multiplicative on v_sd, matching the datasheet
    read-off error mode.
    """
    from scipy.optimize import brentq

    from .diode import vdmos_diode

    rng = np.random.default_rng(seed)
    i_sd, v_sd, temps = [], [], []

    def v_root(tc: float, i_ref: float) -> float:
        def f(v: float) -> float:
            return float(vdmos_diode(v, truth, t_c=tc)) - i_ref

        return brentq(f, -0.5, 2.5, xtol=1e-14, rtol=8.9e-16)

    for t_c in np.atleast_1d(temps_c):
        for i_ref in np.atleast_1d(i_sd_vals):
            v = v_root(float(t_c), float(i_ref))
            if noise:
                v *= 1.0 + noise * rng.standard_normal()
            i_sd.append(float(i_ref))
            v_sd.append(v)
            temps.append(float(t_c))
    return DiodeCurveSet(
        np.array(i_sd), np.array(v_sd), np.array(temps), label="synthetic_diode"
    )
