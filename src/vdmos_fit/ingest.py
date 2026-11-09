"""Input contract for the release pipeline (session 30, Phase A).

One canonical grammar turns any datasheet CSV into a typed table:

* Column headers are ``NAME`` optionally followed by a unit in
  ``[...]`` or ``(...)`` with whitespace anywhere ignored:
  ``V_ds``, ``V_ds[V]``, ``V_ds [V]``, ``V_ds [KV]``, ``T_j (C)``,
  ``R_dson (mOhm)`` all parse. The NAME part is matched
  case-insensitively with ``_``/spaces/``/`` stripped; the unit is
  case-SENSITIVE except that ``k``/``K`` (kilo) and ``u``/``µ``/``μ``
  (micro) are equivalent, per the agreed contract.
* Reverse-polarity names are negated into the canonical column:
  ``V_sd`` -> ``V_ds`` negated, ``I_sd`` -> ``I_ds`` negated,
  ``V_sg`` -> ``V_gs`` negated.
* SI prefixes f p n u(µ,μ) m (none) k(K) M G T P; bases V, A, F,
  OHM/Omega, C (coulomb for charge fields; Celsius/Kelvin for
  temperature -- context decides), s, and composite A/s for di/dt
  (``[MA/s]`` == ``[A/us]`` numerically).
* Per-FILE defaults (not per-row): no ``V_ds`` column -> ``V_ds = V_gs``
  (diode-connected data); no gate column -> ``V_gs = 0`` (gate-off,
  e.g. body-diode data); no temperature column -> the caller's
  reference temperature (default 25 degC). A bare temperature column
  is Celsius (datasheet convention) -- never Kelvin-by-default, which
  would read "25" as 25 K.
* Files are classified by their column set into IV / CAP / QG / RDS
  (Rds(on) chord specs) / VTH (threshold specs) / REC (recovery specs)
  tables. Every applied default, negation and ignored column lands in
  the table's ``notes`` so the GUI/report can show it.

The fitter-facing split (channel rows vs diode rows vs third-quadrant
channel rows) is Phase B -- this module only produces clean tables.
"""

from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# SI prefixes (case-sensitive; k/K and u/mu-signs folded at parse).
_PREFIX = {
    "f": 1e-15, "p": 1e-12, "n": 1e-9, "u": 1e-6, "µ": 1e-6, "μ": 1e-6,
    "m": 1e-3, "": 1.0, "k": 1e3, "K": 1e3, "M": 1e6, "G": 1e9,
    "T": 1e12, "P": 1e15,
}

# canonical field -> (kind, sign): kind is the semantic column, sign is
# the polarity factor applied to the stored values (reverse-polarity
# names negate).  Names are normalized: lower case, non-alphanumerics
# stripped ("V_ds", "Vsd", "I_D" -> "vds", "vsd", "id").
_FIELD_BY_NAME: dict[str, tuple[str, int]] = {
    "vgs": ("vgs", 1), "vsg": ("vgs", -1), "gs": ("vgs", 1),
    "vds": ("vds", 1), "vsd": ("vds", -1), "ds": ("vds", 1),
    "v": ("vds", 1),
    "ids": ("ids", 1), "isd": ("ids", -1), "id": ("ids", 1),
    "i": ("ids", 1),
    "tj": ("tj", 1), "t": ("tj", 1), "temp": ("tj", 1),
    "temperature": ("tj", 1),
    "vth": ("vth", 1), "vth0": ("vth", 1),
    "rds": ("rds", 1), "rdson": ("rds", 1), "rdsone": ("rds", 1),
    "ron": ("rds", 1),
    "cap": ("cap", 1), "c": ("cap", 1),
    "ciss": ("cap_iss", 1), "coss": ("cap_oss", 1), "crss": ("cap_rss", 1),
    "cgs": ("cap_gs", 1), "cgd": ("cap_gd", 1), "cds": ("cap_ds", 1),
    "qg": ("qg", 1), "gatecharge": ("qg", 1),
    "qrr": ("qrr", 1), "qrrc": ("qrr", 1), "trr": ("trr", 1),
    "irm": ("irm", 1), "ir": ("irm", 1),
    "if": ("if_", 1), "idf": ("if_", 1), "idfwd": ("if_", 1),
    "didt": ("didt", 1),
    "condition": ("condition", 1), "index": ("index", 1),
    "type": ("type", 1), "note": ("note", 1),
}

# unit string -> multiplier, validated per field kind (single-char prefix).
# Multi-char bases OHM/omega/sec fold case (no collisions); single-letter
# bases stay strictly case-sensitive per the contract (mV != MV).
_UNIT_RE = re.compile(
    r"^(?P<pre>[fpnuµμmkKMGTPT]?)(?P<base>V|A|F|(?i:OHM)|Ω|C|S|s|(?i:sec))$")
# di/dt composite: [MA/s] == [A/us] numerically
_DIDT_RE = re.compile(
    r"^(?P<npre>[fpnuµμmkKMGTPT]?)(?P<nbase>A)\s*/\s*"
    r"(?P<dpre>[fpnumµμ]?)(?P<dbase>s|sec)$")
# trailing (...)/[...] unit group, whitespace before it allowed
_UNIT_GROUP_RE = re.compile(r"[\[(](?P<unit>[^\]()]*)[\])]\s*$")


class IngestError(ValueError):
    """Raised with file/column context so the GUI can show one clear line."""

    def __init__(self, source: str, column: str, message: str) -> None:
        super().__init__(f"{source}: column {column!r}: {message}")
        self.source, self.column, self.message = source, column, message


def _norm_name(raw: str) -> tuple[str, str | None]:
    """Split a header into (normalized name, raw unit text or None)."""
    m = _UNIT_GROUP_RE.search(raw)
    unit = m.group("unit") if m else None
    name = raw[:m.start()] if m else raw
    return re.sub(r"[^a-z0-9]", "", name.lower()), unit


def _parse_unit(source: str, raw: str, kind: str, unit: str | None) -> float:
    """Multiplier to canonical units for one field kind.

    Canonical: volts, amps, ohms, farads, coulombs, degC, seconds, A/s.
    """
    if unit is None or unit.strip() == "":
        return 1.0  # bare column: canonical unit (Celsius for temp)
    u = re.sub(r"\s+", "", unit)
    if kind == "didt":
        m = _DIDT_RE.match(u)
        if not m:
            raise IngestError(source, raw, f"di/dt unit {unit!r} not A/<time>")
        den = _PREFIX[m.group("dpre")]
        return _PREFIX[m.group("npre")] / den
    if kind == "tj":
        if u == "C":
            return 1.0
        if u == "K":
            return 1.0  # offset handled by caller (kind == "tj" + base K)
        if u == "F":
            raise IngestError(source, raw, "degF not supported")
        raise IngestError(source, raw, f"temperature unit {unit!r} not C/K")
    m = _UNIT_RE.match(u)
    if not m:
        raise IngestError(source, raw, f"unrecognized unit {unit!r}")
    pre, base = m.group("pre"), m.group("base")
    want = {"vgs": "V", "vds": "V", "vth": "V", "ids": "A", "if_": "A",
            "irm": "A", "rds": "OHM", "cap": "F", "qg": "C", "qrr": "C",
            "trr": "s"}[kind]
    if want == "OHM":
        ok = base.upper() == "OHM" or base == "Ω"
    elif want == "s":
        ok = base == "s" or base.lower() == "sec"   # 'S' = siemens: reject
    else:
        ok = base == want
    if not ok:
        raise IngestError(source, raw,
                          f"unit {unit!r} not {want} for a {kind} column")
    return _PREFIX[pre]


@dataclass
class _Col:
    raw: str
    kind: str          # canonical field kind
    sign: int
    scale: float       # multiplier to canonical units
    cap_kind: str = ""  # for cap_* header names: iss/oss/rss/gs/gd/ds


@dataclass
class IVTable:
    vgs: np.ndarray
    vds: np.ndarray
    ids: np.ndarray
    tj: np.ndarray
    gate_known: bool
    source: str
    notes: tuple[str, ...] = ()


@dataclass
class CapTable:
    v: np.ndarray          # reverse voltage (V_ds)
    c: np.ndarray          # farads
    kind: str              # iss/oss/rss/gs/gd/ds
    source: str
    notes: tuple[str, ...] = ()


@dataclass
class QgTable:
    vgs: np.ndarray
    qg: np.ndarray         # coulombs
    condition: np.ndarray  # object array of raw condition strings
    source: str
    notes: tuple[str, ...] = ()


@dataclass
class RdsSpec:
    vgs: np.ndarray
    rds: np.ndarray        # ohms
    i_test: np.ndarray     # amps (0 = unknown, resolve from condition/user)
    tj: np.ndarray
    condition: np.ndarray
    source: str
    notes: tuple[str, ...] = ()


@dataclass
class VthSpec:
    tj: np.ndarray
    vth: np.ndarray        # volts
    i_test: np.ndarray     # amps
    source: str
    notes: tuple[str, ...] = ()


@dataclass
class RecTable:
    qrr: np.ndarray        # coulombs
    trr: np.ndarray        # seconds
    irm: np.ndarray        # amps (peak reverse; OPTIONAL input -- NaN when
                           # the CSV omits it; the fit measures it itself)
    i_f: np.ndarray        # amps (forward current before commutation)
    didt: np.ndarray       # A/s
    source: str
    notes: tuple[str, ...] = ()


@dataclass
class _MultiCap:
    tables: list[CapTable]


Table = IVTable | CapTable | QgTable | RdsSpec | VthSpec | RecTable | _MultiCap

_REC_KINDS = {"qrr", "trr", "irm", "if_", "didt"}


def _read_columns(source: str, header: list[str]) -> list[_Col]:
    cols: list[_Col] = []
    for raw in header:
        name, unit = _norm_name(raw)
        if name == "":
            raise IngestError(source, raw, "empty column name")
        if name not in _FIELD_BY_NAME:
            raise IngestError(source, raw,
                              f"unknown column {raw!r} (rename or delete)")
        kind, sign = _FIELD_BY_NAME[name]
        cap_kind = ""
        if kind.startswith("cap_"):
            kind, cap_kind = "cap", kind[4:]
        scale = _parse_unit(source, raw, kind, unit)
        cols.append(_Col(raw, kind, sign, scale, cap_kind))
    return cols


def _col_data(rows: list[dict], col: _Col) -> np.ndarray:
    out = np.empty(len(rows))
    for i, r in enumerate(rows):
        txt = (r[col.raw] or "").strip()
        try:
            out[i] = float(txt) * col.scale * col.sign
        except ValueError as exc:
            raise IngestError(col.raw, f"row {i + 2}",
                              f"not a number: {txt!r}") from exc
    return out


def load_table(path: str | Path, ref_temp: float = 25.0) -> Table:
    """Parse one CSV into a typed table (see module docstring)."""
    path = Path(path)
    source = path.name
    with path.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        if not reader.fieldnames:
            raise IngestError(source, "-", "empty file")
        cols = _read_columns(source, list(reader.fieldnames))
        rows = list(reader)
    if not rows:
        raise IngestError(source, "-", "no data rows")

    by_kind: dict[str, list[_Col]] = {}
    for c in cols:
        by_kind.setdefault(c.kind, []).append(c)
    notes: list[str] = []
    for c in cols:
        if c.sign == -1:
            notes.append(f"{c.raw!r}: reverse-polarity name, negated into "
                         f"{c.kind.upper()}")
        if c.kind == "tj" and _norm_name(c.raw)[1] is None:
            notes.append(f"{c.raw!r}: bare temperature column, assumed "
                         "Celsius")

    kinds = set(by_kind)
    if kinds & _REC_KINDS:
        return _load_rec(source, by_kind, rows, notes)
    if "vth" in kinds:
        return _load_vth(source, by_kind, rows, notes)
    if "rds" in kinds:
        return _load_rds(source, by_kind, rows, notes)
    if "cap" in kinds:
        return _load_cap(source, by_kind, rows, notes)
    if "qg" in kinds:
        return _load_qg(source, by_kind, rows, notes)
    if "ids" in kinds:
        return _load_iv(source, by_kind, rows, notes, ref_temp)
    raise IngestError(source, "-",
                      "no current/cap/charge/recovery column -- cannot "
                      "classify this CSV")


def _get(source: str, by_kind: dict, rows: list[dict], kind: str,
         required: bool) -> np.ndarray | None:
    if kind not in by_kind:
        if required:
            raise IngestError(source, "-", f"missing required {kind} column")
        return None
    col = by_kind[kind][0]
    return _col_data(rows, col)


def _cond_vals(cond: list[dict[str, float]], key: str,
               n: int) -> np.ndarray | None:
    """Per-row values from parsed condition strings (IMCQ120R004M2H_-
    Transfer.csv style: "V_ds=20V,T_j=25C" per row). None unless EVERY
    row carries the key -- partial condition data must not half-fill a
    column."""
    if cond and len(cond) == n and all(key in c for c in cond):
        return np.array([c[key] for c in cond], dtype=float)
    return None


def _tj_arr(source: str, by_kind: dict, rows: list[dict],
            ref_temp: float,
            cond: list[dict[str, float]] | None = None
            ) -> tuple[np.ndarray, list[str]]:
    notes: list[str] = []
    if "tj" in by_kind:
        col = by_kind["tj"][0]
        unit = _norm_name(col.raw)[1]
        arr = _col_data(rows, col)
        if unit is not None and re.sub(r"\s+", "", unit) == "K":
            arr = arr - 273.15
            notes.append(f"{col.raw!r}: Kelvin converted to Celsius")
    else:
        arr = _cond_vals(cond or [], "tj", len(rows))
        if arr is not None:
            notes.append("T_j taken from condition column")
        else:
            arr = np.full(len(rows), ref_temp)
            notes.append(f"no temperature column: defaulted to "
                         f"{ref_temp:g} degC")
    return arr, notes


def _load_iv(source: str, by_kind: dict, rows: list[dict],
             notes: list[str], ref_temp: float) -> IVTable:
    ids = _get(source, by_kind, rows, "ids", True)
    cond = [parse_condition(str(c)) for c in _cond_arr(source, by_kind, rows)]
    vgs = _get(source, by_kind, rows, "vgs", False)
    gate_known = vgs is not None
    if vgs is None:
        vgs = _cond_vals(cond, "vgs", len(rows))
        if vgs is not None:
            gate_known = True
            notes.append("V_gs taken from condition column")
        else:
            vgs = np.zeros(len(rows))
            notes.append("no gate column: V_gs = 0 (gate-off data, e.g. "
                         "body-diode)")
    vds = _get(source, by_kind, rows, "vds", False)
    if vds is None:
        vds = _cond_vals(cond, "vds", len(rows))
        if vds is not None:
            notes.append("V_ds taken from condition column")
        else:
            vds = vgs.copy()
            notes.append("no V_ds column: defaulted to V_ds = V_gs "
                         "(diode-connected)")
    tj, tn = _tj_arr(source, by_kind, rows, ref_temp, cond=cond)
    return IVTable(vgs, vds, ids, tj, gate_known, source,
                   (*notes, *tn))


def _load_cap(source: str, by_kind: dict, rows: list[dict],
              notes: list[str]) -> Table:
    v = _get(source, by_kind, rows, "vds", True)
    out: list[Table] = []
    typecol = by_kind.get("type")
    types = (np.array([r[typecol[0].raw].strip().lower() for r in rows])
             if typecol else None)
    for col in by_kind["cap"]:
        c = _col_data(rows, col)
        kind = col.cap_kind
        if kind == "" and types is None:
            raise IngestError(source, col.raw,
                              "cap column without ciss/coss/... header "
                              "needs a Type column")
        if kind == "":
            kinds = types
        else:
            kinds = np.full(len(rows), kind)
        for kk in sorted(set(kinds.tolist())):
            sel = kinds == kk
            out.append(CapTable(v[sel], c[sel], str(kk), source,
                                (*notes, f"cap column {col.raw!r} -> "
                                         f"{kk}")))
    return out[0] if len(out) == 1 else _MultiCap(out)


@dataclass
class _MultiCap:
    tables: list[CapTable]


def _load_qg(source: str, by_kind: dict, rows: list[dict],
             notes: list[str]) -> QgTable:
    vgs = _get(source, by_kind, rows, "vgs", True)
    qg = _get(source, by_kind, rows, "qg", True)
    cond = _cond_arr(source, by_kind, rows)
    return QgTable(vgs, qg, cond, source, notes)


def _load_rds(source: str, by_kind: dict, rows: list[dict],
              notes: list[str]) -> RdsSpec:
    rds = _get(source, by_kind, rows, "rds", True)
    cond = [parse_condition(str(c)) for c in _cond_arr(source, by_kind, rows)]
    vgs = _get(source, by_kind, rows, "vgs", False)
    if vgs is None:
        vgs = _cond_vals(cond, "vgs", len(rows))
        if vgs is not None:
            notes.append("V_gs taken from condition column")
    vgs = np.zeros(len(rows)) if vgs is None else vgs
    i_test = _get(source, by_kind, rows, "ids", False)
    if i_test is None:
        i_test = _cond_vals(cond, "ids", len(rows))
        if i_test is not None:
            notes.append("I_test taken from condition column")
    i_test = np.zeros(len(rows)) if i_test is None else i_test
    if i_test.max() == 0:
        notes.append("Rds(on) chord test current not in a column: resolve "
                     "from condition text or user input")
    tj, tn = _tj_arr(source, by_kind, rows, 25.0, cond=cond)
    return RdsSpec(vgs, rds, i_test, tj, np.array(
        [str(c) for c in _cond_arr(source, by_kind, rows)], dtype=object),
        source, (*notes, *tn))


def _load_vth(source: str, by_kind: dict, rows: list[dict],
              notes: list[str]) -> VthSpec:
    vth = _get(source, by_kind, rows, "vth", True)
    cond = [parse_condition(str(c)) for c in _cond_arr(source, by_kind, rows)]
    i_test = _get(source, by_kind, rows, "ids", False)
    if i_test is None:
        i_test = _cond_vals(cond, "ids", len(rows))
        if i_test is not None:
            notes.append("I_test taken from condition column")
    if i_test is None:
        raise IngestError(source, "-", "missing required ids column")
    tj, tn = _tj_arr(source, by_kind, rows, 25.0, cond=cond)
    return VthSpec(tj, vth, i_test, source, (*notes, *tn))


def _load_rec(source: str, by_kind: dict, rows: list[dict],
              notes: list[str]) -> RecTable:
    for req in ("qrr", "trr", "didt"):
        if req not in by_kind:
            raise IngestError(source, "-", f"recovery CSV missing column "
                                           f"{req.upper()}")
    # I_RM (peak reverse recovery current) is OPTIONAL (user, session 30:
    # datasheets often omit it). The extraction measures its own peak
    # reverse current internally -- irm is kept for reporting only.
    if "irm" in by_kind:
        irm = _get(source, by_kind, rows, "irm", True)
    else:
        irm = np.full(len(rows), np.nan)
        notes.append("I_RM column absent (optional): the extraction "
                     "produces the peak reverse current itself")
    if "if_" in by_kind:
        i_f = _get(source, by_kind, rows, "if_", True)
    elif "ids" in by_kind:
        # contract: forward current may arrive as I_sd (negated into I_ds)
        i_f = np.abs(_get(source, by_kind, rows, "ids", True))
        notes.append("I_sd column used as forward current I_F (magnitude)")
    else:
        raise IngestError(source, "-", "recovery CSV missing column "
                                       "I_F (or I_SD)")
    arrs = [_get(source, by_kind, rows, k, True)
            for k in ("qrr", "trr", "didt")]
    return RecTable(arrs[0], arrs[1], irm, i_f, arrs[2], source, notes)


def _cond_arr(source: str, by_kind: dict, rows: list[dict]) -> np.ndarray:
    col = by_kind.get("condition")
    if col is None:
        return np.array([""] * len(rows), dtype=object)
    return np.array([r[col[0].raw] or "" for r in rows], dtype=object)


_COND_PAIR_RE = re.compile(
    r"(?P<name>[A-Za-z_][A-Za-z_0-9]*)\s*=\s*"
    r"(?P<val>-?\d+(?:\.\d*)?(?:[eE][+-]?\d+)?)(?P<unit>[A-Za-zµμ]*)")


def parse_condition(text: str) -> dict[str, float]:
    """Parse a datasheet condition string.

    ``"V_ds=20V,T_j=25C"`` -> ``{"vds": 20.0, "tj": 25.0}``;
    ``"I_d=12.5A, V_ds=20V"`` -> ``{"ids": 12.5, "vds": 20.0}``.
    Temperatures are canonical CELSIUS ("T_j=298.15K" -> 25.0).
    Unknown names keep their normalized spelling.
    """
    out: dict[str, float] = {}
    for m in _COND_PAIR_RE.finditer(text or ""):
        name = re.sub(r"[^a-z0-9]", "", m.group("name").lower())
        kind, sign = _FIELD_BY_NAME.get(name, (name, 1))
        if kind.startswith("cap_"):
            kind = "cap"
        unit = m.group("unit") or None
        scale = _parse_unit("condition", m.group("name") + m.group("unit"),
                            kind, unit)
        if kind == "tj" and unit is not None and unit.strip() == "K":
            out[kind] = float(m.group("val")) - 273.15  # canonical degC
        else:
            out[kind] = float(m.group("val")) * scale * sign
    return out


def load_tables(paths: list[str | Path],
                ref_temp: float = 25.0) -> tuple[list[Table], list[str]]:
    """Parse every file; return (tables, report lines)."""
    tables: list[Table] = []
    report: list[str] = []
    for p in paths:
        t = load_table(p, ref_temp=ref_temp)
        if isinstance(t, _MultiCap):
            tables.extend(t.tables)
            report.extend(f"{x.source}: {n}" for x in t.tables for n in x.notes)
        else:
            tables.append(t)
            report.extend(f"{t.source}: {n}" for n in t.notes)
    return tables, report
