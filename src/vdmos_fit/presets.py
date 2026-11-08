"""Fit presets for the release pipeline (session 30, Phase B).

Three presets, distilled from the three official records (part mapping
confirmed by the user 2026-11-08):

* ``trench`` -- Si power trench (ISC040N10NM8 recipe): SJ charge-cap
  family (v_ovmax/k_ovsharp free, k_lambda = 0), non-Kelvin
  (r_s = r_sk = 0), diode with k_n free and the LINEAR r_diode TC.
* ``sj``     -- superjunction (ST8L65N044M9 recipe): k_lambda + r_s free
  (transfer curves make the gate-loop degeneration identifiable), diode
  with k_n FIXED at 2.0 and all three TCs free.
* ``sic``    -- SiC/GaN (IMCQ120R044M2H recipe): trap-filling gain
  k_beta_trap + k_theta + k_lambda free, cap off (v_ovmax = 0), no
  body-diode data expected (stage skipped, placeholder card + warning).

Free sets are the official records' free sets; starts are generic
datasheet-scale guesses with multi-start variation of the shape cluster
(k_p / v_ovmax / k_ovsharp / k_kneesharp / v_dbsat). The TC stage law
form (k_p power-law, v_th0 linear, r_d quadratic -- NEVER linear, it
goes negative below ~ -150 degC) is common to all presets; the free TC
names are preset-configurable.

Presets are STARTING POINTS: stages.py takes them and the user may
override free sets per stage.
"""

from __future__ import annotations

from dataclasses import dataclass, field

# TC stage: same law form everywhere (see module docstring).
TC_FREE_DEFAULT = ("s_vth", "beta_p", "a_rd", "b_rd")
TC_START = (-6.5e-3, 1.5, 6.0e-3, 2.5e-5)
TC_LO = (-0.03, 0.0, -5e-3, -1e-4)
TC_HI = (0.0, 8.0, 5e-2, 1e-3)


@dataclass
class Preset:
    name: str
    note: str
    # ---- channel, reference-temperature stage
    channel_free25: tuple[str, ...]
    channel_fixed: dict[str, float] = field(default_factory=dict)
    channel_start: dict[str, float] = field(default_factory=dict)
    channel_multistart: tuple[dict[str, float], ...] = ()
    # ---- channel, temperature-coefficient stage
    tc_free: tuple[str, ...] = TC_FREE_DEFAULT
    # ---- diode (joint DC + TC; empty = stage skipped)
    diode_free: tuple[str, ...] = ()
    diode_fixed: dict[str, float] = field(default_factory=dict)
    diode_start: dict[str, float] = field(default_factory=dict)
    diode_multistart: tuple[dict[str, float], ...] = ()
    diode_bounds: dict[str, tuple[float, float]] = field(default_factory=dict)


_CHANNEL_START = dict(k_p=200.0, v_th0=3.0, v_sub=0.25, v_ovmax=5.0,
                      k_ovsharp=2.5, k_kneesharp=1.5, k_eta=0.03,
                      v_dbsat=0.7, r_d=5e-3)

_CAP_MULTISTART = (
    dict(k_p=222, v_ovmax=5.8, k_ovsharp=2.2, k_kneesharp=1.8, v_dbsat=0.53),
    dict(k_p=400, v_ovmax=3.0, k_ovsharp=6.0, k_kneesharp=1.0, v_dbsat=1.0),
    dict(k_p=150, v_ovmax=8.0, k_ovsharp=1.5, k_kneesharp=2.5, v_dbsat=0.3),
    dict(k_p=600, v_ovmax=2.0, k_ovsharp=10.0, k_kneesharp=0.6, v_dbsat=2.0),
    dict(k_p=300, v_ovmax=4.5, k_ovsharp=3.5, k_kneesharp=1.5, v_dbsat=0.8),
)

_DIODE_START = dict(i_s=1e-12, k_n=1.1, i_sr=0.0, k_nr=2.0, i_kf=100.0,
                    a_ikf=0.0, r_diode=1e-3, a_rdiode=0.0, b_rdiode=0.0,
                    t_nom=25.0, k_xti=3.0, v_eg=1.11)

_DIODE_MULTISTART = (
    dict(i_s=1e-12, k_n=1.2, r_diode=1e-3, i_kf=100, v_eg=1.11,
         a_rdiode=0.0),
    dict(i_s=1e-10, k_n=1.05, r_diode=5e-4, i_kf=1000, v_eg=1.0,
         a_rdiode=0.002),
    dict(i_s=1e-14, k_n=1.01, r_diode=2e-3, i_kf=50, v_eg=1.2,
         a_rdiode=0.001),
)

TRENCH = Preset(
    name="trench",
    note="Si power trench (ISC040 recipe): SJ-cap family, non-Kelvin, "
         "diode k_n free + linear r_diode TC",
    channel_free25=("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp",
                    "k_kneesharp", "k_eta", "v_dbsat", "r_d"),
    channel_fixed=dict(k_theta=0.0, v_vsat=1e4, k_lambda=0.0, g_k=0.0,
                       r_s=0.0, r_sk=0.0),
    channel_start=dict(_CHANNEL_START),
    channel_multistart=_CAP_MULTISTART,
    diode_free=("i_s", "k_n", "r_diode", "i_kf", "v_eg", "a_rdiode"),
    diode_start=dict(_DIODE_START),
    diode_multistart=_DIODE_MULTISTART,
    diode_bounds={"r_diode": (0.0, 0.020)},
)

SJ = Preset(
    name="sj",
    note="superjunction (ST8L65 recipe): k_lambda + r_s free (transfer "
         "curves), diode k_n fixed 2.0 + all three TCs free",
    channel_free25=("k_p", "v_th0", "v_sub", "k_lambda", "r_d", "r_s",
                    "v_dbsat", "v_ovmax", "k_ovsharp", "k_kneesharp"),
    channel_fixed=dict(k_theta=0.0, g_k=0.0, r_sk=0.0),
    channel_start=dict(_CHANNEL_START, k_lambda=0.004, r_s=5e-3),
    channel_multistart=_CAP_MULTISTART,
    diode_free=("i_s", "r_diode", "v_eg", "i_kf", "a_ikf", "a_rdiode",
                "b_rdiode"),
    diode_fixed=dict(k_n=2.0),
    diode_start=dict(_DIODE_START, k_n=2.0),
    diode_multistart=_DIODE_MULTISTART,
    diode_bounds={"r_diode": (0.0, 0.020)},
)

SIC = Preset(
    name="sic",
    note="SiC/GaN (IMCQ recipe): k_beta_trap + k_theta + k_lambda free, "
         "cap off; body-diode stage skipped without diode CSVs",
    channel_free25=("k_p", "v_th0", "v_sub", "k_beta_trap", "k_theta",
                    "k_kneesharp", "k_eta", "v_dbsat", "k_lambda", "r_d"),
    channel_fixed=dict(v_ovmax=0.0, g_k=0.0, v_vsat=1e4, r_s=0.0, r_sk=0.0),
    channel_start=dict(_CHANNEL_START, k_beta_trap=0.0, k_theta=0.0,
                       k_lambda=0.004),
    channel_multistart=(
        dict(k_p=200, k_kneesharp=1.5, v_dbsat=0.7, k_lambda=0.004),
        dict(k_p=400, k_kneesharp=0.8, v_dbsat=1.5, k_lambda=0.001),
        dict(k_p=100, k_kneesharp=2.5, v_dbsat=0.3, k_lambda=0.008),
    ),
    diode_free=(),
)

PRESETS: dict[str, Preset] = {p.name: p for p in (TRENCH, SJ, SIC)}
