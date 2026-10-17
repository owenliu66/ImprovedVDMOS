"""The official fitted device records -- single source of truth.

Promoted 2026-10-13: fit A of ``scripts/fit_st8l65_nd.py`` (10-param, drain
knee softness ``k_kneesharp`` free, full 41-point set) is the quoted ST8L65N044M9
model. It supersedes the 9-param fit B of ``scripts/fit_st8l65.py``
(5.54 %, k_kneesharp = 0), which is kept only as the k_kneesharp = 0 baseline.

Why fit A and not the lower-RMSE fit B (1.52 %): fit A drops no measured
curve, so its RMSE is the honest device-quality number, while fit B exists
to isolate how much error the physically unreproducible v_gs = 5 curve
contributes. Fit A is also better conditioned (cond 2.19e3 vs 1.27e4) and
more physical (k_p 579 vs 2354-on-the-bound, r_d/r_s both interior).

Quote the curve shape and r_d+r_s, never the individual knee parameters: the
cluster (k_p, v_ovmax, k_ovsharp, k_kneesharp, v_dbsat, k_lambda) is entangled
at |corr| up to 0.999 and k_lambda/v_dbsat carry sigma = inf because the
dataset stops at v_ds = 20 V and never enters deep saturation. That is
expected, not a fit failure -- see STATUS.md, "Current ST8L65 fit".

``tests/test_st8l65_official.py`` pins this record against the measured CSV,
so a model or bounds change that silently moves the official fit fails the
suite instead of quietly invalidating every quoted number.

(The record was re-labelled with the 2026-11 naming convention; the fitted
values are unchanged and still reproduce the published RMSE exactly.)
"""

from __future__ import annotations

from .params import Params

__all__ = [
    "ST8L65_OFFICIAL",
    "ST8L65_OFFICIAL_RMSE_REL",
    "ST8L65_FIT_A_FREE",
    "ST8L65_DIODE",
    "ST8L65_DIODE_RMSE_REL",
    "ST8L65_DIODE_FREE",
]

# Free set of the promoted fit (order matters only for readability).
ST8L65_FIT_A_FREE: tuple[str, ...] = (
    "k_p", "v_th0", "v_sub", "k_lambda", "r_d", "r_s", "v_dbsat", "v_ovmax",
    "k_ovsharp", "k_kneesharp",
)

# Relative RMSE of the record on data/ST8L65N044M9_IV.csv (41 points,
# relative weighting) as reported by FitResult.rmse_rel.
ST8L65_OFFICIAL_RMSE_REL = 0.0379

ST8L65_OFFICIAL = Params(
    k_p=579.453,
    v_th0=5.06184,
    v_sub=0.243283,
    v_ovmax=1.06399,
    k_ovsharp=7.25141,
    k_theta=0.0,
    v_vsat=1e4,
    k_kneesharp=0.435907,
    k_eta=0.0,
    k_lambda=5.12947e-10,
    g_k=0.0,
    v_thk=1.6,
    v_kovsoft=0.1,
    v_dbsat=339.02,
    r_d=0.0243908,
    r_s=0.00729014,
)

# ---------------------------------------------------------------------------
# Body diode (session 12; TCs promoted session 14; Trs2 promoted session 15):
# fit A of ``scripts/fit_st8l65_diode.py`` -- free i_s, r_diode, v_eg, i_kf,
# a_ikf, a_rdiode, b_rdiode on all 30 points of data/ST8L65N044M9_diode.csv
# (-55/25/150 degC, 5-55 A). k_n FIXED at 2.0 (physical diffusion limit --
# freeing it buys ~1 % RMSE, lands at k_n ~ 1.3-2.1 depending on the TC set
# (the valley direction moves), with corr(i_s, k_n) ~ +1.0, and drifts the
# reverse-blocking floor by orders of magnitude: the ideality/prefactor
# split is a degenerate valley, so k_n comes from physics, not this data),
# k_xti fixed at 3.0 (LTspice default; freeing it gains 0.04 % and pins
# k_xti <-> v_eg at -0.99), i_sr = 0 (recombination branch unsensed above 5 A).
#
# a_ikf/a_rdiode/b_rdiode are the LTspice Tikf/Trs1/Trs2 temperature
# coefficients (Tnom-referenced, verified to ~1e-5 against LTspice 26.0.1 --
# test_diode.py::test_ltspice_tc_laws). Promotion history: TC-off fit A
# (session 12) 4.78 % -> linear TCs (session 14) 4.15 % -> + Trs2 quadratic
# (session 15, user override) 3.91 %. The Trs2 gain is modest and the
# five-knob r_diode/i_kf/a_ikf/a_rdiode/b_rdiode block is STRONGLY correlated
# (|corr| up to 0.985, sigma_rel(i_kf) ~ 15): quote the CURVE, never the
# split. The stubborn 150 degC residual is a ~12 mV digitization error on
# the 5 A point (exponential region -- r_diode carries only ~10 mV of the
# 562 mV drop), not missing series-R physics; excluding that point the
# 150 degC rms is ~3.8 % for both the 6- and 7-param fits. Baselines live on
# as script probes: tc_off_probe 4.78 %, trs2_off_probe 4.15 %.
#
# i_s is weakly determined (sigma_rel ~ 1.1 -- a pure exponential prefactor,
# no low-current tail in the data) and only meaningful as the PAIR with k_n
# that reproduces the curve; v_eg ~ 1.08 eV (silicon bandgap) is genuinely
# sensed (sigma_rel ~ 0.18) and r_diode ~ 1.9 mOhm is the diode-path series
# resistance (the diode bypasses the package r_d/r_s).
# ---------------------------------------------------------------------------

ST8L65_DIODE_FREE: tuple[str, ...] = (
    "i_s", "r_diode", "v_eg", "i_kf", "a_ikf", "a_rdiode", "b_rdiode",
)

ST8L65_DIODE_RMSE_REL = 0.0391

ST8L65_DIODE = Params(
    i_s=3.85056e-06,
    k_n=2.0,
    i_sr=0.0,
    k_nr=2.0,
    i_kf=21.9663,
    a_ikf=0.00832096,
    r_diode=0.00189604,
    a_rdiode=0.00483208,
    b_rdiode=-1.78406e-05,
    t_nom=25.0,
    k_xti=3.0,
    v_eg=1.08196,
)
