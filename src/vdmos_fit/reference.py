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
# Body diode (session 12): fit A of ``scripts/fit_st8l65_diode.py`` -- free
# i_s, r_diode, v_eg, i_kf on all 30 points of data/ST8L65N044M9_diode.csv
# (-55/25/150 degC, 5-55 A). k_n FIXED at 2.0 (physical diffusion limit --
# freeing it buys ~0.1 % RMSE, lands at k_n = 2.07 with corr(i_s, k_n) =
# +1.000, and drifts the reverse-blocking floor: the ideality/prefactor
# split is a degenerate valley, so k_n comes from physics, not this data),
# k_xti fixed at 3.0 (LTspice default; freeing it gains 0.04 % and pins
# k_xti <-> v_eg at -0.99), i_sr = 0 (recombination branch unsensed above 5 A).
#
# i_s is weakly determined (sigma_rel ~ 0.9 -- a pure exponential prefactor,
# no low-current tail in the data) and only meaningful as the PAIR with k_n
# that reproduces the curve; v_eg ~ 1.10 eV (silicon bandgap) is genuinely
# sensed (sigma_rel ~ 0.09) and r_diode ~ 1.6 mOhm is the diode-path series
# resistance (the diode bypasses the package r_d/r_s). r_diode <-> i_kf
# correlation +0.95: quote the curve, not the split.
# ---------------------------------------------------------------------------

ST8L65_DIODE_FREE: tuple[str, ...] = ("i_s", "r_diode", "v_eg", "i_kf")

ST8L65_DIODE_RMSE_REL = 0.0478

ST8L65_DIODE = Params(
    i_s=3.70127e-06,
    k_n=2.0,
    i_sr=0.0,
    k_nr=2.0,
    i_kf=17.4506,
    r_diode=0.00158622,
    t_nom=25.0,
    k_xti=3.0,
    v_eg=1.09674,
)
