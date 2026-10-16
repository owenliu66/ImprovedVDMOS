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

__all__ = ["ST8L65_OFFICIAL", "ST8L65_OFFICIAL_RMSE_REL", "ST8L65_FIT_A_FREE"]

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
