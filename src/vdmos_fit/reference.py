"""The official fitted device records -- single source of truth.

Promoted 2026-10-13: fit A of ``scripts/fit_st8l65_nd.py`` (10-param, drain
knee softness ``nd`` free, full 41-point set) is the quoted ST8L65N044M9
model. It supersedes the 9-param fit B of ``scripts/fit_st8l65.py``
(5.54 %, nd = 0), which is kept only as the nd = 0 baseline.

Why fit A and not the lower-RMSE fit B (1.52 %): fit A drops no measured
curve, so its RMSE is the honest device-quality number, while fit B exists
to isolate how much error the physically unreproducible Vgs = 5 curve
contributes. Fit A is also better conditioned (cond 2.19e3 vs 1.27e4) and
more physical (kp 579 vs 2354-on-the-bound, rd/rs both interior).

Quote the curve shape and rd+rs, never the individual knee parameters: the
cluster (kp, kcap, m, nd, vk, lam) is entangled at |corr| up to 0.999 and
lam/vk carry sigma = inf because the dataset stops at Vds = 20 V and never
enters deep saturation. That is expected, not a fit failure -- see
STATUS.md, "Current ST8L65 fit".

``tests/test_st8l65_official.py`` pins this record against the measured CSV,
so a model or bounds change that silently moves the official fit fails the
suite instead of quietly invalidating every quoted number.
"""

from __future__ import annotations

from .params import Params

__all__ = ["ST8L65_OFFICIAL", "ST8L65_OFFICIAL_RMSE_REL", "ST8L65_FIT_A_FREE"]

# Free set of the promoted fit (order matters only for readability).
ST8L65_FIT_A_FREE: tuple[str, ...] = (
    "kp", "vth0", "ksub", "lam", "rd", "rs", "vk", "kcap", "m", "nd",
)

# Relative RMSE of the record on data/ST8L65N044M9_IV.csv (41 points,
# relative weighting) as reported by FitResult.rmse_rel.
ST8L65_OFFICIAL_RMSE_REL = 0.0379

ST8L65_OFFICIAL = Params(
    kp=579.453,
    vth0=5.06184,
    ksub=0.243283,
    kcap=1.06399,
    m=7.25141,
    theta=0.0,
    vsat=1e4,
    nd=0.435907,
    eta=0.0,
    lam=5.12947e-10,
    gk=0.0,
    vthk=1.6,
    kk=0.1,
    vk=339.02,
    rd=0.0243908,
    rs=0.00729014,
)
