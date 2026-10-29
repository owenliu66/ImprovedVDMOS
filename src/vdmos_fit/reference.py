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

from dataclasses import dataclass

from .caps import CapParams, OriginAnchor
from .params import Params

__all__ = [
    "ST8L65_OFFICIAL",
    "ST8L65_OFFICIAL_RMSE_REL",
    "ST8L65_FIT_A_FREE",
    "ST8L65_DIODE",
    "ST8L65_DIODE_RMSE_REL",
    "ST8L65_DIODE_FREE",
    "IMCQ120_JOINT",
    "IMCQ120_JOINT_RMSE_REL",
    "IMCQ120_JOINT_GROUP_RMSE_REL",
    "IMCQ120_JOINT_FREE",
    "ImcqJointRecord",
    "CAP_ST8L65",
    "CAP_IMCQ120",
    "CAP_ISC040",
    "CAP_ANCHOR_ST8L65",
    "CAP_ANCHOR_IMCQ120",
    "CAP_ANCHOR_ISC040",
    "CAP_PLATEAU_DEFICIT_REL",
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

# ---------------------------------------------------------------------------
# IMCQ120R004M2H (CoolSiC 1200 V, 4 mOhm class, Kelvin source) -- JOINT
# multi-temperature channel record (session 21, user promotion).
#
# Protocol: scripts/fit_imcq120_joint.py -- ONE least-squares vector over
# 200 points / 29 temperatures: IV 25 C (78) + IV 175 C (73)
# + transfer @ V_DS = 20 V, 25/175/200 C (23, data/IMCQ120R004M2H_Transfer.csv)
# + Rds(on)-Tj chords @185.2 A, Vgs 15/18 (22) + Vth-Tj @58.3 mA (4).
# The TC LAWS below are the only extra unknowns -- no independent per-T fits.
# Report: out/IMCQ120_JOINT_FIT_REPORT.md (the session-18 per-T report
# out/IMCQ120_TEMPERATURE_REPORT.md was deleted in session 21; its durable
# findings are carried in that report's section 8).
#
# TC laws (25 C base, evaluated by params_at):
#   k_p(T)         = k_p25 * (T_K/298.15)^(-beta_p)
#   v_th0(T)       = v_th0_25 + s_vth*(T - 25)          [-7.16 mV/K]
#   v_sub          = const (measured invariant -- a shape knob, not a thermal
#                    voltage; ln10*v_sub ~ 1.8 V/dec is a model floor)
#   k_beta_trap(T) = max(k_bt25 + s_bkt*(T - 25), 0)    [linear; the Arrhenius
#                    variant ties at 1.92 vs 1.93 % with Ea sigma = inf --
#                    the joint fit decided: linear]
#   r_d(T)         = r_d25*(1 + a_rd*(T-25) + b_rd*(T-25)^2)  [QUADRATIC,
#                    Trs1/Trs2 law form -- user call, confirmed 7x over the
#                    power law (chord fits: quad 0.067/0.065 % RMSE vs power
#                    9.5/9.1 %, exp 5.9/6.3 %, linear 13.4/10.9 %). NEVER
#                    linear (goes negative below -26 C); vertex ~ -38 C --
#                    valid over the -46..200 C data span, not below.]
#   k_theta, k_kneesharp, k_eta, v_dbsat, k_lambda: T-independent.
# Fixed: r_s = 0 (Kelvin source), g_k = 0, v_ovmax = 0 (cap off -- not sensed
# on this SiC part), v_vsat = 1e4 (off).
#
# QUOTE THE CURVES + TC LAWS, NEVER THE PARAMETER SPLIT:
#  * k_p25 ~ k_bt25 corr -0.95 (9.32 / 0.469 vs the per-T gain-fit anchors
#    25.74 / 0.057): the pair is a model-form split, the composite
#    k_p(T)*v_ov_eff2(T) is the record. beta_p ~ s_bkt +0.97: the channel TC
#    is shared between them. Same doctrine as the SJ cap<->gain record.
#  * k_eta pinned at the 0.05 search cap (sigma inf), v_dbsat 320 (sigma inf),
#    k_lambda ~ 4.4e-10 (sigma inf): the DIBL/CLM block is INERT on this
#    device (k_db(20 V) ~ 0.06) -- bounds, not measurements. Re-freeing
#    k_lambda (session-21 ablation) confirms CLM is unnecessary: knee
#    saturation (k_theta = 0.109, v_knee -> 9.2 V) + the r_d divider + the
#    transfer data carry the high-v_gs shape.
#  * beta_p is interpretable ONLY because the Rds(on)-Tj chords pin r_d(T)
#    inside the vector (IV-only ablation: a_rd sigma_rel 295, corr(a_rd,b_rd)
#    -1.00, out-of-sample Rds(on) fails at 20 % RMSE / 34 % cold).
#  * Gain derivative dv_ov_eff2/dv_ov_eff = 1 + 2*k_bt*v_ov reaches 12.5 at
#    v_gs = 18 -- legal for the package-R bracket ONLY because r_s = 0
#    (monotonicity, not <= 1, is what the proof needs). Re-check before any
#    r_s > 0 use of this record (Kelvin-aware subckt, HANDOVER queue).
#
# Acceptance at promotion (independent specs, never fitted; full table in the
# report): overall RMSE 1.93 % (IV 1.73 / transfer 3.39 / Rds(on) 0.65 /
# Vth 2.06); high-v_gs corner v_gs=18 @20 V: 2118 A at 175 C vs ~2070
# digitized (the old gain-fit overshoot was 3650 -- FIXED); Rds(on)-Tj full
# -46..195 C span max err 1.08 % (was -15.3 % at -46 C -- FIXED); Ron chords
# 3.722/4.652 mOhm @25 C (data 3.70/4.70); Vth spec 57.0/59.1 mA; ZTC 25/175
# at 14.92 V/1703 A (matches the digitized output-curve crossing ~14.9 V).
# KNOWN OPEN GAPS (pinned by tests/test_imcq120_official.py): g_fs 129.5 S vs
# 122 spec (+6.2 %) and the 25 C transfer corner -5..-7 % at v_gs 15-19 --
# both point at the unmodeled Kelvin source-side power path (r_d25 drops
# 7.2 V at 3 kA), not at channel physics.
#
# tests/test_imcq120_official.py pins this record against all four measured
# CSVs (record evaluation, no re-fit), so a model/bounds/law change that
# silently moves the official fit fails the suite.
# ---------------------------------------------------------------------------

IMCQ120_JOINT_FREE: tuple[str, ...] = (
    "k_p25", "beta_p", "v_th0_25", "s_vth", "v_sub", "k_bt25", "s_bkt",
    "k_theta", "k_ksh", "k_eta", "v_dbsat", "k_lambda",
    "r_d25", "a_rd", "b_rd",
)

# Relative RMSE on the 200-point joint set (relative weighting, same residual
# definition as vdmos_fit.fit.residuals: (model - ref)/max(|ref|, 1e-3)).
IMCQ120_JOINT_RMSE_REL = 0.0193
IMCQ120_JOINT_GROUP_RMSE_REL = {
    "iv": 0.017283,
    "transfer": 0.033930,
    "rdson": 0.006482,
    "vth": 0.020617,
}


@dataclass(frozen=True)
class ImcqJointRecord:
    """25 C base parameters + TC laws; ``params_at(t_c)`` builds the Params.

    Field names match the free-set names in scripts/fit_imcq120_joint.py
    (k_ksh = k_kneesharp, k_bt25/s_bkt = k_beta_trap level/slope).
    """

    k_p25: float
    beta_p: float
    v_th0_25: float
    s_vth: float
    v_sub: float
    k_bt25: float
    s_bkt: float
    k_theta: float
    k_ksh: float
    k_eta: float
    v_dbsat: float
    k_lambda: float
    r_d25: float
    a_rd: float
    b_rd: float
    # Fixed (not fitted) -- Kelvin source, no kink, cap off, v_vsat off.
    r_s: float = 0.0
    v_vsat: float = 1e4
    g_k: float = 0.0
    v_ovmax: float = 0.0

    def params_at(self, t_c: float) -> Params:
        """Full Params at temperature ``t_c`` [degC] through the record TC laws.

        Must stay consistent with ``params_at`` in scripts/fit_imcq120_joint.py
        (primary variant: quadratic r_d, linear clamped k_beta_trap); the
        script re-evaluates this record's RMSE through its own machinery and
        prints OK/STALE, and tests/test_imcq120_official.py evaluates the
        record directly -- a divergence between the two fails loudly.
        """
        tk = t_c + 273.15
        dt = t_c - 25.0
        kbt = max(self.k_bt25 + self.s_bkt * dt, 0.0)
        fac = 1.0 + self.a_rd * dt + self.b_rd * dt * dt
        return Params(
            k_p=self.k_p25 * (tk / 298.15) ** (-self.beta_p),
            v_th0=self.v_th0_25 + self.s_vth * dt,
            v_sub=self.v_sub,
            k_beta_trap=kbt,
            k_theta=self.k_theta,
            v_vsat=self.v_vsat,
            k_kneesharp=self.k_ksh,
            k_eta=self.k_eta,
            k_lambda=self.k_lambda,
            g_k=self.g_k,
            v_dbsat=self.v_dbsat,
            r_d=self.r_d25 * max(fac, 1e-6),  # Trs1/Trs2 clamp precedent
            r_s=self.r_s,
            v_ovmax=self.v_ovmax,
        )


IMCQ120_JOINT = ImcqJointRecord(
    k_p25=9.31762,
    beta_p=0.377994,
    v_th0_25=5.74578,
    s_vth=-0.00716141,
    v_sub=0.793309,
    k_bt25=0.468722,
    s_bkt=3.40378e-4,
    k_theta=0.109044,
    k_ksh=1.46481,
    k_eta=0.0499999,      # pinned at the 0.05 search cap -- INERT, sigma inf
    v_dbsat=319.838,      # sigma inf -- inert (k_db(20 V) ~ 0.06)
    k_lambda=4.36994e-10,  # ~ 0: CLM unnecessary (sigma inf)
    r_d25=0.00239852,
    a_rd=0.00672705,
    b_rd=5.36729e-5,
)

# ---------------------------------------------------------------------------
# Capacitance / charge (CAP_SPEC.md; extraction scripts/fit_cmax.py,
# 2026-10-18): origin anchors (first-4-point quadratic convention) and the
# extracted c_gdmax per part. The tables themselves are NOT duplicated here --
# they are rebuilt from data/<part>_Cap_Vds.csv; these records pin the splice
# parameters and the diagnostics that go with them.
#
# HEADLINE FINDING (quote the CURVE, not the split): all three parts land at
# c_gdmax ~ c_gd0 -- the exponential-saturation branch is nearly FLAT, and
# the seed c_gdmax = c_gs(0) over-predicts by 4-8x. The Qg segment-3 slope
# says C_gs + C_gd(end) ~ the c_gs table alone, i.e. these datasheets give no
# evidence for an LTspice-style Cgdmax >> Cgdmin channel cap. The channel
# charge the off-state tables miss shows up instead as the Miller-plateau
# DEFICIT (measured plateau charge minus table prediction, CAP_SPEC
# limitation 1): worst on the super-junction part (+61 %), then trench
# (+22-43 %), smallest on SiC (+10 %). That deficit -- not c_gdmax -- is the
# quantity the deferred v_gs-dependent blend must eventually explain.
#
# IMCQ120 c_gdmax is CLAMPED at the lower bracket (flat branch): the measured
# segment-3 charge is below even the flat-branch prediction under the c_gs
# table convention (spec diagnostic step 4). ISC040 c_gdmax is the median of
# three conditions (647 / 798 / 663 pF raw). ST8L65's raw 2750 pF is within
# digitiser noise of its c_gd0 (the V_on=0 probe clamps it LOWER).
# ---------------------------------------------------------------------------

CAP_ANCHOR_ST8L65 = OriginAnchor(
    c_gd0=2.700257e-09, s_gd0=1.515041e-10, b2=-2.710536e-10,
)
CAP_ANCHOR_IMCQ120 = OriginAnchor(
    c_gd0=3.822921e-09, s_gd0=6.637159e-10, b2=4.088134e-11,
)
CAP_ANCHOR_ISC040 = OriginAnchor(
    c_gd0=4.908948e-10, s_gd0=5.391591e-11, b2=3.350697e-12,
)

CAP_ST8L65 = CapParams(c_gd0=2.700257e-09, s_gd0=1.515041e-10,
                       c_gdmax=2.750175e-09)
CAP_IMCQ120 = CapParams(c_gd0=3.822921e-09, s_gd0=6.637159e-10,
                        c_gdmax=3.822921e-09)   # clamped LOWER (flat branch)
CAP_ISC040 = CapParams(c_gd0=4.908948e-10, s_gd0=5.391591e-11,
                       c_gdmax=6.628331e-10)    # median of 3 conditions

# Plateau deficit / measured plateau charge, per part (worst condition):
#   ST8L65 (400 V): +61.2 %; IMCQ120 (800 V): +9.7 %;
#   ISC040 (20/50/80 V): +42.8 / +24.5 / +22.0 %.
CAP_PLATEAU_DEFICIT_REL = {
    "ST8L65N044M9": 0.612,
    "IMCQ120R004M2H": 0.097,
    "ISC040N10NM8": 0.428,
}
