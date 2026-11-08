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
    "ISC040_RECOV",
    "RecovRecord",
    "ISC040_CHANNEL",
    "IscChannelRecord",
    "ISC040_CHANNEL_FREE25",
    "ISC040_CHANNEL_TC_FREE",
    "ISC040_CHANNEL_RMSE25_REL",
    "ISC040_CHANNEL_RMSE_TC_REL",
    "ISC040_DIODE",
    "ISC040_DIODE_FREE",
    "ISC040_DIODE_RMSE_REL",
    "ISC040_DIODE_GROUP_RMSE_REL",
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
# Fixed: r_s = 0 (Kelvin source), r_sk = 0, g_k = 0, v_ovmax = 0 (cap off --
# not sensed on this SiC part), v_vsat = 1e4 (off).
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
#    (monotonicity, not <= 1, is what the proof needs). RE-CHECKED session
#    22 with r_g > 0 AND the gain active: the bracket survives (monotonicity
#    is all the proof needs) -- tests/test_kelvin_split.py.
#
# Acceptance at promotion (independent specs, never fitted; full table in the
# report): overall RMSE 1.93 % (IV 1.73 / transfer 3.39 / Rds(on) 0.65 /
# Vth 2.06); high-v_gs corner v_gs=18 @20 V: 2118 A at 175 C vs ~2070
# digitized (the old gain-fit overshoot was 3650 -- FIXED); Rds(on)-Tj full
# -46..195 C span max err 1.08 % (was -15.3 % at -46 C -- FIXED); Ron chords
# 3.722/4.652 mOhm @25 C (data 3.70/4.70); Vth spec 57.0/59.1 mA; ZTC 25/175
# at 14.92 V/1703 A (matches the digitized output-curve crossing ~14.9 V).
# KNOWN OPEN GAPS (pinned by tests/test_imcq120_official.py): g_fs 129.5 S vs
# 122 spec (+6.2 %) and the 25 C transfer corner -5..-7 % at v_gs 15-19.
# SESSION-22 UPDATE (out/IMCQ120_KELVIN_SPLIT_REPORT.md): the Kelvin-aware
# r_s split (model.py r_sk) was implemented and the joint vector refit --
# NEITHER gap is a series-R topology problem. Gate-loop degeneration is
# REJECTED by the data (r_s -> 0 boundary optimum from every start; forcing
# r_g = 0.25/0.5/1.0 mO monotonically worsens RMSE 1.93 -> 2.71 % AND g_fs);
# the current-loop split is DC-degenerate (only r_d + r_s sensed, corr -1.00).
# g_fs is pinned by the transfer curve itself (~130 S local slope at 185 A --
# spec-vs-digitized-slope softness); the corner is a channel-shape
# compromise. r_sk stays in the model for the 4-terminal SUBCKT (AC gate
# current through r_sk); the DC record keeps r_s = r_sk = 0.
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
    # r_sk (kelvin tap segment, session 22) is inert in the DC model at
    # r_s = 0; it exists for the 4-terminal subckt (AC gate current).
    r_s: float = 0.0
    r_sk: float = 0.0
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
            r_sk=self.r_sk,
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
# ISC040N10NM8 (100 V trench super-junction, NON-Kelvin package) -- channel
# record (promoted this session; closes the HANDOVER queue item that was
# BLOCKED ON DATA). Protocol scripts/fit_isc040_channel.py.
#
# SEQUENTIAL protocol (user-directed): fit 25 C and LOCK, then fit the TC
# coefficients. No high-temperature I-V data -- the TC comes from Rds(on)-vs-Vgs
# at -55/25/100/175 C (I_test = 25 A) and Vth-vs-Tj (520 uA / 52 uA).
#
# Model family: the SJ charge-cap form (v_ovmax / k_ovsharp free, k_theta = 0,
# v_vsat = 1e4 off, k_lambda = 0, g_k = 0) -- the same family as the ST8L65N044M9
# super-junction record. Beats the velocity-saturation form (1.33 % vs 1.73 %
# 25 C RMSE) and is better conditioned; the cap is genuinely sensed (it ceilings
# the V_gs = 15 channel conductance -- without it the Ron chord is ~13 % low).
#
# NO transfer curve (datasheet V_ds conditions vague) -> the gate-loop
# degeneration r_s is NOT identifiable and the data REJECTS it (r_s -> 0 from
# every start, monotonically worse as r_s grows -- the IMCQ Kelvin finding,
# session 22). The DC record carries the whole current-loop series R in r_d;
# r_s = r_sk = 0. NON-Kelvin package -> r_sk = 0 (no kelvin tap; gate loop
# returns at the power source); the 4-terminal subckt keeps the K pin for the
# driver return, inert at r_s = r_sk = 0.
#
# TC laws (25 C base, evaluated by params_at):
#   k_p(T)   = k_p25 * (T_K/298.15)^(-beta_p)          [beta_p ~ 2.20]
#   v_th0(T) = v_th0_25 + s_vth*(T - 25)               [-7.49 mV/K]
#   r_d(T)   = r_d25*(1 + a_rd*(T-25) + b_rd*(T-25)^2) [QUADRATIC, Trs1/Trs2
#              law form; vertex ~-130 C -- NEVER linear (goes negative); valid
#              over the -55..175 C data span]
#   v_sub, v_ovmax, k_ovsharp, k_kneesharp, k_eta, v_dbsat: T-independent.
#
# QUOTE THE CURVES + TC LAWS, NEVER THE PARAMETER SPLIT:
#  * r_d <-> r_s corr -1.0 (only the sum sensed); gate degeneration rejected.
#  * k_p <-> v_th0 <-> v_ovmax <-> r_d cluster (|corr| up to 0.86): the composite
#    channel conductance k_p*v_ov_eff2(v_gs) is the record.
#  * beta_p weakly determined (sigma_rel ~ 4.4): the k_p(T)/r_d(T)/v_th0(T) TC
#    split is model-form dependent; beta_p is interpretable ONLY because
#    Rds(on)-Vgs(T) pins r_d(T) at high V_gs (the r_d-dominated floor), leaving
#    the low-V_gs channel shape for beta_p -- the IMCQ doctrine.
#  * v_dbsat sigma_rel ~ 77, k_eta ~ 4, k_kneesharp ~ 7, v_ovmax ~ 17: weakly
#    sensed (data stops at V_ds ~ 5 V, no deep saturation); quote the curve.
#
# Acceptance (independent specs, never fitted): 25 C RMSE 1.33 % (132 pts =
# IV 69 + Rds(on)-Id 55 + Rds(on)-Vgs 25 C 8); TC overall 2.97 % (Rds(on) 2.35,
# Vth-voltage 4.04). Ron@25 A 25 C [mO] 3.35/3.56/3.89 (data 3.33/3.52/3.87);
# Ron@Vgs15 vs T [mO] 2.12/3.34/5.20/7.76 at -55/25/100/175 (data 2.05/3.33/
# 5.23/8.00); Vth@520uA vs T [V] 3.52/3.06/2.50/2.07 (data 3.37/3.05/2.58/2.07).
# The stubborn Vth residual is the single-softplus subthreshold model-form limit
# (the two-current swing rises with T; a constant v_sub cannot follow it) -- the
# Vth SLOPE is captured, the absolute current is not (IMCQ precedent).
#
# tests/test_isc040_official.py pins this record against the measured CSVs.
# ---------------------------------------------------------------------------

ISC040_CHANNEL_FREE25: tuple[str, ...] = (
    "k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp", "k_kneesharp",
    "k_eta", "v_dbsat", "r_d",
)
ISC040_CHANNEL_TC_FREE: tuple[str, ...] = ("s_vth", "beta_p", "a_rd", "b_rd")
ISC040_CHANNEL_RMSE25_REL = 0.013288
ISC040_CHANNEL_RMSE_TC_REL = 0.029697


@dataclass(frozen=True)
class IscChannelRecord:
    """ISC040N10NM8 25 C base parameters + TC laws; ``params_at(t_c)`` builds Params.

    Field names match scripts/fit_isc040_channel.py (k_p25/v_th0_25/r_d25 are
    the 25 C base; s_vth/beta_p/a_rd/b_rd are the TC coefficients).
    """

    k_p25: float
    v_th0_25: float
    v_sub: float
    v_ovmax: float
    k_ovsharp: float
    k_kneesharp: float
    k_eta: float
    v_dbsat: float
    r_d25: float
    s_vth: float
    beta_p: float
    a_rd: float
    b_rd: float
    # Fixed (not fitted): SJ cap form, no kink, no CLM, non-Kelvin (r_s=r_sk=0).
    k_theta: float = 0.0
    v_vsat: float = 1e4
    k_lambda: float = 0.0
    g_k: float = 0.0
    r_s: float = 0.0
    r_sk: float = 0.0
    v_kneemin: float = 1e-3

    def params_at(self, t_c: float) -> Params:
        """Full Params at temperature ``t_c`` [degC] through the record TC laws.

        Must stay consistent with ``params_at`` in scripts/fit_isc040_channel.py;
        the script re-evaluates this record's RMSE through its own machinery and
        tests/test_isc040_official.py evaluates the record directly -- a
        divergence fails loudly.
        """
        tk = t_c + 273.15
        dt = t_c - 25.0
        fac = 1.0 + self.a_rd * dt + self.b_rd * dt * dt
        return Params(
            k_p=self.k_p25 * (tk / 298.15) ** (-self.beta_p),
            v_th0=self.v_th0_25 + self.s_vth * dt,
            v_sub=self.v_sub,
            v_ovmax=self.v_ovmax,
            k_ovsharp=self.k_ovsharp,
            k_theta=self.k_theta,
            v_vsat=self.v_vsat,
            k_kneesharp=self.k_kneesharp,
            k_eta=self.k_eta,
            k_lambda=self.k_lambda,
            g_k=self.g_k,
            v_dbsat=self.v_dbsat,
            v_kneemin=self.v_kneemin,
            r_d=self.r_d25 * max(fac, 1e-6),  # Trs1/Trs2 clamp precedent
            r_s=self.r_s,
            r_sk=self.r_sk,
        )


ISC040_CHANNEL = IscChannelRecord(
    k_p25=222.5901694854354,
    v_th0_25=4.286665228917342,
    v_sub=0.2440875052714508,
    v_ovmax=5.774586111568539,
    k_ovsharp=2.1964907755354677,
    k_kneesharp=1.7870065096973613,
    k_eta=0.03130646237241051,
    v_dbsat=0.5252898687616809,
    r_d25=0.002482952377514807,
    s_vth=-0.007486343423130841,
    beta_p=2.2003970923964205,
    a_rd=0.0056205977696523295,
    b_rd=2.0007921690165413e-05,
)

# ---------------------------------------------------------------------------
# ISC040N10NM8 body diode (fit A of ``scripts/fit_isc040_diode.py``, promoted
# session 27): free i_s, k_n, r_diode, i_kf, v_eg, a_rdiode on all 26 points
# of data/ISC040N10NM8_Diode.csv (25/175 degC, 1.2-490 A). k_xti fixed at 3.0
# (LTspice default), i_sr = 0 (recombination branch unsensed -- freeing it
# buys nothing), k_nr = 2.0, t_nom = 25, a_ikf = 0, b_rdiode = 0.
# RMSE 0.514 % (25 C 0.445 %, 175 C 0.575 %).
#
# KEY FINDING: k_n = 2.0 (the ST8L65 fixed value) is WRONG here -- pinned at
# 2.0 the fit lands at 35.2 %; freeing it drops it to 0.514 % and k_n lands
# at 1.053 (near-ideal diffusion, unlike the SJ part). i_s <-> k_n corr
# +0.999, Jacobian rank 5/6: the ideality/prefactor split is a degenerate
# valley (i_s sigma_rel ~12) -- quote the PAIR/curve, never i_s alone.
# v_eg ~ 1.11 eV (silicon) and i_kf ~ 357 A (high-injection knee) are
# genuinely sensed; r_diode <-> i_kf corr +0.994 (quote the curve, not the
# split).
#
# TC: the LINEAR r_diode TC (a_rdiode, LTspice Trs1) is the decisive knob --
# TC-off 3.31 % -> 0.514 %, strongly identified (sigma_rel 0.013). The
# quadratic b_rdiode and the knee TC a_ikf buy NOTHING with only two
# temperatures and go degenerate (a_rdiode <-> b_rdiode corr -1.0) -- DO NOT
# free them (ST8L65 had three temperatures, which is what made Trs2
# identifiable there).
#
# r_diode = 0.961 mOhm is ~10x SMALLER than the recovery bench's old
# R1 = 10 mOhm. The full subckt (scripts/gen_full_subckt.py) reconciles them:
# the DC set goes on D1's model card (Rs = 0 there) and r_diode(T) with the
# Trs1 law goes in the recovery branch's R1, so the forward path carries the
# diode series R exactly ONCE. Recovery is R1-insensitive (measured over
# R1 = 10 u/1m/10m, out/recov_r1_sens.cir; the plateau is carried by B2,
# I(D1) ~ -1 nA -- session 23).
#
# tests/test_isc040_official.py pins this record against the measured CSV.
# ---------------------------------------------------------------------------

ISC040_DIODE_FREE: tuple[str, ...] = (
    "i_s", "k_n", "r_diode", "i_kf", "v_eg", "a_rdiode",
)

ISC040_DIODE_RMSE_REL = 0.00514
ISC040_DIODE_GROUP_RMSE_REL = {25: 0.004450, 175: 0.005751}

ISC040_DIODE = Params(
    i_s=3.767896211320614e-12,
    k_n=1.0528671662498787,
    i_sr=0.0,
    k_nr=2.0,
    i_kf=357.18447376523085,
    a_ikf=0.0,
    r_diode=0.0009613283429766385,
    a_rdiode=0.0010248339388128384,
    b_rdiode=0.0,
    t_nom=25.0,
    k_xti=3.0,
    v_eg=1.1096559091661649,
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


# ---------------------------------------------------------------------------
# Body-diode recovery (DiodeRecov, transit-time pole gauge) -- sessions
# 25a-25c. LTspice-side netlist parameters (dimension-normalized, outside
# the names.xlsx convention): knobs are tau / t_minor / t_pole; q_0 was
# DELETED (pure gauge -- observables depend only on t_minor and
# V_R/(q_0*tau)). Fit doctrine: intrinsic branch current I(x1:R2)
# (datasheet Qrr excludes Q_oss; terminal path is contaminated ~28 nC at
# 100 A/us), t_rr = 10 % of I_RM, no dTmax (the pole rounds the snap),
# no .ic (DC op point IS the pre-event state), solver=alt method=gear
# reltol=1e-4. Extraction: src/vdmos_fit/extract_recov.py (grid ->
# log-log surrogates -> multi-start TRF polish, diff_step=0.03).
# Report: out/recov_extract_isc040_run/extract_report.txt; records
# out/RECOV_POLE_FIT.md + out/RECOV_INTRINSIC_FIT.md.
#
# IDENTIFIABILITY (fix-and-refit penalties, never marginal sigma):
# t_pole sharply identified (x2/0.5 fix costs ~10-15 % RMS); tau is a
# FLAT valley coupled to both other knobs (interior optimum ~3e-4,
# analytic delay law d(Delta t)/d ln tau = -t_pole) -- quote the coupled
# set + RMS, never tau alone. t_minor dominates Qrr.
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class RecovRecord:
    """Official DiodeRecov recovery record (pole gauge).

    targets/achieved: (di/dt [A/us], Qrr [nC], t_rr [ns]) -- intrinsic
    path (q_rri/t_rri), t_rr = 10 % of I_RM. rms_rel: relative RMS over
    ALL datasheet numbers."""
    part: str
    tau: float          # B2 snap threshold (dimension-normalized)
    t_minor: float      # carrier storage lifetime [s] (= C1 farads, Roff=1)
    t_pole: float       # transit/sweep-out pole R3*C3 [s] (R3 = 1 ohm)
    k_chg: float        # forward-recovery branch (not identifiable here)
    r_max: float        # forward-recovery branch (not identifiable here)
    is_diode: float     # D1 DC branch (from the bench; diode.py protocol
    rs_diode: float     #  refines it separately -- do not refit on
    n_diode: float      #  recovery data)
    rms_rel: float
    targets: tuple
    achieved: tuple


ISC040_RECOV = RecovRecord(
    part="ISC040N10NM8",
    tau=3.7892757739719363e-4,
    t_minor=2.2363459877611823e-8,
    t_pole=8.32388859768649e-10,
    k_chg=5.0,
    r_max=0.0,
    is_diode=1e-6,
    rs_diode=10e-3,
    n_diode=1.4,
    rms_rel=0.011132,
    targets=((1000.0, 247.0, 25.0), (100.0, 36.0, 29.0)),
    achieved=((1000.0, 248.20, 24.687), (100.0, 35.710, 29.459)),
)
