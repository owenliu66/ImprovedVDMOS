"""Parameter container, physical bounds and log-space packing for the fitter."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict, dataclass, fields

import numpy as np

__all__ = [
    "Params",
    "DEFAULT",
    "PARAM_NAMES",
    "LOG_PARAMS",
    "BOUNDS",
    "DEFAULT_FREE",
    "pack",
    "unpack",
]


@dataclass(frozen=True)
class Params:
    """Improved VDMOS drain-current parameters.

    Naming convention (README.md, "Naming convention"): ``<type prefix>_<suffix>``,
    prefixes v_ voltage, i_ current, r_ resistance, g_ conductance, k_ factor.
    Greek letters take the k_ prefix so no name collides with a Python keyword
    (k_lambda for lambda); the k_ factors are unitless EXCEPT the industry-
    standard set k_p [A/V^2] and the 1/V set k_lambda, k_eta, k_theta,
    k_beta_trap [1/V]. ``db`` means drain **bias**, not drain-to-body -- this
    model has no bulk terminal.

    Physical rationale for the charge cap (v_ovmax, k_ovsharp) and the package
    resistances (r_d, r_s) is in the project README.md ("Gate-overdrive charge
    cap", "Package resistance"); the notes here are units and gotchas only.

    Attributes
    ----------
    k_p : float
        Process transconductance parameter [A/V^2]. NOTE: not the Level-1 kp
        (28 -> ~34.7); must be re-extracted when k_theta/v_vsat are added.
    v_th0 : float  Zero-bias threshold voltage [V].
    v_sub : float  Subthreshold swing factor [V] (n*k*T/q); the softplus
                   smoothing scale, so the deep-subthreshold slope is 1/v_sub.
    v_ovmax : float  Soft cap (asymptotic ceiling) on the effective gate
                   overdrive [V]; 0 = off. "Off" is reached by leaving v_ovmax
                   out of the free set -- NOT by driving it to the lower bound
                   (v_ovmax -> 0+ means "hard cap at zero volts", the opposite
                   of off).
    k_ovsharp : float  Sharpness of the overdrive cap [-]; >> 1 approaches a
                   hard clip at v_ovmax. Only meaningful when v_ovmax > 0.
    k_beta_trap : float  Interface-trap-filling gain on the effective gate
                   overdrive [1/V] (WBG devices; names.xlsx "beta_trap"):
                   v_ov_eff2 = v_ov_eff*(1 + k_beta_trap*v_ov_eff). Models
                   4H-SiC channel-mobility recovery as the surface potential
                   deepens and interface traps fill -- the counterpart (same
                   knob family, OPPOSITE sign) of the super-junction charge
                   cap. v_ov_eff2 is THE effective overdrive downstream
                   (v_knee and the i_d prefactor both use it; k_theta
                   multiplies v_ov_eff2, not v_ov_eff). 0 = off, exact
                   identity -- "off" is reached by leaving it out of the free
                   set. MUST be >= 0: negative values make v(1 + k*v) non-
                   monotone (derivative 1 + 2*k*v crosses zero), breaking the
                   package-R bracket proof; model.py raises ValueError for
                   k_beta_trap < 0 (physics: trap filling only raises
                   mobility). Linear space (0 is legitimate, never logged);
                   BOUNDS lower bound is exactly 0.0 (user decision) -- safe
                   here because the gain is continuous in k_beta_trap and
                   EXACTLY the identity at 0 (r_d/r_s precedent, unlike the
                   k_kneesharp floor argument). Composition gotcha with the
                   cap active: the v_ov_eff2 asymptote is
                   v_ovmax*(1 + k_beta_trap*v_ovmax), NOT v_ovmax -- expect
                   k_beta_trap <-> v_ovmax anti-correlation if both are free.
    k_theta : float  Mobility-degradation coefficient [1/V] ->
                   v_knee = v_ov_eff2/(1 + k_theta*v_ov_eff2). May be NEGATIVE:
                   a negative k_theta WIDENS the knee (v_knee > v_ov_eff2),
                   which the ISC040 low-voltage fit wants. Safe lower limit is
                   the knee-divisor zero -1/v_ovmax (v_ov_eff2 < v_ovmax keeps
                   the divisor positive inside it); model.py floors the divisor
                   as a backstop. Linear space (not logged) for this reason.
                   NOTE: k_theta's meaning depends on k_kneeform (below).
    k_kneeform : float  Knee-position form selector [-]. 0 (default) = the
                   mobility-degradation form v_knee = v_ov_eff2/(1 + k_theta*v_ov_eff2)
                   with k_theta [1/V] the gate-exponent knob (alpha in (1,2],
                   identity at k_theta = 0). 1 = the linear knee-scale form
                   v_knee = k_theta*v_ov_eff2 with k_theta a dimensionless knee
                   MULTIPLIER (identity 1.0, no divisor/singularity), but then
                   alpha is pinned at 2 (square-law) and all flattening comes
                   from v_vsat. Selector only: 0 or 1, never fitted (kept out of
                   the free set); linear space (0 is a legitimate value).
    v_vsat : float   Velocity-saturation knee [V]; sets i_d,sat ~ v_ov**alpha,
                   alpha < 2. (The doubled v_ is deliberate: v_sat alone would
                   read as the drain-saturation voltage, which is v_knee.)
    k_kneesharp : float  Drain-knee softness exponent [-]; 0 = legacy tanh
                   soft-switch. > 0 replaces v_dseff = v_knee*tanh(v_ds/v_knee)
                   with the SAME power-law soft-min used for the charge cap:
                   v_dseff = v_ds/(1 + (v_ds/v_knee)**k_kneesharp)**(1/k_kneesharp).
                   Shape-only knob: dv_dseff/dv_ds(0) = 1 and v_dseff -> v_knee
                   hold for every value, so Ron and the i_d,sat plateau are
                   untouched while the knee WIDTH scales with it
                   (inf is a hard min, 1 is the harmonic soft-min; ~2.8
                   reproduces the legacy tanh). "Off" is the default 0 --
                   freeing it changes the knee shape, so it stays out of
                   DEFAULT_FREE like v_ovmax/k_ovsharp.
    k_eta : float    DIBL coefficient [1/V].
    v_diblmax : float  Soft cap (asymptotic ceiling) on the DIBL threshold
                    SHIFT k_eta*v_ds*k_db [V]; 0 = off (identity, the pre-cap
                    model bit for bit). The raw shift grows without bound in
                    v_ds, so a fit anchored on low-voltage data extrapolates an
                    aggressive linear DIBL that can drive v_th negative at high
                    v_ds (the ISC040 100 V failure mode); the cap keeps
                    v_th >= v_th0 - v_diblmax. Same power-law soft-min as
                    v_ovmax, same "off = leave it out of the free set" doctrine
                    (driving v_diblmax -> 0+ is a hard clip at zero shift, the
                    opposite of off).
    k_diblsharp : float  Sharpness of the DIBL-shift cap [-]; >> 1 approaches a
                    hard clip at v_diblmax. Only meaningful when v_diblmax > 0.
                    Same role as k_ovsharp for the charge cap.
    k_lambda : float Channel/drift-length modulation [1/V]. (spelled out --
                   ``lambda`` is a Python keyword.)
    g_k : float      Trench kink (parasitic accumulation channel) conductance
                   [S]. One-sided: contributes only for v_gs above v_thk.
    v_thk : float    Kink onset voltage [V].
    v_kovsoft : float  Kink onset softness [V] -- the transition width of the
                   kink sigmoid, the analogue of v_sub for the subthreshold
                   softplus. Bigger = softer.
    v_dbsat : float  Drain-**bias** saturation voltage [V]: the scale at which
                   the drain-bias factor k_db = tanh(v_ds/v_dbsat) saturates.
    v_kneemin : float  Numerical floor on v_knee [V]; keep fixed, not fitted.
    r_d : float      Drain-side series resistance [Ohm] (bond wires, leadframe,
                   pins, plus the undepleted drift). Embedded outside the
                   intrinsic channel; datasheet Rds(on) already contains
                   r_d + r_s, and r_d + r_s is a LOWER bound on the Rds(on)
                   chord (see README, "Package resistance (r_d, r_s)").
    r_s : float      Source series resistance [Ohm]. Also *degenerates* the
                   gate drive: the internal gate-source voltage is
                   v_gs - (r_s - r_sk)*i_d, so
                   g_m,eff = g_m/(1 + (r_s - r_sk)*g_m).
    r_sk : float     Kelvin source package resistance [Ohm] (names.xlsx
                   "sk"): the resistance between the KELVIN source tap and
                   the power source pin. The gate loop returns at the kelvin
                   tap, so it sees only the channel-side segment
                   r_s - r_sk; the FULL r_s stays in the current path
                   ("r_s inside the current path but outside the gate
                   loop"). r_sk = 0 -- no kelvin contact, the gate loop
                   returns at the power source (the legacy r_s semantics,
                   bit-identical). r_sk = r_s -- ideal kelvin tap at the
                   die/source metal: the gate loop is NOT degenerated at
                   all. Constraint 0 <= r_sk <= r_s (model.py raises
                   ValueError otherwise); the kelvin pin carries no channel
                   current in the DC model, so r_sk never appears in the
                   current-loop sum r_d + r_s. Linear space (0 is the
                   legitimate "no kelvin contact" value, r_d/r_s
                   precedent); NOT in DEFAULT_FREE (freeing it silently
                   changes existing fits).

    Body-diode parameters (LTspice diode DC branch; see diode.py and README,
    "Body diode"). Defaults are the OFF sentinel: i_s = i_sr = 0 makes the
    diode current exactly zero and vdmos_id() bit-identical to the pre-diode
    model. Beware near-collisions with channel names: i_s is the diode
    saturation current, NOT the drain current i_d; k_n is the diode emission
    coefficient, unrelated to v_sub; r_diode is the diode's own series
    resistance, separate from the package r_d/r_s (the diode bypasses them).

    i_s : float      Diffusion saturation current [A] (LTspice Is); 0 = off.
                   Only meaningful PAIRED with k_n -- the datasheet diode data
                   has no low-current tail, so i_s alone is weakly determined
                   (sigma_rel ~ 1 in fits).
    k_n : float      Diffusion emission coefficient [-] (LTspice N). Physical
                   diffusion range ~1-2; > 2 absorbs high-injection bending
                   into ideality and inflates the reverse floor.
    i_sr : float     Recombination saturation current [A] (LTspice Isr);
                   0 = branch off (default -- unsensed above 5 A).
    k_nr : float     Recombination emission coefficient [-] (LTspice Nr).
    i_kf : float     High-injection knee current [A] (LTspice Ikf); 0 = off
                   (LTspice Ikf = inf equivalent). Fitted regime is tens of A.
    a_ikf : float    Linear temperature coefficient of i_kf [1/degC] (LTspice
                   Tikf): i_kf(T) = i_kf*(1 + a_ikf*(T - t_nom)), exactly the
                   LTspice law (verified against LTspice 26.0.1 .meas probes).
                   0 = no temperature dependence (LTspice default) -- the
                   pre-TC model, bit-identical. Linear space: can be negative.
    r_diode : float  Diode series resistance [Ohm] (LTspice Rs, renamed to
                   avoid the r_s collision). On the diode path only -- the
                   diode connects to the device terminals and bypasses r_d/r_s.
    a_rdiode : float Linear temperature coefficient of r_diode [1/degC]
                   (LTspice Trs1): r_diode(T) = r_diode*(1 + a_rdiode*(T -
                   t_nom) + b_rdiode*(T - t_nom)**2). 0 = off (LTspice
                   default). NEAR-COLLISION WARNING: scales r_diode (the
                   DIODE's Rs), NOT the channel package resistance r_s.
    b_rdiode : float Quadratic temperature coefficient of r_diode [1/degC^2]
                   (LTspice Trs2), completing the Trs1/Trs2 law above (the
                   exact LTspice form, verified against LTspice 26.0.1
                   .meas probes). 0 = pure linear law -- bit-identical to
                   the pre-Trs2 model. Linear space: can be negative.
    t_nom : float    Temperature [degC] at which i_s/i_sr are defined
                   (LTspice Tnom). Fixed, never fitted.
    k_xti : float    Saturation-current temperature exponent [-] (LTspice Xti).
                   Divided by k_n / k_nr per branch in the temperature law.
    v_eg : float     Activation energy [eV, numerically volts] (LTspice Eg);
                   silicon ~1.11. Sensed by multi-temperature diode data.
    """

    k_p: float = 34.7
    v_th0: float = 3.0
    v_sub: float = 0.05
    v_ovmax: float = 0.0
    k_ovsharp: float = 1.0
    k_beta_trap: float = 0.0
    k_theta: float = 0.015
    k_kneeform: float = 0.0
    v_vsat: float = 1.5
    k_kneesharp: float = 0.0
    k_eta: float = 0.02
    v_diblmax: float = 0.0
    k_diblsharp: float = 1.0
    k_lambda: float = 0.02
    g_k: float = 0.4
    v_thk: float = 1.6
    v_kovsoft: float = 0.1
    v_dbsat: float = 10.0
    v_kneemin: float = 1e-3
    r_d: float = 0.0
    r_s: float = 0.0
    r_sk: float = 0.0
    # Body diode (all defaults = off sentinel; see diode.py docstring).
    i_s: float = 0.0
    k_n: float = 1.0
    i_sr: float = 0.0
    k_nr: float = 2.0
    i_kf: float = 0.0
    a_ikf: float = 0.0
    r_diode: float = 0.0
    a_rdiode: float = 0.0
    b_rdiode: float = 0.0
    t_nom: float = 25.0
    k_xti: float = 3.0
    v_eg: float = 1.11

    def as_dict(self) -> dict:
        return asdict(self)

    def replace(self, **kw) -> Params:
        return Params(**{**asdict(self), **kw})


DEFAULT = Params()

PARAM_NAMES: tuple[str, ...] = tuple(f.name for f in fields(Params))

# Parameters fitted in log space (strictly positive, span decades).
# v_ovmax is NOT here: 0 is a legitimate value ("cap off") and cannot be logged.
# v_diblmax likewise: 0 is the legitimate "DIBL cap off" identity (same doctrine
# as v_ovmax). k_kneesharp likewise: 0 is the legitimate "legacy tanh" value. k_beta_trap
# likewise: 0 is the legitimate "gain off" identity (and negative values are
# rejected outright by model.py, so the parameter is not strictly positive).
# k_theta likewise: it may be NEGATIVE (a negative k_theta WIDENS the knee,
# v_knee > v_ov_eff2 -- the ISC040 low-voltage fit wants a wider knee than
# k_theta >= 0 allows), so it is not strictly positive and cannot be logged.
# Its search lower bound is the knee-divisor zero -1/v_ovmax (see BOUNDS); the
# model floors the divisor as a numerical backstop (model._KNEE_DENOM_FLOOR).
# Diode saturation currents are in LOG_PARAMS even though 0 is the off
# sentinel: "off" is reached by leaving them out of the free set (same
# doctrine as v_ovmax/k_kneesharp), and a FREE i_s must never wander near 0
# -- its lower bound (1e-22) is far from the sentinel's neighbourhood.
LOG_PARAMS: frozenset[str] = frozenset(
    {
        "k_p", "v_sub", "v_vsat", "k_eta", "k_lambda",
        "g_k", "v_kovsoft", "v_dbsat", "k_ovsharp", "k_diblsharp",
        "i_s", "k_n", "i_sr", "k_nr", "i_kf", "k_xti", "v_eg",
    }
)

# Physically sensible search bounds, same units as Params.
# Rationale for the non-obvious entries: README, "Bounds and fitting-space
# decisions".
BOUNDS: dict[str, tuple[float, float]] = {
    "k_p": (1e-2, 1e4),
    "v_th0": (0.0, 10.0),
    "v_sub": (5e-3, 1.0),
    # Linear space; lower bound deliberately well above 0 (a free v_ovmax must
    # never enter the "hard clip at ~0 V" region).
    "v_ovmax": (0.05, 100.0),
    "k_ovsharp": (0.5, 40.0),
    # Linear space; lower bound EXACTLY 0.0 (user decision: k_beta_trap >= 0).
    # No sentinel-mimicry hazard unlike the k_kneesharp floor argument: the
    # gain is continuous in k_beta_trap and exactly the identity at 0 (same
    # precedent as r_d/r_s), and model.py rejects k_beta_trap < 0 outright
    # (non-monotone v(1 + k*v) kills the package-R bracket proof). Upper
    # bound 2.0 is a search guard -- the measured physical scale is ~0.06 1/V
    # (IMCQ120R004M2H, 25/175 degC).
    "k_beta_trap": (0.0, 2.0),
    # Linear space (k_theta may be negative). The PHYSICALLY safe lower bound is
    # the knee-divisor zero -1/v_ovmax (below it, 1 + k_theta*v_ov_eff2 <= 0 at
    # the maximum overdrive and v_knee blows up / flips sign). That bound is
    # dynamic (depends on the free v_ovmax), so the static table carries a
    # permissive -1.0 guard; the GUI and any fit that frees k_theta should tighten
    # the lower bound to -1/v_ovmax, and model.py floors the divisor as a
    # numerical backstop. Upper bound 5.0 unchanged.
    "k_theta": (-1.0, 5.0),
    # Knee-form selector: 0 = mobility form (default), 1 = linear knee-scale
    # form. Never fitted (kept out of the free set); a discrete switch, so the
    # box is just [0, 1] and it is linear space (0 is a legitimate value).
    "k_kneeform": (0.0, 1.0),
    "v_vsat": (0.05, 10000.0),
    # Linear space; lower bound above 0 for the same reason as v_ovmax -- a
    # free k_kneesharp must never enter the region near the sentinel 0 ("legacy
    # tanh") where a fitted value and the sentinel are numerically
    # indistinguishable. Relaxed 0.5 -> 0.1 (session 10): the fit-B relaxed-
    # bound probe landed at 0.47, just under the old guard, so the data wants a
    # knee softer than the old floor allowed.
    "k_kneesharp": (0.1, 40.0),
    "k_eta": (1e-5, 1.0),
    # Linear space; lower bound deliberately well above 0 (a free v_diblmax must
    # never enter the "hard clip at ~0 shift" region, which is the opposite of
    # off -- same v_ovmax doctrine). Upper bound 20 V is a loose search guard:
    # any real threshold shift is a few volts, so a cap this high never binds and
    # the fit collapses to the uncapped model.
    "v_diblmax": (0.05, 20.0),
    "k_diblsharp": (0.5, 40.0),
    "k_lambda": (1e-12, 1.0),
    "g_k": (1e-12, 100.0),
    "v_thk": (-5.0, 10.0),
    "v_kovsoft": (1e-3, 5.0),
    "v_dbsat": (0.1, 1300.0),
    "v_kneemin": (1e-6, 1e-1),
    # Linear space (the physical answer can be exactly zero); 0.5 Ohm is a
    # loose global ceiling only -- per device cap at datasheet Rds(on).
    # r_sk: linear space, 0 is the legitimate "no kelvin contact" value
    # (r_d/r_s precedent). The physical constraint r_sk <= r_s is enforced
    # by model.py (ValueError), not by the bound table -- a bound cannot
    # reference another parameter.
    "r_d": (0.0, 0.5),
    "r_s": (0.0, 0.5),
    "r_sk": (0.0, 0.5),
    # Body diode (LTspice D-model DC branch). i_s/i_sr lower bounds sit far
    # from the 0 sentinel so a free value can never mimic "off"; i_kf lower
    # bound 1e-3 A keeps a free i_kf away from the 0 sentinel ("off" = leave
    # it out of the free set). r_diode linear like r_d/r_s (can be exactly 0);
    # 0.5 Ohm global ceiling -- per device cap (scripts/fit_st8l65_diode.py
    # DEVICE_BOUNDS: the diode dV/dI slope bounds it at ~20 mOhm on ST8L65).
    "i_s": (1e-22, 1e-3),
    "k_n": (1.0, 4.0),
    "i_sr": (1e-22, 1e-3),
    "k_nr": (1.0, 4.0),
    "i_kf": (1e-3, 1e6),
    # Linear TCs (LTspice Tikf/Trs1/Trs2): can be negative, 0 is the legitimate
    # "no temperature dependence" default -- linear space, never logged. The
    # bounds keep the TC factor 1 + a*dT + b*dT^2 POSITIVE across the device
    # range (-55..150 C => |dT| <= 125 K) for the LINEAR part at the lower
    # bound (1 - 5e-3*125 = 0.375); a near the UPPER bound or any b can still
    # drive the factor negative at the range edges -- diode.py clamps it at a
    # tiny floor in that case (never hit by physical fits, only by absurd
    # probe excursions). Physical scale: silicon/metal resistance TCs are
    # ~+0.4..1 %/K linear, quadratic terms ~1e-6..1e-5 /K^2; the ceilings are
    # deliberately generous search guards, not physics.
    "a_ikf": (-5e-3, 5e-2),
    "r_diode": (0.0, 0.5),
    "a_rdiode": (-5e-3, 5e-2),
    "b_rdiode": (-1e-4, 1e-3),
    "k_xti": (0.5, 10.0),
    "v_eg": (0.1, 2.0),
}

# Sensible starting free set: everything that shapes the static I-V curves.
# v_ovmax/k_ovsharp are deliberately absent -- freeing them would silently
# change every existing fit and anchor (see README, "Bounds and fitting-space
# decisions"). k_beta_trap likewise: default 0.0 is the identity, freeing it
# silently changes existing fits (and on silicon SJ parts it anti-correlates
# with v_ovmax -- README, "Interface-trap-filling gain").
DEFAULT_FREE: tuple[str, ...] = (
    "k_p",
    "v_th0",
    "k_theta",
    "v_vsat",
    "k_eta",
    "k_lambda",
    "g_k",
    "v_thk",
    "r_d",
    "r_s",
)


def pack(params: Params, free: Sequence[str]) -> np.ndarray:
    """Free parameters -> flat optimisation vector (log space where appropriate)."""
    out = []
    for name in free:
        val = getattr(params, name)
        out.append(np.log(val) if name in LOG_PARAMS else val)
    return np.asarray(out, dtype=float)


def unpack(x: np.ndarray, free: Sequence[str], base: Params = DEFAULT) -> Params:
    """Flat optimisation vector -> full :class:`Params` (non-free values from `base`)."""
    if len(x) != len(free):
        raise ValueError(f"expected {len(free)} values, got {len(x)}")
    values = asdict(base)
    for name, xi in zip(free, x, strict=True):
        values[name] = float(np.exp(xi)) if name in LOG_PARAMS else float(xi)
    return Params(**values)


def vector_bounds(
    free: Sequence[str],
    overrides: dict[str, tuple[float, float]] | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Bounds for the flat vector, mapped through the log transform.

    `overrides` replaces individual entries of :data:`BOUNDS` for one call
    without mutating the module-level table (e.g. a per-device cap of r_d/r_s
    at the datasheet Rds(on)).
    """
    table = {**BOUNDS, **(overrides or {})}
    lo, hi = [], []
    for name in free:
        a, b = table[name]
        if name in LOG_PARAMS:
            lo.append(np.log(a))
            hi.append(np.log(b))
        else:
            lo.append(a)
            hi.append(b)
    return np.asarray(lo), np.asarray(hi)
