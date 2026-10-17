"""Body-diode current: LTspice diode equation, verified against LTspice 26.0.1.

    i_di = i_diffuse + i_recomb                    [A, positive anode -> cathode]

    v_di      = v_sd - r_diode * i_di              implicit (diode's own series R)
    i_diffuse = i_s(T) * (exp(v_di/(k_n*v_therm)) - 1)
                / sqrt(1 + i_s(T)/i_kf * (exp(v_di/(k_n*v_therm)) - 1))
    i_recomb  = i_sr(T) * (exp(v_di/(k_nr*v_therm)) - 1)
    v_therm   = k_B * T / q                        T = t_c + 273.15 [K]

    i_s(T)  = i_s  * (T/T_nom)**(k_xti/k_n)  * exp(v_eg*q/(k_B*k_n)  * (1/T_nom - 1/T))
    i_sr(T) = i_sr * (T/T_nom)**(k_xti/k_nr) * exp(v_eg*q/(k_B*k_nr) * (1/T_nom - 1/T))
    i_kf(T)    = i_kf    * (1 + a_ikf    * (T - T_nom))     [LTspice Tikf law]
    r_diode(T) = r_diode * (1 + a_rdiode * (T - T_nom)
                                  + b_rdiode * (T - T_nom)**2) [Trs1 + Trs2]

LTspice model-card mapping (DC branch only -- capacitance, transit time,
breakdown BV/Ibv and noise are out of scope):

    Is -> i_s   N -> k_n   Isr -> i_sr   Nr -> k_nr   Ikf -> i_kf
    Rs -> r_diode   Tnom -> t_nom   Xti -> k_xti   Eg -> v_eg
    Tikf -> a_ikf   Trs1 -> a_rdiode   Trs2 -> b_rdiode   (TCs per-degC /
    per-degC^2, Tnom-referenced; the Trs1+Trs2 pair is the complete LTspice
    Rs law)

LTspice's instance multiplicity (``m=``) is deliberately NOT modelled: it is a
wiring-level parallel-device count, not a model parameter (and the model-card
``M`` is the junction grading coefficient, a capacitance term). ``M = 1``.

Verified against LTspice 26.0.1 ``.dc``/``.step temp``/``.meas`` probes
(tests/test_diode.py::test_ltspice_cross_validation):

* The high-injection divisor above matches LTspice to ~1e-5 relative across
  3.4 A .. 375 kA.
* The saturation-current temperature law uses **N** (not Berkeley's NTT) as
  the divisor of both k_xti and v_eg; the recombination branch uses **Nr**.
* LTspice's recombination branch deviates from its own documented
  ``Isr*(exp(V/NrVt)-1)`` by ~1 % at 0.02 V growing to ~17 % at 0.3 V
  (V-dependent, exactly linear in Isr -- an implementation quirk). We
  implement the *documented* form; the cross-validation test therefore
  exercises the diffusion + Ikf + temperature terms only. Our data lives at
  0.56-1.05 V where i_sr is unsensed anyway.
* The TC laws are exactly LTspice's, referenced to **Tnom** (not to 27 degC):
  Tikf linear ``i_kf*(1 + TC*(T - Tnom))`` and the Trs1+Trs2 quadratic
  ``Rs*(1 + Trs1*(T - Tnom) + Trs2*(T - Tnom)**2)``. Probes (Tikf=0.006
  Trs1=0.008 Trs2=2e-5) match to ~1e-5 relative at -55/25/100/150 degC and
  0.6/0.9 V; a 27 degC-referenced variant misses by 0.1-0.8 %, and dropping
  the Trs2 term misses by up to 10.6 % at 150 degC (linear-only) -- the
  curvature is real in the law. (tests/test_diode.py::test_ltspice_tc_laws)

Sentinels (same doctrine as v_ovmax = 0 / k_kneesharp = 0):

* ``i_s = 0 and i_sr = 0`` -> diode current is exactly zero: the model is the
  exact pre-diode channel model, bit for bit. "Off" is reached by leaving
  i_s/i_sr out of the free set, never by driving them to the lower bound.
* ``i_kf = 0`` -> high-injection divisor off (equivalent to LTspice
  ``Ikf = inf``). The divisor is applied to the forward branch only; reverse
  saturation is physically unaffected by high-level injection.
* ``a_ikf = a_rdiode = b_rdiode = 0`` (the defaults, and LTspice's own
  defaults) -> no temperature dependence: bit-identical to the pre-TC model
  (b_rdiode = 0 alone reduces the Rs law to pure Trs1). The TCs act on
  i_kf/r_diode only -- the i_s(T)/i_sr(T) Xti/Eg law is a separate mechanism.

Temperature handling: ``t_nom`` (Params) is the temperature i_s/i_sr were
extracted at; ``t_c`` (function argument) is the device operating temperature.
Both in degC (LTspice ``Tnom``/``.temp`` convention). The future ``Tj`` node
supersedes ``t_c``.
"""

from __future__ import annotations

import numpy as np

from .params import DEFAULT, Params

__all__ = ["vdmos_diode"]

# Physical constants (SI). q/k_B = 11604.518 K is the activation-law scale.
_Q = 1.602176634e-19
_K_B = 1.380649e-23
_Q_OVER_K = _Q / _K_B  # [K/V]

# Exponent clip shared with the channel softplus: exp(700) ~ 1e304 stays in
# float64 and the clipped value is far beyond any physical diode current.
_EXP_CLIP = 700.0

# Floor for the TC polynomial 1 + a*(T - Tnom) + b*(T - Tnom)**2: it
# extrapolates below zero far outside any physical range (e.g. Trs1 = 0.008
# at -150 C -> factor -0.4); the clamp keeps r_diode/i_kf strictly positive so
# the series-R bracket proof (dg/di <= -1) and the log-space i_kf evaluation
# survive absurd probe temperatures. Never active for physical TC values over
# the device range (-55..150 C).
_TC_FLOOR = 1e-6


def _thermal_voltage(t_k: np.ndarray) -> np.ndarray:
    """v_therm = k_B*T/q [V], T in kelvin."""
    return _K_B * t_k / _Q


def _tc_factor(a: float, t_k: np.ndarray, t_nom_k: float,
               b: float = 0.0) -> np.ndarray:
    """LTspice TC multiplier 1 + a*(T - Tnom) + b*(T - Tnom)^2.

    a is per-degC (LTspice Tikf/Trs1 units), b per-degC^2 (Trs2; 0 = pure
    linear law). The kelvin DIFFERENCE equals the degC difference, so
    t_k - t_nom_k is the right argument. Clamped at _TC_FLOOR -- see the
    constant's comment.
    """
    dt = t_k - t_nom_k
    return np.maximum(1.0 + a * dt + b * dt * dt, _TC_FLOOR)


def _scaled_saturation(i_0: float, n: float, p: Params, t_k: np.ndarray) -> np.ndarray:
    """LTspice saturation-current temperature law for one branch.

    i_0 * (T/T_nom)**(k_xti/n) * exp(v_eg*q/k_B/n * (1/T_nom - 1/T)), with the
    divisor n = k_n for the diffusion branch and n = k_nr for recombination
    (verified against LTspice 26.0.1 -- NOT Berkeley's NTT).
    """
    t_nom_k = p.t_nom + 273.15
    return i_0 * (t_k / t_nom_k) ** (p.k_xti / n) * np.exp(
        p.v_eg * _Q_OVER_K / n * (1.0 / t_nom_k - 1.0 / t_k)
    )


def _diode_current(v_di, p: Params, t_k: np.ndarray) -> np.ndarray:
    """Explicit diode current [A] at the *internal* diode voltage (no r_diode).

    Overflow-safe: the forward diffusion branch is evaluated through log space
    so the high-injection divisor cannot overflow even at extreme bound
    combinations (v_eg/k_xti at their ceilings). The reverse branch (v_di <= 0)
    is evaluated directly -- |i| <= i_s(T), no overflow possible -- and is
    exactly 0 at v_di = 0, which is what the series-R bracket proof needs.
    """
    v_therm = _thermal_voltage(t_k)

    # --- Diffusion branch (LTspice Is/N/Ikf) -------------------------------
    if p.i_s > 0.0:
        i_s_t = _scaled_saturation(p.i_s, p.k_n, p, t_k)
        x = np.clip(v_di / (p.k_n * v_therm), -_EXP_CLIP, _EXP_CLIP)
        pos = x > 0.0
        with np.errstate(divide="ignore", invalid="ignore"):
            # Reverse / zero side: direct, |i| <= i_s(T) so no overflow, and
            # expm1(0) = 0 exactly (the series-R bracket needs f(0) = 0).
            i_neg = i_s_t * np.expm1(np.minimum(x, 0.0))
            # Forward side in log space: log(expm1(x)) = x + log(1 - exp(-x))
            # for x > 0 (nan for x <= 0, discarded by the where below).
            log_j = np.log(i_s_t) + x + np.log1p(-np.exp(-x))
            if p.i_kf > 0.0:
                # LTspice Tikf law: the knee current itself moves with T.
                i_kf_t = p.i_kf * _tc_factor(p.a_ikf, t_k, p.t_nom + 273.15)
                # I = J / sqrt(1 + J/Ikf) = exp(logJ - log1p(exp(logJ - logIkf))/2)
                i_pos = np.exp(log_j - 0.5 * np.logaddexp(0.0, log_j - np.log(i_kf_t)))
            else:
                i_pos = np.exp(log_j)
        i_diffuse = np.where(pos, i_pos, i_neg)
    else:
        i_diffuse = 0.0

    # --- Recombination branch (documented LTspice form; see module docstring
    #     for the measured deviation of LTspice's own implementation at low V)
    if p.i_sr > 0.0:
        i_sr_t = _scaled_saturation(p.i_sr, p.k_nr, p, t_k)
        x_r = np.clip(v_di / (p.k_nr * v_therm), -_EXP_CLIP, _EXP_CLIP)
        i_recomb = i_sr_t * np.expm1(x_r)
    else:
        i_recomb = 0.0

    return i_diffuse + i_recomb


def _solve_diode_series_r(
    v_sd: np.ndarray,
    r_diode: float,
    p: Params,
    t_k: np.ndarray,
    *,
    atol: float = 1e-12,
    rtol: float = 1e-10,
    max_iter: int = 60,
) -> np.ndarray:
    """Solve i_di = f(v_sd - r(T)*i_di) for the terminal diode current.

    r(T) = r_diode * (1 + a_rdiode*(T - T_nom) + b_rdiode*(T - T_nom)**2) is
    the LTspice Trs1+Trs2 law, evaluated per point from the broadcast t_k (so
    multi-temperature datasets get their own r_diode at each point). Bracketed
    safeguarded
    Newton on i_di in [lo, hi] with
    lo = min(0, v_sd/r), hi = max(0, v_sd/r): the terminal
    conductance is never below 1/r, and g(i) = f(v_sd - r*i) - i is
    strictly decreasing (df/dv >= 0, monotone diode) with f(0) = 0 exactly --
    the same bracket proof as the package-R solver. r_diode = 0 takes the
    explicit fast path. Raises RuntimeError on non-convergence; never returns
    an unconverged current.
    """
    if r_diode < 0.0:
        raise ValueError(f"r_diode must be >= 0 (got {r_diode})")
    if r_diode == 0.0:
        return _diode_current(v_sd, p, t_k)

    r = r_diode * _tc_factor(p.a_rdiode, t_k, p.t_nom + 273.15, p.b_rdiode)
    lo = np.minimum(0.0, v_sd / r)
    hi = np.maximum(0.0, v_sd / r)
    i_di = 0.5 * (lo + hi)
    stall = np.zeros(v_sd.shape, dtype=int)
    prev_width = None

    g = None
    for _ in range(max_iter):
        v_di = v_sd - r * i_di
        g = _diode_current(v_di, p, t_k) - i_di

        # g > 0 -> root above i_di; g < 0 -> root below (g strictly decreasing).
        lo = np.where(g > 0.0, i_di, lo)
        hi = np.where(g < 0.0, i_di, hi)

        width = hi - lo
        if prev_width is not None:
            stall = np.where(width > 0.9 * prev_width, stall + 1, 0)
        prev_width = width

        done = np.abs(g) <= atol + rtol * np.abs(i_di)
        if done.all():
            return i_di

        # Central-difference df/dv_di (the diode equation is C-infinity).
        hv = 1e-6 * np.maximum(1.0, np.abs(v_di))
        df = (
            _diode_current(v_di + hv, p, t_k) - _diode_current(v_di - hv, p, t_k)
        ) / (2.0 * hv)
        dg = -r * df - 1.0  # <= -1 always (r >= _TC_FLOOR*r_diode > 0)

        cand = i_di - g / dg
        inside = (cand > lo) & (cand < hi)
        moved = np.abs(cand - i_di) > 1e-12 * (hi - lo)
        cand = np.where(inside & moved & (stall < 3), cand, 0.5 * (lo + hi))
        i_di = np.where(done, i_di, cand)

    raise RuntimeError(
        f"diode r_diode solve failed after {max_iter} iterations: "
        f"max |residual| = {np.abs(g).max():.3g} A"
    )


def vdmos_diode(v_sd, p: Params = DEFAULT, *, t_c: float = 25.0) -> np.ndarray:
    """Body-diode forward current [A] at terminal source-drain voltage v_sd.

    Sign convention: ``v_sd = -v_ds`` (anode = source, cathode = drain for an
    N-channel FET); positive current flows source -> drain. Includes the
    diode's own series resistance r_diode (implicit solve; its LTspice Trs1
    + Trs2 TCs ``a_rdiode``/``b_rdiode`` are evaluated at ``t_c`` -- likewise
    ``a_ikf`` moves
    the high-injection knee i_kf(T)); the package
    r_d/r_s are NOT part of the diode path (they carry channel current only).

    Parameters
    ----------
    v_sd : float or array-like
        Source-drain terminal voltage [V] (forward-biased when positive).
    p : Params
        Model parameter set; diode parameters i_s, k_n, i_sr, k_nr, i_kf,
        r_diode, t_nom, k_xti, v_eg (see :mod:`vdmos_fit.params`).
    t_c : float or array-like
        Device temperature [degC] (LTspice ``.temp`` analogue); broadcast to
        v_sd, so per-point temperatures (multi-temperature datasets) work.
    """
    v_sd = np.asarray(v_sd, dtype=float)
    t_k = np.asarray(t_c, dtype=float) + 273.15
    v_sd, t_k = np.broadcast_arrays(v_sd, t_k)
    return _solve_diode_series_r(v_sd, p.r_diode, p, t_k)
