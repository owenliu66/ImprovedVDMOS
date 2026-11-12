"""Parameterized FULL/FAST 4-terminal subckt export (session 30, Phase C).

Moved out of scripts/gen_full_subckt.py so the release orchestrator
(release.py) can export ANY part from fitted records, not just ISC040.
The ISC040 script is now a thin wrapper that supplies the ISC040 records
and its hand-written header text; the generated files stay byte-identical
(pinned by tests).

Variants (same subckt NAME -> drop-in .include swap; never both):

* FULL : body diode + DiodeRecov recovery network (pole gauge). Only when
  recovery knobs were extracted from datasheet Qrr/t_rr data.
* FAST : plain body diode. With recovery knobs available but the network
  dropped for speed, tt = t_minor + t_pole/0.69 maps the pole-gauge
  storage (user mapping, session 29). WITHOUT recovery data (release
  rule, user 2026-11-08): plain diode, NO tt -- standard SPICE.

B3 doctrine (user, session 30): the forward-recovery branch modulator B3
is only emitted when r_max != 0. At r_max = 0 the branch is R1 tied
straight to K (R2's 1 uOhm alone still breaks the bridge loop).

Channel expression: collapsed (ISC040 SJ-cap family: k_theta =
k_beta_trap = k_lambda = g_k = v_diblmax = k_kneeform = 0 -- byte-identical
to the pinned files) vs full family (those terms emitted). The full form
transcribes model.py _channel_current including the DIBL-shift soft cap
(v_diblmax/k_diblsharp) and the k_kneeform knee selector (mobility vs linear
scale) and the negative-k_theta divisor guard -- it has NO LTspice parity
bench yet for the non-collapsed family; the export report says so.
channel_expr_auto picks collapsed only when every one of those is off.
"""

from __future__ import annotations

from dataclasses import dataclass

from .params import Params


@dataclass
class RecovParams:
    """The five recovery knobs (ISC040_RECOV record is duck-type compatible)."""
    tau: float
    t_minor: float
    t_pole: float
    k_chg: float = 5.0
    r_max: float = 0.0


def _n(v: float) -> str:
    """SPICE engineering suffix for time-like values (s -> n)."""
    return f"{v * 1e9:.6g}n"


# Channel temperature coefficients (the TC law in stages.params_at):
#   k_p(T)   = k_p25 * ((T+273.15)/298.15)^(-beta_p)
#   v_th0(T) = v_th0_25 + s_vth*(T-25)
#   r_d(T)   = r_d25*(1 + a_rd*(T-25) + b_rd*(T-25)^2)
# When any is non-zero the export emits them as .params referencing LTspice's
# built-in TEMP (verified: TEMP tracks `.temp` in .params AND in B-source
# expressions), so the subckt responds to `.temp`/`.step temp`. With all zero
# the export stays byte-identical to the 25 C-only pinned files.
_TC_NAMES = ("s_vth", "beta_p", "a_rd", "b_rd")


def _tc_active(tc: dict | None) -> bool:
    return bool(tc) and any(float(tc.get(n, 0.0)) for n in _TC_NAMES)


def channel_expr(p: Params, kp_str: str | None = None,
                 vth0_str: str | None = None) -> str:
    """Collapsed (SJ-cap) family: k_theta = k_beta_trap = g_k = k_lambda = 0.

    ONE flat B-source expression, mirror-symmetric third quadrant (the
    forward expression at |V(dx,sx)| times a sign factor -- LTspice has
    no sign() function, nested if() is the exact mirror). Parity with
    model.py vdmos_id (session 29).

    DO NOT buffer the shared sub-expressions on helper nodes: measured
    negative (out/turnon_flat.cir vs turnon_buf.cir, identical physics) --
    the 8 extra voltage-source nodes cost the solver more than the
    duplicated evaluation: 6517 vs 9161 iterations, 0.29 vs 0.38 s.
    LTspice's expression evaluator is cheap; matrix/convergence overhead
    dominates. Softplus overflow-safe: log(1+e^z) = max(z,0)+log(1+e^-|z|).
    """
    vd, vg = _vd_vg()
    kp = kp_str if kp_str is not None else f"{p.k_p:.6g}"
    vth = _vth(vd, p, vth0_str)
    voveff = _vov(vd, vg, p, vth)
    vknee = f"max({voveff},{p.v_kneemin:.6g})"          # k_theta = 0
    ich = (f"{kp}*{voveff}*({_vdseff(vd, vknee, p)}"
           f"-0.5*{_vdseff(vd, vknee, p)}*{_vdseff(vd, vknee, p)}/{vknee})"
           f"/(1+{_vdseff(vd, vknee, p)}/{p.v_vsat:.6g})")
    return f"{_SGN}*{ich}"


def channel_expr_full(p: Params, kp_str: str | None = None,
                      vth0_str: str | None = None) -> str:
    """Full family: k_beta_trap gain, k_theta knee, k_lambda CLM, g_k kink.

    Direct transcription of model.py _channel_current (magnitudes inside,
    mirror sign outside). NOT yet parity-checked against LTspice for a
    non-collapsed fit -- the release report flags this.
    """
    vd, vg = _vd_vg()
    kp = kp_str if kp_str is not None else f"{p.k_p:.6g}"
    kdb = f"tanh({vd}/{p.v_dbsat:.6g})"
    vth = _vth(vd, p, vth0_str)
    vov = _vov(vd, vg, p, vth)
    if p.k_beta_trap:
        vov = f"({vov}*(1+{p.k_beta_trap:.6g}*{vov}))"   # v_ov_eff2
    # knee, mirroring model.py's k_kneeform selector
    if p.k_kneeform >= 0.5:
        vknee = f"max({p.k_theta:.6g}*{vov},{p.v_kneemin:.6g})"
    elif p.k_theta < 0.0:
        vknee = (f"max({vov}/max(1+{p.k_theta:.6g}*{vov},{_KNEE_FLOOR}),"
                 f"{p.v_kneemin:.6g})")
    elif p.k_theta:
        vknee = f"max({vov}/(1+{p.k_theta:.6g}*{vov}),{p.v_kneemin:.6g})"
    else:
        vknee = f"max({vov},{p.v_kneemin:.6g})"
    vds = _vdseff(vd, vknee, p)
    ich = (f"{kp}*{vov}*({vds}-0.5*{vds}*{vds}/{vknee})"
           f"/(1+{vds}/{p.v_vsat:.6g})")
    if p.k_lambda:
        ich = f"({ich}*(1+{p.k_lambda:.6g}*{vd}*{kdb}))"
    if p.g_k:
        kkink = (f"0.5*(1+tanh((V(G,sx)-{p.v_thk:.6g})/{p.v_kovsoft:.6g}))")
        ich = (f"({ich}+{p.g_k:.6g}*{vd}*{kkink}*{kdb})")
    return f"{_SGN}*{ich}"


_SGN = "if(V(dx,sx)<0,-1,if(V(dx,sx)>0,1,0))"

# Mirror of model.py _KNEE_DENOM_FLOOR: guards the mobility-knee divisor
# 1 + k_theta*v_ov_eff2 against a NEGATIVE k_theta (wider knee) driving it
# through zero. A no-op for k_theta >= 0 (divisor >= 1).
_KNEE_FLOOR = "1e-9"


def _vth(vd: str, p: Params, vth0_str: str | None = None) -> str:
    """DIBL-shifted threshold, mirroring model.py: the raw shift
    k_eta*v_ds*k_db is soft-capped at v_diblmax (the same power-law soft-min
    as the charge cap) so v_th >= v_th0 - v_diblmax for every v_ds. Without
    the cap a low-voltage fit extrapolates an unbounded linear DIBL that drives
    v_th negative at high v_ds, which pins the overdrive at the charge cap and
    freezes the saturation current (the exported-subckt 950 A failure).
    v_diblmax = 0 -> the raw shift, byte-identical to the pre-cap export.

    vth0_str is the threshold BASE expression: a literal when TC is off
    (byte-identical) or the TC-aware param name 'vth0t' when TC is emitted."""
    vth0 = vth0_str if vth0_str is not None else f"{p.v_th0:.6g}"
    shift = f"{p.k_eta:.6g}*{vd}*tanh({vd}/{p.v_dbsat:.6g})"
    if p.v_diblmax > 0.0:
        shift = (f"{shift}/pow(1+pow({shift}/{p.v_diblmax:.6g},"
                 f"{p.k_diblsharp:.6g}),1/{p.k_diblsharp:.6g})")
    return f"({vth0}-{shift})"


def channel_expr_auto(p: Params, tc: dict | None = None) -> str:
    collapsed = (p.k_theta == 0.0 and p.k_beta_trap == 0.0
                 and p.k_lambda == 0.0 and p.g_k == 0.0
                 and p.v_diblmax == 0.0 and p.k_kneeform == 0.0)
    if _tc_active(tc):
        kp_str, vth0_str = "kp", "vth0t"
    else:
        kp_str = vth0_str = None
    return (channel_expr(p, kp_str, vth0_str) if collapsed
            else channel_expr_full(p, kp_str, vth0_str))


def _vd_vg() -> tuple[str, str]:
    return ("abs(V(dx,sx))",   # mirror: forward expr at |v_ds|
            "V(G,sx)")         # internal gate-source (degeneration via sx)


def _vov(vd: str, vg: str, p: Params, vth: str) -> str:
    z = f"({vg}-{vth})/{p.v_sub:.6g}"
    vov = f"{p.v_sub:.6g}*(max({z},0)+log(1+exp(-abs({z}))))"
    if p.v_ovmax and p.k_ovsharp:
        vov = (f"{vov}/pow(1+pow({vov}/{p.v_ovmax:.6g},{p.k_ovsharp:.6g}),"
               f"1/{p.k_ovsharp:.6g})")
    return vov


def _vdseff(vd: str, vknee: str, p: Params) -> str:
    # mirror model.py: k_kneesharp = 0 is the legacy tanh soft-switch (the
    # power-law form below divides by k_kneesharp in the exponent, so it must
    # not be emitted at 0); k_kneesharp > 0 is the power-law soft-min.
    if p.k_kneesharp == 0.0:
        return f"({vknee}*tanh({vd}/{vknee}))"
    return (f"{vd}/pow(1+pow({vd}/{vknee},{p.k_kneesharp:.6g}),"
            f"1/{p.k_kneesharp:.6g})")


def device_common(p25: Params, caps_name: str,
                  tc: dict | None = None) -> str:
    """Channel + package R + caps delta block (shared by both variants).

    tc carries the channel temperature coefficients (s_vth/beta_p/a_rd/b_rd).
    When any is non-zero, k_p(T)/v_th0(T)/r_d(T) are emitted as .params that
    reference LTspice's built-in TEMP, so the subckt responds to `.temp`. With
    tc off/zero the block is byte-identical to the 25 C-only pinned files."""
    # LTspice rejects R=0 and a 0 V short in series with the channel current
    # source makes the DC matrix singular (node sx has no conductance to
    # ground). Stand in for the exact 0 with a 1 uOhm floor.
    rs_line = (f"Rs  sx S   {p25.r_s:.6g}" if p25.r_s > 0
               else "Rs  sx S   1e-6")
    rsk_line = (f"Rsk sx K   {p25.r_sk:.6g}" if p25.r_sk > 0
                else "Rsk S K    1e-6")
    if _tc_active(tc):
        dt = "(TEMP-25)"
        tc_block = (
            "* ---- channel temperature coefficients: k_p(T), v_th0(T), "
            "r_d(T) ----\n"
            "* reference LTspice's built-in TEMP (tracks `.temp`); 25 C base.\n"
            f".param kp={{{p25.k_p:.6g}*pow((TEMP+273.15)/298.15,"
            f"-{tc['beta_p']:.6g})}}\n"
            f".param vth0t={{{p25.v_th0:.6g}+({tc['s_vth']:.6g})*{dt}}}\n"
            f".param rdt={{{p25.r_d:.6g}*(1+{tc['a_rd']:.6g}*{dt}"
            f"+{tc['b_rd']:.6g}*{dt}*{dt})}}\n")
        chan = f"Bch dx sx I={{{channel_expr_auto(p25, tc)}}}"
        rd = "{rdt}"
    else:
        tc_block = ""
        chan = f"Bch dx sx I={channel_expr_auto(p25)}"
        rd = f"{p25.r_d:.6g}"
    return f"""* ---- intrinsic channel (model.py expression) ----
* ONE flat B-source expression; mirror-symmetric third quadrant (session 29).
* Buffered helper nodes were tried and measured WORSE (see channel_expr).
{tc_block}{chan}
* high-Z leak across the channel: silences the "node dx floating on a current
* source" initial-check warning; 100 MOhm changes i_d by <1e-6 relative.
Rleak dx sx 100Meg
* ---- package series R (channel path only; diode bypasses) ----
Rd  D  dx  {rd}
{rs_line}
{rsk_line}
* ---- capacitance delta: Cgs/Cgd at terminals ----
Xcaps G D S {caps_name}
"""


def build_full(part: str, p25: Params, dio: Params, rec: RecovParams,
               cds_q_expr: str, caps_text: str, header: str,
               tc: dict | None = None) -> str:
    """FULL: body diode + DiodeRecov. cds_q_expr = the bare 'Q=...' clause
    (the Cds charge-master lives INSIDE the recovery branch at Ki-A)."""
    xrec_params = " ".join([
        f"Is={dio.i_s:.6g}", f"Rs={dio.r_diode:.6g}", f"N={dio.k_n:.6g}",
        f"Ikf={dio.i_kf:.6g}", f"Eg={dio.v_eg:.6g}", f"Xti={dio.k_xti:g}",
        f"Tnom={dio.t_nom:g}", f"Trs1={dio.a_rdiode:.6g}",
        f"tau={rec.tau:.6g}", f"t_minor={_n(rec.t_minor)}",
        f"t_pole={_n(rec.t_pole)}", f"k_chg={rec.k_chg:g}",
        f"r_max={rec.r_max:g}",
    ])
    xrec = f"Xrec s d DiodeRecov params: {xrec_params}"

    diode_note = ("* ---- body diode + recovery: anode = source, cathode = "
                  "drain, at terminals ----")
    device = f""".subckt {part} G D S K
{device_common(p25, f"{part.lower()}_caps", tc)}{diode_note}
{xrec}
.ends {part}

"""

    if rec.r_max != 0.0:
        r_block = ("R1 N001 Ki {Rs*(1+Trs1*(temp-Tnom))}\n"
                   "R2 Ai A 1u\n"
                   "B3 N001 K V=if(r_max>0,-I(R2)*(r_max/(1+k_chg*V(q_j))"
                   "+1u),0)\n")
    else:
        # session 30 (user): B3 only with a non-zero r_max; at r_max = 0 the
        # branch is R1 tied straight to K (R2's 1 uOhm alone still breaks
        # the near-singular B/R bridge loop).
        r_block = ("R1 K Ki {Rs*(1+Trs1*(temp-Tnom))}\n"
                   "R2 Ai A 1u\n")

    recov = f""".subckt DiodeRecov A K
B2 Ki Ai I=max(0,V(q_j)*(V(Ki)-V(Ai))/tau)
C1 q_j0 0 {_n(rec.t_minor)} Rser=0
B1 0 q_j0 I=-I(R2)
D2 0 q_j0 ClampD
* R1 carries the diode's series R (instance param Rs = record r_diode)
* with its LTspice Trs1 law, Tnom-referenced -- verified against diode.py to
* ~1e-5 (tests/test_diode.py::test_ltspice_tc_laws). D1 itself has Rs=0 so
* the forward path carries the series R exactly ONCE.
{r_block}C2 Ki A {cds_q_expr}
D1 Ai Ki MyDiode
E1 N002 0 q_j0 0 1
R3 q_j N002 1
C3 q_j 0 {_n(rec.t_pole)}
.model ClampD D(Ron=1u,Roff=1,Vfwd=0)
.model MyDiode D(Is=Is,Rs=0,N=N,Ikf=Ikf,Eg=Eg,Xti=Xti,Tnom=Tnom)
.ends DiodeRecov
"""
    return header + device + caps_text + "\n" + recov


def build_fast(part: str, p25: Params, dio: Params, cds_line: str,
               caps_text: str, header: str,
               tt: float | None = None, tc: dict | None = None) -> str:
    """FAST: plain body diode; tt only when recovery knobs justify it.

    cds_line is the full instance line 'Cds D S Q=...' (moved to the D-S
    terminals -- no recovery branch here; the Q= expression is in the
    element voltage x, placement-agnostic).
    """
    tt_part = f",tt={_n(tt)}" if tt is not None else ""
    dio_card = (f".model MyDiode D(Is={dio.i_s:.6g},Rs={dio.r_diode:.6g},"
                f"N={dio.k_n:.6g},Ikf={dio.i_kf:.6g},Eg={dio.v_eg:.6g},"
                f"Xti={dio.k_xti:g},Tnom={dio.t_nom:g},"
                f"Trs1={dio.a_rdiode:.6g}{tt_part})")

    diode_note = ("* ---- body diode: plain D, series R in the card "
                  "(single count) ----")
    device = f""".subckt {part} G D S K
{device_common(p25, f"{part.lower()}_caps", tc)}{diode_note}
D1 s d MyDiode
{dio_card}
* ---- Cds charge-master at the D-S terminals (x = V(d,s) = v_ds) ----
{cds_line}
.ends {part}

"""
    return header + device + caps_text
