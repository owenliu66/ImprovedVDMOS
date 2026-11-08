"""Generate the FULL and FAST 4-terminal ISC040N10NM8 subckts.

Both variants assemble the same intrinsic channel, package series R and
capacitance delta, and share the subckt NAME ISC040N10NM8 -- the fast file
is a drop-in .include swap. Never .include both files in one circuit
(duplicate subckt name).

FULL (out/ISC040N10NM8_full.subckt): body-diode + minority-carrier recovery
network (DiodeRecov, pole gauge) at the terminals.

FAST (out/ISC040N10NM8_fast.subckt, session 29): NO recovery network. The
body diode is a plain D whose minority storage is set by the SPICE diode's
own charge control, tt = t_minor + t_pole/0.69 (user-specified mapping of
the pole-gauge knobs onto tt); its series R (r_diode + Trs1 law) lives in
the model card (single count). The Cds charge-master (which lives inside
the recovery branch in the FULL variant) moves to the D-S terminals.
Rationale: in a synchronous buck the xrec branch dominated LTspice's
convergence report (~10x slower than native VDMOS); the fast variant
removes that cost at the price of the exponential storage tail instead of
the pole-gauge snap -- first-order recovery, datasheet Qrr/t_rr not
guaranteed (the FULL variant owns those numbers).

CHANNEL (both variants, session 29):
  * mirror-symmetric third quadrant: i_ch = sign(V(dx,sx)) x the forward
    expression evaluated at |V(dx,sx)| -- synchronous-rectifier channel
    conduction, parity with the package model (model.py vdmos_id).
  * the channel stays ONE flat B-source expression: buffering the shared
    sub-expressions on helper nodes was tried and measured WORSE (8 extra
    voltage-source nodes: 6517 -> 9161 iterations on an identical-physics
    turn-on bench, out/turnon_flat.cir vs turnon_buf.cir). LTspice's
    expression evaluator is cheap; matrix/convergence overhead dominates.

Assembly:

  * the intrinsic channel as an LTspice behavioral current source (B-source)
    generated from the reference.py ISC040_CHANNEL record -- the model.py
    _channel_current expression with the fitted constants substituted, written
    in terms of the INTERNAL node voltages (the package R is wired explicitly
    as Rd/Rs/Rsk, so the B-source sees the degenerated v_ds'/v_gs' for free);
  * the package series R (r_d drain leg, r_s source leg, r_sk kelvin tap);
  * the capacitance delta (Cgs/Cgd only -- Cds lives in the recovery branch
    in FULL, at the D-S terminals in FAST);
  * FULL: the body-diode + recovery network (DiodeRecov, pole gauge);
    FAST: the plain body diode D1 with tt.

Terminals: D G S K. The channel branch carries the channel current through
the package R; the diode branch connects at the terminals and bypasses
r_d/r_s. For ISC040 r_s = r_sk = 0 (the static data rejects the gate-loop
degeneration -- see the record notes), so sx = S = K and the gate loop is
not degenerated; the K pin is kept for the driver return.

The channel B-source control is V(G,sx): with Rs between sx and S carrying the
channel current, V(G)-V(sx) = v_gs - r_s*i_d = v_gs' automatically (the
non-Kelvin degeneration, r_sk = 0). At r_s = 0 this is just v_gs.

    python scripts/gen_full_subckt.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import (  # noqa: E402
    ISC040_CHANNEL,
    ISC040_DIODE,
    ISC040_RECOV,
    Params,
)

OUT = ROOT / "out"


def _n(v: float) -> str:
    """SPICE engineering suffix for time-like values (s -> n)."""
    return f"{v * 1e9:.6g}n"


def channel_expr(p: Params) -> str:
    """LTspice B-source current expression for the intrinsic channel: ONE
    flat expression, mirror-symmetric third quadrant (the forward
    expression at |V(dx,sx)| times a sign factor -- LTspice has no sign()
    function, nested if() is the exact mirror). Parity with model.py
    vdmos_id (session 29).

    DO NOT buffer the shared sub-expressions on helper nodes: measured
    negative (out/turnon_flat.cir vs turnon_buf.cir, identical physics) --
    the 8 extra voltage-source nodes cost the solver more than the
    duplicated evaluation: 6517 vs 9161 iterations, 0.29 vs 0.38 s.
    LTspice's expression evaluator is cheap; matrix/convergence overhead
    dominates.

    Built from model.py _channel_current with the fitted constants
    substituted. k_theta = 0, k_beta_trap = 0, g_k = 0, k_lambda = 0 are
    the record's fixed values, so those terms collapse. Softplus is written
    overflow-safe: log(1+e^z) = max(z,0) + log(1+e^-|z|).
    """
    vd = "abs(V(dx,sx))"             # mirror: forward expr at |v_ds|
    vg = "V(G,sx)"                   # internal gate-source (degeneration via sx)
    kdb = f"tanh({vd}/{p.v_dbsat:.6g})"
    vth = f"({p.v_th0:.6g}-{p.k_eta:.6g}*{vd}*{kdb})"
    z = f"({vg}-{vth})/{p.v_sub:.6g}"
    vov = f"{p.v_sub:.6g}*(max({z},0)+log(1+exp(-abs({z}))))"
    voveff = (f"{vov}/pow(1+pow({vov}/{p.v_ovmax:.6g},{p.k_ovsharp:.6g}),"
              f"1/{p.k_ovsharp:.6g})")
    vknee = f"max({voveff},{p.v_kneemin:.6g})"          # k_theta = 0
    vdseff = (f"{vd}/pow(1+pow({vd}/{vknee},{p.k_kneesharp:.6g}),"
              f"1/{p.k_kneesharp:.6g})")
    ich = (f"{p.k_p:.6g}*{voveff}*({vdseff}-0.5*{vdseff}*{vdseff}/{vknee})"
           f"/(1+{vdseff}/{p.v_vsat:.6g})")
    sgn = "if(V(dx,sx)<0,-1,if(V(dx,sx)>0,1,0))"
    return f"{sgn}*{ich}"


def _device_common(p25: Params) -> str:
    """Channel + package R + caps delta block (shared by both variants)."""
    # LTspice rejects R=0 and a 0 V short in series with the channel current
    # source makes the DC matrix singular (node sx has no conductance to
    # ground). Stand in for the exact 0 with a 1 uOhm floor -- negligible vs
    # r_d (2.48 mO) and it gives sx a conductive path (clean DC solve).
    rs_line = (f"Rs  sx S   {p25.r_s:.6g}" if p25.r_s > 0
               else "Rs  sx S   1e-6")
    rsk_line = (f"Rsk sx K   {p25.r_sk:.6g}" if p25.r_sk > 0
                else "Rsk S K    1e-6")
    chan = f"Bch dx sx I={channel_expr(p25)}"
    return f"""* ---- intrinsic channel (model.py expression, ISC040_CHANNEL record) ----
* ONE flat B-source expression; mirror-symmetric third quadrant (session 29).
* Buffered helper nodes were tried and measured WORSE (see channel_expr).
{chan}
* high-Z leak across the channel: silences the "node dx floating on a current
* source" initial-check warning; 100 MOhm changes i_d by <1e-6 relative.
Rleak dx sx 100Meg
* ---- package series R (channel path only; diode bypasses) ----
Rd  D  dx  {p25.r_d:.6g}
{rs_line}
{rsk_line}
* ---- capacitance delta: Cgs/Cgd at terminals ----
Xcaps G D S isc040n10nm8_caps
"""


def build_full() -> str:
    p25 = ISC040_CHANNEL.params_at(25.0)
    rec = ISC040_RECOV
    dio = ISC040_DIODE
    cds = (ROOT / "data" / "ISC040N10NM8_cds_expr.txt").read_text(
        encoding="utf-8").strip()
    cds = "Q=" + cds.split("Q=", 1)[1]
    caps = (OUT / "ISC040N10NM8_caps.subckt").read_text(encoding="utf-8")

    header = f"""* ISC040N10NM8 -- FULL ImprovedVDMOS device subckt (4 terminals D G S K).
* GENERATED by scripts/gen_full_subckt.py from the reference.py records
* ISC040_CHANNEL (channel) + ISC040_DIODE (body diode DC) + ISC040_RECOV
* (recovery) + CAP_ISC040 (caps). Edit the records (or rerun the
* fit/extraction scripts), never hand-edit.
*
* ASSEMBLY:
*   D o--[Rd]--(dx)=== channel B-source ===(sx)--[Rs]--o S
*   |                                                  |
*   +---------------- DiodeRecov (D1 + recovery) ------+   anode = source
*   G o--[Cgs/Cgd via Xcaps to D,S terminals]
*   K = kelvin source tap (driver return); inert at r_s = r_sk = 0.
*
* CHANNEL (ISC040_CHANNEL, SJ charge-cap family, 25 C base + TC laws):
*   fit 25 C first (locked) then TC; 25 C RMSE 1.33 %, TC overall 2.97 %.
*   r_d = {p25.r_d:.6g} (whole current-loop series R; r_s = r_sk = 0 -- the
*   static data rejects the gate-loop degeneration, IMCQ session-22 finding).
*   TC laws (apply by re-evaluating the record at T; the CHANNEL branch is
*   25 C):
*     k_p(T)=k_p25*(T/298.15)^-2.2004, v_th0(T)=v_th0_25-7.486m*(T-25),
*     r_d(T)=r_d25*(1+5.621m*(T-25)+2.001u*(T-25)^2).
*   The channel is MIRROR-SYMMETRIC (session 29): sign(V(dx,sx)) times the
*   forward expression at |V(dx,sx)| -- third-quadrant synchronous-rectifier
*   conduction in PARALLEL with DiodeRecov (real-device topology). The
*   mirror cannot capture V_gs modulation of the diode forward drop
*   (accepted limitation; the channel takes over before it matters in Si).
*
* BODY DIODE (ISC040_DIODE record, scripts/fit_isc040_diode.py -- DC fit RMSE
*   0.51 % on 26 pts at 25/175 C): D1's model card carries Is/N/Ikf/Eg/Xti/
*   Tnom from the record (Rs=0 there); the diode's own series R lives in R1
*   with its LTspice Trs1 law, r_diode(T) = 0.961m*(1+1.0248m*(T-25)), so
*   the forward path carries the series R exactly ONCE and the branch is
*   temperature-complete (the diode.py TC laws are LTspice-verified to ~1e-5,
*   tests/test_diode.py). Recovery knobs (tau/t_minor/t_pole/k_chg/r_max)
*   come from ISC040_RECOV; recovery is insensitive to R1 (measured over
*   10u..10m, out/recov_r1_sens.cir) and to the D1 DC-set swap (verified,
*   out/recov_d1_check.cir).
*
* FAST VARIANT: ISC040N10NM8_fast.subckt (same subckt name, drop-in swap)
*   drops the recovery network for speed -- see that file's header.
*
* SOLVER (recovery snap): .options method=gear reltol=1e-4 solver=alt
"""

    xrec = "Xrec s d DiodeRecov params: " + " ".join([
        f"Is={dio.i_s:.6g}", f"Rs={dio.r_diode:.6g}", f"N={dio.k_n:.6g}",
        f"Ikf={dio.i_kf:.6g}", f"Eg={dio.v_eg:.6g}", f"Xti={dio.k_xti:g}",
        f"Tnom={dio.t_nom:g}", f"Trs1={dio.a_rdiode:.6g}",
        f"tau={rec.tau:.6g}", f"t_minor={_n(rec.t_minor)}",
        f"t_pole={_n(rec.t_pole)}", f"k_chg={rec.k_chg:g}",
        f"r_max={rec.r_max:g}",
    ])

    diode_note = ("* ---- body diode + recovery: anode = source, cathode = "
                  "drain, at terminals ----")
    device = f""".subckt ISC040N10NM8 G D S K
{_device_common(p25)}{diode_note}
{xrec}
.ends ISC040N10NM8

"""

    recov = f""".subckt DiodeRecov A K
B2 Ki Ai I=max(0,V(q_j)*(V(Ki)-V(Ai))/tau)
C1 q_j0 0 {_n(rec.t_minor)} Rser=0
B1 0 q_j0 I=-I(R2)
D2 0 q_j0 ClampD
* R1 carries the diode's series R (instance param Rs = ISC040_DIODE.r_diode)
* with its LTspice Trs1 law, Tnom-referenced -- verified against diode.py to
* ~1e-5 (tests/test_diode.py::test_ltspice_tc_laws). D1 itself has Rs=0 so
* the forward path carries the series R exactly ONCE (the old bench set put
* Rs=10m in R1 with a bench D1; the DC fit says 0.961 m -- recovery is
* R1-insensitive over 10u..10m, out/recov_r1_sens.cir).
R1 N001 Ki {{Rs*(1+Trs1*(temp-Tnom))}}
R2 Ai A 1u
B3 N001 K V=if(r_max>0,-I(R2)*(r_max/(1+k_chg*V(q_j))+1u),0)
C2 Ki A {cds}
D1 Ai Ki MyDiode
E1 N002 0 q_j0 0 1
R3 q_j N002 1
C3 q_j 0 {_n(rec.t_pole)}
.model ClampD D(Ron=1u,Roff=1,Vfwd=0)
.model MyDiode D(Is=Is,Rs=0,N=N,Ikf=Ikf,Eg=Eg,Xti=Xti,Tnom=Tnom)
.ends DiodeRecov
"""

    return header + device + caps + "\n" + recov


def build_fast() -> str:
    p25 = ISC040_CHANNEL.params_at(25.0)
    rec = ISC040_RECOV
    dio = ISC040_DIODE
    cds_line = (ROOT / "data" / "ISC040N10NM8_cds_expr.txt").read_text(
        encoding="utf-8").strip()
    caps = (OUT / "ISC040N10NM8_caps.subckt").read_text(encoding="utf-8")
    # User-specified mapping (session 29): pole-gauge storage knobs onto the
    # plain SPICE diode's storage lifetime. 0.69 ~ ln 2 (the pole delay).
    tt = rec.t_minor + rec.t_pole / 0.69
    dio_card = (f".model MyDiode D(Is={dio.i_s:.6g},Rs={dio.r_diode:.6g},"
                f"N={dio.k_n:.6g},Ikf={dio.i_kf:.6g},Eg={dio.v_eg:.6g},"
                f"Xti={dio.k_xti:g},Tnom={dio.t_nom:g},"
                f"Trs1={dio.a_rdiode:.6g},tt={_n(tt)})")

    header = f"""* ISC040N10NM8 -- FAST ImprovedVDMOS device subckt (4 terminals D G S K).
* GENERATED by scripts/gen_full_subckt.py (same subckt NAME as the FULL
* file -- drop-in .include swap; never .include both in one circuit).
*
* WHY: in a synchronous buck the DiodeRecov branch dominated LTspice's
* convergence report (device ~10x slower than native VDMOS, session 29).
* This variant deletes the recovery network and puts minority storage on
* the plain body diode's own charge control: tt = t_minor + t_pole/0.69
* (0.69 ~ ln 2) = {_n(tt)}, from the ISC040_RECOV pole-gauge knobs
* (t_minor {_n(rec.t_minor)}, t_pole {_n(rec.t_pole)}).
*
* WHAT YOU GET / LOSE vs FULL:
*   + no xrec branch: storage is the SPICE diode's exponential decay
*     (Qs = tt*I), solved by the junction itself -- no B2/B3 snap, no
*     storage-node network, far fewer convergence events.
*   - recovery SHAPE is the exponential storage tail, not the pole-gauge
*     snap: I_RM/t_rr/Qrr-vs-di/dt follow tt, NOT the fitted datasheet
*     match. The FULL variant owns the four datasheet recovery numbers
*     (RMS 1.12 %); treat this one as first-order storage for speed.
*   + forward diode DC is IDENTICAL physics: same Is/N/Ikf/Eg/Xti/Tnom,
*     series R counted ONCE in the model card (Rs + Trs1 law, no R1).
*   - Cds charge-master moves from the recovery branch (Ki-A) to the
*     D-S terminals (no recovery branch here).
*
* CHANNEL/package-R/caps blocks are IDENTICAL to FULL (mirror-symmetric
* third quadrant, single flat channel expression, single-count series R) --
* see ISC040N10NM8_full.subckt header for the record/TC documentation.
*
* SOLVER: .options method=gear reltol=1e-4 recommended (trapezoidal may
* ring on the diode turn-off); solver=alt no longer needed for the snap.
"""

    diode_note = ("* ---- body diode: plain D, storage via tt, series R in "
                  "the card (single count) ----")
    device = f""".subckt ISC040N10NM8 G D S K
{_device_common(p25)}{diode_note}
D1 s d MyDiode
{dio_card}
* ---- Cds charge-master at the D-S terminals (x = V(d,s) = v_ds) ----
{cds_line}
.ends ISC040N10NM8

"""

    return header + device + caps


def main() -> None:
    OUT.mkdir(exist_ok=True)
    for name, text in (("ISC040N10NM8_full.subckt", build_full()),
                       ("ISC040N10NM8_fast.subckt", build_fast())):
        path = OUT / name
        path.write_text(text, encoding="utf-8")
        print(f"wrote {path}  ({path.stat().st_size} bytes)")


if __name__ == "__main__":
    main()
