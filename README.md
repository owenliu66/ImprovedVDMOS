# Improved VDMOS — continuous drain-current model + parameter fitter

A compact power-MOSFET model that replaces the built-in LTspice Level-1 VDMOS
two-branch (triode/saturation) drain current with **one C¹ continuous
expression**, plus least-squares parameter extraction from static I-V data,
an LTspice-exact body-diode DC branch, a charge-master capacitance model, and
a behavioral body-diode recovery network for LTspice export.

Session history, official fitted records, and open problems live in
[`src/vdmos_fit/STATUS.md`](src/vdmos_fit/STATUS.md) and
[`HANDOVER.md`](HANDOVER.md). This file is the model reference.

## Why

LTspice's Level-1 VDMOS ties everything together through one shape knob:

| mtriode | knee = v_ov/m | Ron conductance = m·Kp·v_ov | i_d,sat |
|---|---|---|---|
| 1 | 2.0 V | 4.00 S | 4.000 A |
| 2 | 1.0 V | 8.00 S | 4.000 A |
| 4 | 0.5 V | 16.0 S | 4.000 A |

- **Ron, knee shape and saturation level are not independently settable**:
  fitting Rds(on) forces Ron, which forces the knee, which forces
  $i_{d,sat} = \tfrac{1}{2}K_p V_{ov}^2$ — the square law then over-predicts
  saturation current by 4.3–5.3× over v_gs = 10–12 V at v_ds = 20 V, the
  regime that matters for buck/bridge inrush.
- **The saturation exponent is locked at α = 2.00.** Real power devices sit at
  α ≈ 1.2–1.6 because of drift-velocity saturation.
- **No subthreshold conduction, no mobility degradation, no DIBL, no kink, no
  defined off-state slope.**

| Quantity | Level-1 | Improved |
|---|---|---|
| C¹ continuity at knee | yes (by construction) | yes |
| i_d,sat ∝ v_ov^α | α = 2.00, fixed | α tunable (see below) |
| Ron vs i_d,sat coupling | locked together | independent |
| Subthreshold | none (i_d = 0 below Vto) | smooth softplus |
| Off-state slope | undefined | gate-dependent softplus floor |
| Trench kink | none | one-sided sigmoid |
| Gate-charge saturation (SJ parts) | none | soft cap `v_ovmax`, `k_ovsharp` |
| Interface-trap-filling gain (WBG) | none | gain `k_beta_trap` |
| Knee width vs knee position | locked (mtriode moves both) | independent knob `k_kneesharp` |

`vdmos_id_level1()` implements Level-1 **as actually implemented** (knee at
$v_{ov}/m_{tri}$, verified against LTspice 26.0.1 `.meas`) with the same
package-R embedding, for a fair A/B.

## Channel model

One C¹ continuous expression for the intrinsic drain current
($v_{gs}, v_{ds}$ source-referenced; all names are code identifiers, Greek
parameters carry the `k_` prefix):

$$
i_d = k_p\, v_{ov}^{eff2}\,
\frac{v_{dseff} - \dfrac{v_{dseff}^2}{2\,v_{knee}}}{1 + \dfrac{v_{dseff}}{v_{vsat}}}
\bigl(1 + k_\lambda\, v_{ds}\, k_{db}\bigr)
\;+\; g_k\, v_{ds}\, k_{kink}\, k_{db}
$$

built from the chain

$$
k_{db} = \tanh\!\left(\frac{v_{ds}}{v_{dbsat}}\right)
\qquad\text{(drain-bias factor: 0 ohmic, ~1 above the knee)}
$$

$$
v_{th} = v_{th0} - k_\eta\, v_{ds}\, k_{db}
\qquad\text{(DIBL, bounded)}
$$

$$
v_{ov} = v_{sub}\,\ln\!\Bigl(1 + e^{(v_{gs} - v_{th})/v_{sub}}\Bigr)
\qquad\text{(softplus: subthreshold → strong inversion)}
$$

$$
v_{ov}^{eff} = \frac{v_{ov}}{\Bigl[1 + \bigl(v_{ov}/v_{ovmax}\bigr)^{k_{ovsharp}}\Bigr]^{1/k_{ovsharp}}}
\qquad\text{(soft gate-charge cap, } v_{ovmax}=0 \Rightarrow \text{identity)}
$$

$$
v_{ov}^{eff2} = v_{ov}^{eff}\,\bigl(1 + k_{\beta\text{-trap}}\, v_{ov}^{eff}\bigr)
\qquad\text{(interface-trap-filling gain, off at } k_{\beta\text{-trap}}=0)
$$

$$
v_{knee} = \max\!\left(\frac{v_{ov}^{eff2}}{1 + k_\theta\, v_{ov}^{eff2}},\; v_{kneemin}\right)
\qquad\text{(vertical-field mobility degradation)}
$$

$$
v_{dseff} =
\begin{cases}
v_{knee}\,\tanh\!\bigl(v_{ds}/v_{knee}\bigr) & k_{kneesharp}=0 \text{ (default)}\\[6pt]
\dfrac{v_{ds}}{\Bigl[1 + \bigl(v_{ds}/v_{knee}\bigr)^{k_{kneesharp}}\Bigr]^{1/k_{kneesharp}}} & \text{power-law soft-min}
\end{cases}
$$

$$
k_{kink} = \tfrac{1}{2}\Bigl[1 + \tanh\!\Bigl(\frac{v_{gs} - v_{thk}}{v_{kovsoft}}\Bigr)\Bigr]
\qquad\text{(one-sided trench-kink gate)}
$$

### Limits and the saturation exponent

For $v_{ds} \ll v_{knee}$: $i_d \to k_p\, v_{ov}^{eff2}\, v_{ds}$ (ohmic; the
Ron conductance is $k_p v_{ov}^{eff2}$, independent of any knee knob).
For $v_{ds} \gg v_{knee}$: $v_{dseff} \to v_{knee}$ and

$$
i_{d,sat} = \frac{k_p\, v_{ov}^{eff2}\, v_{knee}}{2\,\bigl(1 + v_{knee}/v_{vsat}\bigr)}\,\bigl(1 + k_\lambda v_{ds}\bigr)
$$

With $v_{vsat} \to \infty$ the gate exponent is
$\alpha \equiv \partial \ln i_{d,sat}/\partial \ln v_{ov}^{eff2}
= 1 + \bigl(1 + k_\theta v_{ov}^{eff2}\bigr)^{-1} \in (1, 2]$:
$k_\theta = 0$ recovers the square law, large $k_\theta$ approaches
$\alpha \to 1$. Velocity saturation ($v_{vsat}$) lowers α further; the
DEFAULT set gives α = 1.08 over v_gs = 6–12 V, and α ≈ 1.2–1.6 — the physical
WBG/SJ regime — is reachable with $v_{vsat} \approx 3\text{–}12$ V
($k_\theta$ secondary).

Setting $k_\theta = v_{vsat} = k_\eta = g_k = 0$ makes the model a smooth
approximation of Level-1 with $m_{tri}=1$: exact in the ohmic limit and in
deep saturation, up to ~6 % low around the knee (the tanh soft-switch rounds
the hard corner).

### Mechanism notes

- **Gate-charge cap** (super-junction parts): the measured gate exponent
  collapses from ~1.7 to ~0.1 above v_gs ≈ 8 V — the lightly doped pillars
  hold finite free charge, so effective overdrive saturates. The soft
  power-law cap is identity for $v_{ov} \ll v_{ovmax}$, asymptote $v_{ovmax}$
  beyond, and
  $\dfrac{d v_{ov}^{eff}}{d v_{ov}} = \bigl(1 + u^{k_{ovsharp}}\bigr)^{-1/k_{ovsharp} - 1}$
  with $u = v_{ov}/v_{ovmax}$, hence $0 \le dv_{ov}^{eff}/dv_{ov} \le 1$
  everywhere. Measured collapse is abrupt → $k_{ovsharp} \approx 10$, not 1.
  Applied **before** $v_{knee}$: the pinch-off voltage is set by the final
  channel charge. Evaluated in log space (`logaddexp`) — $k_{ovsharp}$ ~10–40
  is a fitted regime, not a hypothetical one.
- **Trap-filling gain** (4H-SiC): implied channel overdrive slope
  $dv_{ov}/dv_{gs}$ rises past 1 as interface traps fill — impossible for the
  pre-gain model (softplus slope saturates at 1, the cap only reduces). The
  gain gives $dv_{ov}^{eff2}/dv_{ov}^{eff} = 1 + 2k_{\beta\text{-trap}} v_{ov}^{eff} > 1$.
  **$k_{\beta\text{-trap}} \ge 0$ is a hard constraint** (ValueError below):
  negative values make $v(1+k v)$ non-monotone and $i_d$ multi-valued in
  $v_{gs}$. With the cap active the composite asymptote is
  $v_{ovmax}(1 + k_{\beta\text{-trap}} v_{ovmax})$, so the two knobs
  anti-correlate on silicon — never mix records that free one vs the other.
- **Knee softness** $k_{kneesharp}$: $k_\theta$ and $v_{vsat}$ move knee
  position *and* width together (~84 % correlated). The power-law soft-min
  breaks the coupling: $dv_{dseff}/dv_{ds}(0) = 1$ and $v_{dseff} \to v_{knee}$
  hold for every $k_{kneesharp}$ — Ron and the saturation plateau are
  untouched while knee width shrinks monotonically
  ($k \to \infty$ hard min, $k = 1$ harmonic, $k \approx 2.8$ reproduces the
  tanh). Same $0 \le dv_{dseff}/dv_{ds} \le 1$ derivative bound as the cap.
- **Trench kink**: the gate is the one-sided sigmoid above, not bare tanh —
  tanh is odd and injects *negative* current below onset (non-physical, a
  convergence hazard in the off state).
- **Off state**: below threshold the softplus channel is the only channel
  current, falling gate-dependently with slope $1/v_{sub}$. With the body
  diode on, its reverse saturation $i_s(T) + i_{sr}(T)$ is the floor instead.
  Real HV leakage needs a dedicated BV/Ibv element, not a fitted resistor.
- **Polarity**: the channel formulation is forward-only ($v_{ds} \ge 0$) —
  the velocity divisor is singular near $v_{ds} = -v_{vsat}$ and $k_{db}$ is
  odd. All reverse current is diode current; sign-symmetric synchronous-
  rectifier conduction is a queue item, not this model.

### Sentinels (defaults = legacy behavior, bit-identical)

| off switch | effect |
|---|---|
| `v_ovmax = 0` | cap is the exact identity |
| `k_beta_trap = 0` | gain is the exact identity |
| `k_kneesharp = 0` | exact legacy tanh knee (fitted floor 0.1 — k→0 collapses i_d) |
| `r_d = r_s = 0` | explicit fast path, no implicit solve |
| `i_s = i_sr = 0` | diode current exactly 0 |
| `i_kf = 0` | high-injection divisor off (LTspice `Ikf=inf`) |
| TCs `= 0` | no diode temperature dependence |

"Off" is reached by leaving a parameter out of the free set, never by driving
it to the lower bound.

## Package resistance (r_d, r_s, r_sk)

Bond wires/leadframe/DBC are embedded **outside** the intrinsic channel,
exactly as LTspice's Level-1 `Rd`/`Rs` do:

$$
v_{ds}' = v_{ds} - (r_d + r_s)\, i_d \qquad\text{(current loop: all of } i_d)
$$
$$
v_{gs}' = v_{gs} - (r_s - r_{sk})\, i_d \qquad\text{(gate loop: returns at the kelvin tap)}
$$
$$
i_d = f_{int}(v_{ds}', v_{gs}') \qquad\text{(implicit, solved per point)}
$$

$r_{sk}$ is the resistance between the kelvin tap and the power source pin
(carrying only gate current): $r_{sk}=0$ is no kelvin contact (full
degeneration $g_{m,eff} = g_m/(1 + r_s g_m)$), $r_{sk} = r_s$ is an ideal tap
at the die. Guard $0 \le r_{sk} \le r_s$ (ValueError). DC data senses only
$r_d + r_s$ and $r_s - r_{sk}$; the split itself needs 4-terminal AC data.

The implicit equation is solved with bracketed safeguarded Newton (fixed
point diverges: the contraction factor is $(r_d+r_s) k_p v_{ov} \approx 1.6$
at DEFAULT). The bracket $i_d \in [0, v_{ds}/(r_d+r_s)]$ is invariant because

$$
\frac{dg}{di_d} = -(r_d + r_s)\frac{\partial f}{\partial v_{ds}} - (r_s - r_{sk})\frac{\partial f}{\partial v_{gs}} - 1 \le -1
$$

whenever the partials are non-negative — monotonicity is all the proof needs
(the trap gain inflates the bound but preserves monotonicity for
$k_{\beta\text{-trap}} \ge 0$). $f_{int}(0, \cdot) = 0$ exactly, which makes
clamping internal $v_{ds}$ at 0 harmless on overshoot. Non-convergence raises
RuntimeError; it never silently returns an unconverged current.

## Body diode

The **LTspice diode DC branch** (`diode.py`), wired at the device TERMINALS so
it **bypasses r_d/r_s** — only channel current flows through the package
resistances, keeping the two implicit solves independent 1-D problems:

```
   drain o──[r_d]──(D')══════ channel ══════(S')──[r_s]──┬─o source
         │                                               │
         └───────── diode + r_diode ─────────────────────┘   anode = source
```

Terminal current ($v_{sd} = -v_{ds}$, positive $i_{di}$ flows source → drain):

$$
i_d^{term} = i_{ch}(v_{ds} > 0) - i_{di}(v_{sd})
$$

$$
i_{di} =
\underbrace{\frac{i_s(T)\bigl[e^{v_{di}/(k_n v_T)} - 1\bigr]}{\sqrt{1 + \dfrac{i_s(T)}{i_{kf}}\bigl[e^{v_{di}/(k_n v_T)} - 1\bigr]}}}_{\text{diffusion + high-injection compression}}
+ \underbrace{i_{sr}(T)\bigl[e^{v_{di}/(k_{nr} v_T)} - 1\bigr]}_{\text{recombination}},
\qquad
v_{di} = v_{sd} - r_{diode}(T)\, i_{di}\ \text{(implicit)}
$$

Temperature laws (LTspice-exact, all referenced to $T_{nom}$):

$$
i_s(T) = i_s \left(\frac{T}{T_{nom}}\right)^{k_{xti}/k_n}
\exp\!\left[\frac{q\, v_{eg}}{k_B k_n}\left(\frac{1}{T_{nom}} - \frac{1}{T}\right)\right]
\qquad\text{(same law for } i_{sr}, \text{ divisor } k_{nr})
$$

$$
i_{kf}(T) = i_{kf}\bigl[1 + a_{ikf}(T - T_{nom})\bigr]
\qquad
r_{diode}(T) = r_{diode}\bigl[1 + a_{rd}(T - T_{nom}) + b_{rd}(T - T_{nom})^2\bigr]
$$

(code spellings `a_rdiode`, `b_rdiode`; the $a$/$b$ pair is the complete
LTspice Trs1+Trs2 Rs law; the TC polynomial factor is clamped at a tiny floor
where extrapolation would go negative, protecting the solver bracket.)

LTspice model-card mapping:

| LTspice | `Is` | `N` | `Isr` | `Nr` | `Ikf` | `Rs` | `Tnom` | `Xti` | `Eg` | `Tikf` | `Trs1` | `Trs2` |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| package | `i_s` | `k_n` | `i_sr` | `k_nr` | `i_kf` | `r_diode` | `t_nom` | `k_xti` | `v_eg` | `a_ikf` | `a_rdiode` | `b_rdiode` |

Verified against LTspice 26.0.1 to ~1e-5 relative (high-injection divisor
across 3.4 A … 375 kA; TC laws at −55/25/100/150 °C) — details and the known
LTspice `Isr`-branch quirk in STATUS/tests. Capacitance (Cjo/Vj/M), transit
time, breakdown and noise are out of scope here: junction charge is the caps
delta (below), recovery is the network (below), BV/Ibv is a queue item.

## Body-diode recovery (minority-carrier storage) — LTspice-side

Recovery is a **behavioral LTspice network** (`DiodeRecov` subckt), not
package Python — it lives in the export layer
(`out/ISC040N10NM8_recov.subckt`; records `out/RECOV_POLE_FIT.md`,
`out/RECOV_INTRINSIC_FIT.md`, `out/RECOV_REPORT.md`; extraction driver
`src/vdmos_fit/extract_recov.py`). D1 is the diode DC branch above; the
storage node `q_j0` carries $C_1 = t_{\text{minor}}$ farads to ground
(leak $R_{off} = 1$), fed by the injection mirror $B1 = -i_{R2}$ and
clamped $V_{q0} \ge 0$ by ClampD. A unity buffer $E1$ + $R3 = 1\,\Omega$ +
$C3 = t_{\text{pole}}$ inserts **one pole** — $v_{q_j} = v_{q0}$ filtered
with time constant $t_{\text{pole}} = R3\,C3$ (no C1 loading) — and the
conduction drive reads the FILTERED node:

$$
i_{rec} = \max\!\left[\frac{V_{q_j}\,(v_{Ki} - v_{Ai})}{\tau},\; 0\right]
$$

**Gauge (session 25b).** With $R_{off} = 1$, $C_1 = t_{\text{minor}}$:
storage time constant $= t_{\text{minor}}$, $V_{q0,ss} = I_F$ volts,
stored charge $= I_F\,t_{\text{minor}}$. $q_0$ is DELETED — observables
depended only on $t_{\text{minor}}$ and $V_R/(q_0\tau)$; it was a pure
gauge.

**Three knobs:** $t_{\text{minor}}$ (carrier lifetime — stored charge and
Qrr), $t_{\text{pole}}$ (transit/sweep-out memory — the pole makes
$I_{RM}$ occur 1–2 ns **after** charge exhaustion and makes the Qrr
di/dt-ratio tunable), and $\tau$ (B2 snap threshold). $\tau$ has a
genuine **interior optimum ≈ 3e-4**: the clamp→$I_{RM}$ delay obeys

$$
\Delta t = t_{\text{pole}}\,\ln\!\frac{V_{q_j}|_{\text{clamp}}}{V_{th}},
\qquad V_{th} = \frac{I_{RM}\,\tau}{v_{Ki}-v_{Ai}},
\qquad \frac{d\,\Delta t}{d\ln\tau} = -\,t_{\text{pole}}
$$

(verified to ~1 %), so shrinking $\tau$ adds tail time that must be paid
by shortening $t_{\text{pole}}$ (sharply identified, wall below ~0.6 ns)
and $t_{\text{minor}}$ — the minimum is where those costs balance.

**Measurement doctrine (session 25a).** The datasheet Qrr excludes the
$C_{oss}$ excursion, so fit the **intrinsic** branch current $i_{R2}$
(which carries no C2 displacement current); terminal-path numbers differ
by ≈ 28 nC at 100 A/µs and are reported for reference only.

**Fitted on ISC040N10NM8** ($V_R = 50$ V, $I_F = 12.5$ A, 25 °C, both
datasheet di/dt rows): $t_{\text{minor}} = 22.363$ n,
$t_{\text{pole}} = 0.8324$ n, $\tau = 3.79\times 10^{-4}$ →
**RMS 1.11 % over all four datasheet numbers** (1000 A/µs: Qrr 248.2 nC
+0.5 %, t_rr 24.69 ns −1.3 %; 100 A/µs: Qrr 35.71 nC −0.8 %, t_rr 29.46 ns
+1.6 %; t_rr = 10 % of $I_{RM}$, intrinsic path; record `ISC040_RECOV`). The two-knob model's
Qrr-ratio floor (0.175 vs datasheet 0.146) is DEAD — the pole reproduces
both di/dt rows from one physical state; the session-23 "Qrr rises with
di/dt is unreproducible" verdict is superseded. $t_{\text{pole}}$ is
**sharply identified** (fix-and-refit: 0.3 n → 15.6 %, 2 n → 11.6 %,
8 n → 52.8 %, → 0 → 17.45 % RMS). Physical anchor:
$L = \sqrt{2 D_a t_{\text{pole}}} = 1.4\text{–}2.7\ \mu$m for
$D_a = 12\text{–}36$ cm²/s — a diffusive scale, right for the
depletion-adjacent region of a 100 V SJ device, NOT a drift transit time.
Falsifiable out-of-sample prediction at 300 A/µs: Qrr ≈ 104 nC,
t_rr ≈ 29 ns, $I_{RM}$ ≈ 6.7 A.

**Extraction protocol (session 25c, `extract_recov.py`).** One batched
`.step param k` index sweep per di/dt point (mixed-radix decomposition —
deterministic, no `.step` nesting ambiguity), log-log multilinear
surrogates for (Qrr, t_rr), multi-start surrogate inverse solve, then TRF
polish on real runs in log space with `diff_step = 0.03` (smaller steps
put the cost difference below the simulator's own ~1e-4 relative
measurement noise — the Jacobian goes noise-dominated and the fit stalls;
measured), a generous step budget (tau is a flat valley coupled to both
other knobs), tau free (≥ 1e-6), **no dTmax** (the pole rounds the snap:
none/2.5n/1n/0.5n agree to ≤0.01 % and none is faster — the session-24
"dTmax mandatory" finding was pole-free-era), **no `.ic`** (the DC
operating point IS the pre-event state: no DC current in R3 ⇒
$V_{q_j} = V_{q0}$), `solver=alt method=gear reltol=1e-4`.
Identifiability by fix-and-refit penalties, never marginal σ.

**Integration requirements** (LTspice 26.x, all verified):
`method=gear reltol=1e-4` (trapezoidal rings on the snap); the measured
C_ds charge-master sits **inside** the recovery branch at Ki–A
(electrically D–S through B3+R1/R2, fc ≈ 3 GHz) — at the D–S terminals
hard forward re-turn-on crashes (D1 has no Cjo), at Ki–Ai the network is
singular; the caps delta therefore exports **Cgs/Cgd only**
(`ltspice.py` `cgs_cgd_only=True`) — never export Cds in both places
(double-counts C_oss); `R2 = 1 µΩ` + B3 `+1 µΩ` floor break
near-singular B3/R1/R2/switch bridge loops (< 0.1 % waveform impact).

## Capacitance / charge (summary — full spec in CAP_SPEC.md)

Charge-master LUTs $q_k(v) = \int C\,dv$ (PCHIP on ln C vs linear v) built
from digitized Ciss/Coss/Crss, with the $c_{gd}(v_{gd})$ splice (off-state
table for $v_{gd} < 0$, exponential-saturation branch for $v_{gd} \ge 0$,
$c_{gdmax}$ extracted from Qg curves). Exported to LTspice as the general
nonlinear capacitor in **incremental ramp-sum** form,

$$
q(x) = q_0 + m_1 (x - x_0) + \sum_k \Delta m_k \max(0, x - x_k)
$$

so LTspice's symbolic $dq/dx$ is exactly the piecewise-constant C staircase
and the end-segment slopes carry the clamp-C policy for free. The delta
(Cgs G–S, Cgd G–D, Cds D–S) is terminal-wired; charge conservation is
structural. Live parity ≤ 0.3 % on all three catalog parts.

## Naming convention

Every model quantity is `<type prefix>_<suffix>` (source of truth:
`names.xlsx`):

| prefix | type | examples |
|---|---|---|
| `v_` | voltage [V] | `v_th0`, `v_th`, `v_ov`, `v_ov_eff`, `v_ov_eff2`, `v_ovmax`, `v_sub`, `v_knee`, `v_kneemin`, `v_dseff`, `v_vsat`, `v_dbsat`, `v_thk`, `v_kovsoft` |
| `i_` | current [A] | `i_d` |
| `r_` | resistance [Ω] | `r_d`, `r_s`, `r_sk` |
| `g_` | conductance [S] | `g_k` — trench-kink conductance, *not* transconductance |
| `k_` | factor | `k_db`, `k_kink`, `k_ovsharp`, `k_kneesharp`, `k_p`, `k_theta`, `k_eta`, `k_lambda`, `k_beta_trap`, `k_agd` |
| `c_` | capacitance [F] | `c_gd`, `c_gs`, `c_ds`, `c_gd0`, `c_gdmax` |
| `q_` | charge [C] | `q_gd`, `q_gs`, `q_ds` — the charge-master state variables |

- Greek parameters carry `k_` so no name collides with a Python keyword
  (θ → `k_theta`, λ → `k_lambda`). The `k_` factors are unitless EXCEPT the
  industry-standard set: `k_p` [A/V²] and the 1/V quartet `k_lambda`,
  `k_eta`, `k_theta`, `k_beta_trap`.
- `db` = drain **bias**, not drain-to-body — this model has no bulk terminal.
- Sharp/soft duality: `*_sharp` names are power-law **exponents** (bigger =
  sharper: `k_ovsharp`, `k_kneesharp`); `*_soft` names are voltage **scales**
  (bigger = softer: `v_kovsoft`).
- `m_triode` is the one carve-out: `vdmos_id_level1()` is a foreign-model
  adapter and keeping the model-card spelling `mtriode` *is* the mapping.
- The recovery parameters (`t_minor`, `t_pole`, `tau`, `k_chg`, `r_max`) are
  deliberately outside the convention — LTspice-side netlist values,
  dimension-normalized (see the recovery section).

## Extraction

1. **Relative weighting** first — subthreshold and high-current points weigh
   equally, so the low-current shape (`v_th0`, `v_sub`) is not swamped;
   **ampere weighting** second — makes the Ron/Rds(on) corner dominate, which
   is what datasheet specs pin down.
2. Positive parameters are optimised in **log space** with physical bounds
   (zero-able linear parameters — `r_d`, `r_s`, `v_ovmax`, `k_kneesharp` —
   stay linear; their 0 value is the "off" sentinel and cannot be logged).
3. The **free set IS the protocol**: check `identifiability()` — dead
   Jacobian columns report σ = inf, never 0; order-unity σ means weakly
   determined. `k_theta`/`v_vsat` are ~84 % correlated (both flatten the
   knee); quote the curve, never the split.
4. Static data pins $r_d + r_s$ tightly (±0.3 % synthetic) but the split is
   −0.999 correlated and leans on r_s's gm-degeneration signature (±15 %).
   Quote the sum; split from geometry or 4-terminal data.
5. $r_d + r_s$ is always **below** the measured Rds(on) chord — the rest is
   the intrinsic channel the model computes separately.

## Layout

```
ImprovedVDMOS/
├── pyproject.toml          metadata, deps, pytest + ruff config
├── src/vdmos_fit/          the importable package
│   ├── model.py            vdmos_id(), vdmos_id_level1(), package-R solver
│   ├── diode.py            vdmos_diode(), LTspice diode DC branch + r_diode solver
│   ├── caps.py             charge-master C–V LUTs + c_gd splice + c_gdmax extraction
│   ├── params.py           Params dataclass, bounds, log-space packing
│   ├── data.py             CurveSet / DiodeCurveSet / CapCurveSet / QgCurveSet, CSV I/O
│   ├── fit.py              least_squares extraction + identifiability (+ fit_diode)
│   ├── extract_recov.py    datasheet → DiodeRecov extraction (LTspice-driven)
│   ├── reference.py        official fitted device records
│   ├── ltspice.py          LTspice export: charge-master nonlinear-cap subckt
│   ├── plotting.py         output/transfer/gds/overlay plots
│   ├── cli.py              `vdmos-fit` command
│   ├── STATUS.md           session log, official records, open problems
│   └── CAP_SPEC.md         settled capacitance/charge-model design spec
├── src/curve_digitizer/    PDF curve digitizer — `curve-digitizer` console script
├── tests/                  pytest suite (records pinned; LTspice cross-validation
│                           and live cap parity run when LTspice is installed)
├── scripts/                thin entry points (no logic here)
├── data/                   measured/digitised CSVs (see data/README.md)
├── notebooks/              exploratory work
└── out/                    generated plots, fitted params, subckt templates (gitignored)
```

## Setup

```bash
cd Projects/ImprovedVDMOS
python -m venv .venv
.venv\Scripts\activate            # Windows
pip install -e ".[dev]"
pytest                            # suite green (headless Tk tests self-skip; totals drift)
vdmos-fit check                   # reference numbers + plots
python scripts/roundtrip_demo.py  # synthetic recovery of known parameters
```

Without installing, the `scripts/` entry points add `src/` to `sys.path`
themselves.

## Usage as a library

```python
from vdmos_fit import Params, vdmos_id, fit, load_csv

p = Params(k_p=30.0, v_th0=3.1, k_theta=0.02, v_vsat=1.2)
id_amps = vdmos_id(v_ds=20.0, v_gs=10.0, p=p)

res = fit(load_csv("data/IPW065N04N.csv"), weight="relative")
print(res.summary())
```

## Scope / not yet implemented

In scope: channel current + package series R + body-diode DC branch +
charge-master capacitance + behavioral recovery (LTspice-side).
Open: $v_{gs}$-dependent $C_{gd}$ blend (the Miller-plateau deficit
quantifies it: SJ +61 %, trench +22–43 %, SiC +10 %), full-device subckt
channel B-source + kelvin wiring, HV leakage/BV element, $T_j$ thermal node,
third-quadrant channel, recovery-parameter datasheet extraction protocol.
Validation target: double-pulse Eon/Eoff within ~10 %, which static fitting
alone cannot deliver.

## Prior art

Hefner's physics-based MOSFET model, Infineon Level-100, Wolfspeed SiC
models, and the Vayser / "Modelling Farmers" LTspice lineage.
