# Capacitance / charge model — settled design spec

Status: **design settled with the user (2026-10-18); code not started.**
This file is the single source of truth for the C–V / Q(V) implementation.
Digitized data on disk: `data/<part>_Cap_Vds.csv` (three parts). Qg curves:
to be digitized (contract in §Qg data).

## Scope

Three-capacitor delta (G–S, G–D, D–S) per device, **charge-based**, built from
digitized datasheet C–V curves measured at **V_GS = 0** (off-state).
No closed-form C–V equations per device (rejected — see §Rejected).
No temperature dependence (all curves are 25 °C typical). The channel-charge
v_gs dependence of C_gd is deferred (see §Limitations) — the first suspect if
simulated turn-on is too fast.

## Decomposition (datasheet Ciss/Coss/Crss → delta)

```
c_gd = Crss                       (table + splice at v_gd = 0, below)
c_gs = Ciss − Crss                (plain table)
c_ds = Coss − Crss                (plain table)
```

Positivity verified on all three digitized sets (min Ciss−Crss = 3023 / 4592 /
12 800 pF, min Coss−Crss = 254 / 72 / 496 pF — all > 0). Fitting the three
branches independently would double-count; the decomposition is the contract.

## Storage: charge is master

```
q_k(v) = ∫₀^v C_k(x) dx        stored / interpolated
C_k(v) = dq_k/dv               solver differentiates (FTC: exact by construction)
```

- **Interpolation: PCHIP on (v, ln C)** — linear v axis, log C. Positive by
  construction, monotonicity-preserving through the SJ cliff, C¹ in v up to the
  splice. NOT log-log PCHIP: the splice sits at v_gd = 0 and log 0 = −∞.
- q_k(v): per-segment Gauss–Legendre quadrature of exp(PCHIP), cumulative sums
  precomputed at load; the exp branch of c_gd has a closed-form integral (§Splice).
- **Out-of-range policy:** clamp C at the last digitized node (q grows linearly
  beyond); record the digitized v range in the provenance record. Closed forms
  extrapolate wildly for SJ curves; a clamped table is safer, not worse.

## c_gd(v_gd) composite — the only splice

v_gd = v_g − v_d. The datasheet Crss(V_DS) at V_GS = 0 is exactly the
v_gd < 0 branch with v_ds = −v_gd.

```
C_gd(v_gd) = PCHIP_lnC[Crss](−v_gd)                            v_gd < 0
           = Cmax − (Cmax − C0)·exp(−A·v_gd)                   v_gd ≥ 0
A = −s / (Cmax − C0)          (C1 at the splice is then automatic for ANY Cmax > C0)
```

- **Origin anchor convention:** C0 and s = dC_rss/dv_ds|₀ come from a quadratic
  fit `C = C0 + s·v + b2·v²` to the **first 4 digitized points** of each Crss
  branch. The window is part of the convention (SJ s is window-sensitive:
  −42 pF/V at n=3, −152 at n=4, −468 at n=5 — the SJ plateau→cliff inflection
  is real, not noise).
- C0 continuity: exp branch at v_gd = 0 equals C0. C1: slope at 0⁺ is
  (Cmax−C0)·A = −s = slope of the table branch at 0⁻ (dC_gd/dv_gd = −dC_rss/dv_ds).
- **Cmax is the single tunable knob**; A is tied to it by C1.

### Measured anchors (first-4-point convention, 2026-10-18)

| part | C0 [pF] | s [pF/V] | b2 [pF/V²] | Ciss(0) [pF] | c_gs(0) [pF] | seed Cmax [pF] | A [1/V] | τ=1/A [V] |
|---|---|---|---|---|---|---|---|---|
| ST8L65N044M9 (SJ 650 V) | 2700 | −151.5 | −271 | 7658 | 4958 | 4958 | 0.0671 | 14.9 |
| IMCQ120R004M2H (SiC 1200 V) | 3823 | −663.7 | +41 | 16208 | 12385 | 12385 | 0.0775 | 12.9 |
| ISC040N10NM8 (trench 100 V) | 491 | −53.9 | +3.4 | 3448 | 2957 | 2957 | 0.0219 | 45.7 |

Caveats: the SiC plot's x-axis starts at 1 V, so its C0/s are extrapolations
(±10–15 %). The origin slopes are NOT flat — the log-log datasheet axes made
them look flat (log-log slope ≈ −0.03…−0.6, but linear slope is tens–hundreds
of pF/V). Independent cross-check on the trench: graded-junction fit
C₀(1+v/V_bi)^−M gives M = 0.721, V_bi = 6.41 V, max residual 10.6 %, implying
s = −MC₀/V_bi = −55 pF/V vs the direct fit −54 — the splice slope is real.
The same single junction law fails on the SiC part (≈100 % residual): the
table approach is confirmed for SiC and SJ.

## Cmax protocol

1. **Seed:** `Cmax = c_gs(0) = Ciss(0) − C0`. Physical basis: when the channel
   forms (v_gd > 0) the gate's counter-electrode is the channel, a resistive
   sheet splitting its capacitance roughly evenly to source and drain, so
   C_gd(triode) ≈ ½C_gc ≈ C_gs. Consistent with vendor VDMOS cards (LTspice
   help example: Cgs = 3.1 nF, Cgdmax = 1.9 nF).
2. **Bracket:** `C0 ≤ Cmax ≤ 2·c_gs(0)` (fitter bounds).
3. **Extraction from the Qg curve (the real answer — a datasheet measurement,
   no hardware needed).** CORRECTED SPLIT (2026-10-18, implementation day):
   the Miller plateau ENDS at v_gd = 0 (v_ds falls to v_gs), so the plateau
   trajectory never enters the exp branch — the plateau charge is a PURE
   table-branch prediction and carries NO c_gdmax information. It VALIDATES
   the tables:

   ```
   Q_plateau_meas  ?=  ∫_0^{V_DS−V1} Crss(v) dv + ∫_{V1}^{V2} c_gs(v) dv
   ```

   (V1, V2 = plateau start/end voltages; the c_gs term is the sloped-plateau
   recharge — a flat-plateau behavioral model would miss it.) c_gdmax is
   extracted from the POST-PLATEAU segment (v_gd > 0, v_ds → i_d·Rds(on)):

   ```
   remainder = (Q_end − Q2) − ∫_{V2}^{V_end} c_gs(v) dv
   remainder = Cmax·V − (Cmax−C0)(1−e^{−A·V})/A,   V = V_end − i_d·Rds(on)
   ```

   Everything except Cmax is known from the digitized C–V tables; solve the
   scalar equation (nonlinear through A(Cmax); 1-D bracketed brentq).
4. **Diagnostic:** if the remainder is BELOW the flat-branch value
   q_exp(V; Cmax = C0) (negative channel contribution under the c_gs table
   convention), clamp Cmax = C0 (flat branch) and record the discrepancy —
   do not chase a negative charge. Symmetrically, remainder above
   q_exp(V; 2·c_gs(0)) clamps UPPER. Both paths implemented and tested.
5. **Quote the curve, not the split:** Cmax and s are correlated in this
   extraction (bigger measured slope → smaller needed Cmax), and s is the
   noisiest number in the model. Same rule as the diode fit.

RESULTS (2026-10-18, scripts/fit_cmax.py; records in reference.py):

| part | plateau Qgs err | plateau deficit (meas − table) | c_gdmax extracted | vs seed |
|---|---|---|---|---|
| ST8L65 (400 V, 29 A) | +4.1 % | **+61.2 %** | 2750 pF (≈ c_gd0) | 0.56× seed |
| IMCQ120 (800 V, 185 A) | +16.2 % | +9.7 % | 3823 pF = c_gd0 (LOWER clamp) | 0.29× seed |
| ISC040 (20/50/80 V, 12.5 A) | −4.5…−1.0 % | +42.8/+24.5/+22.0 % | 663 pF (median of 647/798/663) | 0.22× seed |

HEADLINE: **all three parts land at c_gdmax ≈ c_gd0 — the exp branch is
nearly FLAT; the seed c_gs(0) over-predicts 4–8×.** These Qg curves give no
evidence for an LTspice-style Cgdmax ≫ Cgdmin channel cap. The channel charge
the off-state tables miss shows up instead as the Miller-PLATEAU DEFICIT
(limitation 1): worst on the SJ part (+61 %, its plateau is nearly flat =
huge effective Miller capacitance), then trench, smallest on SiC. That
deficit — not c_gdmax — is the quantity the deferred v_gs-dependent blend
must eventually explain. r_ds_on for the ohmic tail is a documented
datasheet-class assumption (36/4/40 mΩ part-code classes); the V_on = 0
probe moves c_gdmax by < 10 % (ISC040) to a LOWER clamp (ST8L65).

Sensitivity is bounded: the exp branch only acts for v_ds < v_gs (last slice
of the plateau + ohmic), and the extra charge scales ~linearly with Cmax —
a 2× Cmax error moves turn-on time modestly. Double-pulse turn-on time is the
final arbiter.

## Rejected alternatives (on the record)

- **Closed-form C–V per device.** Not one equation but three different sets
  (junction law for trench; constant+graded-tail for SiC; sigmoid-blend for the
  SJ plateau→cliff + non-monotonic Crss bump), each with 6–8 parameters fit to
  one typical curve at one frequency/temperature — a `kappa`-grade mechanism by
  this project's own standard. Vendor SPICE models for SJ/SiC parts use C–V
  tables internally.
- **Zero-slope tanh² positive branch** (earlier proposal): wrong — the table
  slope at the origin is not flat (log-scale illusion); C1 requires a branch
  that leaves the splice with slope −s > 0.
- **LTspice's tanh branch verbatim.** Solving C0+C1 gives
  Cgdmin = 2·C0 − Cmax, so Cgdmin > 0 requires Cmax < 2·C0 — **violated** for
  SiC (seed 12 385 vs 2·3823) and trench (2957 vs 982); only the SJ barely
  passes. The exponential-saturation form has no such ceiling. (LTspice's own
  tanh/atan splice is not even C1 — π/2 slope mismatch.)
- **Escape hatch kept in reserve:** tanh with inflection shifted to v₀ > 0 lifts
  the ceiling at the cost of one extra shape knob — only if double-pulse says
  the Cmax rise must be convex (abrupt channel formation).

## Limitations (agreed, on the record)

1. At v_gd < 0 with the channel ON (Miller region), the off-state table
   underestimates C_gd. A v_gs-dependent blend is the cure; deferred. First
   suspect if turn-on simulates too fast.
2. No temperature dependence (25 °C typical curves only).
3. Frequency mismatch: C–V at 250 kHz (ST) / 100 kHz (Infineon) / 1 MHz
   (trench) vs quasi-static Qg test — expect a few % real discrepancy in the
   plateau cross-check, not a bug.
4. c_gs(v_gs) uses the (Ciss−Crss)(v_ds) table shape against v_gs. Justified
   because all three Ciss curves are nearly flat vs v_ds (≈ constant), matching
   LTspice VDMOS's constant-Cgs approximation; the SJ iss branch is flat by
   construction (accumulation plateau).

## Qg data contract (digitized 2026-10-18)

V_GS vs Q_g at specified V_DS, I_D (standard quasi-static gate-charge test).
CSV: `index,Q_g(nC),V_gs,condition` (condition = `I_d=..A, V_ds=..V`; header
match ignores case and spaces). Extract: V_gp (plateau voltage), Qgs (charge
to plateau start), Qgd (plateau span). Uses: c_gdmax extraction (§Cmax, via
the post-plateau segment), plateau validation (deficit = limitation-1
quantifier), and (later) V_gp as the anchor for any channel-charge model.

## Validation plan

1. Trench junction-law anchor: pipeline reproduces M ≈ 0.72, V_bi ≈ 6.4 V.
2. Qgd plateau charge vs ∫Crss dv over the plateau span (quasi-static caveat).
3. Qoss = ∫₀^V Coss dv vs any datasheet Qoss spec.
4. Double-pulse Eon/Eoff within ~10 % (standing project target); turn-on time
   arbitrates Cmax.
5. LTspice A/B: export as `.subckt` with `table()` behavioral capacitors.

## Naming (convention extension)

New prefixes for `names.xlsx`: `c_` capacitance [F] (`c_gd`, `c_gs`, `c_ds`),
`q_` charge [C] (`q_gd`, `q_gs`, `q_ds`). Splice parameters follow the existing
style: `c_gdmax` [F], `c_gd0` [F], `k_agd` [1/V] (= A above; unitless-except-1/V
like the k_lambda trio). `db`-style ambiguity does not arise: `gd`/`gs`/`ds`
are terminal-pair names here.

## Code TODO (implementation checklist)

**LANDED 2026-10-18: Phases 1–3 + tests + records. LANDED 2026-10-21:
Phase 4 (LTspice export + parity).**
Phase 1: `CapCurveSet`/`QgCurveSet`, `load_cap_csv`/`load_qg_csv` (+savers)
in `data.py`. Phase 2: `caps.py` — `_BranchTable` (PCHIP on ln C, GL-5
quadrature, clamp policy), `OriginAnchor` (first-4-point convention),
`CapParams` (c_gd0, s_gd0, c_gdmax; k_agd derived; flat-branch limit),
`CapModel` (c_gd/q_gd composite with origin-quadratic head, c_gs/c_ds union-
node difference tables). Phase 3: `extract_cmax` (plateau detection by max
dQ/dV, segment decomposition, bracketed brentq, LOWER/UPPER clamp
diagnostics) + `scripts/fit_cmax.py`. Records: `CAP_*` in `reference.py`;
tests: `tests/test_caps.py` (CSV contracts, decomposition positivity, dq/dv
≡ C identity, splice C0/C1, monotonicity through the SJ cliff, clamp policy,
trench junction anchor, pinned extraction records, synthetic clamp
diagnostics).

Phase 4 (`ltspice.py` + `scripts/export_caps.py`): the delta exports as a
`.subckt <name> G D S` of three **charge-master nonlinear capacitors**
(`Cxxx n1 n2 Q=<expr>` — LTspice's one primitive that IS a charge master;
`x` = device voltage, symbolic dQ/dx = capacitance). The Q= expression is
the **incremental ramp-sum** form of the sampled charge table:

```
q(x) = q0 + m1·(x − x0) + Σ_k (m_{k+1} − m_k)·max(0, x − x_k)
```

(m_k = charge-table segment slopes, table sampled at N_SUB=12 subdivisions
per digitised node interval). LTspice's symbolic derivative is then exactly
the piecewise-CONSTANT C staircase, and the end-segment slopes extend the
clamp-C policy (flat C, linearly growing q) for free — no extension nodes,
no B-source d().

LTspice 26.0.1 gotchas (measured 2026-10-21, on the record):
- `table()` in a B source: rejected — "Function not supported in this
  context".
- `table()` inside Q=: compiles, but its symbolic derivative is GARBAGE
  (measured currents ~1e11 A). Never use.
- Absolute-slope ramp sums (no increments): wrong — they SUM segment slopes
  instead of stepping them (caught by parity at exactly 5×/n× factors).
- `.MODEL name C (C=...)` nonlinear-capacitor model form: rejected as
  "Unknown parameter" on the instance.

Parity (validation step 5, `test_caps.py::test_ltspice_caps_parity` +
live runs on all three parts): PWL-ramp drive, three delta configurations
(A drives G with D=S=0 → c_gs+c_gd; B drives G=D with S=0 → c_gs+c_ds;
C drives D with G=S=0 → c_ds+c_gd(−v) — the delta cannot isolate one branch
by DC wiring), 14 `.meas tran FIND` points per part at segment midpoints
(FIND must land strictly INSIDE staircase segments — a target on a knot
measures the solver's blend; and maxstep must keep the Gear-2 companion's
q-history window ≪ the segment width). Results: measured vs exported-
staircase ≤0.3 % (most ≤0.05 %); staircase vs PCHIP model C ≤0.4 %
(incl. the SJ cliff). Sign convention verified: I(Vdrive) = −dq/dt.

Still open:
- Full-device subckt (channel B source + diode + r_d/r_s + this delta)
  for A/B switching simulation.
- Clean the ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF
  points) on the next re-digitization pass.
- v_gs-dependent blend for the plateau deficit (limitation 1) — now
  QUANTIFIED per part by the deficit numbers above; it is the next mechanism
  to build if double-pulse turn-on misses.
