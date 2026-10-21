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
   no hardware needed).** During the Miller plateau v_gs ≈ V_gp is constant,
   C_gs charges nothing, and the whole plateau charge flows through C_gd along
   v_gd: V_gp − V_DS → V_gp. Split at v_gd = 0 (v_ds = V_gp):

   ```
   Q_plateau_meas = ∫_{V_gp}^{V_DS} Crss(v) dv        (table branch — from the LUT)
                  + Cmax·V_gp − (Cmax−C0)(1−e^{−A·V_gp})/A    (exp branch — solve for Cmax)
   ```

   Everything except Cmax is known from the digitized C–V tables; solve the
   scalar equation (nonlinear through A(Cmax); 1-D bisection/Newton).
4. **Diagnostic:** if Q_plateau_meas < table integral (negative difference),
   the off-state table over-counts (quasi-static mismatch: 250 kHz C–V vs the
   slower Qg test; real on SJ parts) or V_gp was misread. Then clamp
   Cmax = C0 (flat branch) and record the discrepancy — do not chase a
   negative charge.
5. **Quote the curve, not the split:** Cmax and s are correlated in this
   extraction (bigger measured slope → smaller needed Cmax), and s is the
   noisiest number in the model. Same rule as the diode fit.

Predicted exp-branch plateau excess charge (beyond the table integral) at the
seed, for sanity-checking when Qg data arrives:

| part | V_gp = 6 V | 8 V | 10 V |
|---|---|---|---|
| ST8L65 | 18.6 nC | 25.7 nC | 33.1 nC |
| IMCQ120 | 33.2 nC | 48.0 nC | 64.3 nC |
| ISC040 | 3.9 nC | 5.6 nC | 7.4 nC |

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

## Qg data contract (to digitize)

V_GS vs Q_g at specified V_DS, I_D (standard quasi-static gate-charge test).
CSV: `index,Q_g,V_gs,condition` (condition = `Vds=..V,Id=..A`). Extract:
V_gp (plateau voltage), Qgs (charge to plateau start), Qgd (plateau span).
Uses: Cmax extraction (§Cmax), Qgd vs ∫Crss dv cross-check, and (later) V_gp as
the anchor for any channel-charge model.

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

## Code TODO (implementation checklist — not yet done)

Phase 0 — data (not code):
- Digitize Qg curves for all three parts; write `.meta.txt` sidecars (figure
  number, frequency, test condition) for the three Cap CSVs.
- Clean the ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF
  points) on the next re-digitization pass.

Phase 1 — data plumbing (`data.py`):
- `CapCurveSet` (v_ds, C, branch ∈ {iss,oss,rss}) + `load_cap_csv()` — header
  contract `index,V_ds,Cap (pF),Type`, case-insensitive, group by Type, sort by
  v_ds.
- `load_qg_csv()` for `index,Q_g,V_gs,condition`.

Phase 2 — core (`caps.py`, new module):
- PCHIP-on-ln-C builder + per-segment GL quadrature → q(v), C(v)=dq/dv,
  clamp policy.
- Origin-anchor fit (quadratic, first 4 points of Crss) → c_gd0, s, b2.
- `c_gd(v_gd)` composite + closed-form q on the exp branch:
  `q_gd(v) = q_gd(0) + c_gdmax·v + (c_gdmax−c_gd0)/k_agd·(1 − e^{−k_agd·v})`.
- `CapParams` dataclass per part (c_gd0, s, c_gdmax; k_agd derived), bounds
  c_gdmax ∈ [c_gd0, 2·c_gs(0)]; log-space packing like `params.py`.

Phase 3 — extraction (`fit.py` or `fit_cmax.py`):
- `extract_cmax(qg_curve, cap_tables)`: plateau detection (V_gp, Qgd span),
  table-integral subtraction, 1-D solve for c_gdmax; negative-difference
  diagnostic → clamp + warning.
- Cross-check report: Qgd vs ∫Crss dv, Qoss vs spec.

Phase 4 — integration:
- Wire the delta (q_gd(v_gs−v_d), q_gs(v_gs), q_ds(v_ds)) into the device
  subckt / behavioral current sources i = dq/dv · dv/dt.
- LTspice export: `.subckt` with `table()` caps; parity vs the Python model.
- `reference.py`: `CAP_ST8L65` / `CAP_IMCQ120` / `CAP_ISC040` records + drift
  checks (diode protocol style).

Phase 5 — tests:
- Decomposition positivity; splice C0/C1 (numerical derivative jump < tol);
  dq/dv ≡ C identity; monotonicity through the SJ cliff; clamp behavior;
  trench junction anchor (M ≈ 0.72, V_bi ≈ 6.4, residual ≤ 10.6 %);
  c_gdmax bracket; origin-fit determinism (first-4-point window);
  pinned RMSE-style bounds on the Qgd cross-check once Qg data exists.
