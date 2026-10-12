# Project Status — ImprovedVDMOS

*Snapshot date: 2026-10-11 (blocking problems 1–2 resolved); gate-overdrive-cap
section appended after sessions 4–5. Companion to [README.md](README.md).*

## What exists

```
Projects/ImprovedVDMOS/
├── pyproject.toml          deps, pytest, ruff config, `vdmos-fit` console script
├── README.md               model + corrected Level-1 analysis
├── .gitignore
├── src/vdmos_fit/          model, params, data, fit, plotting, cli, __main__
├── tests/                  test_model, test_fit, test_package_r, test_gate_cap,
│                           test_cli, test_digitizer
├── scripts/                fit_datasheet.py, fit_st8l65.py, roundtrip_demo.py
├── data/                   README.md + ST8L65N044M9_IV.csv (41 digitised points)
└── notebooks/              README.md only (empty)
```

## What it currently does

- **Model** — `vdmos_id()` evaluates the unified continuous drain current
  (softplus subthreshold, gate-overdrive charge cap kcap/m, mobility
  degradation via Vk, velocity saturation, DIBL, gated CLM, drift conductance,
  one-sided kink sigmoid), embedded with package series resistance rd/rs
  (implicit external-node solve with stall guard, LTspice Rd/Rs semantics;
  Vds ≥ 0 only). `vdmos_id_level1()` is now the **correct** LTspice reference
  with knee at `Vov/mtriode` and takes the same rd/rs embedding for fair A/B.
- **Extraction** — `fit()` runs bounded `least_squares` with positive
  parameters in log space; two-stage relative→ampere weighting is wired into
  `scripts/fit_datasheet.py`.
- **Diagnostics** — `sensitivity_report()` (pre-fit, mean |d ln Id/d ln p|) and
  `identifiability()` (post-fit: Jacobian rank, condition number,
  per-parameter σ, correlated pairs). Zero-sensitivity columns report `inf`,
  not 0.
- **Data I/O** — flat `CurveSet` (vgs, vds, ids), CSV load/save, synthetic
  generator for round-trip validation.
- **Plots** — output curves, transfer curves, gds across knee, fit-vs-data
  overlay, and `mtriode_shape()` reproducing the LTspice mtriode sweep figure.
- **CLI** — `vdmos-fit check | fit <csv> | demo`, plus `python -m vdmos_fit`.
  Shared flags (`--show`, `--out`/`-o`) parse before *or* after the subcommand
  (fixed 2026-10-11 via a SUPPRESS-default parent parser; `-o` alias added to
  match the documented usage). Pinned by `tests/test_cli.py`.

## Verified state

- **86 tests passed, 1 skipped** (skip = digitizer test needing a Tk display;
  verified run before the cap commit: 78 pre-cap + 7 `test_gate_cap.py` +
  dead-column identifiability test + Newton two-cycle stall test).
  `ruff check` clean.
- All pre-package-R anchors unchanged (rd = rs = 0 fast path is bit-identical);
  all pre-cap anchors unchanged (kcap = 0 is the exact identity, pinned by
  `tests/test_gate_cap.py`):
  Ron 4.27 mΩ, saturation table 31.5/91.0/157.1/225.9/297.7 A, α = 1.08,
  Level-1 reference 960.4 A at Vgs = 10/Vds = 20.
- Round-trip fit recovers kp/vth0/vsat/eta/lam/gk within ~1 %, relative RMSE
  0.86 %.
- Corrected Level-1 matches LTspice 26.0.1 `.meas` exactly: Id = 4.000000 A
  either side of the knee for m = 1, 2, 4; min gds ≥ 0; Ron conductance =
  m·Kp·Vov; Id,sat independent of mtriode.

## Problems to fix

### Blocking / correctness — both resolved 2026-10-11

1. **"4–5× over-prediction" claim re-measured against the corrected Level-1.
   RESOLVED — the number survived.** Measured: **4.25× at Vgs = 10 V, 4.80× at
   11 V, 5.33× at 12 V** (Vds = 20 V; Level-1 defaults kp = 28, mtriode = 1.2,
   λ = 0.02 vs improved DEFAULT). "Matched Ron" needs no kp rescaling: the two
   Ron values agree within 0.15% at Vds = 50 mV across 10–12 V. The knee bug
   never entered this ratio — at Vds = 20 V (≫ Vov/m) both the old (knee = Vov)
   and corrected (knee = Vov/m) references evaluate the identical saturation
   branch ½·Kp·Vov²(1+λVds), and Ron = 1/(m·Kp·Vov) is the same in both. The
   earlier suspicion that the number would change was itself wrong; it was the
   gds/continuity claims that the bug corrupted. README updated to the precise
   "4.3× → 5.3×"; pinned by `test_matched_ron_square_law_overprediction`.
2. **DEFAULT re-derived. RESOLVED — parameters verified, not changed.**
   α = **1.078** over Vgs 6–12 V (1.098 over 10–12 V) for the current DEFAULT —
   README's "α ≈ 1.08" confirmed and now pinned to ±0.01 by
   `test_saturation_exponent_below_two`. Saturation table (31.5/91.0/157.1/
   225.9/297.7 A) and the 4.27 mΩ Ron anchor re-checked and current. DEFAULT is
   safe to keep because its only external anchor (LTspice A/B Ron at Vgs = 10)
   lives in the ohmic limit, where the knee bug had no effect. Tunability
   verified: α = 1.21 (vsat = 3 V), 1.32 (vsat = 5 V), 1.50 (vsat = 12 V), θ
   secondary. **Flag:** DEFAULT's α = 1.08 sits *below* the physical 1.2–1.6
   band; recalibrate vsat ≈ 3–5 V when real datasheet data arrives (problem 6).

### Known weak spots (measured, not guessed)

3. `theta` and `vsat` are **84 % correlated** (σ ≈ 9.7 and 0.69) — output
   curves alone cannot separate them.
4. `vthk` is barely sensed (σ ≈ 477 %, sensitivity ~2000× below kp) because
   `kk = 0.1 V` saturates the kink sigmoid above Vgs = 2 V. Needs dense Vgs
   points through the kink, or fix it from measurement.
5. `gds` σ ≈ 455 % — only the off-state constrains it, and the synthetic grid
   has no off-state points.

### Missing

6. **Real datasheet data. RESOLVED (session 4)** — `data/ST8L65N044M9_IV.csv`:
   41 digitised points from the ST8L65N044M9 (650 V super-junction, 36–44 mΩ)
   datasheet, 25 °C pulsed, Vds 0–20 V, Vgs 3.7–10 V, including the Vgs(th)
   and Rds(on) spec anchors. Provenance and the removed 650 V off-state point
   documented in `data/README.md`. First fits: see the cap section below.
7. No Q(V) charge model, no body diode, no thermal node, no `Tj` — the rest of
   the agreed architecture. The body diode is now the **biggest model hole**:
   the package-R path is Vds ≥ 0 only and reverse operation is undefined.
8. No LTspice export path: fitted parameters cannot yet be written back as a
   `.subckt`/`.model` for A/B simulation.
9. No CI. (The git repo now exists — `main` tracking `origin/main`.)
10. `notebooks/` is empty scaffolding.

## Key corrections made this session

| Claim | Before | After (LTspice-verified) |
|---|---|---|
| Level-1 knee location | `Vds = Vov` | `Vds = Vov / mtriode` |
| Level-1 continuity | "gds jump ~+400 %, −40 S NDR dip" | **C1 continuous for any mtriode**; artifacts were a bug in the project's own reference formula |
| Kink gate | bare `tanh((Vgs−Vthk)/kk)` (odd → negative off-state current) | one-sided `0.5·(1+tanh(…))` |
| Level-1 limit of improved model | "strict superset, exact" | smooth approximation: exact in ohmic + deep saturation, ~6 % low at the knee |
| `parameter_uncertainty` | zero Jacobian column → σ = 0 ("perfectly determined") | `identifiability()` → σ = `inf` ("NOT SENSED"), rank/cond/correlations |

## Package resistance added (rd, rs) — 2026 session 2

- **Embedding:** `Vds' = Vds − (rd+rs)·Id`, `Vgs' = Vgs − rs·Id`, implicit per
  point; shared bracketed safeguarded-Newton solver in `model.py`
  (`_solve_package_r`), used by both `vdmos_id` and `vdmos_id_level1`.
  Fixed-point iteration deliberately avoided: contraction factor
  (rd+rs)·kp·Vov ≈ 1.6 > 1 in the ohmic region for DEFAULT.
- **Fitted by default** (`DEFAULT_FREE` + 2), linear space, bounds [0, 0.5] Ω —
  datasheet Rds(on) is 2-terminal and already contains rd+rs.
- **Measured properties (pinned by `tests/test_package_r.py`):**
  - Solver matches independent `brentq` roots to 3e-16 relative; grid
    self-consistency residual < 1e-8.
  - Synthetic round trip (truth rd = 2 mΩ, rs = 1 mΩ): rd+rs recovered to
    0.3 %, individually −6.5 % / +14 %; identifiability reports
    **rd ↔ rs = −0.999** — only the sum is strongly pinned by static data.
  - Off state: Id = gds·Vds/(1 + gds·(rd+rs)) exactly.
  - Terminal Ron = intrinsic + rd + rs to within ~70 µΩ (chord-resistance
    curvature at the shifted internal Vds accounts for the residual).
  - When truth rd = rs = 0, the free rd/rs columns are near-parallel →
    identifiability cond rises 1e4 → 1.4e4 (test threshold raised to 1e5);
    the theta↔vsat > 0.8 pair drops out of `worst_pairs` with 11 free params
    (test now evaluates it on the legacy 9-param set).
- **Scope guard:** package-R path raises `ValueError` for Vds < 0. Discovered
  while testing: the *intrinsic* channel model is itself undefined in reverse
  (velocity-saturation divisor 1 + v/vsat singular near Vds ≈ −vsat; e.g.
  f(−20, 10) = +1038 A garbage). Body diode (problem 7) must define reverse
  operation; the fast path's finiteness sweep (`test_no_nan_or_inf_anywhere`)
  does not imply a valid reverse characteristic.
- `sensitivity_report` now probes bound-anchored parameters (rd/rs at 0)
  one-sided; the log-space formula is meaningless at p = 0 and a two-sided
  probe would request negative resistance.

## Gate-overdrive cap (kcap, m) — sessions 4+5

**Why.** The ST8L65N044M9 datasheet curves show a bimodal gate exponent
α = d ln Id/d ln Vov: **1.68 / 1.61 / 1.00 / 0.16 / 0.07** across
Vgs 5→6→7→8→9→10 at the top of each curve. No setting of θ/vsat/λ drops α by
an order of magnitude. Mechanism: finite free charge in the lightly doped SJ
pillars → soft cap on the effective gate overdrive,
`Vov_eff = Vov/(1 + (Vov/kcap)^m)^(1/m)`, applied **before Vk** (variant A:
charge and knee saturate together). `kcap = 0` is the exact identity (every
pre-cap anchor bit-identical); `0 ≤ dVov_eff/dVov ≤ 1` preserves the
package-R bracket invariant `dg/dId ≤ −1`.

**`kappa` deleted.** The previous mechanism (drain-bias-dependent series
resistance) was a voltage divider, not a charge limit: the fit made rd carry
75 % of the terminal voltage already at the Ron anchor (Vds = 0.483 V,
Id = 14 A). Cap beats kappa: **6.88 % vs 7.45 %** RMSE (session 4,
unconstrained lam), full Jacobian rank, ~1000× better conditioning.

**Solver stall guard.** With the cap, `rs·Id` can drag the internal gate
across the cap knee and make `g(Id)` flat-then-cliff: bracketed Newton
two-cycles between the flat shoulders (measured at Vds ≈ 178 V with fitted
params, rd = 26 mΩ) and failed at 60 iterations. A per-point stall counter
vetoes Newton after three consecutive <10 % bracket shrinks and forces
bisection. Changes no converged result; `test_no_ndr_with_fitted_cap` fails
without it, and `test_solver_survives_newton_two_cycle` pins it directly.

**Identifiability fix (`fit.py`).** The dead-column test is now RELATIVE to
the strongest Jacobian column (1e-5 of max), not absolute 1e-10×max, and
correlations/cond are computed over live columns only. On the ST8L65 fit the
ksub column measured 7.8e-7 of max — provably dead (pinning ksub at
0.007–0.09 V leaves RMSE and all other parameters identical) — yet the old
test called it live, reported a confident small σ, and manufactured the
ksub↔vk = +1.000 pair plus ~3 decades of cond inflation. Pinned by
`test_dead_parameter_reports_inf_not_a_small_sigma`.

**Fits (`scripts/fit_st8l65.py`, rewritten; outputs `out/st8l65_capfit.png`,
`out/st8l65_capfit_params.txt`).** Free = kp, vth0, ksub, lam, gds, rd, rs,
vk, kcap, m; fixed θ = 0, vsat = 1e4, η = 0, gk = 0, vthk = 1.6, kk = 0.1
(measured dead on this data); device bounds rd, rs ≤ 36 mΩ (datasheet
Rds(on)) and lam ≤ 0.004 1/V. Four multi-starts, all converge to the same
basin (the cap + capped lam made the landscape tame).

| | fit A (full, 41 pts) | fit B (no Vgs=5, 36 pts) |
|---|---|---|
| RMSE_rel | **9.17 %** | **5.54 %** |
| kp | 63.59 | 50.60 |
| vth0 | 4.6656 | 4.4460 |
| ksub | 0.02336 | 0.02645 |
| lam | **0.004 PINNED UPPER** | **0.004 PINNED UPPER** |
| gds | 6.766e-5 S | 6.757e-5 S |
| rd | 31.92 mΩ | 29.85 mΩ |
| rs | 3.8e-9 **PINNED LOWER** | 0.66 mΩ |
| vk | **4.245** | **5.562** |
| kcap | 2.610 V | 2.904 V |
| m | 6.203 | 10.241 |
| rd+rs | 31.92 mΩ | 30.51 mΩ (chord 34.5 mΩ) |
| plateau Vgs=10 @19.7 V | **226.4 A** (meas 240.2, −5.7 %) | **223.9 A** (−6.8 %) |
| rank / cond | 10/10 / 3.1e9 | 10/10 / 2.2e9 |
| Vgs≥6 only | 5.83 % | 5.62 % |

(rank/cond measured with the OLD identifiability tolerance; see the fix
above — re-measure after the next fit run.)

**Open problems from these fits (replace the session-4 list):**

1. **lam pinned at its 0.004 bound in both fits** — the data wants more output
   slope than the cap allows, so the plateau is now 5–7 % LOW (226/224 A vs
   240 A measured). A bias-dependent drain element is genuinely needed — but
   NOT kappa-as-is (pole + divider physics). Candidate: a bounded bias-
   dependent rd (e.g. `rd·(1 + kd·Vds)` with the same `dV/dI ≤ bound`
   discipline, or JFET-style depletion of rd). Not yet tried.
2. **ksub is dead on this dataset** (one off-state point only; see the
   identifiability fix above). Fix ksub from a subthreshold spec or tie
   vk↔Vk; re-run; cond should drop from ~2e9 back to ~1e4–1e5. The
   three-gate coupling of vk (DIBL/CLM/kink all use `tanh(Vds/vk)`) is still
   unexamined.
3. **vk rose 0.65 → 4.2–5.6 V** — much softer knee, consistent with the
   "knee too sharp" complaint; σ(vk) ≈ 112–150, still barely sensed.
4. **Vgs = 5 curve costs 21.3 % on fit A** — unreproducible by any C¹ model
   (identical Id at Vds = 1.127 and 5.153 V requires β = 0 over a 4.6×
   interval, then β ≈ 0.12). Fit B (5.54 %) is the fair model-quality number.
5. **rd/rs split still unidentifiable** (rs pinned at 0 in fit A; 0.66 mΩ in
   fit B) — quote only rd+rs ≈ 30–32 mΩ.

Note: the root `STATUS.md` deletion (the move of this file to
`src/vdmos_fit/`, committed as `d399939`) lands in the same cap commit.

