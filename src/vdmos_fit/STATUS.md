# Project Status — ImprovedVDMOS

*Snapshot: 2026-10-11 (session 8). Companion to [README.md](../../README.md).*

## Current state

Eight sessions, all committed to `main`. The package provides one C¹ continuous
drain-current model — softplus subthreshold, gate-overdrive charge cap (kcap, m),
mobility degradation via Vk, velocity saturation, bounded DIBL, gated CLM,
one-sided trench kink — embedded with package series resistance rd/rs (implicit
external-node solve with a Newton stall guard, LTspice Rd/Rs semantics,
Vds ≥ 0 only), plus bounded least-squares extraction with identifiability
diagnostics. Along the way the project's Level-1 reference was corrected (knee
at `Vov/mtriode`, verified against LTspice 26.0.1 — the earlier NDR/continuity
alarms were a bug in that reference, not in Level-1), and two mechanisms were
added then deleted as unsupported by the data: `kappa` (a drain-bias-dependent
series resistance that was a voltage divider, not a charge limit; replaced by
the cap) and a `gds·Vds` drift leak (a one-point parameter that flattened the
subthreshold curve to a gate-independent floor; the off state is now the
gate-dependent softplus floor). Diagnostics were hardened to match what the
data can actually support: dead Jacobian columns report σ = inf rather than a
confident small σ, and `on_bounds` flags pinned parameters relative to the
bound's magnitude, not the search-box width. Current fit on ST8L65N044M9:
**5.54 % relative RMSE (fit B), rank 9/9, cond 3.75e4, ksub live**, the
Vgs(th) spec point reproduced exactly by the channel. 89 tests green, ruff
clean. Datasheet Figure 3 (typical output curves) is on disk in `data/`.

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
│                           + ST8L65N044M9_Id_Vds_Vgs.png (datasheet Figure 3)
└── notebooks/              README.md only (empty)
```

## What it currently does

- **Model** — `vdmos_id()` evaluates the unified continuous drain current
  (softplus subthreshold, gate-overdrive charge cap kcap/m, mobility
  degradation via Vk, velocity saturation, DIBL, gated CLM, one-sided kink
  sigmoid), embedded with package series resistance rd/rs. No drift-leakage
  term: the off state is the gate-dependent softplus floor.
  `vdmos_id_level1()` is the **correct** LTspice reference with knee at
  `Vov/mtriode` and takes the same rd/rs embedding for fair A/B.
- **Extraction** — `fit()` runs bounded `least_squares` with positive
  parameters in log space; two-stage relative→ampere weighting is wired into
  `scripts/fit_datasheet.py`.
- **Diagnostics** — `sensitivity_report()` (pre-fit, mean |d ln Id/d ln p|) and
  `identifiability()` (post-fit: Jacobian rank, condition number,
  per-parameter σ, correlated pairs). Zero-sensitivity columns report `inf`,
  not 0. `FitResult.on_bounds()` flags bound-pinned parameters.
- **Data I/O** — flat `CurveSet` (vgs, vds, ids), CSV load/save, synthetic
  generator for round-trip validation.
- **Plots** — output curves, transfer curves, gds across knee, fit-vs-data
  overlay, and `mtriode_shape()` reproducing the LTspice mtriode sweep figure.
- **CLI** — `vdmos-fit check | fit <csv> | demo`, plus `python -m vdmos_fit`.
  Shared flags (`--show`, `--out`/`-o`) parse before *or* after the subcommand.
  Pinned by `tests/test_cli.py`.

## Verified state

- **89 tests, all passing**; `ruff check` clean. Caveat: 8 of the digitizer
  window tests drive Tk and self-skip when the process has no display — in
  headless runs 0–2 of them skip (flaky per-process), so "N passed, M skipped"
  totals vary; don't pin the count.
- All model anchors hold (rd = rs = 0 and kcap = 0 fast paths are
  bit-identical to the pre-package-R / pre-cap model):
  Ron 4.27 mΩ, saturation table 31.5/91.0/157.1/225.9/297.7 A, α = 1.08,
  Level-1 reference 960.4 A at Vgs = 10/Vds = 20, Level-1 C¹ continuity,
  matched-Ron over-prediction 4.25×/4.80×/5.33× at Vgs 10/11/12.
- Corrected Level-1 matches LTspice 26.0.1 `.meas` exactly: Id = 4.000000 A
  either side of the knee for m = 1, 2, 4; min gds ≥ 0; Ron conductance =
  m·Kp·Vov; Id,sat independent of mtriode.
- Synthetic round trips recover kp/vth0/vsat/eta/lam/gk within ~1 %
  (relative RMSE 0.86 %) and rd+rs to 0.3 % (individually −6.5 %/+14 %,
  rd ↔ rs correlation −0.999 — only the sum is pinned by static data).
- Package-R solver matches independent `brentq` roots to 3e-16 relative;
  grid self-consistency residual < 1e-8 relative, or < 1e-12 A absolute at
  the subthreshold floor (solver atol).

## Current ST8L65 fit (official, 9-param, gds-free)

`scripts/fit_st8l65.py`; outputs `out/st8l65_capfit.png`,
`out/st8l65_capfit_params.txt`. Free = kp, vth0, ksub, lam, rd, rs, vk,
kcap, m; fixed θ = 0, vsat = 1e4, η = 0, gk = 0, vthk = 1.6, kk = 0.1
(measured dead on this data); device bounds rd, rs ≤ 36 mΩ (datasheet
Rds(on)) and lam ≤ 0.004 1/V. Four multi-starts, all converge to the same
basin in both fits.

| | fit A (full, 41 pts) | fit B (no Vgs=5, 36 pts) |
|---|---|---|
| RMSE_rel | **9.97 %** | **5.54 %** |
| kp | 63.59 | 50.58 |
| vth0 | 4.6656 | 4.4459 |
| ksub | 0.0234 — dead (σ = inf) | **0.1836 — live, σ = 3.5** |
| lam | **0.004 PINNED UPPER** | **0.004 PINNED UPPER** |
| rd | 31.92 mΩ | 29.85 mΩ |
| rs | 3.6e-9 **PINNED LOWER** | 0.66 mΩ |
| vk | **4.244** | **5.562** |
| kcap | 2.610 V | 2.904 V |
| m | 6.203 | 10.241 |
| rd+rs | 31.92 mΩ | 30.50 mΩ (chord 34.5 mΩ) |
| plateau Vgs=10 @19.7 V | **226.4 A** (meas 240.2, −5.7 %) | **223.9 A** (−6.8 %) |
| Vgs(th) point (3.7, 3.7) | 0 µA — abandoned | **250.0 µA — exact, via channel** |
| rank / cond | **8/9** / 2.9e4 | **9/9** / 3.75e4 |
| Vgs≥6 only | 5.83 % | 5.62 % |

Fit B is the fair model-quality number: the Vgs = 5 curve is genuine
datasheet data but unreproducible by any C¹ model (identical Id at
Vds = 1.127 and 5.153 V requires β = 0 over a 4.6× interval, then
β ≈ 0.12), and it costs 21.3 % on fit A.

## Open problems

1. **lam pinned at its 0.004 bound in both fits** — the data wants more
   output slope than the cap allows, so the plateau is 5–7 % LOW (226/224 A
   vs 240 A measured). A bias-dependent drain element is genuinely needed —
   but NOT kappa-as-is (pole + divider physics). Candidate: a bounded
   bias-dependent rd (e.g. `rd·(1 + kd·Vds)` with the same `dV/dI ≤ bound`
   discipline, or JFET-style depletion of rd). Not yet tried.
2. **ksub one-point-anchored** — live in fit B (σ = 3.5) but ksub = 0.184 V
   implies ~420 mV/dec swing, not physical (~60–100 mV/dec); the single
   Vgs(th) point pins the *level*, not the *slope*. Dead in fit A (the Vgs=5
   curve fights the anchor). Fix the swing from a digitised near-threshold
   transfer curve; do NOT hard-fix ksub (pushes the anchor into vth0,
   corr +0.99).
3. **Knee shape vs datasheet Figure 3** — the figure is now on disk
   (`data/ST8L65N044M9_Id_Vds_Vgs.png`); the visual check of fit-B knee
   softness (vk = 5.56 V, σ ≈ 150, barely sensed) against it is still due.
4. **theta ↔ vsat 84 % correlated** (synthetic grid) — output curves alone
   cannot separate them; and `vthk` is barely sensed (σ ≈ 477 %, sensitivity
   ~2000× below kp) because kk = 0.1 V saturates the kink sigmoid above
   Vgs = 2 V. Fix from dense kink-region data or measurement.
5. **rd/rs split unidentifiable** (rs pinned at 0 in fit A; 0.66 mΩ in
   fit B) — quote only rd+rs ≈ 30–32 mΩ.
6. The three-gate coupling of vk (DIBL/CLM/kink all use `tanh(Vds/vk)`) is
   unexamined.

## Missing (agreed architecture, not yet built)

- **Body diode + HV leakage** — the biggest model hole: the package-R path
  is Vds ≥ 0 only and reverse operation is undefined; real high-voltage
  leakage belongs here, not in a fitted resistor.
- **Q(V) charge model** (charge-based capacitors via symbolic differentiation)
  and the `Tj` thermal node for Cauer-ladder co-simulation.
- **LTspice export**: fitted parameters → `.subckt`/`.model` for A/B.
- **CI**. (`notebooks/` is empty scaffolding.)
- Validation target remains double-pulse Eon/Eoff within ~10 %, which static
  fitting alone cannot deliver.
