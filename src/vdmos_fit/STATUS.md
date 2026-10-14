# Project Status — ImprovedVDMOS

*Snapshot: 2026-10-13 (session 10). Companion to [README.md](../../README.md).*

## Current state

Ten sessions, all committed to `main` (session 9 added the drain-knee softness
exponent `nd`; session 10 re-fitted ST8L65 with it — see the nd section below).
The package provides one C¹ continuous
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
bound's magnitude, not the search-box width. Current OFFICIAL fit on
ST8L65N044M9: **3.79 % relative RMSE (fit A, 10-param with `nd` free, full
41-point set, rank 8/10, cond 2.19e3)** — promoted this session; the 9-param
fit B (5.54 %, `scripts/fit_st8l65.py`) is kept below as the superseded
reference. 103 tests green, ruff clean. Datasheet Figure 3 (typical output
curves) is on disk in `data/`.

## What exists

```
Projects/ImprovedVDMOS/
├── pyproject.toml          deps, pytest, ruff config, `vdmos-fit` console script
├── README.md               model + corrected Level-1 analysis
├── .gitignore
├── src/vdmos_fit/          model, params, data, fit, plotting, reference
│                           (official device records), cli, __main__
├── tests/                  test_model, test_fit, test_package_r, test_gate_cap,
│                           test_drain_softness, test_st8l65_official, test_cli,
│                           test_digitizer
├── scripts/                fit_datasheet.py, fit_st8l65.py (superseded),
│                           fit_st8l65_nd.py (official recipe), roundtrip_demo.py
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

- **103 tests, all passing**; `ruff check` clean. Caveat: 8 of the digitizer
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

## Current ST8L65 fit — OFFICIAL: 10-param with `nd`, fit A (full data)

**Promoted 2026-10-13.** `scripts/fit_st8l65_nd.py` fit A; outputs
`out/st8l65_nd_params.txt`, `out/st8l65_nd_probe.txt`, `out/st8l65_nd_fit.png`,
`out/st8l65_nd_percurve.png`, `out/st8l65_nd_overlay.png`. The parameter record
lives in that script's `OFFICIAL` dict, which now re-exports
`vdmos_fit.reference.ST8L65_OFFICIAL` (single source of truth), and `main()`
checks each run reproduces it to 1e-3 relative (`OFFICIAL_TOL`) — if a future
refit drifts, the run says `[STALE]` and the record must be re-promoted
deliberately. `tests/test_st8l65_official.py` pins the same record against the
measured CSV (RMSE, Ron chord, plateau, Vgs(th), per-curve, nd interior, no
negative gds) so a model/bounds change that moves the official fit fails the
suite instead of silently invalidating these numbers.

Free = kp, vth0, ksub, lam, rd, rs, vk, kcap, m, **nd**; fixed θ = 0,
vsat = 1e4, η = 0, gk = 0, vthk = 1.6, kk = 0.1 (all measured dead on this
data); device bounds rd, rs ≤ 36 mΩ (datasheet Rds(on)). Global bounds
widened this session: kp ≤ 1e4, vsat ≤ 1e4, lam 1e-12–1.0, vk ≤ 1300,
nd ≥ 0.1. Four multi-starts, all four converge to the same basin.

| | **fit A (OFFICIAL, 41 pts)** | fit B (no Vgs=5, 36 pts) |
|---|---|---|
| RMSE_rel | **3.79 %** | 1.52 % |
| kp | 579.45 | 2354.2 (pinned UPPER 1e4) |
| vth0 | 5.0618 | 5.4469 |
| ksub | 0.2433 (σ_rel 0.95) | 0.2475 (σ_rel 2.07) |
| kcap | 1.0640 V | 0.5135 V |
| m | 7.2514 | 2.8434 |
| **nd** | **0.4359 — interior** | **0.4081 — interior** |
| lam | 5.1e-10 **PINNED LOWER 1e-12** | 1.9e-12 **PINNED LOWER** |
| rd / rs | 24.39 / 7.29 mΩ (both interior) | 27.15 / 6.68 mΩ |
| vk | 339.0 (σ = inf) | 0.2487 (σ = inf) |
| rd+rs | **31.68 mΩ** | 33.83 mΩ (chord 34.5 mΩ) |
| plateau Vgs=10 @19.7 V | **246.9 A (+2.8 %)** | 244.7 A (+1.9 %) |
| Vgs(th) point (3.7, 3.7) | 259.8 µA (meas 250) | 249.3 µA |
| rank / cond | **8/10** / 2.19e3 | 8/10 / 1.27e4 |
| Vgs≥6 only | 2.60 % | — |

**Why fit A and not the lower-RMSE fit B** (the promotion decision, made
deliberately): fit A uses **all 41 measured points** — no curve is dropped —
so its 3.79 % is the honest device-quality number, while fit B's 1.52 % is a
diagnostic that isolates how much error the physically unreproducible Vgs = 5
curve contributes. Fit A is also the better-conditioned and more physical
record: cond 2.19e3 vs 1.27e4, kp 579 vs 2354 (fit B's kp sits on the bound),
and rd/rs are both interior rather than leaning on the bound. Fit B stays in
the script and in the tables for comparison.

**lam and vk reporting σ = inf is expected, not a fit failure.** The dataset
stops at Vds = 20 V and never enters deep saturation, so it carries no CLM or
drain-gate information for those two to act on; `lam → 0` and `vk` either ≫ or
≪ the measured range both give tanh(Vds/vk) ≈ 1. What the data *does* pin is
the knee SHAPE (nd + kcap + m) and rd+rs. Confirmed by the lam-cap probe:
capping lam at the old (1e-6, 0.004) or freeing it to 1e-12 gives identical
RMSE and identical nd/kp/kcap/m to 4 significant figures — only vk moves
(1300 ↔ 0.249) to compensate.

**Model sanity at the official set:** finite everywhere on Vds 0–60 V ×
Vgs −2…16 V (max 293.7 A); min gds over Vds 0–200 V × Vgs −2…16 V is
−3.1e-12 S at Vgs = −0.8 / Vds = 95.7 V, i.e. finite-difference noise at the
1e-14 A subthreshold floor, not a physical NDR region. The vsat-divisor gds dip
does not apply (vsat fixed at 1e4).

## Superseded reference: 9-param fit (gds-free, nd = 0)

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

## ST8L65 fit with `nd` free (session 10) — now the OFFICIAL recipe

`scripts/fit_st8l65_nd.py`: same fit-A/fit-B recipe, `nd` added to the free
tuple (start 2.8 = legacy tanh; linear space, bound **0.1–40** — relaxed from
0.5 this session because the optimum sits below the old guard). Outputs
`out/st8l65_nd_params.txt`, `out/st8l65_nd_probe.txt` (pre-fit sensitivity +
identifiability + probes), `out/st8l65_nd_fit.png`,
`out/st8l65_nd_percurve.png` (Id vs Vds at each measured Vgs — the Figure 3
knee check), `out/st8l65_nd_overlay.png`. The fit-A column below **is** the
official record; the table is kept for the A-vs-B comparison.

| | nd fit A (41 pts, OFFICIAL) | nd fit B (36 pts) | superseded 9-param fit B (nd = 0) |
|---|---|---|---|
| RMSE_rel | **3.79 %** (was 9.97) | **1.52 %** (was 5.54) | 5.54 % |
| nd | **0.4359 — interior** | **0.4081 — interior** | 0 |
| lam | 5.1e-10 — pinned LOWER 1e-12 | 1.9e-12 — pinned LOWER | 0.004 pinned UPPER |
| plateau Vgs=10 @19.7 V | 246.9 A (+2.8 %) | 244.7 A (+1.9 %) | 223.9 A (−6.8 %) |
| Vgs(th) point (3.7, 3.7) | 259.8 µA | 249.3 µA | 250.0 µA |
| rd+rs | 31.68 mΩ | 33.83 mΩ | 30.50 mΩ |
| ksub | 0.2433 V | 0.2475 V | 0.1836 V |
| rank / cond | 8/10 / 2.19e3 | 8/10 / 1.27e4 | 9/9 / 3.75e4 |

Findings:

- **The data wants a knee SOFTER than the legacy tanh, and softer than the old
  0.5 guard allowed**: nd lands at 0.436 (fit A) / 0.408 (fit B), interior to
  the relaxed (0.1, 40) bound. The nd lower-bound probe reaches the same
  1.52 % from starts at nd = 0.5 and 0.3 — a stable interior optimum, so the
  0.1 floor is no longer load-bearing. nd → 0 collapses the current, so lower
  is not physical anyway; do not lower the floor further.
- **Open problem 1 was a knee-SHAPE problem, not a missing element**: with nd
  free the plateau lands +2.8 % instead of −6.8 % low. The bounded
  bias-dependent drain element (Queue 2) is NOT needed on this data. lam now
  pins at the LOWER 1e-12 — the data wants *zero* CLM once the knee shape is
  free, which is the expected answer for a 650 V part measured only to 20 V.
- gds ≥ 0 everywhere for the official fit (min −3e-12 S = FD noise at the
  deep-off floor); the vsat-divisor dip caveat does not apply (vsat = 1e4).
- **theta↔vsat probe** (theta, vsat freed on fit-B data): both push to their
  bounds (theta → 1e-4 LOWER, vsat → 1e4 UPPER) — the data rejects that
  mechanism outright. Freeing vsat up to 1e4 improved the no-nd fit from
  5.75 % to 4.05 % only because vsat ran to ~2800 V, i.e. divisor ≈ off.
  With nd free: 1.52 %, and nd is the correlation hub (kcap↔nd +0.949,
  kp↔nd −0.962, ksub↔nd +0.920) — the theta↔vsat pair is moot, not loosened.
- **Knee cluster (kp, kcap, m, nd, vk, lam) strongly entangled**
  (kp↔kcap −0.999, vth0↔kcap −0.985, kp↔vth0 +0.982); rank 8/10 with lam and
  vk dead (σ = inf). Quote the curve shape and rd+rs, never the individual
  knee parameters. The probes reach the same RMSE with vk = 0.25, 2.6, 339 or
  1300 — shape pinned, mechanism split not.
- ksub got worse as a one-point artifact (0.243 V ≈ 560 mV/dec). Queue 1
  (digitise a near-threshold transfer curve) is now the top data priority.

## Open problems

1. **lam / vk unsensed (σ = inf) in the official fit** — *accepted, not a
   defect*: the dataset stops at Vds = 20 V and never enters deep saturation,
   so there is no output-slope information for CLM or the drain gate to act on.
   lam pins at the LOWER 1e-12 (data wants zero CLM) and vk is free to be ≫ or
   ≪ the measured range with identical RMSE. Historical note: the old 9-param
   fit pinned lam at its 0.004 cap and that is why the plateau was 5–7 % LOW;
   **superseded by session 10** — with `nd` free the knee SHAPE supplies the
   slope and the plateau lands +2.8 %. Keep a bounded bias-dependent element
   (e.g. `rd·(1 + kd·Vds)` with the same `dV/dI ≤ bound` discipline) in
   reserve for datasets that *do* reach deep saturation; NOT kappa-as-is
   (pole + divider physics). If high-Vds data ever arrives, un-pin lam/vk and
   re-check before trusting any split of the knee cluster.
2. **ksub one-point-anchored** — live in the official fit (σ_rel 0.95) but
   ksub = 0.243 V implies ~560 mV/dec swing, not physical (~60–100 mV/dec);
   the single Vgs(th) point pins the *level*, not the *slope*. Fix the swing
   from a digitised near-threshold transfer curve; do NOT hard-fix ksub
   (pushes the anchor into vth0, corr +0.99).
3. **Knee shape vs datasheet Figure 3** — session 10 overlay plots
   (`out/st8l65_nd_percurve.png`) show the nd fit matching the knee at every
   measured Vgs; a formal digitisation of Figure 3 would quantify it, but the
   visual check is done.
4. **theta ↔ vsat 84 % correlated** (synthetic grid) — output curves alone
   cannot separate them; and `vthk` is barely sensed (σ ≈ 477 %, sensitivity
   ~2000× below kp) because kk = 0.1 V saturates the kink sigmoid above
   Vgs = 2 V. Fix from dense kink-region data or measurement.
5. **rd/rs split unidentifiable** — quote only rd+rs ≈ 32 mΩ. In the official
   fit both are interior (rd 24.4 / rs 7.3 mΩ) rather than pinned, but that
   split still leans entirely on rs's gm-degeneration signature, not on
   independent evidence.
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
