# Project Status — ImprovedVDMOS

*Snapshot: 2026-10-24 (session 20 — CSV header contract: `load_csv()` now
accepts BOTH the legacy `vgs,vds,ids` and the digitizer `index,V_ds,I_d,V_gs`
spelling, HANDOVER queue item 6 CLOSED; session 19 landed the
interface-trap-filling gain `v_ov_eff2`/`k_beta_trap` in the channel model,
off by default — committed as `8574ce2`; session 18 investigated the
IMCQ120R004M2H two-temperature IV fits and the gain probe — committed as
`3385ee2`, findings in `out/IMCQ120_TEMPERATURE_REPORT.md` + HANDOVER.md;
session 17 landed the capacitance/charge model — see CAP_SPEC.md).
Companion to [README.md](../../README.md).*

## Current state

Twenty sessions, committed to `main` (session 18 added the IMCQ120R004M2H
SiC IV data + two-temperature fits + the gain probe (`3385ee2`); session 19
landed the interface-trap-filling gain `v_ov_eff2`/`k_beta_trap` in the channel
model (`8574ce2`); session 20 closed the CSV header contract — `load_csv()`
reads both header spellings natively, see the session-20 section below;
session 11 renamed every parameter
to the `<type prefix>_<suffix>` convention — names.xlsx + README "Naming
convention" are the mapping, tag `pre-rename-2026-10` is the escape hatch, and
the rename is bit-identical; session 12 added the **body diode** —
`diode.py`, LTspice diode DC branch verified against LTspice 26.0.1,
terminal-wired bypassing r_d/r_s, multi-temperature extraction, official
record `ST8L65_DIODE`; session 13 closed session 12's pending list — test
fixes, headless LTspice batch-mode verification, docs). Session 14 added the
diode **linear temperature coefficients** `a_ikf`/`a_rdiode` (LTspice
Tikf/Trs1, Tnom-referenced, verified to ~1e-5 against LTspice 26.0.1;
TC=0 default is bit-identical to session 12) and **promoted the 6-param TC
fit (4.15 %) to the official record**. Session 15 added the **quadratic
r_diode TC `b_rdiode`** (LTspice Trs2 — completing the exact LTspice Rs law,
probed to ~1e-5 at −55/25/100/150 °C) and **promoted the 7-param TC+Trs2 fit
(3.91 %) to the official record (user override)**. The promotion caveats are
recorded in the Trs2 section below: the gain over the linear-TC record is
modest (4.15 → 3.91 %), the five-knob TC/knee block is strongly correlated
(|corr| up to 0.985), and the stubborn 150 °C residual is a ~12 mV
digitization error on the 5 A point, not missing series-R physics.
Session 16 settled the **capacitance/charge design** (CAP_SPEC.md) and
digitized the C–V + Qg data for three parts; session 17 **landed the model**
(`caps.py`: PCHIP-on-lnC charge-master LUTs, first-4-point origin anchor,
exp-saturation splice, `extract_cmax()` from the Qg curves — headline:
all three parts extract at c_gdmax ≈ c_gd0, the seed over-predicts 4–8×,
and the Miller-plateau deficit quantifies limitation 1) plus the
**LTspice export** (`ltspice.py`: charge-master nonlinear capacitors,
Q= ramp-sum of sampled q(v); live parity ≤0.3 % on all three parts).
The package provides one C¹ continuous
drain-current model — softplus subthreshold, gate-overdrive charge cap (v_ovmax, k_ovsharp),
interface-trap-filling gain (k_beta_trap) on v_ov_eff2, mobility degradation via v_knee,
velocity saturation, bounded DIBL, gated CLM,
one-sided trench kink — embedded with package series resistance r_d/r_s (implicit
external-node solve with a Newton stall guard, LTspice Rd/Rs semantics,
channel v_ds ≥ 0 only; reverse polarity is carried by the body diode), plus bounded least-squares extraction with identifiability
diagnostics. Along the way the project's Level-1 reference was corrected (knee
at `v_ov/mtriode`, verified against LTspice 26.0.1 — the earlier NDR/continuity
alarms were a bug in that reference, not in Level-1), and two mechanisms were
added then deleted as unsupported by the data: `kappa` (a drain-bias-dependent
series resistance that was a voltage divider, not a charge limit; replaced by
the cap) and a `gds·v_ds` drift leak (a one-point parameter that flattened the
subthreshold curve to a gate-independent floor; the off state is now the
gate-dependent softplus floor, plus the diode's reverse saturation when the
diode is on). Diagnostics were hardened to match what the data can actually
support: dead Jacobian columns report σ = inf rather than a confident small σ,
and `on_bounds` flags pinned parameters relative to the bound's magnitude, not
the search-box width. Current OFFICIAL fits on
ST8L65N044M9: channel **3.79 % relative RMSE (fit A, 10-param with `k_kneesharp` free, full
41-point set, rank 8/10, cond 2.19e3)**; diode **4.78 % on 30 points
(−55/25/150 °C), rank 4/4, cond 6.1e3, v_eg = 1.097 eV**. 154 tests green,
ruff clean. Datasheet Figure 3 (typical output curves), the body-diode
forward characteristics CSV, and the three-part C–V + Qg CSVs are on disk
in `data/`.

## What exists

```
Projects/ImprovedVDMOS/
├── pyproject.toml          deps, pytest, ruff config, `vdmos-fit` console script
├── README.md               model + corrected Level-1 analysis
├── .gitignore
├── src/vdmos_fit/          model, diode, caps, params, data, fit, plotting,
│                           reference (official device records), ltspice export,
│                           cli, __main__
├── tests/                  test_model, test_fit, test_package_r, test_gate_cap,
│                           test_drain_softness, test_st8l65_official, test_diode,
│                           test_caps (incl. live LTspice parity), test_cli,
│                           test_digitizer
├── scripts/                fit_datasheet.py, fit_st8l65.py (superseded),
│                           fit_st8l65_nd.py (official channel recipe),
│                           fit_st8l65_diode.py (official diode recipe),
│                           fit_cmax.py (official c_gdmax extraction),
│                           export_caps.py (LTspice subckt export),
│                           roundtrip_demo.py
├── data/                   README.md + ST8L65N044M9_IV.csv (41 digitised points)
│                           + ST8L65N044M9_diode.csv (30 pts, −55/25/150 °C)
│                           + ST8L65N044M9_Id_Vds_Vgs.png (datasheet Figure 3)
│                           + <part>_Cap_Vds.csv ×3 + <part>_Vgs_Qg.csv ×3
└── notebooks/              README.md only (empty)
```

## What it currently does

- **Model** — `vdmos_id()` evaluates the unified continuous drain current
  (softplus subthreshold, gate-overdrive charge cap v_ovmax/k_ovsharp, mobility
  degradation via v_knee, velocity saturation, DIBL, gated CLM, one-sided kink
  sigmoid), embedded with package series resistance r_d/r_s, minus the body
  diode `vdmos_diode()` (LTspice diode DC branch: diffusion + high-injection
  compression + recombination + r_diode + Xti/Eg temperature law + LTspice
  Tikf/Trs1/Trs2 temperature coefficients on i_kf/r_diode), wired at
  the terminals and bypassing r_d/r_s. No drift-leakage term: the off state
  is the gate-dependent softplus floor (channel) plus the diode reverse
  saturation when the diode is on.
  `vdmos_id_level1()` is the **correct** LTspice reference with knee at
  `v_ov/mtriode` and takes the same r_d/r_s embedding for fair A/B.
- **Extraction** — `fit()` runs bounded `least_squares` with positive
  parameters in log space; two-stage relative→ampere weighting is wired into
  `scripts/fit_datasheet.py`.
- **Diagnostics** — `sensitivity_report()` (pre-fit, mean |d ln i_d/d ln p|) and
  `identifiability()` (post-fit: Jacobian rank, condition number,
  per-parameter σ, correlated pairs). Zero-sensitivity columns report `inf`,
  not 0. `FitResult.on_bounds()` flags bound-pinned parameters.
- **Data I/O** — flat `CurveSet` (v_gs, v_ds, i_d) and `DiodeCurveSet`
  (i_sd, v_sd, temp_c), CSV load/save, synthetic generators for round-trip
  validation.
- **Plots** — output curves, transfer curves, gds across knee, fit-vs-data
  overlay, and `mtriode_shape()` reproducing the LTspice mtriode sweep figure.
- **CLI** — `vdmos-fit check | fit <csv> | demo`, plus `python -m vdmos_fit`.
  Shared flags (`--show`, `--out`/`-o`) parse before *or* after the subcommand.
  Pinned by `tests/test_cli.py`.

## Verified state

- **145 tests, all passing** (0–2 Tk digitizer window tests self-skip when
  the process has no display); `ruff check` clean. Caveat: 8 of the digitizer
  window tests drive Tk and self-skip when the process has no display — in
  headless runs 0–2 of them skip (flaky per-process), so "N passed, M skipped"
  totals vary; don't pin the count.
- **Body diode verified against LTspice 26.0.1**
  (`test_diode.py::test_ltspice_cross_validation`, runs when LTspice is
  installed): diffusion +
  Ikf + Rs match to ~1e-3 at 0.3–0.9 V (the equation itself matches to ~1e-5
  across 3.4 A … 375 kA), and the Xti/Eg temperature law (divisor N per
  branch, Nr for Isr) matches at −55 °C to ~1 %. LTspice's own Isr branch
  deviates from its documented form by 1–17 % below 0.3 V (measured;
  documented in diode.py/README) — excluded from the cross-validation scope.
- **Diode off-sentinel bit-identity**: with i_s = i_sr = 0 (default),
  `vdmos_id` is bit-identical to the pre-diode model on the v_ds ≥ 0 grid
  (fast path and package-R path both pinned).
- **Diode TC laws verified against LTspice 26.0.1**
  (`test_diode.py::test_ltspice_tc_laws`): Tikf linear and Trs1+Trs2
  quadratic, all **Tnom-referenced** — agreement ~1e-5 relative at
  −55/25/100/150 °C × 0.6/0.9 V; a 27 °C-referenced variant is excluded
  (0.1–0.6 % off), TC-off excluded (up to 25 %), linear-only excluded
  (up to 10.6 % at 150 °C). TC=0 default is bit-identical to the session-12
  model; TC round trip recovers a_ikf/a_rdiode/b_rdiode exactly (unit-column
  SVD rank 7/7; identifiability's raw mixed-unit report says 6/7 -- a
  tolerance artifact, see the Trs2 section); the
  extrapolation clamp keeps r_diode/i_kf positive at absurd temperatures.
- All model anchors hold (r_d = r_s = 0 and v_ovmax = 0 fast paths are
  bit-identical to the pre-package-R / pre-cap model):
  Ron 4.27 mΩ, saturation table 31.5/91.0/157.1/225.9/297.7 A, α = 1.08,
  Level-1 reference 960.4 A at v_gs = 10/v_ds = 20, Level-1 C¹ continuity,
  matched-Ron over-prediction 4.25×/4.80×/5.33× at v_gs 10/11/12.
- Corrected Level-1 matches LTspice 26.0.1 `.meas` exactly: i_d = 4.000000 A
  either side of the knee for m = 1, 2, 4; min gds ≥ 0; Ron conductance =
  m·Kp·v_ov; i_d,sat independent of mtriode.
- Synthetic round trips recover k_p/v_th0/v_vsat/k_eta/k_lambda/g_k within ~1 %
  (relative RMSE 0.86 %) and r_d+r_s to 0.3 % (individually −6.5 %/+14 %,
  r_d ↔ r_s correlation −0.999 — only the sum is pinned by static data).
- Package-R solver matches independent `brentq` roots to 3e-16 relative;
  grid self-consistency residual < 1e-8 relative, or < 1e-12 A absolute at
  the subthreshold floor (solver atol).

## Current ST8L65 fit — OFFICIAL: 10-param with `k_kneesharp`, fit A (full data)

**Promoted 2026-10-13.** `scripts/fit_st8l65_nd.py` fit A; outputs
`out/st8l65_nd_params.txt`, `out/st8l65_nd_probe.txt`, `out/st8l65_nd_fit.png`,
`out/st8l65_nd_percurve.png`, `out/st8l65_nd_overlay.png`. The parameter record
lives in that script's `OFFICIAL` dict, which now re-exports
`vdmos_fit.reference.ST8L65_OFFICIAL` (single source of truth), and `main()`
checks each run reproduces it to 1e-3 relative (`OFFICIAL_TOL`) — if a future
refit drifts, the run says `[STALE]` and the record must be re-promoted
deliberately. `tests/test_st8l65_official.py` pins the same record against the
measured CSV (RMSE, Ron chord, plateau, v_gs(th), per-curve, k_kneesharp interior, no
negative gds) so a model/bounds change that moves the official fit fails the
suite instead of silently invalidating these numbers.

Free = k_p, v_th0, v_sub, k_lambda, r_d, r_s, v_dbsat, v_ovmax, k_ovsharp, **k_kneesharp**; fixed k_theta = 0,
v_vsat = 1e4, k_eta = 0, g_k = 0, v_thk = 1.6, v_kovsoft = 0.1 (all measured dead on this
data); device bounds r_d, r_s ≤ 36 mΩ (datasheet Rds(on)). Global bounds
widened this session: k_p ≤ 1e4, v_vsat ≤ 1e4, k_lambda 1e-12–1.0, v_dbsat ≤ 1300,
k_kneesharp ≥ 0.1. Four multi-starts, all four converge to the same basin.

| | **fit A (OFFICIAL, 41 pts)** | fit B (no v_gs=5, 36 pts) |
|---|---|---|
| RMSE_rel | **3.79 %** | 1.52 % |
| k_p | 579.45 | 2354.2 (pinned UPPER 1e4) |
| v_th0 | 5.0618 | 5.4469 |
| v_sub | 0.2433 (σ_rel 0.95) | 0.2475 (σ_rel 2.07) |
| v_ovmax | 1.0640 V | 0.5135 V |
| k_ovsharp | 7.2514 | 2.8434 |
| **k_kneesharp** | **0.4359 — interior** | **0.4081 — interior** |
| k_lambda | 5.1e-10 **PINNED LOWER 1e-12** | 1.9e-12 **PINNED LOWER** |
| r_d / r_s | 24.39 / 7.29 mΩ (both interior) | 27.15 / 6.68 mΩ |
| v_dbsat | 339.0 (σ = inf) | 0.2487 (σ = inf) |
| r_d+r_s | **31.68 mΩ** | 33.83 mΩ (chord 34.5 mΩ) |
| plateau v_gs=10 @19.7 V | **246.9 A (+2.8 %)** | 244.7 A (+1.9 %) |
| v_gs(th) point (3.7, 3.7) | 259.8 µA (meas 250) | 249.3 µA |
| rank / cond | **8/10** / 2.19e3 | 8/10 / 1.27e4 |
| v_gs≥6 only | 2.60 % | — |

**Why fit A and not the lower-RMSE fit B** (the promotion decision, made
deliberately): fit A uses **all 41 measured points** — no curve is dropped —
so its 3.79 % is the honest device-quality number, while fit B's 1.52 % is a
diagnostic that isolates how much error the physically unreproducible v_gs = 5
curve contributes. Fit A is also the better-conditioned and more physical
record: cond 2.19e3 vs 1.27e4, k_p 579 vs 2354 (fit B's k_p sits on the bound),
and r_d/r_s are both interior rather than leaning on the bound. Fit B stays in
the script and in the tables for comparison.

**k_lambda and v_dbsat reporting σ = inf is expected, not a fit failure.** The dataset
stops at v_ds = 20 V and never enters deep saturation, so it carries no CLM or
drain-gate information for those two to act on; `k_lambda → 0` and `v_dbsat` either ≫ or
≪ the measured range both give tanh(v_ds/v_dbsat) ≈ 1. What the data *does* pin is
the knee SHAPE (k_kneesharp + v_ovmax + k_ovsharp) and r_d+r_s. Confirmed by the k_lambda-cap probe:
capping k_lambda at the old (1e-6, 0.004) or freeing it to 1e-12 gives identical
RMSE and identical k_kneesharp/k_p/v_ovmax/k_ovsharp to 4 significant figures — only v_dbsat moves
(1300 ↔ 0.249) to compensate.

**Model sanity at the official set:** finite everywhere on v_ds 0–60 V ×
v_gs −2…16 V (max 293.7 A); min gds over v_ds 0–200 V × v_gs −2…16 V is
−3.1e-12 S at v_gs = −0.8 / v_ds = 95.7 V, i.e. finite-difference noise at the
1e-14 A subthreshold floor, not a physical NDR region. The v_vsat-divisor gds dip
does not apply (v_vsat fixed at 1e4).

## Superseded reference: 9-param fit (gds-free, k_kneesharp = 0)

`scripts/fit_st8l65.py`; outputs `out/st8l65_capfit.png`,
`out/st8l65_capfit_params.txt`. Free = k_p, v_th0, v_sub, k_lambda, r_d, r_s, v_dbsat,
v_ovmax, k_ovsharp; fixed k_theta = 0, v_vsat = 1e4, k_eta = 0, g_k = 0, v_thk = 1.6, v_kovsoft = 0.1
(measured dead on this data); device bounds r_d, r_s ≤ 36 mΩ (datasheet
Rds(on)) and k_lambda ≤ 0.004 1/V. Four multi-starts, all converge to the same
basin in both fits.

| | fit A (full, 41 pts) | fit B (no v_gs=5, 36 pts) |
|---|---|---|
| RMSE_rel | **9.97 %** | **5.54 %** |
| k_p | 63.59 | 50.58 |
| v_th0 | 4.6656 | 4.4459 |
| v_sub | 0.0234 — dead (σ = inf) | **0.1836 — live, σ = 3.5** |
| k_lambda | **0.004 PINNED UPPER** | **0.004 PINNED UPPER** |
| r_d | 31.92 mΩ | 29.85 mΩ |
| r_s | 3.6e-9 **PINNED LOWER** | 0.66 mΩ |
| v_dbsat | **4.244** | **5.562** |
| v_ovmax | 2.610 V | 2.904 V |
| k_ovsharp | 6.203 | 10.241 |
| r_d+r_s | 31.92 mΩ | 30.50 mΩ (chord 34.5 mΩ) |
| plateau v_gs=10 @19.7 V | **226.4 A** (meas 240.2, −5.7 %) | **223.9 A** (−6.8 %) |
| v_gs(th) point (3.7, 3.7) | 0 µA — abandoned | **250.0 µA — exact, via channel** |
| rank / cond | **8/9** / 2.9e4 | **9/9** / 3.75e4 |
| v_gs≥6 only | 5.83 % | 5.62 % |

Fit B is the fair model-quality number: the v_gs = 5 curve is genuine
datasheet data but unreproducible by any C¹ model (identical i_d at
v_ds = 1.127 and 5.153 V requires β = 0 over a 4.6× interval, then
β ≈ 0.12), and it costs 21.3 % on fit A.

## ST8L65 fit with `k_kneesharp` free (session 10) — now the OFFICIAL recipe

`scripts/fit_st8l65_nd.py`: same fit-A/fit-B recipe, `k_kneesharp` added to the free
tuple (start 2.8 = legacy tanh; linear space, bound **0.1–40** — relaxed from
0.5 this session because the optimum sits below the old guard). Outputs
`out/st8l65_nd_params.txt`, `out/st8l65_nd_probe.txt` (pre-fit sensitivity +
identifiability + probes), `out/st8l65_nd_fit.png`,
`out/st8l65_nd_percurve.png` (i_d vs v_ds at each measured v_gs — the Figure 3
knee check), `out/st8l65_nd_overlay.png`. The fit-A column below **is** the
official record; the table is kept for the A-vs-B comparison.

| | k_kneesharp fit A (41 pts, OFFICIAL) | k_kneesharp fit B (36 pts) | superseded 9-param fit B (k_kneesharp = 0) |
|---|---|---|---|
| RMSE_rel | **3.79 %** (was 9.97) | **1.52 %** (was 5.54) | 5.54 % |
| k_kneesharp | **0.4359 — interior** | **0.4081 — interior** | 0 |
| k_lambda | 5.1e-10 — pinned LOWER 1e-12 | 1.9e-12 — pinned LOWER | 0.004 pinned UPPER |
| plateau v_gs=10 @19.7 V | 246.9 A (+2.8 %) | 244.7 A (+1.9 %) | 223.9 A (−6.8 %) |
| v_gs(th) point (3.7, 3.7) | 259.8 µA | 249.3 µA | 250.0 µA |
| r_d+r_s | 31.68 mΩ | 33.83 mΩ | 30.50 mΩ |
| v_sub | 0.2433 V | 0.2475 V | 0.1836 V |
| rank / cond | 8/10 / 2.19e3 | 8/10 / 1.27e4 | 9/9 / 3.75e4 |

Findings:

- **The data wants a knee SOFTER than the legacy tanh, and softer than the old
  0.5 guard allowed**: k_kneesharp lands at 0.436 (fit A) / 0.408 (fit B), interior to
  the relaxed (0.1, 40) bound. The k_kneesharp lower-bound probe reaches the same
  1.52 % from starts at k_kneesharp = 0.5 and 0.3 — a stable interior optimum, so the
  0.1 floor is no longer load-bearing. k_kneesharp → 0 collapses the current, so lower
  is not physical anyway; do not lower the floor further.
- **Open problem 1 was a knee-SHAPE problem, not a missing element**: with k_kneesharp
  free the plateau lands +2.8 % instead of −6.8 % low. The bounded
  bias-dependent drain element (Queue 2) is NOT needed on this data. k_lambda now
  pins at the LOWER 1e-12 — the data wants *zero* CLM once the knee shape is
  free, which is the expected answer for a 650 V part measured only to 20 V.
- gds ≥ 0 everywhere for the official fit (min −3e-12 S = FD noise at the
  deep-off floor); the v_vsat-divisor dip caveat does not apply (v_vsat = 1e4).
- **k_theta↔v_vsat probe** (k_theta, v_vsat freed on fit-B data): both push to their
  bounds (k_theta → 1e-4 LOWER, v_vsat → 1e4 UPPER) — the data rejects that
  mechanism outright. Freeing v_vsat up to 1e4 improved the no-k_kneesharp fit from
  5.75 % to 4.05 % only because v_vsat ran to ~2800 V, i.e. divisor ≈ off.
  With k_kneesharp free: 1.52 %, and k_kneesharp is the correlation hub (v_ovmax↔k_kneesharp +0.949,
  k_p↔k_kneesharp −0.962, v_sub↔k_kneesharp +0.920) — the k_theta↔v_vsat pair is moot, not loosened.
- **Knee cluster (k_p, v_ovmax, k_ovsharp, k_kneesharp, v_dbsat, k_lambda) strongly entangled**
  (k_p↔v_ovmax −0.999, v_th0↔v_ovmax −0.985, k_p↔v_th0 +0.982); rank 8/10 with k_lambda and
  v_dbsat dead (σ = inf). Quote the curve shape and r_d+r_s, never the individual
  knee parameters. The probes reach the same RMSE with v_dbsat = 0.25, 2.6, 339 or
  1300 — shape pinned, mechanism split not.
- v_sub got worse as a one-point artifact (0.243 V ≈ 560 mV/dec). Queue 1
  (digitise a near-threshold transfer curve) is now the top data priority.

## Current ST8L65 diode fit — OFFICIAL: fit A (7 free params incl. all three TCs, session 15)

`scripts/fit_st8l65_diode.py` fit A; outputs `out/st8l65_diode_params.txt`,
`out/st8l65_diode_probe.txt`, `out/st8l65_diode_fit.png`. Record in
`vdmos_fit.reference.ST8L65_DIODE`, drift-checked by the script
(`OFFICIAL_TOL`, 3.8e-5 [OK] at promotion) and pinned by
`tests/test_diode.py::test_st8l65_diode_record_reproduces`
(RMSE < 4.0 %, per-temperature < 7 %, r_diode interior, reverse floor < 10 µA,
v_eg in (0.9, 1.3) eV, TCs interior — linear ones positive, b_rdiode negative).

Data: `data/ST8L65N044M9_diode.csv` — 30 points, I_sd 5–55 A at −55/25/150 °C.
Free = i_s, r_diode, v_eg, i_kf, a_ikf, a_rdiode, **b_rdiode** (Trs2 promoted
session 15 — see the TC section below); **k_n fixed 2.0** (physical diffusion
limit),
k_xti fixed 3.0 (LTspice default), i_sr = 0, t_nom = 25 °C; device bound
r_diode ≤ 20 mΩ. Four multi-starts, all converge to the same basin.

| | fit A (OFFICIAL, 30 pts) |
|---|---|
| RMSE_rel | **3.91 %** (−55 °C 2.37, 25 °C 3.16, 150 °C 5.49) |
| i_s | 3.851e-6 A (σ_rel 1.11 — weak prefactor, pair with k_n) |
| r_diode | 1.896 mΩ (σ_rel 0.013) |
| v_eg | 1.0820 eV (σ_rel 0.18 — silicon bandgap, sensed) |
| i_kf | 21.97 A (σ_rel 15.4; r_diode↔i_kf corr +0.985) |
| a_ikf | +0.00832 /°C (σ_rel 0.059) |
| a_rdiode | +0.00483 /°C (σ_rel 0.022) |
| b_rdiode | −1.78e-5 /°C² (σ_rel 0.0003 — NOT quotable, see TC section) |
| rank / cond (as reported) | 6/7* / 1.64e4 |
| reverse floor (−1 V, 25 °C) | 3.9 µA |

Findings / protocol decisions:

- **The temperature law is mandatory**: a 25 °C-only fit predicts −55 °C
  out-of-sample at ~260 % error; V_T-only scaling (no Xti/Eg) predicts the
  forward drop going UP with T — the data goes DOWN. v_eg ≈ 1.09 eV is
  genuinely sensed (σ_rel 0.16) — the multi-temperature data pins it.
- **k_n fixed at 2.0 from physics**: the k_n sweep (i_kf off) shows RMSE
  improving monotonically toward k_n ≈ 2.5–3 — ideality absorbing the
  high-injection bending, with the reverse floor scaling with i_s (3 pA at
  k_n = 1 → 0.4 mA at k_n = 3). Freeing k_n on top of the TC fit lands at
  k_n = 1.31 (3.00 %, corr(i_s, k_n) = +0.996, reverse floor drifts to
  3.9 nA — TC-off it was 2.07 / 4.66 % / 0.1–0.4 mA): the valley's DIRECTION
  moves with the TCs, the unconstrained floor does not — k_n stays fixed
  from physics, only the (k_n, i_s) pair and the curve are pinned.
- **i_kf belongs in the primary set**: pinned off (sentinel 0.0 — NOT the
  1e-3 lower bound, which is an active knee that crushes the current) the
  fit is 6.63 %; free, 4.15 %. The high-injection knee carries real
  curvature in the 5–55 A range.
- **i_sr/k_nr unsensed**: freeing them gives identical RMSE with
  σ(i_sr) = inf — recombination lives below the 5 A data floor. Keep off.
- **i_s is a weak prefactor** (σ_rel ≈ 1.1): quote the curve or the pair
  (k_n, i_s), never i_s alone — same doctrine as the channel knee cluster.
- **Log-space packing is mandatory** for the diode fit: linear-space i_s
  gave rank-deficient Jacobians and a runaway basin at r_diode ≈ 0.34 Ω;
  log-space gives rank 4/4 and r_diode stable across starts. The per-device
  r_diode ≤ 20 mΩ bound (measured dV/dI is 3.4–11.7 mΩ) excludes the
  runaway even so.

## Diode temperature coefficients (sessions 14–15) — all three TCs promoted to the official record

`a_ikf` (LTspice **Tikf**), `a_rdiode` (**Trs1**) and `b_rdiode` (**Trs2**,
session 15): the exact LTspice laws `i_kf(T) = i_kf·(1 + a_ikf·(T − t_nom))`
and `r_diode(T) = r_diode·(1 + a_rdiode·(T − t_nom) + b_rdiode·(T − t_nom)²)`,
Tnom-referenced — pinned by `test_diode.py::test_ltspice_tc_laws` (~1e-5 at
−55/25/100/150 °C; linear-only misses up to 10.6 % at 150 °C, proving the
Trs2 term is real LTspice behavior). Defaults 0 = off =
bit-identical to the session-12 model; bounds (−0.005, +0.05) /°C and
(−1e-4, +1e-3) /°C², with a tiny clamp floor so the solver bracket proof
survives absurd extrapolation.

TC probe on the measured CSV (`scripts/fit_st8l65_diode.py`):

| | fit A (OFFICIAL, TC off) | fit A + a_ikf + a_rdiode |
|---|---|---|
| RMSE_rel | 4.78 % | **4.15 %** |
| −55 / 25 / 150 °C rms | 3.40 / 3.85 / 6.49 % | 2.26 / 3.41 / 5.90 % |
| a_ikf | 0 (fixed) | +0.00545 /°C (σ_rel 0.028) |
| a_rdiode | 0 (fixed) | +0.00341 /°C (σ_rel 0.019) |
| i_kf / r_diode | 17.45 A / 1.586 mΩ | 15.64 A / 1.497 mΩ |
| v_eg | 1.097 eV | 1.088 eV |
| rank / cond | 4/4 / 6.13e3 | 6/6 / 7.12e3 |

Findings:

- **The user-flagged effect is real**: a temperature-independent r_diode/i_kf
  bends the multi-T curves; the TCs fix most of the −55 °C residual
  (3.40→2.26 %) and cut the overall RMSE 13 % relative.
- **TCs are the best-sensed diode parameters after r_diode** (σ_rel 2–3 %),
  interior to bounds, one basin across three multi-starts. Physical values:
  a_rdiode = +0.34 %/°C (silicon-ish series R), a_ikf = +0.54 %/°C (the
  high-injection knee rises with T).
- **Promoted 2026-10-17 (user go-ahead)**: `ST8L65_DIODE` is now the 6-param
  TC fit — script drift check 2.3e-6 [OK], pinned-test bound tightened to
  < 4.3 % (sits between the 4.15 % record and the 4.78 % TC-off level),
  per-T < 7 %, TCs asserted interior and positive. The TC-off protocol lives
  on as the script's `tc_off_probe` (reproduces 4.78 %) and as the
  session-12 baseline quoted above.

### Trs2 promotion (session 15) — user override, with the caveats on record

User flag: "the improvement isn't as large as expected, especially the
150 °C branch" — probing `b_rdiode` free on top of the 6-param record
(four starts one basin), then **promoted on the user's explicit decision
(2026-10-17)** despite the modest gain:

| | record (6-param, b=0) | + b_rdiode (7-param) |
|---|---|---|
| RMSE_rel | 4.15 % | 3.91 % |
| −55 / 25 / 150 °C rms | 2.26 / 3.41 / 5.90 % | 2.37 / 3.16 / 5.49 % |
| rank / cond (as reported) | 6/6 / 7.12e3 | 6/7* / 1.64e4 |
| b_rdiode | 0 (fixed) | −1.78e-5 /°C² (a_rdiode +0.48 %/°C) |

\* **rank 6/7 is a tolerance artifact**: `b_rdiode`'s Jacobian column carries
the dT² scale (norm ~1.9e4 vs ~2 for i_kf), inflating identifiability's
relative dead-column cutoff; with unit-column scaling the SVD is full rank
7/7 (smallest singular value 3.6e-2 of max).

Caveats that came with the promotion (why the record must be quoted as a
curve, not as parameters):

- **Correlated five-knob block**: r_diode/i_kf/a_ikf/a_rdiode/b_rdiode are
  strongly entangled on three temperatures (|corr| up to 0.985; σ(i_kf)
  blows from 5.9 to 15.4) — the record would quote an untrustworthy split
  (b_rdiode's apparent σ_rel = 0.0003 is pinv along the near-degenerate
  direction, not real information).
- **The 150 °C residual is not series-R curvature**: it is dominated by the
  5 A point (+13…15 %). The digitized v_sd there (0.562 V) sits ~12 mV above
  what any smooth fit wants — in the exponential region, where r_diode·i_di
  is only ~10 mV of the 562 mV drop, so NO Rs TC (Trs1 or Trs2) can move it.
  Excluding that one point, 150 °C rms is 3.85 % (record) / 3.76 % (Trs2
  fit) — the curves agree once the digitization edge point is out. The 5 A
  point at 25 °C is off the other way (−6 mV): the leftmost digitized point
  of each curve is the data's weak spot, not the model.
- The recombination branch was probed too (it is the physical candidate for
  low-current/high-T bending): 3.50 % overall but σ(i_sr) = 645 (unsensed),
  rank 7/9, and the 5 A/150 °C point stays +11.7 % — same verdict: keep off.

Promotion bookkeeping: drift check 3.8e-5 [OK] (four multi-starts, all 3.91 %
one basin); pinned-test bound 0.043 → **0.040** (between the 3.91 % record
and the 4.15 % linear-TC baseline); `trs2_off_probe` (b_rdiode pinned 0)
reproduces the session-14 4.15 % protocol as the standing baseline, next to
`tc_off_probe` (4.78 %, session 12).

## Interface-trap-filling gain `k_beta_trap` (session 19) — LANDED, off by default

The session-18 mechanism is now in the package model:

    v_ov_eff2 = v_ov_eff·(1 + k_beta_trap·v_ov_eff)     k_beta_trap [1/V], 0 = identity

`v_ov_eff2` is THE effective overdrive downstream: `v_knee =
max(v_ov_eff2/(1 + k_theta·v_ov_eff2), v_kneemin)` and the i_d prefactor both
switch to it. Naming per names.xlsx (source of truth): `k_beta_trap` (type
factor, suffix `beta_trap`) and `v_ov_eff2` — the user renamed the session-18
probe's `k_beta` to `k_beta_trap` before implementation.

Decisions pinned (README, "Interface-trap-filling gain"):

- **k_theta multiplies v_ov_eff2, not v_ov_eff** — v_knee is set by the final
  channel charge (same argument the cap section makes).
- **k_beta_trap >= 0 is a hard constraint (user-directed):**
  `BOUNDS["k_beta_trap"] = (0.0, 2.0)` with the lower bound EXACTLY 0.0, and
  `model.py::_require_nonneg_gain` raises ValueError for negative values at
  both `vdmos_id` and `_channel_current`. The clamp alternative (diode.py
  TC-clamp precedent) was rejected — it would silently evaluate the
  non-monotone branch. Rationale: v(1 + k·v) with k < 0 has derivative
  1 + 2·k·v crossing zero → i_d multi-valued in v_gs, bracket proof dead;
  physics (trap filling only raises mobility) makes the constraint free.
- Linear space (0 legitimate, never logged), NOT in DEFAULT_FREE, default 0.0
  = exact identity → every pre-gain anchor bit-identical
  (`tests/test_overdrive_gain.py::test_k_beta_trap_zero_bit_identical`; the
  official ST8L65 record keeps k_beta_trap out of the free set and
  `test_st8l65_official.py` stays green).
- **Bracket-invariant amendment (the docstring claim was imprecise):** the
  package-R proof needs the partials NON-NEGATIVE (monotonicity), not ≤ 1.
  The gain lifts dv_ov_eff2/dv_ov_eff = 1 + 2·k_beta_trap·v_ov_eff to 2.36
  (25 °C) / 3.00 (175 °C) at v_gs = 18 — breaks the ≤ 1 half, keeps
  monotonicity. The r_s > 0 case (untested per HANDOVER 1.4 — silicon fits
  live there) is now tested: dense-grid convergence + brentq parity + the
  SJ-cap set + gain through the gate-degeneration stall band
  (`test_package_r_bracket_with_gain_cap_off`,
  `test_package_r_stall_guard_with_cap_and_gain`).
- **Composition gotcha pinned:** with the cap active the composite asymptote
  is v_ovmax·(1 + k_beta_trap·v_ovmax), NOT v_ovmax
  (`test_gain_raises_cap_ceiling`); synthetic cap+gain round trip gives
  corr(k_beta_trap, v_ovmax) = −0.85 — freeing both on silicon SJ parts
  requires re-quoting the record (never mix records).

**Anchor verification (session 18 HANDOVER 1.7 table):** the probe protocol
re-run with the REAL package model (scratch
`tooling/workspace/outputs/verify_gain_anchors.py`; cap off, k_lambda = 0,
r_s = 0, 8 free incl. k_beta_trap, 4 seeds, r_d capped at datasheet Rds(on)):
RMSE 4.61 / 3.30 % (25/175 °C), k_beta_trap 0.05700 → 0.07584 1/V, k_p 25.74,
v_th0 6.060, v_sub 0.7644, r_d 1.731 → 7.032 mΩ, Ron chords 3.701/4.729 mΩ
(25 °C), Vth spec 57.3/58.2 mA, gain derivative 2.36/3.01 — all match the
standalone-probe anchors to solver tolerance. Note: at 175 °C the DIBL pair
(k_eta, v_dbsat) is unsensed (σ = inf; k_eta pins at the 1e-5 floor, v_dbsat
at the 1300 ceiling — the package BOUNDS are wider than the probe's 0.05/10
caps); the anchors above are unaffected. The high-v_gs transfer overshoot
(model 3650 A vs datasheet figure ~2300 A at v_gs = 18, V_DS = 20 V) is
UNCHANGED by the implementation — fixing it needs the transfer-curve
digitization + joint k_lambda refit (HANDOVER Part 3 queue item 2).

10 new tests in `tests/test_overdrive_gain.py`; suite green (154–155 passed
depending on how many digitizer Tk tests self-skip headless — don't pin
counts). ST8L65 record re-run (`scripts/fit_st8l65_nd.py`): fit A reproduces
the record **bit-for-bit** (drift 1.8e-06 [OK]) under the session interpreter
(scipy 1.18.1). The first re-run (system Python, scipy 1.16.0) landed at a
different point of the flat σ=inf `k_lambda`/`v_dbsat` valley (v_dbsat 1287.03
vs 339.02, k_lambda 1.63e-9 vs 5.13e-10) and the script's naive
max-relative-drift metric flagged STALE 2.80 — A/B-verified as a pure
scipy-version solver-trajectory artifact (`tooling/workspace/outputs/
_s19_drift_ab*.log`): identical code lands at 339.02 under 1.18.1 and at
1287.03 under 1.16.0, every SENSED parameter identical to ~1e-5 relative
either way. Quote the curve, never the valley point.

## CSV header contract (session 20) — queue item 6 CLOSED

`load_csv()` now accepts **both** external I-V header spellings, matched per
column, case-insensitively, order-free; extra columns (the digitizer `index`)
are ignored and the spellings may be mixed in one file:

    legacy SPICE-style   vgs,vds,ids            (what save_csv() writes)
    digitizer/model      index,V_ds,I_d,V_gs    (what curve_digitizer exports)

Decision (user-directed): extend `data.py` rather than rename the CSVs — the
on-disk files are external artifacts and both spellings are legitimate
contracts. `save_csv()` keeps writing the legacy spelling (bit-identical
round trip pinned by `tests/test_data.py::test_save_load_roundtrip_bit_identical`).
The per-script `load_imcq()` workarounds are deleted from
`scripts/fit_imcq120_temperature.py` and `scripts/fit_imcq120_gain_probe.py`
(both now call `load_csv` directly; the gain probe re-run reproduces the
session-18 anchors, scratch `tooling/workspace/outputs/_s20_gainprobe_loadcsv.log`).
New `tests/test_data.py` (8 tests): both spellings load identically, mixed
spelling/case/spacing, `index` ignored, the error message lists both
accepted spellings, real IMCQ CSVs (78/73 pts) and the ST8L65 official CSV
(41 pts) load natively. `data/README.md` documents both contracts and the
IMCQ IV datasets; the ST8L65 header description was corrected to the actual
`index,Vds,Ids,Vgs` (it always loaded — matching is case-insensitive).

## Open problems

1. **k_lambda / v_dbsat unsensed (σ = inf) in the official fit** — *accepted, not a
   defect*: the dataset stops at v_ds = 20 V and never enters deep saturation,
   so there is no output-slope information for CLM or the drain gate to act on.
   k_lambda pins at the LOWER 1e-12 (data wants zero CLM) and v_dbsat is free to be ≫ or
   ≪ the measured range with identical RMSE. Historical note: the old 9-param
   fit pinned k_lambda at its 0.004 cap and that is why the plateau was 5–7 % LOW;
   **superseded by session 10** — with `k_kneesharp` free the knee SHAPE supplies the
   slope and the plateau lands +2.8 %. Keep a bounded bias-dependent element
   (e.g. `r_d·(1 + kd·v_ds)` with the same `dV/dI ≤ bound` discipline) in
   reserve for datasets that *do* reach deep saturation; NOT kappa-as-is
   (pole + divider physics). If high-v_ds data ever arrives, un-pin k_lambda/v_dbsat and
   re-check before trusting any split of the knee cluster.
2. **v_sub one-point-anchored** — live in the official fit (σ_rel 0.95) but
   v_sub = 0.243 V implies ~560 mV/dec swing, not physical (~60–100 mV/dec);
   the single v_gs(th) point pins the *level*, not the *slope*. Fix the swing
   from a digitised near-threshold transfer curve; do NOT hard-fix v_sub
   (pushes the anchor into v_th0, corr +0.99).
3. **Knee shape vs datasheet Figure 3** — session 10 overlay plots
   (`out/st8l65_nd_percurve.png`) show the k_kneesharp fit matching the knee at every
   measured v_gs; a formal digitisation of Figure 3 would quantify it, but the
   visual check is done.
4. **k_theta ↔ v_vsat 84 % correlated** (synthetic grid) — output curves alone
   cannot separate them; and `v_thk` is barely sensed (σ ≈ 477 %, sensitivity
   ~2000× below k_p) because v_kovsoft = 0.1 V saturates the kink sigmoid above
   v_gs = 2 V. Fix from dense kink-region data or measurement.
5. **r_d/r_s split unidentifiable** — quote only r_d+r_s ≈ 32 mΩ. In the official
   fit both are interior (r_d 24.4 / r_s 7.3 mΩ) rather than pinned, but that
   split still leans entirely on r_s's gm-degeneration signature, not on
   independent evidence.
6. The three-gate coupling of v_dbsat (DIBL/CLM/kink all use `tanh(v_ds/v_dbsat)`) is
   unexamined.

## Missing (agreed architecture, not yet built)

- **HV leakage / reverse breakdown (BV/Ibv)** — the diode now sets the
  off-state floor at i_s(T) + i_sr(T) (µA scale, Shockley-consistent), but
  the 650 V blocking spec and its leakage need a dedicated element; the
  removed (0 V, 650 V, 1 µA) IV point waits on it.
- **Third-quadrant channel conduction** — channel is 0 for v_ds ≤ 0; a
  synchronous rectifier needs `sign(v_ds)·i_d(|v_ds|, v_gs)` (Queue item,
  out of scope session 12).
- **Body-diode Qrr / capacitance** (Cjo/Vj/M, Tt) — static diode I–V is in;
  recovery and junction charge are not.
- **Q(V) charge model** — design settled AND LANDED 2026-10-18/21
  (Phases 1–4): see [`CAP_SPEC.md`](CAP_SPEC.md). `caps.py` builds
  charge-master LUTs
  (PCHIP on ln C vs linear v, GL-5 quadrature, clamp policy) from the
  digitized Ciss/Coss/Crss of three parts, with the `c_gd(v_gd)` splice
  (off-state table for v_gd < 0 + exponential-saturation branch for
  v_gd ≥ 0, C0/C1 at the splice, k_agd tied to c_gdmax) and
  `extract_cmax()` on the digitized Qg curves (`scripts/fit_cmax.py`,
  `CAP_*` records in reference.py, 17 tests in test_caps.py).
  **Headline result: all three parts extract at c_gdmax ≈ c_gd0 (flat exp
  branch); the seed c_gs(0) over-predicts 4–8×.** The channel charge the
  off-state tables miss shows up as the Miller-plateau DEFICIT (measured
  minus table prediction): SJ +61 %, trench +22–43 %, SiC +10 % — the
  quantified target for the deferred v_gs-dependent blend (first suspect if
  double-pulse turn-on simulates too fast). Phase 4 (delta subckt export +
  parity) landed 2026-10-21 via `ltspice.py` — see the LTspice-export
  bullet below. The `Tj`
  thermal node for Cauer-ladder co-simulation remains open (the diode `t_c`
  argument is the placeholder interface for Tj).
- **LTspice export**: the capacitance delta is exported and verified —
  `ltspice.py` writes a `.subckt <name> G D S` of three charge-master
  **nonlinear capacitors** (`Cxxx n1 n2 Q=<expr>`, LTspice's one charge
  primitive; `scripts/export_caps.py` emits `out/<part>_caps.subckt` at
  the extracted c_gdmax records). The Q= expression is the **incremental
  ramp-sum** form of the sampled charge table, `q0 + m₁·(x−x₀) +
  Σ Δm·max(0,x−x_k)`: LTspice's symbolic dQ/dx is then exactly the
  piecewise-constant C staircase, and the end-segment slopes carry the
  clamp-C policy for free. Live parity (`test_caps.py::test_ltspice_caps_parity`,
  PWL-ramp drive, three delta configurations): 14 meas points per part
  match the exported staircase to ≤0.3 % (most ≤0.05 %) for all three
  parts; staircase-vs-PCHIP ≤0.4 %. LTspice 26.0.1 gotchas on record:
  `table()` in a B source is rejected outright, and `table()` inside Q=
  compiles but its symbolic derivative is garbage (~1e11 A) — hence the
  ramp-sum form. Still open: the full-device subckt (channel B source +
  diode + r_d/r_s wiring) for A/B simulation — the diode part is already
  LTspice-exact: `D(Is=… N=… Isr=… Nr=… Ikf=… Rs=…
  Tnom=… Xti=… Eg=…)`.
- **CI**. (`notebooks/` is empty scaffolding.)
- Validation target remains double-pulse Eon/Eoff within ~10 %, which static
  fitting alone cannot deliver.
