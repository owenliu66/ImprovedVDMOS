# Project Status — ImprovedVDMOS

*Snapshot: 2026-10-17 (session 12 — body diode: LTspice diode DC branch,
multi-temperature extraction, ST8L65 diode record). Companion to
[README.md](../../README.md).*

## Current state

Twelve sessions, all committed to `main` (session 11 renamed every parameter
to the `<type prefix>_<suffix>` convention — names.xlsx + README "Naming
convention" are the mapping, tag `pre-rename-2026-10` is the escape hatch, and
the rename is bit-identical; session 12 added the **body diode** —
`diode.py`, LTspice diode DC branch verified against LTspice 26.0.1,
terminal-wired bypassing r_d/r_s, multi-temperature extraction, official
record `ST8L65_DIODE`).
The package provides one C¹ continuous
drain-current model — softplus subthreshold, gate-overdrive charge cap (v_ovmax, k_ovsharp),
mobility degradation via v_knee, velocity saturation, bounded DIBL, gated CLM,
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
(−55/25/150 °C), rank 4/4, cond 6.1e3, v_eg = 1.097 eV**. 120 tests green,
ruff clean. Datasheet Figure 3 (typical output curves) and the body-diode
forward characteristics CSV are on disk in `data/`.

## What exists

```
Projects/ImprovedVDMOS/
├── pyproject.toml          deps, pytest, ruff config, `vdmos-fit` console script
├── README.md               model + corrected Level-1 analysis
├── .gitignore
├── src/vdmos_fit/          model, diode, params, data, fit, plotting, reference
│                           (official device records), cli, __main__
├── tests/                  test_model, test_fit, test_package_r, test_gate_cap,
│                           test_drain_softness, test_st8l65_official, test_diode,
│                           test_cli, test_digitizer
├── scripts/                fit_datasheet.py, fit_st8l65.py (superseded),
│                           fit_st8l65_nd.py (official channel recipe),
│                           fit_st8l65_diode.py (official diode recipe),
│                           roundtrip_demo.py
├── data/                   README.md + ST8L65N044M9_IV.csv (41 digitised points)
│                           + ST8L65N044M9_diode.csv (30 pts, −55/25/150 °C)
│                           + ST8L65N044M9_Id_Vds_Vgs.png (datasheet Figure 3)
└── notebooks/              README.md only (empty)
```

## What it currently does

- **Model** — `vdmos_id()` evaluates the unified continuous drain current
  (softplus subthreshold, gate-overdrive charge cap v_ovmax/k_ovsharp, mobility
  degradation via v_knee, velocity saturation, DIBL, gated CLM, one-sided kink
  sigmoid), embedded with package series resistance r_d/r_s, minus the body
  diode `vdmos_diode()` (LTspice diode DC branch: diffusion + high-injection
  compression + recombination + r_diode + Xti/Eg temperature law), wired at
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

- **120 tests, all passing**; `ruff check` clean. Caveat: 8 of the digitizer
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

## Current ST8L65 diode fit — OFFICIAL: fit A (4 free params, session 12)

`scripts/fit_st8l65_diode.py` fit A; outputs `out/st8l65_diode_params.txt`,
`out/st8l65_diode_probe.txt`, `out/st8l65_diode_fit.png`. Record in
`vdmos_fit.reference.ST8L65_DIODE`, drift-checked by the script
(`OFFICIAL_TOL`) and pinned by `tests/test_diode.py::test_st8l65_diode_record_reproduces`
(RMSE < 5 %, per-temperature < 8 %, r_diode interior, reverse floor < 10 µA,
v_eg in (0.9, 1.3) eV).

Data: `data/ST8L65N044M9_diode.csv` — 30 points, I_sd 5–55 A at −55/25/150 °C.
Free = i_s, r_diode, v_eg, i_kf; **k_n fixed 2.0** (physical diffusion limit),
k_xti fixed 3.0 (LTspice default), i_sr = 0, t_nom = 25 °C; device bound
r_diode ≤ 20 mΩ. Four multi-starts, all converge to the same basin.

| | fit A (OFFICIAL, 30 pts) |
|---|---|
| RMSE_rel | **4.78 %** (−55 °C 3.40, 25 °C 3.85, 150 °C 6.49) |
| i_s | 3.701e-6 A (σ_rel 0.88 — weak prefactor, pair with k_n) |
| r_diode | 1.586 mΩ (σ_rel 0.004) |
| v_eg | 1.0967 eV (σ_rel 0.093 — silicon bandgap, sensed) |
| i_kf | 17.45 A (σ_rel 5.2; r_diode↔i_kf corr +0.948) |
| rank / cond | **4/4** / 6.13e3 |
| reverse floor (−1 V, 25 °C) | 3.7 µA |

Findings / protocol decisions:

- **The temperature law is mandatory**: a 25 °C-only fit predicts −55 °C
  out-of-sample at ~260 % error; V_T-only scaling (no Xti/Eg) predicts the
  forward drop going UP with T — the data goes DOWN. v_eg ≈ 1.10 eV is
  genuinely sensed (σ_rel 0.09) — the multi-temperature data pins it.
- **k_n fixed at 2.0 from physics**: the k_n sweep (i_kf off) shows RMSE
  improving monotonically toward k_n ≈ 2.5–3 — ideality absorbing the
  high-injection bending, with the reverse floor scaling with i_s (3 pA at
  k_n = 1 → 0.4 mA at k_n = 3). Freeing k_n on top of fit A buys 0.12 %
  (4.66 %), lands at k_n = 2.07 with corr(i_s, k_n) = +1.000 — the
  ideality/prefactor split is a degenerate valley; only the PAIR is pinned.
- **i_kf belongs in the primary set**: pinned off (sentinel 0.0 — NOT the
  1e-3 lower bound, which is an active knee that crushes the current) the
  fit is 6.48 %; free, 4.78 %. The high-injection knee carries real
  curvature in the 5–55 A range.
- **i_sr/k_nr unsensed**: freeing them gives identical RMSE with
  σ(i_sr) = inf — recombination lives below the 5 A data floor. Keep off.
- **i_s is a weak prefactor** (σ_rel ≈ 0.9): quote the curve or the pair
  (k_n, i_s), never i_s alone — same doctrine as the channel knee cluster.
- **Log-space packing is mandatory** for the diode fit: linear-space i_s
  gave rank-deficient Jacobians and a runaway basin at r_diode ≈ 0.34 Ω;
  log-space gives rank 4/4 and r_diode stable across starts. The per-device
  r_diode ≤ 20 mΩ bound (measured dV/dI is 3.4–11.7 mΩ) excludes the
  runaway even so.

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
- **Q(V) charge model** (charge-based capacitors via symbolic differentiation)
  and the `Tj` thermal node for Cauer-ladder co-simulation (the diode `t_c`
  argument is the placeholder interface for Tj).
- **LTspice export**: fitted parameters → `.subckt`/`.model` for A/B (the
  diode part is already LTspice-exact: `D(Is=… N=… Isr=… Nr=… Ikf=… Rs=…
  Tnom=… Xti=… Eg=…)`).
- **CI**. (`notebooks/` is empty scaffolding.)
- Validation target remains double-pulse Eon/Eoff within ~10 %, which static
  fitting alone cannot deliver.
