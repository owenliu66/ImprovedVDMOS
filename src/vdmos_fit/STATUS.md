# Project Status — ImprovedVDMOS

*Snapshot: 2026-11-07 (sessions 27–28 — ISC040N10NM8 FULL DEVICE LANDED:
official channel record `ISC040_CHANNEL` (SJ charge-cap family, 25 C RMSE
1.33 %, TC 2.97 %), official body-diode DC record `ISC040_DIODE` (0.51 % —
**k_n ≈ 1.05, the k_n = 2.0 ST8L65 value is WRONG here**, linear `a_rdiode`
is the decisive TC), full 4-terminal subckt generator
`scripts/gen_full_subckt.py` → `out/ISC040N10NM8_full.subckt` (D G S K;
D1 = DC-fit set with Rs=0 + `R1 = {Rs*(1+Trs1*(temp-Tnom))}` — the diode
series R counted EXACTLY ONCE; LTspice parity 0.05 ppm channel / ≤27 ppm
diode at 25 AND 175 °C; recovery insensitive to the D1 swap), pinned
`tests/test_isc040_official.py` (11 tests); caps parity pin
**RE-BASELINED on LTspice 26.1.1** (1 % — stair-breakpoint knife-edge,
worst measured 0.56 %); HANDOVER.md retired, its findings folded into this
STATUS. Before that, session 26 (2026-11-02) completed the
`extract_recov.py` rewrite end-to-end (template verified against the bench
before use; brace bugs fixed), and sessions 25a–25c — RECOVERY REWORKED, LTspice-side:
intrinsic (C_oss-free) `q_rri`/`t_rri` measurement doctrine adopted (datasheet
Qrr excludes Q_oss); the `DiodeRecov` model gained the E1+R3+C3 TRANSIT-TIME
POLE — q_0 deleted as a pure gauge, knobs are now t_minor / t_pole / tau;
all four ISC040 datasheet recovery numbers fit to **RMS 1.11 %**
(tau 3.79e-4, t_minor 22.363 n, t_pole 0.8324 n) — the session-23 "Qrr rises
with di/dt is unreproducible" verdict is DEAD; the datasheet-driven extraction
protocol LANDED as `src/vdmos_fit/extract_recov.py` (pole gauge:
index-decomposed single `.step`, log-log surrogates, multi-start TRF polish
with `diff_step=0.03`, no dTmax, tau free with interior optimum ~3e-4);
durable records `out/RECOV_INTRINSIC_FIT.md` + `out/RECOV_POLE_FIT.md`.
Before that, session 23 (2026-10-30) — BODY-DIODE RECOVERY INTEGRATED +
FITTED (ISC040N10NM8, LTspice-side only, NO package code touched): the user's
behavioral `DiodeRecov` network (charge-storage node q_j, charge-controlled
B2 reverse-recovery source, B3 forward-recovery resistor) is now part of the
full-device subckt template `out/ISC040N10NM8_recov.subckt` — the measured
C_ds charge-master relocated INSIDE the recovery branch at Ki–A (D–S
placement crashes hard forward re-turn-on; Ki–Ai is singular), B3 +1 µΩ
floor and R2 = 1 µΩ break near-singular bridge loops (<0.1 % waveform
impact), `method=gear reltol=1e-4` required; two-knob fit **q_0 = 35 n,
tau = 5** hits the datasheet 1000 A/µs recovery (Qrr 249.6 vs typ 247 nC,
t_rr 27.6 vs 25 ns — two-knob fit SUPERSEDED by the pole gauge above);
report `out/RECOV_REPORT.md`. Before that, session 22
(2026-10-29): KELVIN-AWARE r_s SPLIT LANDED + queue
item 2 CLOSED BY EXPERIMENT: `Params.r_sk` (names.xlsx "sk") gives the
Kelvin-aware package-R split `v_gs' = v_gs − (r_s − r_sk)·i_d` (gate loop
returns at the kelvin tap; full `r_s` stays in the current path; guard
`0 ≤ r_sk ≤ r_s`; `r_sk = 0` bit-identical legacy), pinned by
`tests/test_kelvin_split.py`; the joint-vector refit
(`scripts/fit_imcq120_kelvin.py` + `fit_imcq120_kelvin_probe.py`, report
`out/IMCQ120_KELVIN_SPLIT_REPORT.md`) shows the two open IMCQ gaps are NOT
series-R topology problems: gate-loop degeneration is REJECTED by the data
(`r_s → 0` boundary optimum from every start; forced `r_g` worsens RMSE
1.93 → 2.71 % AND g_fs), the current-loop split is DC-degenerate (only
`r_d + r_s` sensed), and terminal g_fs is pinned by the transfer curve
itself (~130 S). Record UNCHANGED. Before that: session 21 closed HANDOVER
PART 1 and PROMOTED the
IMCQ120R004M2H transfer characteristic (V_DS = 20 V, 25/175/200 °C) was
digitized and a JOINT multi-temperature fit (both IV sets + transfer +
Rds(on)-Tj + Vth-Tj, TC laws as the only extra unknowns) landed as
`scripts/fit_imcq120_joint.py` and is now the **OFFICIAL IMCQ record**
(`reference.py::IMCQ120_JOINT`, pinned by `tests/test_imcq120_official.py`);
the session-18 `out/IMCQ120_TEMPERATURE_REPORT.md` was **deleted** (durable
findings carried into `out/IMCQ120_JOINT_FIT_REPORT.md` §8 — the joint-fit
report is now the IMCQ document). Headline results in that report and the
session-21 section below; session 20 closed the CSV header contract —
`load_csv()` reads both header spellings natively, see the session-20 section
below; session 19 landed the interface-trap-filling gain `v_ov_eff2`/
`k_beta_trap` in the channel model, off by default — committed as `8574ce2`;
session 18 investigated the IMCQ120R004M2H two-temperature IV fits and the
gain probe — committed as `3385ee2` (its report was deleted session 21;
durable findings in the joint report §8); session 17 landed the
capacitance/charge model — see CAP_SPEC.md).
Companion to [README.md](../../README.md).*

## Current state

Twenty-eight sessions, committed to `main` (sessions 27–28 landed the
ISC040N10NM8 FULL DEVICE — official channel record `ISC040_CHANNEL`, official
diode DC record `ISC040_DIODE`, the 4-terminal subckt generator
`scripts/gen_full_subckt.py` and pinned `tests/test_isc040_official.py`; see
the sessions-27–28 section below; session 26 completed the
`extract_recov.py` rewrite end-to-end; sessions 25a–25c reworked the
body-diode recovery into the TRANSIT-TIME POLE gauge — intrinsic doctrine,
q_0 deleted, RMS 1.12 % on all four ISC040 datasheet numbers — and landed
the extraction protocol `src/vdmos_fit/extract_recov.py`, the record
`ISC040_RECOV`, and the regenerated pole-gauge template
`out/ISC040N10NM8_recov.subckt` — all LTspice-side; the session-23 section
below keeps the integration decisions (its fit is superseded); session 22 landed the Kelvin-aware
`r_sk` package-R split and closed HANDOVER queue item 2 by experiment — the
IMCQ joint vector refit rejects gate-loop degeneration and proves the
current-loop split DC-degenerate (session-22 section below); session 21
closed HANDOVER PART 1:
the IMCQ120R004M2H joint multi-temperature fit — transfer-curve digitization +
one residual vector over IV/transfer/Rds(on)-Tj/Vth-Tj with TC laws as the only
extra unknowns, quadratic `r_d(T)` and freed `k_theta` per user direction,
then **promoted to the official IMCQ record** `IMCQ120_JOINT` with the
session-18 temperature report deleted (uncommitted at snapshot time:
`scripts/fit_imcq120_joint.py`, `reference.py` record,
`tests/test_imcq120_official.py`, data-README entry, README, this section); session 18 added the IMCQ120R004M2H
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
channel v_ds ≥ 0 only; reverse polarity is carried by the body diode) with the
Kelvin-aware gate-loop split r_g = r_s − r_sk (session 22,
`tests/test_kelvin_split.py`), plus bounded least-squares extraction with identifiability
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
41-point set, rank 8/10, cond 2.19e3)**; diode **3.91 % on 30 points
(−55/25/150 °C), v_eg = 1.082 eV ≈ silicon, all three TCs interior**
(session-15 record; the 4.15 % linear-TC and 4.78 % TC-off fits live on as
script baselines/probes). Official record on IMCQ120R004M2H (SiC):
**`IMCQ120_JOINT` — joint multi-temperature fit, 1.93 % relative RMSE on 200
points / 29 temperatures** (IV 25/175 °C + transfer 25/175/200 °C +
Rds(on)-Tj + Vth-Tj; TC laws as the only extra unknowns; quadratic `r_d(T)`;
see the session-21 section below and `out/IMCQ120_JOINT_FIT_REPORT.md`).
Official records on ISC040N10NM8 (trench SJ, non-Kelvin; sessions 27–28):
**`ISC040_CHANNEL`** — 25 C fit 1.33 % (132 pts) + TC fit 2.97 % (Rds(on)
chords at −55/25/100/175 °C + Vth-vs-Tj voltage residual), pinned by
`tests/test_isc040_official.py`; **`ISC040_DIODE`** — 0.51 % on 26 pts at
25/175 °C; **`ISC040_RECOV`** — pole-gauge recovery (sessions 25a–25c).
The full-device deliverable is `out/ISC040N10NM8_full.subckt`.
Suite green, ruff clean — test totals drift with headless Tk skips, don't pin
counts (250 tests at the session-30 snapshot; the caps-parity pin was
re-baselined on LTspice 26.1.1, see the sessions-27–28 section). Datasheet
Figure 3 (typical output curves), the body-diode forward
characteristics CSV, the IMCQ120R004M2H
IV + Rds(on)/V_th-vs-T_j CSVs, and the three-part C–V + Qg CSVs are on disk
in `data/`.

## What exists

```
Projects/ImprovedVDMOS/
├── pyproject.toml          deps, pytest, ruff config, `vdmos-fit` console script
├── README.md               model + corrected Level-1 analysis
├── .gitignore
├── src/vdmos_fit/          model, diode, caps, params, data, fit, plotting,
│                           reference (official device records), ltspice export,
│                           extract_recov (recovery extraction protocol),
│                           cli, __main__
├── tests/                  test_model, test_fit, test_package_r, test_gate_cap,
│                           test_drain_softness, test_overdrive_gain, test_data
│                           (CSV header contract), test_st8l65_official,
│                           test_imcq120_official (joint multi-T record),
│                           test_isc040_official (channel + diode records +
│                           generated-subckt structure pin),
│                           test_kelvin_split, test_recov_record, test_diode,
│                           test_caps (incl. live LTspice parity), test_cli,
│                           test_digitizer
├── scripts/                fit_datasheet.py, fit_st8l65.py (legacy baseline),
│                           fit_st8l65_nd.py (official channel recipe),
│                           fit_st8l65_diode.py (official diode recipe),
│                           fit_imcq120_temperature.py + fit_imcq120_gain_probe.py
│                           (SiC two-temperature protocol + k_beta_trap probe),
│                           fit_imcq120_joint.py (session-21 joint multi-T fit,
│                           TC laws + quadratic r_d(T) — PART 1 recipe),
│                           fit_isc040_channel.py (25 C-locked then TC —
│                           sessions 27), fit_isc040_diode.py (DC fit, session 28),
│                           gen_full_subckt.py (4-term device subckt generator),
│                           fit_cmax.py (official c_gdmax extraction),
│                           export_caps.py (LTspice subckt export),
│                           roundtrip_demo.py
├── data/                   README.md + ST8L65N044M9_IV.csv (41 digitised points)
│                           + ST8L65N044M9_diode.csv (30 pts, −55/25/150 °C)
│                           + ST8L65N044M9_Id_Vds_Vgs.png (datasheet Figure 3)
│                           + IMCQ120R004M2H IV ×2 + Transfer + Rdson_Tj +
│                             Vth_Tj CSVs
│                           + ISC040N10NM8 IV_25C + Rdson_Id_25C + Rdson_Vgs +
│                             Vth_Tj + Diode CSVs (sessions 27–28)
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

- **Suite green, `ruff check` clean** (161 passed / 2 skipped at the
  session-20 snapshot). Caveat: 8 of the digitizer
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
- **ISC040 full-device subckt verified against LTspice 26.1.1** (sessions
  27–28): channel B-source parity 0.05 ppm vs the netlist-constant model,
  forward body-diode ≤27 ppm at 25 AND 175 °C, recovery insensitive to the
  D1 DC-set swap — full numbers in the sessions-27–28 section below.

## Current ST8L65 fit — OFFICIAL: 10-param with `k_kneesharp`, fit A (full data)

**Promoted 2026-10-13.** `scripts/fit_st8l65_nd.py` fit A; outputs
`out/st8l65_nd_params.txt`, `out/st8l65_nd_probe.txt`, `out/st8l65_nd_fit.png`,
`out/st8l65_nd_percurve.png`, `out/st8l65_nd_overlay.png`. The parameter record
lives in that script's `OFFICIAL` dict, which now re-exports
`vdmos_fit.reference.ST8L65_OFFICIAL` (single source of truth), and `main()`
checks each run reproduces it to 1e-3 relative (`OFFICIAL_TOL`) — if a future
refit drifts, the run says `[STALE]` and the record must be re-promoted
deliberately. Caveat (session 19): a `[STALE]` flag can also be an
environment artifact — the fit lands at different points of the flat
σ = inf `k_lambda`/`v_dbsat` valley across scipy versions while every
SENSED parameter stays put; A/B the interpreter before re-promoting (see
the session-19 section). `tests/test_st8l65_official.py` pins the same record against the
measured CSV (RMSE, Ron chord, plateau, v_gs(th), per-curve, k_kneesharp interior, no
negative gds) so a model/bounds change that moves the official fit fails the
suite instead of silently invalidating these numbers.

Free = k_p, v_th0, v_sub, k_lambda, r_d, r_s, v_dbsat, v_ovmax, k_ovsharp, **k_kneesharp**; fixed k_theta = 0,
v_vsat = 1e4, k_eta = 0, g_k = 0, v_thk = 1.6, v_kovsoft = 0.1 (all measured dead on this
data); device bounds r_d, r_s ≤ 36 mΩ (datasheet Rds(on)). Global bounds
widened at promotion (session 10): k_p ≤ 1e4, v_vsat ≤ 1e4, k_lambda 1e-12–1.0, v_dbsat ≤ 1300,
k_kneesharp ≥ 0.1. Four multi-starts (k_kneesharp starts 2.8 = legacy tanh,
1.5, 6), all four converge to the same basin.

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
the script and in the table above for comparison.

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

**Promotion history (session 10):** the retired 9-param `k_kneesharp = 0`
recipe (`scripts/fit_st8l65.py`, kept as the legacy baseline probe) gave
9.97 % (fit A) / 5.54 % (fit B); freeing `k_kneesharp` roughly halved it to
the current 3.79 / 1.52 % — the missing physics was the knee SHAPE, not a
missing circuit element. The v_gs = 5 curve is genuine datasheet data but
unreproducible by any C¹ model (identical i_d at v_ds = 1.127 and 5.153 V
requires β = 0 over a 4.6× interval, then β ≈ 0.12); it costs ~21 % on
fit A — that is what fit B isolates.

Findings (promotion record; probe details in `out/st8l65_nd_probe.txt`):

- **The data wants a knee SOFTER than the legacy tanh, and softer than the old
  0.5 guard allowed**: k_kneesharp lands at 0.436 (fit A) / 0.408 (fit B), interior to
  the relaxed (0.1, 40) bound (relaxed from 0.5 at promotion because the optimum
  sits below the old guard). The k_kneesharp lower-bound probe reaches the same
  1.52 % from starts at k_kneesharp = 0.5 and 0.3 — a stable interior optimum, so the
  0.1 floor is no longer load-bearing. k_kneesharp → 0 collapses the current, so lower
  is not physical anyway; do not lower the floor further.
- **Open problem 1 was a knee-SHAPE problem, not a missing element**: with k_kneesharp
  free the plateau lands +2.8 % instead of −6.8 % low. The bounded
  bias-dependent drain element (Queue 2) is NOT needed on this data. k_lambda now
  pins at the LOWER 1e-12 — the data wants *zero* CLM once the knee shape is
  free, which is the expected answer for a 650 V part measured only to 20 V.
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

## Joint multi-temperature fit — IMCQ120R004M2H (session 21) — HANDOVER PART 1 CLOSED

One residual vector, 200 points / 29 temperatures: IV 25 °C (78) + IV 175 °C (73)
+ digitized transfer @ V_DS = 20 V, 25/175/200 °C (23, `data/IMCQ120R004M2H_Transfer.csv`)
+ Rds(on)-Tj chords @185.2 A, Vgs 15/18 (22) + Vth-Tj @58.3 mA (4). TC LAWS are the
only extra unknowns (no more independent per-T fits). Script
`scripts/fit_imcq120_joint.py`, report `out/IMCQ120_JOINT_FIT_REPORT.md`,
log `out/_s21_joint_run.log`.

User-directed protocol changes vs the HANDOVER spec:

- **`r_d(T)` QUADRATIC in T(°C)** (Trs1/Trs2 law form), replacing the power law /
  exponential of the temperature report. The digitized chords fit the quadratic to
  **0.067 % / 0.065 % RMSE** over −46…196 °C (power law 9.5/9.1 %, exponential
  5.9/6.3 %, linear 13.4/10.9 %) — the quadratic curvature IS the cold-end
  flattening. Joint-fit Rds(on) check: **RMSE 0.66 %, max 1.08 %, cold end 0.70 %**
  (was −15.3 % at −46 °C). The `rd_powerlaw` ablation degrades the series path 7×
  (4.59 %, cold 6.98 %) — the user call is confirmed in the joint fit too.
- **`k_theta` FREED** — lands at 0.109 1/V (v_knee → 9.2 V) and is the mechanism
  that bends the high-v_gs transfer over. Ablation `no_ktheta`: overall RMSE
  1.93 → 4.14 % (25 °C v_gs=8 curve 2.7 → 9.3 %).
- **`k_lambda` re-freed — and the data says NO**: it lands at 4.4e−10, σ = inf;
  the `no_klambda` ablation is identical. The HANDOVER hypothesis "bending the
  transfer is what k_lambda does" is disconfirmed — knee saturation (k_theta) +
  the package-r_d divider + the transfer data itself do it.

Headline (primary basin, 7 starts, costs 3.7105–3.7327e−2): **overall RMSE
1.93 %** (IV 1.73 / transfer 3.39 / Rds(on) 0.65 / Vth 2.06 %) — the session-18
8.53/6.13 % floor is gone. High-v_gs corner v_gs=18 @20 V: 175 °C **2118 A**
(digitized ~2070; the gain fit predicted 3650) — the "don't quote above 1500 A"
gate is LIFTED. ZTC 25/175 at 14.92 V/1703 A (matches the digitized output-curve
crossing ~14.9 V; the datasheet-figure ~14.2 V read stays +5 % off). Ron chords
3.722/4.652 mΩ @25 °C (data 3.70/4.70). Vth spec 57.0/59.1 mA.

TC laws (25 °C base): k_p25 9.318 A/V² with beta_p 0.378 (k_p ~ T^−β);
v_th0 5.746 V with slope **−7.16 mV/K** (σ 0.5 %; datasheet −6.67; fitted offset
≈ +1.5 V constant — never substitute the datasheet V_th); v_sub 0.7933 invariant;
k_beta_trap 0.4687 with +3.40e−4 1/V/K (linear law; the Arrhenius variant ties at
1.92 vs 1.93 % with Ea σ = inf → linear kept, "the joint fit decides"); r_d 2.3985 mΩ
with a 6.727e−3 /K, b 5.367e−5 /K² (vertex ≈ −38 °C — law valid over the data span,
not below); k_kneesharp 1.465; k_eta pinned at the 0.05 cap (σ inf), v_dbsat 320
(σ inf), k_lambda ≈ 0 — the DIBL/CLM block is inert (bounds, not measurements).

Identifiability doctrine CONFIRMED (HANDOVER 1.3): the `no_transfer` ablation
(IV-only, same TC laws) fits IV *better* (1.59 %) but its r_d(T) TC is degenerate
(a_rd σ_rel 295, corr(a_rd,b_rd) −1.00) and its OUT-OF-SAMPLE Rds(on)-Tj fails at
20.2 % RMSE / 34.1 % cold — the chord data is what pins r_d(T); only then is
beta_p interpretable. Composition gotchas on record: k_p25~k_bt25 −0.95 (the pair
is a model-form split — the composite k_p·v_ov_eff2 is the record; quote curves,
never the pair), beta_p~s_bkt +0.97 (the channel TC is shared). The gain
derivative dv_ov_eff2/dv_ov_eff reaches **12.5 at v_gs=18** — legal for the
package-R bracket only because r_s = 0; re-check before any r_s > 0 use of this
record (Kelvin-aware subckt, queue 4).

Open gaps from this fit (in the report): g_fs 129.5 S vs 122 spec (+6.2 %), the
25 °C transfer corner −5…−7 % at v_gs 15–19 (r_d25 2.40 mΩ drops 7.2 V at 3 kA —
the Kelvin source-side power path is the expected fix), ZTC 175/200 crossing
12.97 V vs digitized 12.03 V (hot-curve separation).

### Promotion to OFFICIAL (session 21, user go-ahead) — `IMCQ120_JOINT`

The joint fit is now the quoted IMCQ120R004M2H channel model, single source of
truth in `reference.py`:

- **Record:** `IMCQ120_JOINT` — an `ImcqJointRecord` (25 °C base + TC-law
  coefficients, `params_at(t_c)` builds the `Params`). First record that is not
  a bare `Params`: the multi-T laws ARE the model, so the record carries them
  (quadratic `r_d(T)`, linear-clamped `k_beta_trap(T)`, power-law `k_p(T)`,
  linear `v_th0(T)`). `IMCQ120_JOINT_RMSE_REL = 0.0193` + per-group constants.
- **Pinned test:** `tests/test_imcq120_official.py` (11 tests, ~1 s) evaluates
  the RECORDED set against all four measured CSVs (no re-fit): overall +
  per-group RMSE reproduction, full-span Rds(on) chords (max < 1.5 %, cold rmse
  < 1.2 % — pins the −15.3 %-cold-end fix), Ron chords, the high-v_gs corner
  (< 2400 A ceiling — regression guard against the old 3650 A gain-fit
  overshoot), ZTC crossings incl. the digitized 175/200 crossing point (±2 %),
  V_th spec, TC-law physical invariants (k_beta_trap ≥ 0, r_d > 0 over the
  span, v_th0 slope, v_sub invariant), and the two KNOWN OPEN GAPS pinned as
  tolerances (g_fs +6.2 %, 25 °C corner −7 %) so they cannot silently worsen
  while the Kelvin-aware subckt is built.
- **Drift check:** `scripts/fit_imcq120_joint.py` re-checks the record every
  run — sensed-parameter drift only (session-19 lesson: unsensed directions
  are solver-trajectory noise, scipy-version sensitive) + record-RMSE
  re-evaluation through the script's own machinery (cross-checks the
  `reference.py` `params_at` law implementation against the script's).
- **Bound-pin exceptions documented in the record:** `k_eta` at the 0.05 cap,
  `k_lambda` ≈ 0, `v_dbsat` σ = inf — the inert DIBL/CLM block; everything
  that matters is interior (`r_d25 < 3.7 mΩ` chord guard asserted by the test).
- **Legacy demoted:** `scripts/fit_imcq120_temperature.py` (per-T protocol) and
  `scripts/fit_imcq120_gain_probe.py` (mechanism probe) stay as historical
  baselines with supersession notes in their docstrings.
- **Temperature report deleted (user direction):**
  `out/IMCQ120_TEMPERATURE_REPORT.md` removed; its durable SiC findings
  (cap inertness evidence, datasheet frame-truncation corrections, linear-r_d
  −26 °C argument, v_th0 offset history, v_sub invariance, k_p-TC
  model-form-dependence, the 4 suspect digitized points, digitized ZTC
  crossing, Kelvin caveat) are carried in `out/IMCQ120_JOINT_FIT_REPORT.md`
  §8 — that report is now the IMCQ document.

## Kelvin-aware r_s split — IMCQ120R004M2H (session 22) — queue item 2 CLOSED

**Model support landed + experiment done; the record is UNCHANGED.**
`Params.r_sk` (names.xlsx "sk", linear space, bounds (0, 0.5), not in
DEFAULT_FREE) implements the asymmetric split the queue asked for:

```
v_ds' = v_ds − (r_d + r_s)·i_d      current loop: ALL of i_d
v_gs' = v_gs − (r_s − r_sk)·i_d     gate loop: returns at the KELVIN tap
```

`r_sk` = resistance between the kelvin tap and the power source pin (the
downstream segment of `r_s`, carrying only gate current — zero in the DC
model). `r_sk = 0`: no kelvin contact, legacy semantics (bit-identical,
pinned); `r_sk = r_s`: ideal tap at the die, gate loop clean. Guard
`0 ≤ r_sk ≤ r_s` raises ValueError (negative degeneration is not a passive
network and kills the bracket proof). The solver takes the gate-loop
coefficient `r_g` explicitly; the HANDOVER "re-check before r_s > 0" note is
DONE: the bracket survives with `r_g > 0` AND the trap-filling gain active
(gate derivative ~12.5) — monotonicity is all the proof needs. Pinned by
`tests/test_kelvin_split.py` (9 tests: identity, ideal-kelvin equivalence,
brentq cross-checks, r_sk-monotone current lift, guards, fitting space).

**The experiment** (`scripts/fit_imcq120_kelvin.py`, same 200-pt joint
vector/TC laws/acceptance as the official fit; starts from the record
optimum; log `out/_s22_kelvin_run.log`, report
`out/IMCQ120_KELVIN_SPLIT_REPORT.md`):

| variant | landed | RMSE % | g_fs S | corner25 A |
|---|---|---|---|---|
| control (record) | — | 1.93 (drift 4e-6 OK) | 129.5 | 2990 |
| kelvin_ideal (r_sk = r_s, gate clean) | r_s25 σ_rel 75, corr(r_d25,r_s25) −1.00, r_tot25 invariant | 1.93 (identical) | 129.5 | 2990 |
| degenerated (r_sk = 0, degeneration = r_s) | **r_s25 → 0 pinned LOWER** (3 starts) | 1.93 | 129.5 | 2990 |
| split_free (r_s25, f_sk = r_sk/r_s) | r_s25 → 0 pinned; f_sk σ inf | 1.93 | 129.5 | 2990 |

Forced-degeneration probe (`fit_imcq120_kelvin_probe.py`): r_g = 0.25/0.5/1.0
mΩ monotonically WORSENS RMSE (1.93 → 2.11 → 2.28 → 2.71 %), deepens the
corner deficit (−6.2 → −10 %) and does NOT lower g_fs (the channel re-tunes
to re-pass through the same transfer curve).

**Verdicts (on record):**
- The current-loop split is DC-invisible (only `r_d + r_s` and the total TC
  are sensed) — the doctrine's −0.999/±1.00 correlation result, proven.
- Gate-loop degeneration ("the datasheet static measurement did not use the
  kelvin pin") is REJECTED by the data: boundary optimum r_g = 0 from every
  start; the corner needs MORE channel drive, degeneration removes it.
- g_fs +6.2 % is a spec-vs-digitized-slope softness: the transfer curve
  itself pins the terminal slope at 185 A to ~130 S; no series-R topology
  can move it without breaking the curve fit.
- The 25 °C corner −6 % is a channel-shape compromise (k_p25~k_bt25~k_theta
  squeezed by the hot curves + chords), NOT a series-path gap.
- `r_sk` stays in the model for the 4-terminal SUBCKT (queue 3): the kelvin
  pin carries only gate CURRENT — it matters during switching (Miller /
  displacement current through r_sk), which is exactly the full-device
  subckt. DC record keeps r_s = r_sk = 0.

## Body-diode recovery network — ISC040N10NM8 (session 23) — integrated, LTspice-side only
**(the q_0/tau fit and the "unreproducible Qrr-vs-di/dt" verdict below are SUPERSEDED by
sessions 25a–25c — see the next section; kept for the integration decisions, which stand)**

The user's behavioral recovery model (`DiodeRecov` subckt; benches in the
Strata workspace `outputs/vdmos_recov/`) is integrated into the full-device
subckt template `out/ISC040N10NM8_recov.subckt` (report:
`out/RECOV_REPORT.md`). **No package Python touched** — recovery lives in
the LTspice export layer.

Network: D1 (diode DC branch, `Rs` externalized to R1 = 5 mΩ) in parallel
with B2 (charge-controlled reverse current `V(q_j)·(V(Ki)−V(Ai))/tau`) and
the forward-recovery path R1 + B3 (`V = −I(R2)·(r_max/(1+k_chg·V(q_j)) + 1µ)`);
the storage node q_j = C1 (`q_0/tau` farads to ground) fed by B1 (`−I(R2)/tau`)
and clamped ≥ 0 by ClampD (`Roff=tau`).

**Decoded semantics (simulation-verified — the model is more physical than it
looks):**
- ClampD `Roff=tau` closes the storage loop: forward steady state
  **V(q_j) = I_F volts, tau-independent** (verified 12.5 V at 12.5 A); q_j
  time constant = Roff·C1 = **q_0** (35 ns at q_0 = 35 n — the 1 Ω
  normalization makes the user's "stored charge at 1 A" and "carrier
  lifetime" readings the same number); stored charge = C1·I_F = q_0·I_F/tau
  ∝ I_F — quasi-steady, not an unbounded integrator.
- The −I_RM plateau is carried by **B2, not D1** (I(D1) ~ −1 nA during the
  plateau); B2 is a charge-controlled conductance g = V(q_j)/tau; the plateau
  ends when C1 (drained at I_rec/tau by B1) can no longer sustain g·V = I_RM:
  snap at V(q_j) ≈ I_RM·tau/V_R = 1.25 V predicted, 1.33 V observed.
- **q_0 sets timing (Qrr, t_rr ≈ linear in q_0); tau sets shape**
  (plateau→tail departure, snap hardness; Qrr nearly tau-independent).

**Fit (bench = user's DiodeTest: 500 settled periods, ±12.5 A PULSE current
source, di/dt = 1000 A/µs, 50 µF/8 Ω self-biased V_R ≈ 49.3 V; k_chg = 5,
r_max = 0.1 fixed; datasheet typ @1000 A/µs: t_rr 25 ns, Qrr 247 nC):**

| variant | q_0 | tau | Qrr [nC] | t_rr(25 % I_RM) [ns] |
|---|---|---|---|---|
| **recommended** | **35 n** | **5** | **249.6 (+1.1 %)** | 27.6 (+10 %) |
| datasheet-strict | 34 n | 3 | ≈248 | ≈26.8 |
| waveform-image-strict | 37 n | 5 | 266 | 29.0 |

No (q_0, tau) hits both targets: at Qrr = 247 the t_rr floor is ≈26 ns. At
di/dt = 100 A/µs the fit gives Qrr 77 nC vs datasheet 36 — the datasheet's
Qrr *rises* with di/dt (36 → 247 nC), unreproducible by any single
stored-charge state (those low-di/dt numbers are "defined by design"); the
model is tuned to the hard-switching (1000 A/µs) condition. **Next session:
q_0/tau extraction protocol from datasheet info** (t_rr/Qrr/I_RM pairs,
di/dt scaling, two-knob limitation quantified first; candidate extra knobs
only if the residual demands them).

**Integration decisions (all simulation-verified, LTspice 26.0.1):**
- Caps delta exports **Cgs/Cgd only**; the measured C_ds charge-master
  (Cds branch of `out/ISC040N10NM8_caps.subckt`, ramp-sum inline in the
  template) moves **inside the recovery branch at Ki–A**. At the D–S
  terminals: convergence failure at hard forward re-turn-on (D1 has no Cjo;
  B3+R1 isolate the branch — timestep 1.25e-19 crash). At Ki–Ai: singular
  (user's original finding, reproduced). Ki–A is electrically D–S (B3 is a
  voltage source; R1 5 mΩ + R2 1 µΩ → fc ≈ 3 GHz) so the caps Δ stays
  complete. The 500 p placeholder was a good approximation: q_ds(50 V) =
  27 nC ⇒ avg C(0→50 V) = 541 pF; with the real Cds, Qrr moves <1 %.
- **Half-bridge stability (user fix, adopted as template default):**
  `R2 = 1 µΩ` (was 1 nΩ) and B3 carries a **+1 µΩ floor**:
  `V=-I(R2)*(r_max/(1+k_chg*V(q_j))+1u)`. In a bridge, B3 (transimpedance
  source) + R1/R2 + switch Ron + the opposite device form near-singular
  loops. A/B over r2 × floor: waveform identity < 0.1 % (t_rr 30.369 →
  30.382 ns, Qrr 248.55 → 248.73 nC, V_SD +12 µV @ 12.5 A).
- `method=gear reltol=1e-4` required (trapezoidal rings on the snap).
- Bench gotcha: `.tran` Tstart shifts the raw/meas time axis to 0 — `.meas
  TD` must be relative (TD=0 happened to work in the original bench).
- Forward recovery: B3 = 0.1 Ω cold → 1.6 mΩ hot (V(q_j) = I_F), τ_decay = q_0.
- DC diode params Is=1n N=1.4 Rs=5m → V_SD(25 A) = 0.91 V (datasheet typ
  0.82, max 1.0; Is ≈ 2.5 n would hit typ; q_0/tau insensitive).

## Body-diode recovery — transit-time pole, intrinsic doctrine, extraction protocol (sessions 25a–25c)

**Model reworked to the pole gauge; all four ISC040 datasheet recovery
numbers fit to RMS 1.12 %; the datasheet-driven extraction protocol now
exists as `src/vdmos_fit/extract_recov.py`. Durable records:
`out/RECOV_INTRINSIC_FIT.md` (25a), `out/RECOV_POLE_FIT.md` (25b/25c),
figure `out/RECOV_pole_fit.png`; bench `outputs/vdmos_recov/DiodeTest_pole.net`.**

1. **Intrinsic doctrine (25a).** The datasheet Qrr does NOT include the
   C_oss excursion: fit `q_rri`/`t_rri` measured on the INTRINSIC branch
   current `I(x1:R2)` (C2 displacement bypasses R2 at node A). The terminal
   path is contaminated by ≈28 nC / ≈15 ns at 100 A/µs — the intrinsic
   measurement is what makes the low-di/dt datasheet point comparable at all.
2. **Transit-time pole (25b, user's idea B).** `E1` (unity buffer) +
   `R3 = 1 Ω` + `C3 = {t_pole}` between the storage node and the
   conduction-driving node — an EXACT pole (no C1 loading) modelling finite
   carrier transit/sweep-out memory. Gauge cleaned: `C1={t_minor}` +
   `ClampD Roff=1` → storage τ = t_minor, V(q_j0)_ss = I_F volts, stored
   charge = I_F·t_minor; **q_0 deleted (pure gauge)**. Mechanism verified on
   raw: I_RM occurs 1–2 ns AFTER charge exhaustion (V(q_j) holds +0.29 V at
   the q_j0 clamp instant, decays with exactly R3·C3); snap condition
   `V(q_j) = I_RM·tau/ΔV` exact to all printed digits.
   **The session-23 verdict "no single stored-charge state follows both
   di/dt points" is DEAD**: the pole-free Qrr_lo/Qrr_hi floor (0.175
   intrinsic; 0.30 terminal) was a model artifact — the datasheet ratio
   0.146 is reachable through the pole.
3. **Protocol + fitter (25c).** dTmax DROPPED (the pole rounds the snap:
   none/2.5n/1n/0.5n agree to ≤0.01 % on q_rri/t_rri and none is FASTER —
   4088 vs 8368 iterations at 100 A/µs; the "dTmax mandatory" finding was
   pole-free-era). tau stays FREE: genuine interior minimum at ~3e-4
   (analytic delay law Δt = t_pole·ln(V(q_j)|clamp/V_th),
   V_th = I_RM·tau/ΔV, d(Δt)/d ln tau = −t_pole, verified to ~1 %).
   Polish TRF in log space with `diff_step=0.03` (1e-4 puts cost
   differences below the simulator's own ~1e-4 relative noise →
   noise-dominated Jacobian, tau stalls at 0.119, RMS 5.57 %), generous
   step budget (tau is a flat valley coupled to both other knobs), bounds
   tau ≥ 1e-6. Single `.step param k` mixed-radix index sweep per di/dt
   point (deterministic; kills the `.step` nesting gotcha). Chained
   `{i_fi_min}` stepped-table gotcha RETRACTED (bookkeeping error; verified
   working in `outputs/vdmos_recov/decomp.log`).

**Official fit (ISC040N10NM8, V_R = 50 V, I_F = 12.5 A, 25 °C, intrinsic
path, t_rr = 10 % of I_RM, record `ISC040_RECOV`):** tau = 3.789e-4,
t_minor = 22.363 n, t_pole = 0.8324 n → Qrr 248.2 nC (+0.5 %) / t_rri
24.69 ns (−1.3 %) @1000 A/µs; Qrr 35.71 nC (−0.8 %) / t_rri 29.46 ns
(+1.6 %) @100 A/µs; **RMS 1.11 % over all four numbers** (17.45 %
pole-free two-knob). Fix-and-refit penalties at the optimum
(x0.5/x2 perturbations): tau +0.13/+0.08 pp (flat valley), t_minor
+9.6/+86 pp (dominant), t_pole +8.1/+5.4 pp (sharp).
t_pole SHARPLY identified (fix-and-refit: 0.3 n → 15.6 %, 0.6 n → 10.8 %,
2 n → 11.6 %, 4 n → 27.5 %, 8 n → 52.8 %, →0 → 17.45 %). Physical anchor:
L = sqrt(2·D_a·t_pole) = 1.4–2.7 µm (diffusive scale, right for the
depletion-adjacent region of a 100 V SJ device; NOT a drift transit time).
Falsifiable out-of-sample prediction at 300 A/µs: Qrr ≈ 104 nC,
t_rr ≈ 29 ns, I_RM ≈ 6.7 A (a real DPT discriminates against the ~90 nC
power-law interpolation).

**Extraction driver** `src/vdmos_fit/extract_recov.py` (pole gauge,
`RecovSpec` per part): batched grid runs → log-log multilinear surrogates
(Qrr, t_rr) per di/dt point → multi-start surrogate inverse solve → TRF
polish on real runs → verification at every datasheet point (intrinsic +
terminal-path reference numbers) → fix-and-refit identifiability →
intrinsic-vs-terminal waveform plot. Report
`out/recov_extract_isc040_run/extract_report.txt`, figure
`out/RECOV_pole_extract_waveforms.png`. Official record `ISC040_RECOV` in
`reference.py`, pinned by `tests/test_recov_record.py`. The subckt template
`out/ISC040N10NM8_recov.subckt` is regenerated in the pole gauge — do NOT
mix the old q_0/tau gauge with this one.

## ISC040N10NM8 full device — channel record, diode DC record, 4-terminal subckt (sessions 27–28)

**Official records `ISC040_CHANNEL` + `ISC040_DIODE` in `reference.py`,
pinned by `tests/test_isc040_official.py` (11 tests); reproducible protocols
`scripts/fit_isc040_channel.py` + `scripts/fit_isc040_diode.py`; the
deliverable `out/ISC040N10NM8_full.subckt` generated by
`scripts/gen_full_subckt.py` (edit the records, never hand-edit).**

1. **Channel (session 27).** SJ charge-cap family (v_ovmax/k_ovsharp free —
   genuinely sensed: without the cap the Ron chord is ~13 % low; k_theta=0,
   v_vsat=1e4 off, k_lambda=0, g_k=0) — same family as ST8L65, beats the
   velocity-sat form (1.33 % vs 1.73 % at 25 C). User protocol (three notes):
   (i) NO transfer curve (datasheet V_ds conditions vague) → output curves +
   Rds(on) chords only; (ii) NOT Kelvin → r_s = r_sk = 0 and the data REJECT
   gate-loop degeneration (r_s → 0 from every start, monotonically worse as
   r_s grows — the IMCQ session-22 finding), whole current-loop series R in
   r_d; (iii) fit 25 C FIRST and LOCK, then TC (no high-T IV — TC comes from
   Rds(on)-vs-Vgs at −55/25/100/175 °C (I_test = 25 A chord, residual as
   CURRENT at R_data·25 A) + Vth-vs-Tj at 520 µA/52 µA with a Vth-VOLTAGE
   residual (solve the model gate V at V_ds = V_gs giving I_d, vs datasheet
   Vth) — fits the Vth(T) slope without the exponential-current sensitivity
   the single-softplus subthreshold cannot reproduce). Stage 1: 132 pts,
   RMSE 1.33 %, r_d ≤ 3.4 mΩ (datasheet 25 C chord). Stage 2: s_vth/beta_p/
   a_rd/b_rd free; TC laws k_p(T) = k_p25·(T_K/298.15)^−β_p,
   v_th0(T) = v_th0_25 + s_vth·ΔT, r_d(T) quadratic Trs1/Trs2 — NEVER linear
   (goes negative below ~−150 °C; vertex ~−130 °C, valid over the data span);
   shape knobs T-invariant. TC overall 2.97 % (Rds(on) 2.35 %, Vth 4.04 %).
   β_p essential (without it Rdson RMSE → 13 %) and identifiable ONLY because
   Rds(on)-Vgs(T) pins r_d(T) at high V_gs (the r_d-dominated floor) leaving
   the low-V_gs channel shape for β_p — IMCQ doctrine. Acceptance specs
   (never fitted; pinned): Ron@25A 3.35/3.56/3.89 mΩ @Vgs 15/10/8 (data
   3.33/3.52/3.87); Ron@Vgs15.7 vs T 2.12/3.34/5.20/7.76 mΩ @−55/25/100/175
   (data 2.05/3.33/5.23/8.00); Vth@520µA 3.52/3.06/2.50/2.07 V. s_vth
   −7.49 mV/K vs the datasheet Vth_Tj slope ~−6.3 mV/K (within model form).
   Identifiability (quote curves, never the split): r_d↔r_s corr −1.0 (only
   the sum sensed); k_p↔v_th0↔v_ovmax↔r_d cluster |corr| ≤ 0.86 (the
   composite channel conductance is the record); β_p σ ≈ 4; v_dbsat σ ≈ 77
   (inert — data stops at V_ds ≈ 5 V).

2. **Diode DC (session 28).** `ISC040_DIODE` (plain `Params`): free i_s,
   k_n, r_diode, i_kf, v_eg, a_rdiode; fixed k_xti = 3, i_sr = 0, k_nr = 2,
   t_nom = 25. **KEY FINDING: k_n = 2.0 (the ST8L65 fixed value) is WRONG
   here** — with k_n free RMSE drops 35 % → 3.31 %, k_n lands ≈1.05
   (near-ideal diffusion, unlike the SJ channel). i_s↔k_n corr +0.999
   (degenerate prefactor pair — quote the pair/curve, never i_s alone; the
   test pins the pair through the forward curve: vdmos_diode(1 V, 25 C) =
   152.6 A). i_kf (≈357 A knee) and v_eg (≈1.11 eV, silicon) well sensed;
   recombination branch NOT needed; **linear a_rdiode is the decisive TC
   (3.31 % → 0.51 %, σ 0.013)**; quadratic b_rdiode and a_ikf buy nothing
   with only 2 temperatures and go degenerate (a_rdiode↔b_rdiode corr −1.0)
   → stay pinned at 0. RMSE 0.51 % (25 °C 0.45 %, 175 °C 0.58 %). r_diode
   ≈ 0.96 mΩ is ~10× SMALLER than the old recovery-bench R1 = 10 mΩ.

3. **Full 4-terminal subckt (D G S K).** Assembly:
   `D–[Rd]–(dx)≡channel B-source≡(sx)–[Rs]–S` plus Xcaps (Cgs/Cgd delta at
   terminals) plus DiodeRecov (anode = source, cathode = drain, at terminals;
   Cds lives inside the branch); K = kelvin tap, inert at r_s = r_sk = 0.
   **Series-R single-count reconciliation:** the D1 model card
   `D(Is Rs=0 N Ikf Eg Xti Tnom)` carries the DC set, and the diode's own
   r_diode(T) lives in R1 as `{Rs*(1+Trs1*(temp-Tnom))}` — putting Rs on
   BOTH D1 and R1 would double-count (~1.92 mΩ → ~0.7 V error at 367 A).
   The `temp` reserved variable works in R value expressions (verified), so
   the diode branch is temperature-complete while the channel branch stays
   25 C (TC laws documented in the subckt header, apply by re-evaluating the
   record). Recovery knobs from `ISC040_RECOV` (record UNCHANGED — it
   documents the extraction-run bench conditions; D1 params come from
   `ISC040_DIODE`). LTspice gotchas (do not re-derive): R=0 is rejected →
   emit 1 µΩ floors for Rs/Rsk when r_s = r_sk = 0; `Rleak dx sx 100Meg`
   silences the benign floating-node warning (<1e−6 effect); the Xrec
   instance line must be ONE line; `.include` in a run dir needs an ABSOLUTE
   path; tie K→S by node in benches (a 3rd 0 V source makes a V-loop);
   chained `.meas` references DO work in stepped runs (25c correction).

4. **LTspice verification (26.1.1 batch, `-b -w`):**
   - Channel parity (`out/parity_isc040_full_meas.cir`, reltol 1e-8): worst
     333 ppm vs the RAW record — fully explained by the 1 µΩ Rs floor (it
     carries the full channel current: I·1 µΩ ≈ 1 mV gate degeneration at
     1 kA) + .6g netlist formatting; vs the NETLIST-CONSTANT Python model
     (floor + formatted values): **0.05 ppm**. The session-27 "<0.001 %"
     claim holds only against a floor-aware comparison.
   - Forward body-diode V_SD (`out/diode_fwd_25.cir` / `_175.cir`, `.dc
     LIST` at the exact CSV V_sd points — no FIND interpolation error):
     LTspice vs netlist-constant model **≤27 ppm at BOTH 25 and 175 °C**
     (residual = LTspice's internal physical constants, same level as the
     test_diode cross-validation); vs the CSV the RMSE equals the fit's own
     residual — the branch adds nothing. Confirms the Eg/Xti/Tnom card TC +
     R1 Trs1 law end-to-end.
   - Recovery D1 A/B (`out/recov_d1_check.cir`, clean-ASCII rewrite of
     recov_r1_sens.cir which carries stray 0xA7/0xB5 bytes; `.step d1 0/1`
     old bench D1 vs DC set, R1 = 0.961 m both): q_rri 248.17 → 248.06 nC,
     t_rri 24.684 → 24.702 ns (≤0.075 %) — recovery INSENSITIVE to the D1
     Is/N/Ikf swap (plateau carried by B2, I(D1) ≈ −1 nA carries nothing),
     and the record's achieved 248.20/24.687 reproduces at the new R1
     (consistent with the measured R1 insensitivity over 10 µ–10 m).

5. **Caps parity pin RE-BASELINED on LTspice 26.1.1.**
   `test_caps.py::test_ltspice_caps_parity` moved from rel=2e-3 (written on
   26.0.1, whose adaptive step sequence happened to land the FIND samples
   cleanly inside the fine stair segments) to **rel=1e-2**: the sample
   targets (digitized-interval midpoints) coincide with fine stair
   breakpoints (12 samples per digitized interval), so which side of a
   breakpoint the sample straddles is solver/version sensitive — measured
   worst on 26.1.1 is 0.56 % (c2; c1 0.32 %, all others ≤0.11 %). A
   fixed-step harness only MOVES which point straddles (c1 → a2): the
   knife-edge is intrinsic to the sample placement. A future redesign should
   move the sample points off the fine breakpoints rather than tightening
   this pin back.

6. **Data/loader gotchas:** the ISC040 CSVs carry a `T_j (C)` temperature
   column → package `load_diode_csv()` (expects `T`) FAILS; use script-local
   loaders (like the IMCQ spec CSVs). The Rdson_Vgs figure is a **25 A
   chord** (verified: Vgs=15.7 → 3.33 mΩ matches the IV_25C row
   Vgs=15, Vds=0.08325, I=25) — convert to V_ds = R·I for the fit. The two
   interpreters on this PC (scipy 1.16 vs 1.18) can drift fits in flat
   valleys — pinned tests run under `tooling/engineering313`.

## Third-quadrant channel, expression cost, fast variant — buck-converter field report (session 29)

User ran the subckt in a synchronous buck and reported three problems; all
three landed.

1. **Third-quadrant channel conduction — LANDED.** `vdmos_id` is now
   mirror-symmetric: `i_ch = sign(v_ds)·i_ch_fwd(|v_ds|, v_gs)` (the
   package-R solver runs on |v_ds| magnitudes, so the forward bracket
   proof is untouched; the solver keeps its v_ds ≥ 0 guard). The channel
   sits in parallel with the diode at v_ds < 0 — real-device topology.
   Multiplying by ±1.0 is exact, so every forward-polarity anchor is
   bit-identical. **Known limitation (user-accepted):** the mirror cannot
   capture V_gs modulating the diode's *forward* voltage — matters for
   SiC with negative drive; GaN (no body diode) relies on the channel
   anyway. In the LTspice subckt the mirror is `abs(V(dx,sx))` inside the
   forward expression times a nested-`if()` sign factor (LTspice has no
   `sign()` function). Verified: new bench `out/parity_isc040_rev_meas.cir`
   (v_gs=10, v_ds −0.01…−0.5 V, up to −140 A): worst **7.5 ppm** vs the
   netlist-constant model (281 ppm vs the raw record — the 1 µΩ Rs floor +
   .6g formatting, same as forward). Buck effect: V(sw) low
   −0.751 V → **−0.022 V** (channel takes over the freewheel; the diode
   only carries dead time), V_out 3.55 → 4.00 V, iterations 45088 → 39554
   (the recovery branch is no longer engaged all low-side phase).
2. **Duplicate calculations in the channel expression — measured, kept
   flat.** The flat B-source tree re-expands the drain factor ~90×, the
   DIBL tanh and the softplus 28× each, the power-law cap/min 10× per
   evaluation. Buffering each shared sub-expression on its own internal
   voltage-source node was tried and is a **measured regression**: 8 extra
   matrix/FD nodes cost more than the duplicated evaluation —
   `out/turnon_flat.cir` vs `out/turnon_buf.cir` (identical physics):
   6517 vs 9161 iterations, 0.29 vs 0.38 s. LTspice's expression evaluator
   is cheap; matrix/convergence overhead dominates. Flat expression kept;
   the negative result is documented in `scripts/gen_full_subckt.py`.
   (The 10835-iteration re-run with mirror enabled is *physics*: gear
   ringing undershoots v_ds below 0 during turn-on and the mirror now
   conducts there, as a real MOSFET does.)
3. **Fast variant — `out/ISC040N10NM8_fast.subckt`.** Same subckt name,
   drop-in `.include` swap (never include both): NO `DiodeRecov`/Xrec
   branch; plain `D1 s d MyDiode` with **tt = t_minor + t_pole/0.69 =
   23.5698n** (0.69 ≈ ln 2: the pole-gauge storage maps to the SPICE
   diode's exponential decay time constant); series R single-counted in
   the model card (`Rs=` + `Trs1=`, no R1); Cds charge-master relocated
   verbatim to the D-S terminals (its Q= expression is placement-agnostic).
   Buck A/B (same body, 60 µs): full 1.73–1.91 s wall / 39554 iterations /
   1.43 s sim; **fast 1.58–1.60 s / 31643 iterations / 1.12 s sim**
   (−20 % iterations) and the DC operating point now solves directly —
   the Gmin-stepping struggle came from the recovery network's floating
   storage nodes. **Tradeoff, documented in the file header:** fast gives
   first-order storage (Q ≈ tt·I_F) but a hard snap — measured recovery
   spike −152 A vs −3.9 A soft snap for the full model. FULL owns the
   datasheet recovery numbers (RMS 1.12 %); fast is for speed.
4. **User's large-buck field data (2026-11-08, their board netlist, not
   ours):** ImprovedVDMOS 22.266 s / 1,936,021 iterations; native VDMOS +
   external FRD diode 6.323 s / 2,388,137; ImprovedVDMOS + FRD 28.068 s /
   3,149,122. Per-iteration cost: **11.5 / 2.6 / 8.9 µs/iter** — on their
   bench per-evaluation cost dominates (big matrix; the B-source tree is
   evaluated as bytecode with FD Jacobian perturbations), and the external
   FRD interacting with the model's own diode+channel branch makes it
   *worse* on both axes (stiff current division near the zero crossing).
   The fast variant is the aimed answer; the third-quadrant mirror also
   removes the always-on diode/recovery conduction that made xrec the
   convergence dominator in their report.
5. **Parity re-verified on the final subckt pair:** forward channel 40 pts
   worst 332.8 ppm (floor effect, unchanged); diode fwd 25/175 °C worst
   26.5/25.6 ppm (a transient 953 ppm scare was the *comparison script*
   double-counting the Trs1 law and forgetting the 1 µΩ R2 — the subckt
   was fine); recovery A/B reproduces the pinned q_rri/t_rri/i_fi_min
   exactly.
6. **LTspice gotchas learned (batch, 26.1.1):** no `sign()` function
   (nested `if()`); `.meas` window keyword is `INTEG`; netlist source is
   `PULSE(...)`; `.tran Tstart` shifts the RAW axis to 0 but `.meas
   FROM/TO` stay ABSOLUTE (an absolute window with Tstart silently FAILs).

## Release pipeline — input contract, presets, staged fit, one-command build (session 30, 2026-11-08)

Six problems raised by the user (no unified driver; no input contract;
protocol choices hardcoded; recovery doesn't eat CSVs; environment not
pinned; forward-recovery params in the subckt unfitted) — Phases A–D
landed; Phase E (GUI) pending.

- **Phase A — input contract (`ingest.py`, 27 tests in `test_ingest.py`):**
  4-column IV contract `NAME[unit]`. Grammar: unit in `[...]` or `(...)`;
  name = case-insensitive, non-alnum stripped; unit case-SENSITIVE except
  k/K and u/µ/μ fold; prefixes f p n u m (none) k M G T P; bases
  V A F OHM/Ω C s sec; composite di/dt `A/<time>` ([MA/s] ≡ [A/µs]);
  reverse-biased spellings NEGATE (V_sd→−v_ds, I_sd→−i_ds, V_sg→−v_gs);
  per-FILE defaults (no V_ds column → V_ds = V_gs only for diode-connected
  data, warned; no gate column → V_gs = 0; no temp column → ref_temp;
  bare `T` is CELSIUS, never Kelvin-by-default); typed classification
  IV/CAP/QG/RDS/VTH/REC; condition strings ("V_ds=20V,T_j=25C") parseable.
  Regression test: all 18 repo data CSVs parse.
- **Row-split doctrine (diode vs channel):** i_ds > 0 → channel;
  i_ds < 0 with gate OFF → diode; i_ds < 0 with gate ON → third-quadrant
  mirror data (reported, excluded from the diode fit).
- **Phase B — presets + stages (`presets.py`, `stages.py`, 6 tests in
  `test_stages.py`):** fit order channel@ref → diode@ref → TC. Presets
  `trench | sj | sic` = the three official recipes (user mapping:
  ISC040 = trench, ST8L65 = SJ, IMCQ120 = SiC; the cap-family choice
  lives in the preset). TC residuals: Rds(on) chords as CURRENT at the
  per-row test current (column > condition string > `rdson_i_test` arg),
  Vth as VOLTAGE via brentq at v_ds = v_gs. Equivalence pin: the staged
  refit reproduces the ISC040 record's s_vth/beta_p/a_rd/b_rd BIT-FOR-BIT
  (rel 1e-6).
- **Phase C — export engine (`export_subckt.py`):** parameterized
  FULL/FAST builders (part name, channel p25 + TC, diode, recovery knobs,
  cds/caps text, r_max); `scripts/gen_full_subckt.py` is now a thin
  wrapper supplying the ISC040 records + its header docs.
  `channel_expr_auto`: collapsed form (byte-identical to the pinned
  ISC040 files) when k_theta = k_beta_trap = k_lambda = g_k = 0, else the
  full-family transcription of `model.py` — NOT yet LTspice-parity-benched
  for a non-collapsed fit; the export says so.
- **B3 rule (user, session 30):** B3 is emitted ONLY with a non-zero
  `r_max`; at r_max = 0 (the ISC040 record) the branch is R1 tied
  straight to K and orphan node N001 is gone. Physics re-verified on the
  regenerated pair — every bench number BIT-IDENTICAL to the session-29
  pins (parity 332.8 ppm fwd / 7.5 ppm rev; diode fwd 26.5/25.6 ppm;
  recovery q_rri 248.171/248.059 nC, t_rri 24.684/24.702 ns, i_fi_min
  −18.436/−18.445 A; buck −0.022 V switch undershoot): the 1 µΩ B3 floor
  was below the noise, exactly as the R1-insensitivity study predicted.
  No re-pinning needed. Structure pins updated
  (`test_isc040_official.py`: `R1 K Ki {Rs*(1+Trs1*(temp-Tnom))}`, no
  B3/N001 in the netlist).
- **Phase C — orchestrator (`release.py` + `scripts/build_part.py`):**
  one command, CSVs → ingest → channel@ref → diode@ref → TC → caps
  (+ c_gdmax from Qg via `extract_cmax`) → recovery (OPT-IN `--recovery`:
  needs a REC CSV — Q_rr, t_rr, I_RM, I_sd, di/dt — and LTspice; drives
  the `extract_recov.py` grid) → export. Export rule (user): recovery
  extracted → FULL subckt; no recovery data → FAST subckt with a plain
  diode and NO tt. Outputs: `out/<PART>_full|fast.subckt`,
  `out/<PART>_caps.subckt`, `data/<PART>_cds_expr.txt`,
  `out/<PART>_record.json` (params + per-stage rmse + warnings + library
  versions). No diode data → generic placeholder card with a LOUD warning.
- **Phase D — e2e golden test (`test_release_e2e.py`, 3 tests):** ISC040
  CSVs → orchestrator → FAST subckt (no recovery CSV in repo → no tt,
  no recovery network) + record JSON quality bands. `requirements.txt`:
  floors (numpy ≥ 2.0, scipy ≥ 1.13, matplotlib ≥ 3.8) + verified stack
  (Python 3.13.5, numpy 2.5.3, scipy 1.18.1, matplotlib 3.11.2, LTspice
  26.1.1) + the scipy flat-valley drift note.
- **Phase E — GUI (`scripts/release_gui.py`, Tkinter):** five columns
  over the package API (nothing re-implemented in the shell):
  (1) file list + preset dropdown + ref temp + Parse with the ingest
  report; (2) parsed tables + rows, displayNum-formatted; (3) channel
  and diode parameter trees with click-to-toggle FREE checkboxes
  (`free_override` added to `fit_channel_ref`/`fit_diode_stage`), fit
  buttons running on a worker thread (window stays live); (4) caps
  points + ramp-sum resolution knob `n_sub`; (5) plots (IV
  measured-vs-model, caps) + optional recovery extraction (LTspice
  grid, thread) + Save subckt through `release.export_part` — same
  export rule (recovery -> FULL, none -> FAST, no tt). Reusable
  release helpers `caps_model_from_tables` / `recov_from_table` /
  `export_part` extracted so CLI and GUI share one path.
  `displayNum` lives in `vdmos_fit.display`, FIXED for negatives/zero
  (the log10 domain-error crash) — 14 tests in `test_display.py`.
  Headless smoke: full 5-column build + ISC040 parse verified.
- **Condition-column ingest (session 30, user request):** the IV/RDS/VTH
  loaders now pull per-row V_ds, V_gs, T_j and I_test from a `condition`
  column in the IMCQ120R004M2H_Transfer.csv style ("V_ds=20V,T_j=25C")
  when no real column carries them — precedence column > condition >
  documented default, ALL-or-nothing per file (partial condition rows
  never half-fill a column); condition temperatures are canonical Celsius
  ("T_j=298.15K" -> 25). This makes the Transfer file's three T_j curves
  (25/175/200 °C at V_ds = 20 V, one CSV) parse correctly instead of
  defaulting to ref-temp/diode-connected. User data edits land the same
  doctrine as explicit columns: I_d=25 A into ISC040_Rdson_Vgs.csv,
  I_ds=185.2 A into IMCQ_Rdson_Tj, T_j[C] into IMCQ_IV_175C.csv — the
  resolve_test_current pin updated to match (the raise path is now pinned
  on a synthetic no-current file). GUI Save gained an overwrite
  confirmation after a part-name collision clobbered the official ISC040
  artifacts once (regenerated via gen_full_subckt.py).
- Kelvin fit crash fix (GUI field report, IMCQ with r_s/r_sk unlocked):
  the model guard 0 <= r_sk <= r_s (model.py) fired MID-FIT because the
  solver explored r_sk and r_s independently. fit() now adopts the
  session-22 doctrine: when r_sk is free the optimisation column is the
  RATIO f_sk = r_sk/r_s in [0, 1], r_sk = f_sk * r_s -- the guard is
  structurally unfireable during a fit, FitResult still reports r_sk,
  and f_sk pinned at 1 means the ideal kelvin tap. Numerically neutral
  for every non-kelvin fit (all official records reproduce bit-for-bit);
  GUI warns when r_sk is freed without r_s.
- Session-30 GUI round 3 (five user changes): (1) SIC preset diode
  support ENABLED (trench-style free set, v_eg start 1.4 for the wider
  SiC gap, r_diode bounds to 50 mOhm) -- the stage still self-skips when
  the IV tables carry no gate-off reverse rows; (2) parameters are now
  hand-editable: double-click a value cell, engineering suffixes
  accepted (display.parseEng, inverse of displayNum, pinned round-trip),
  the trees are the live parameter state and fits write back into them;
  (3) user rule v_ovmax = 0 whenever v_ovmax is NOT free -- enforced in
  fit_channel_ref (start build, after multistart merge) and again in GUI
  _save, so a locked preset start value can never leak a charge-cap term
  into an export; (4) Plot diode button (measured diode rows vs the
  card in the tree, per T_j) and (5) Plot TC button (measured
  Rds(on)/Vth vs T_j and the TC law applied to the tree params; flat TC
  when unfitted, brentq failures plotted as gaps). All plots now rebuild
  their axes via fig.clear (no stale self.ax).
- Session-30 GUI round 4 (user list): (1) TC plot now splits VERTICALLY
  (Rds(on) top, Vth bottom); (2) the parameter trees gained editable
  min/max columns (param|min|value|max|free) -- defaults from
  params.BOUNDS (diode tree shows BOUNDS widened by the preset's
  diode_bounds, i.e. what fit_diode_stage actually uses), edits flow
  into bounds_overrides at fit time (stages fit_channel_ref/
  fit_diode_stage gained bounds_override; the r_d chord ceiling wins
  over any user max -- data beats guess; r_sk ratio bounds still win in
  fit()); inverted boxes rejected with a clear error. NOTE: extreme
  boxes CAN make the model unsolvable (a k_p ceiling starved the
  package-R solve at the 25 A chords -> RuntimeError, not a graceful
  pin) -- the pin test uses v_dbsat instead; (3) caps plot overlays the
  EXPORTED model curves (c_gs/c_gd/c_ds branch functions over each
  measured branch's range, labelled "model <kind> (exported)");
  (4) fit progress: fit()/fit_diode() take progress(msg) called every
  25 evals ("start i/n: N evals, cost X"), stages thread it through
  all three fits, the GUI shows it live in the fit label; (5) I_RM is
  now an OPTIONAL REC column (datasheets often omit the peak reverse
  recovery current; the extraction measures i_fi_min itself) -- missing
  column -> NaN + note; (6) recovery-fit zip crash FIXED: extract()/
  polish()/polish_fix_one() had a fixed 4-tuple weights default that
  strict-zip rejected for any other point count (the GUI hit it with a
  1-2-row datasheet table) -- _weights() now broadcasts (None -> ones,
  scalar -> replicated, exact length -> as-is, else clear ValueError).
 - Recovery I_F doctrine (user, session 30): the extraction NEVER
   hardcodes the forward current -- recov_spec_from_table (release.py,
   new pure helper) builds the RecovSpec from the REC table, taking
   I_F per row from the I_sd/I_F column (sign-normalized). RecovSpec
   gained i_f_points (per-point, aligned with points; empty = single
   i_f); make_netlist takes a per-bench i_f override, and extract()
   passes it to the grid runs, the polish benches AND the waveform
   runs; the report header prints the per-point list. The 12.5 A on
   RecovSpec is only the ISC040 RECORD default (pinned netlist tests
   unchanged).
 - Session-30 GUI round 5: (1) CAPS PLOT DOCTRINE FIXED — plotting raw
   branch functions against the measured V_ds axis is meaningless (the
   gd branch is indexed by v_gd ~= -V_ds; gs/ds tables by V_ds; the
   first overlay showed flat "model iss/rss" lines). The overlay must
   RECOMBINE branches into the terminal capacitance the subckt
   realizes: Ciss = c_gs(v) + c_gd(-v), Coss = c_ds(v) + c_gd(-v),
   Crss = c_gd(-v); verified model-vs-measured max deviation 0.00 /
   0.98 / 0.31 % (iss/oss/rss, ISC040 Cap_Vds). (2) RECOVERY RMS
   DISPLAY BUG: extract() returns rms_rel ALREADY in percent; the GUI
   multiplied by 100 again -> "103.55 %" for a healthy 1.04 % fit.
   Fixed; the fitted tau/t_minor/t_pole are now also logged to the
   report panel after extraction. (3) recovery extraction gained a
   progress callback (grid sweeps, per-polish RMS, LTspice-run tally
   every 10 runs, identifiability steps) wired to the GUI status line;
   (4) the "Extract recovery (slow)" button moved to the bottom of
   column 4 (with the caps data it operates on).
- Suite: 278 tests green, ruff clean.

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
- **Third-quadrant channel conduction — CLOSED (session 29)**:
  `vdmos_id` and both subckt variants are mirror-symmetric
  `sign(v_ds)·i_ch_fwd(|v_ds|, v_gs)`; see the session-29 section.
  Residual (accepted): the mirror cannot capture V_gs modulating the
  diode's forward voltage.
- **Unified driver / release pipeline — CLOSED (session 30)**:
  `scripts/build_part.py PART data/*.csv --preset trench|sj|sic
  [--recovery]` runs the whole protocol (ingest → channel@ref →
  diode@ref → TC → caps/cmax → optional recovery → export) and writes
  the subckt + record JSON; package API `vdmos_fit.release.build_part`.
  See the session-30 section. The Phase-E GUI (`scripts/release_gui.py`,
  five columns, worker-thread fits) landed on top of it (session 30).
- **Body-diode recovery — package-side state not built**: recovery exists
  as the behavioral LTspice network in the pole gauge (sessions 25a–25c:
  `DiodeRecov` + E1/R3/C3 transit-time pole, record `ISC040_RECOV` in
  `reference.py`, template `out/ISC040N10NM8_recov.subckt`, extraction
  protocol `src/vdmos_fit/extract_recov.py` — grid → log-log surrogates →
  multi-start polish → fix-and-refit identifiability); junction charge is
  covered by the caps delta (the Cds branch relocates into the recovery
  branch). Open: other parts (need recovery datasheet specs), T_j
  dependence, and whether it earns a `diode.py` counterpart (the
  `reference.py` *_RECOV record landed). In reserve, only with data:
  two-section transit ladder, a 300 A/µs DPT point, T_j recovery TC —
  do NOT invent. The recovery branch is integrated into the FULL-DEVICE
  subckt `out/ISC040N10NM8_full.subckt` (sessions 27–28) with the DC-fit
  D1 + R1 = r_diode(T); recovery verified insensitive to that swap.
  Session 29 adds the speed-first sibling `out/ISC040N10NM8_fast.subckt`
  (no recovery network; plain D with tt = t_minor + t_pole/0.69), same
  subckt name, drop-in `.include` swap — see the session-29 section.
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
  match the exported staircase to ≤0.3 % (most ≤0.05 %) on LTspice 26.0.1;
  on 26.1.1 the pin is re-baselined to 1 % (stair-breakpoint knife-edge,
  worst measured 0.56 % — see the sessions-27–28 section);
  staircase-vs-PCHIP ≤0.4 %. LTspice 26.0.1 gotchas on record:
  `table()` in a B source is rejected outright, and `table()` inside Q=
  compiles but its symbolic derivative is garbage (~1e11 A) — hence the
  ramp-sum form. The full-device subckt (channel B source + diode +
  r_d/r_s wiring) LANDED sessions 27–28 for ISC040N10NM8 —
  `scripts/gen_full_subckt.py` (thin wrapper over the generic engine `vdmos_fit/export_subckt.py`
 since session 30) → `out/ISC040N10NM8_full.subckt`; the
  diode part is LTspice-exact: `D(Is=… N=… Isr=… Nr=… Ikf=… Rs=…
  Tnom=… Xti=… Eg=…)`.
- **CI**. (`notebooks/` is empty scaffolding.)
- Validation target remains double-pulse Eon/Eoff within ~10 %, which static
  fitting alone cannot deliver.
