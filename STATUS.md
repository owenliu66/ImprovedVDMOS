# Project Status — ImprovedVDMOS

*Snapshot date: 2026-10-11 (blocking problems 1–2 resolved). Companion to [README.md](README.md).*

## What exists

```
Projects/ImprovedVDMOS/
├── pyproject.toml          deps, pytest, ruff config, `vdmos-fit` console script
├── README.md               model + corrected Level-1 analysis
├── .gitignore
├── src/vdmos_fit/          model, params, data, fit, plotting, cli, __main__
├── tests/                  test_model.py, test_fit.py
├── scripts/                fit_datasheet.py, roundtrip_demo.py
├── data/                   README.md only (no datasets yet)
└── notebooks/              README.md only (empty)
```

## What it currently does

- **Model** — `vdmos_id()` evaluates the unified continuous drain current
  (softplus subthreshold, mobility degradation via Vk, velocity saturation,
  DIBL, gated CLM, drift conductance, one-sided kink sigmoid).
  `vdmos_id_level1()` is now the **correct** LTspice reference with knee at
  `Vov/mtriode`.
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

## Verified state

- **31 tests passed** (verified 2026-10-11 run, including the new
  `test_matched_ron_square_law_overprediction` regression and the tightened
  alpha pin).
- **`ruff check` still reports 1 error**: unused `import math` in
  `tests/test_fit.py:5`. The `E501` and `E731` issues are resolved.
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

6. No real datasheet data in `data/` — everything is synthetic, so the fitter
   has never seen a real device.
7. No Q(V) charge model, no body diode, no thermal node, no `Tj` — the rest of
   the agreed architecture.
8. No LTspice export path: fitted parameters cannot yet be written back as a
   `.subckt`/`.model` for A/B simulation.
9. No git repo initialised; no CI.
10. `notebooks/` is empty scaffolding.

## Key corrections made this session

| Claim | Before | After (LTspice-verified) |
|---|---|---|
| Level-1 knee location | `Vds = Vov` | `Vds = Vov / mtriode` |
| Level-1 continuity | "gds jump ~+400 %, −40 S NDR dip" | **C1 continuous for any mtriode**; artifacts were a bug in the project's own reference formula |
| Kink gate | bare `tanh((Vgs−Vthk)/kk)` (odd → negative off-state current) | one-sided `0.5·(1+tanh(…))` |
| Level-1 limit of improved model | "strict superset, exact" | smooth approximation: exact in ohmic + deep saturation, ~6 % low at the knee |
| `parameter_uncertainty` | zero Jacobian column → σ = 0 ("perfectly determined") | `identifiability()` → σ = `inf` ("NOT SENSED"), rank/cond/correlations |
