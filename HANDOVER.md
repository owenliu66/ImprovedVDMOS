# Handover — ImprovedVDMOS session 18: IMCQ120R004M2H two-temperature IV fits DONE
## NEXT SESSION = implement the interface-trap-filling gain `v_ov_eff2` / `k_beta`

Session 17 (caps/charge model, LTspice export + parity) is COMPLETE and committed;
its durable findings live in `STATUS.md` + `CAP_SPEC.md` + `reference.py` (the
LTspice 26.0.1 cap gotchas are recorded there — do not re-derive them).
This file replaces the session-17 handover.

Session 18 added the **SiC IV data + two-temperature fits** (uncommitted: two new
`scripts/`, `out/` artifacts) and diagnosed *why* the fit quality is capped. The
fix is a new mechanism the user has now specified and named. **Implement it first**
— it unblocks the C_gd blend item in the queue.

---

# PART 1 — THE TASK: `v_ov_eff2` and `k_beta`

## 1.1 Exact spec (user-directed; do not change the cap)

Keep the existing soft cap **exactly as is** — it is the silicon super-junction
knob and must stay tunable for silicon FETs:

```
v_ov      = v_sub·ln(1 + exp((v_gs − v_th)/v_sub))                 unchanged
v_ov_eff  = v_ov/(1 + (v_ov/v_ovmax)^k_ovsharp)^(1/k_ovsharp)      UNCHANGED
v_ov_eff2 = v_ov_eff·(1 + k_beta·v_ov_eff)                         NEW
```

`v_ov_eff2` becomes **the** effective overdrive for everything downstream:

```
v_knee  = max(v_ov_eff2/(1 + k_theta·v_ov_eff2), v_kneemin)   <- switches to v_ov_eff2
v_dseff = v_knee·tanh(v_ds/v_knee)  or  softmin(v_ds, v_knee, k_kneesharp)   unchanged
i_d     = k_p·v_ov_eff2·(v_dseff − v_dseff²/(2·v_knee))/(1 + v_dseff/v_vsat)
          ·(1 + k_lambda·v_ds·k_db) + g_k·v_ds·k_kink·k_db     <- prefactor switches
```

Decision to make explicitly (state it in the docstring): `k_theta` multiplies
`v_ov_eff2`, not `v_ov_eff`. Rationale: `v_knee` is the pinch-off/pinball knee set
by the *final* channel charge, same argument the cap section already makes.

## 1.2 Edit sites (verified against the current `model.py`)

| line | what |
|---|---|
| 3, 10, 12 | module docstring equation block — add the `v_ov_eff2` line |
| 130 | `v_ov_eff = _soft_cap(v_ov, p.v_ovmax, p.k_ovsharp)` → **insert `v_ov_eff2 = v_ov_eff*(1 + p.k_beta*v_ov_eff)` immediately after** |
| 134 | `v_knee = np.maximum(v_ov_eff/(1 + p.k_theta*v_ov_eff), p.v_kneemin)` → both `v_ov_eff` → `v_ov_eff2` |
| 152 | `i_ch = p.k_p*v_ov_eff*(v_dseff − 0.5*v_dseff**2/v_knee)` → `v_ov_eff2` |
| ~189 | `_soft_cap` docstring: it currently claims `0 <= dv_ov_eff/dv_ov <= 1` "keeps the package-R bracket invariant intact" — **amend** (see 1.4) |

Those are the only code uses of `v_ov_eff` in the file (grep-verified). Do **not**
touch `vdmos_id_level1()` — the Level-1 A/B adapter must stay pure.

## 1.3 `params.py` + naming

- Add `k_beta: float = 0.0` to `Params`. **Default 0.0 = exact identity**
  (`v_ov_eff2 ≡ v_ov_eff`), so every pre-gain anchor stays bit-identical — same
  doctrine as `v_ovmax = 0` and `k_kneesharp = 0`. "Off" is reached by leaving it
  out of the free set, never by driving it to the lower bound.
- **NOT in `LOG_PARAMS`** (0 is a legitimate value and cannot be logged → linear
  space, like `v_ovmax`/`k_kneesharp`).
- `BOUNDS["k_beta"] = (1e-3, 2.0)` — lower bound strictly above 0 so a *free*
  value can never mimic the sentinel (the `k_kneesharp` 0.1-floor precedent;
  measured physical scale is ~0.06 1/V, so 1e-3 is far below any real answer).
- **NOT in `DEFAULT_FREE`** — freeing it silently changes existing fits.
- Units: **1/V**. It joins the dimensional `k_` trio `k_lambda`/`k_eta`/`k_theta`;
  update the README naming-table note ("the 1/V trio" → quartet).
- `names.xlsx`: the **`v_ov_eff2` row already exists** (user's description:
  "Effective OV considering drift region saturation and interface trap fill (WBG
  devices)", suffix `ov_eff2`). **`k_beta` is MISSING — add it**: description
  "Interface-trap-filling coefficient on the effective gate overdrive (WBG
  devices)", Type `factor`, Suffix `beta`, New name `k_beta`.

## 1.4 Monotonicity and the package-R bracket — the real hazard

The cap section's claim (`0 <= dv_ov_eff/dv_ov <= 1` protects the solver) needs
re-derivation, because the gain makes that derivative **larger**:
`dv_ov_eff2/dv_ov_eff = 1 + 2·k_beta·v_ov_eff` → measured **2.36 (25 °C) and
3.00 (175 °C)** at v_gs = 18 on the IMCQ fit.

Analysis (verify, don't trust): the bracket sign only needs `f` monotone
increasing in both arguments. With
`g(i) = f(v_ds − (r_d+r_s)i, v_gs − r_s i) − i`,
`dg/di = −(r_d+r_s)·∂f/∂v_ds − r_s·∂f/∂v_gs − 1 <= −1` whenever `∂f >= 0` —
so a **positive** gain preserves the bracket. Session 18's probe used `r_s = 0`
(Kelvin source), where this is trivially true; **the `r_s > 0` case is untested**
and is exactly where silicon fits live. Real consequences to test:

1. The bracket upper bound `i_max = f(v_ds, v_gs)` is inflated by the gain →
   bisection needs more iterations; check `maxiter`/tolerance.
2. Steeper gate dependence → Newton overshoots more; the existing per-point
   **stall guard** (see `test_package_r.py::test_solver_survives_newton_two_cycle`)
   will fire more often. Add a gain analogue of that test.
3. The ohmic fixed-point contraction factor `(r_d+r_s)·∂f/∂v_gs` grows —
   fixed-point iteration stays forbidden (solver is already Newton+bisection).
4. **`k_beta < 0` is the genuine danger:** `v(1 + k_beta·v)` is non-monotone,
   derivative hits 0 at `v_ov_eff = −1/(2 k_beta)` then goes negative → `f`
   non-monotone in v_gs → multi-valued `i_d`, bracket proof dies. **Recommend
   bounding `k_beta >= 0`** (physics: trap filling only *raises* channel
   mobility). If negative values are ever wanted, clamp
   `1 + 2·k_beta·v_ov_eff >= tiny_floor` — the `diode.py` TC-clamp precedent
   exists for exactly this kind of absurd-extrapolation guard. Decide and pin.
5. `test_improved_model_has_no_negative_gds` may need its tolerance revisited:
   the known `v_vsat`-induced gds dip near the knee deepens with a sharper knee
   and a gain (README already flags this regime).

## 1.5 Composition gotcha for silicon SJ fits (important)

With the cap active, the **asymptote of `v_ov_eff2` is
`v_ovmax·(1 + k_beta·v_ovmax)`, not `v_ovmax`.** On ST8L65 (`v_ovmax ≈ 1.9`,
`k_ovsharp ≈ 10`) a freed `k_beta ≈ 0.06` lifts the plateau ~11 %, so expect a
strong **`k_beta` ↔ `v_ovmax` anti-correlation** and a compensating drop in
`v_ovmax`. Therefore:

- The official ST8L65 record (fit A, 10-param, 3.79 %, `reference.py`) must stay
  **bit-identical with `k_beta` out of the free set** — verify by re-running
  `scripts/fit_st8l65_nd.py`.
- If `k_beta` is ever freed on ST8L65, re-quote the record and report the
  correlation; do not mix records.

## 1.6 Test list (`tests/test_overdrive_gain.py` + regressions)

1. `test_k_beta_zero_bit_identical` — `k_beta = 0` reproduces the pre-gain model
   bit-for-bit on a grid (pattern of `test_gate_cap.py` / `test_drain_softness.py`).
2. `test_gain_monotone_in_vgs` — for `k_beta >= 0`, `i_d` strictly increasing in
   v_gs at fixed v_ds over the whole sampled grid (forward polarity only).
3. `test_gain_raises_cap_ceiling` — assert the composite asymptote is
   `v_ovmax(1 + k_beta·v_ovmax)`; this pins the composition semantic of 1.5.
4. `test_negative_k_beta_guard` — pin whichever decision 1.4(4) lands on.
5. `test_package_r_bracket_with_gain` — **`r_s > 0` and `k_beta > 0`**: solver
   converges on a dense grid, matches a high-accuracy `brentq` reference, no
   NaN/negative; assert the stall-guard path.
6. `test_synthetic_roundtrip_gain` — `make_synthetic` with `k_beta` in the truth
   (cap off) on a grid with dense low-v_ds transfer points → `fit` recovers it
   with finite `sigma_rel`; then the joint cap+gain case (expect the 1.5
   correlation, report it, do not fail on it).
7. Existing anchors must not move: `test_level1_limit_agreement`,
   `test_matched_ron_square_law_overprediction`, `reference.py` drift checks,
   `test_drain_softness.py`, `test_gate_cap.py`.
8. Docs: README equation block + a new section "Interface-trap-filling gain
   (`k_beta`)" mirroring the cap section (mechanism, why WBG devices want the
   opposite sign of the SJ cap, sentinel, bounds); `STATUS.md` session entry.

## 1.7 Anchor numbers to validate the implementation against (session 18, measured)

IMCQ120R004M2H, cap **off**, `k_lambda = 0`, `r_s = 0`, `k_theta = 0`,
`v_vsat = 1e4`, free = `k_p, v_th0, v_sub, r_d, k_kneesharp, k_eta, v_dbsat`
(+ gain `k_beta` in the right column):

| | baseline | + `k_beta` |
|---|---|---|
| RMSE 25 °C / 175 °C | 8.53 % / 6.13 % | **4.61 % / 3.31 %** |
| `k_beta` [1/V] | — | 0.0569 → 0.0757 (**+33 %, a real TC**) |
| `k_p` / `v_th0` / `v_sub` (25 °C) | 50.90 / 6.615 / 0.8364 | 25.75 / 6.060 / 0.7644 |
| `r_d` (25 °C → 175 °C) | 1.893 → 7.251 mΩ | 1.731 → 7.031 mΩ |
| Ron chord Vgs=18 / 15 @25 °C | 3.700 / 4.446 mΩ | **3.701 / 4.729** (data 3.70 / 4.70) |
| V_th spec point (58.3 mA) | 59.6 / 58.5 mA | 57.3 / 58.2 mA |
| worst 25 °C curve | v_gs=8: 13.6 % | v_gs=8: 9.3 % (v_gs=18: **0.37 %**) |
| `dv_ov_eff2/dv_ov_eff` @v_gs=18 | 1.000 | 2.36 / 3.00 |

The **diagnostic that motivated the mechanism** — implied channel overdrive from
the near-origin ohmic points of the 25 °C curves (subtract `r_d·i_d`, divide by
`k_p`):

| v_gs | 7 | 8 | 10 | 12 | 15 | 18 |
|---|---|---|---|---|---|---|
| implied v_ov [V] | 0.24 | 0.72 | 2.09 | 3.83 | 6.99 | 10.85 |
| dv_ov/dv_gs | 0.48 | 0.69 | 0.87 | 1.05 | **1.29** | — |

Slope rising **past 1** is impossible for the pre-gain model (softplus slope
saturates at 1, cap only reduces). A synthetic set generated with `k_beta > 0`
should reproduce this signature — good regression check.

**Do not trust the high-v_gs corner yet:** with the gain the model overshoots the
datasheet transfer figure at V_DS = 20 V (3650 A vs ~2300 A at v_gs = 18). The
output-curve data cannot see that region (curves exit the 1600 A plot frame).
Bending it over is what `k_lambda` does → **digitize the transfer characteristic
(V_DS = 20 V; 25/175/200 °C, datasheet p. 9) and re-free `k_lambda` in a joint
fit** before quoting anything above ~1500 A.

---

# PART 2 — WHAT SESSION 18 DID (verified; NOT committed)

- **New scripts** (ruff clean, compile clean; package tests still green):
  - `scripts/fit_imcq120_temperature.py` — 9-param protocol
    (`k_p, v_th0, v_sub, r_d, k_kneesharp, k_eta, v_dbsat, v_ovmax, k_ovsharp`),
    `k_lambda = 0` fixed, `r_s = 0` fixed (Kelvin source), 6-start multistart,
    prints signed residuals + a cap-activity check + the `g_fs` spec check +
    the Rds(on)-T_j forward prediction. Takes ~4 min → **exceeds the 120 s tool
    cap; run detached and poll the log.**
  - `scripts/fit_imcq120_gain_probe.py` — the `k_beta` mechanism as a standalone
    probe (own channel function + vectorized bisection for `r_d`), with all
    caveats documented in its docstring. This is the seed for PART 1.
- **Artifacts** (`out/`, gitignored): `IMCQ120_TEMPERATURE_REPORT.md` (rev 2 —
  read it first), `imcq120_fit_params.json`, `imcq120_gain_params.json`,
  `imcq120_iv_{25,175}C_fit.png`, `imcq120_gain_fit.png`,
  `imcq120_transfer_ztc.png`, `imcq120_ron_Tj.png`.
- **Datasheet downloaded** for cross-checking:
  `tooling/workspace/outputs/IMCQ120R004M2H_datasheet.pdf` (Infineon rev 1.10,
  2025-10-13).

## Durable findings (SiC temperature dependence)

- **Fitted `v_th0` is NOT the datasheet V_th.** 6.615 → 5.653 V, slope
  **−6.41 mV/K** ≈ datasheet −6.67 mV/K, but a **constant +2.4 V offset**
  (58.3 mA / v_gs = v_ds spec vs strong-inversion apparent threshold). Never
  substitute the datasheet value; the softplus bridges the offset.
- **`v_sub` is temperature-invariant** (0.8364 → 0.8390 V, +0.3 %). It is a
  *shape knob* here, not a thermal voltage (a real swing would rise ~50 %
  298→448 K). Keep it T-independent in a two-temperature model.
- **`k_p`'s TC is not identifiable** from independent per-T fits: +8.5 %
  (β = −0.20) with `r_d` free vs −65 % (β = +2.54, RMSE 17.6 %) with `r_d`
  pinned; corr(`k_p`,`r_d`) = +0.82…+0.88. And it **flips sign** under the gain
  (β = +0.44) → `k_p(T)` is model-form dependent. Pin `r_d(T)` externally first.
- **`r_d` carries the whole Ron(T) TC:** 1.893 → 7.251 mΩ (×3.83, γ = 3.30 for
  `r_d ∝ T^γ`). **Linear `r_d(T)` goes NEGATIVE below −26 °C** — use a power law
  `(T/298 K)^3.3` or exponential (≈0.90 %/K at 25 °C).
- **Two-temperature interpolation predicts Rds(on)-T_j to max |err| 0.95 %
  (RMSE 0.49 %) over 25…195 °C but −15.3 % at −46 °C** — TC extrapolation fails
  cold; needs a third temperature or a joint fit.
- **Charge cap is NOT sensed on SiC** (as designed — it is a silicon SJ knob):
  `v_ovmax` ran to 67–77 V vs a max sampled `v_ov` of 13.4 V, `k_ovsharp` to
  39.4 (at the 40 ceiling), `sigma_rel = inf`, Jacobian rank 7/7 → **6/9**,
  `dv_ov_eff/dv_ov = 1.00000` at every sampled bias. It also *cannot* help: the
  model is already **below** the data at high v_gs (−3.9 → −11.4 % residuals)
  and the cap only reduces current.
- **`k_lambda = 0` is the right call** (user-directed): it pinned at the
  2e-3 1/V search guard in every fit with σ_rel 35/64 (a bound, not a
  measurement). Cost: +0.16 pp / +0.06 pp RMSE, and the high-v_gs output slope
  is lost (only DIBL remains, worth ~0.6 % at `v_sub = 0.84 V`).
- **Independent spec check (never fitted):** `g_fs` = 122 S at I_D = 185.2 A,
  V_DS = 20 V (Table 4) → model **128.5 S (+5.4 %)**. ZTC crossing model 14.40 V
  vs datasheet figure ≈14.2 V / 1500 A.
- **Datasheet corrections (rev 1 of the report was wrong, rev 2 fixed it):** the
  ~1550–1600 A "plateau" is **plot-frame truncation at 1600 A**, not a pulsed
  current limit (v_gs = 20 exits the top at ≈5.5 V, 18 at ≈6.4 V, 15 at
  ≈11.7 V); the transfer characteristic is **convex to ≈4200 A at V_GS = 20 V**
  — no gate-charge saturation anywhere in this device. I_DM = 1433 A.
- **Fit-quality floor is partly artifacts:** the 25 °C v_gs = 7 curve has two
  identical i_d = 19.9638 A points at v_ds = 6.90 and 9.36 V (no C¹ model can
  fit that) plus 2 near-origin suspect points; removing those 4 gives 7.03 %
  and **barely moves the parameters** (k_p 50.9→49.5, v_th0 6.615→6.572,
  r_d 1.89→1.87 mΩ) — parameters are robust, the headline number is not.

---

# PART 3 — Queue (updated)

1. **Implement `v_ov_eff2`/`k_beta`** — PART 1. Highest priority.
2. **Digitize the IMCQ transfer characteristic** (V_DS = 20 V, 25/175/200 °C,
   datasheet p. 9) → **joint multi-temperature fit** (both IV sets + Rds(on)-T_j
   + V_th-T_j + transfer, with the TC laws as the only extra unknowns). Fixes
   the −15 % cold-end error and the high-v_gs overshoot. A third temperature
   would pin `r_d(T)` curvature.
3. **`v_gs`-dependent C_gd blend (channel charge, CAP_SPEC limitation 1)** —
   session 17's blocker ("IMCQ120 output curves NOT digitized") is **RESOLVED**;
   this item is now unblocked. Targets: plateau deficits SJ +61 %, SiC +9.7 %,
   trench +42.8/24.5/22.0 %.
4. **Full-device LTspice subckt** (channel B source + diode + r_d/r_s + the caps
   subckt) then double-pulse Eon/Eoff within ~10 %. If PART 1 lands, the B-source
   expression gains the `k_beta` factor — trivial to express.
5. **Kelvin-aware subckt:** `r_s` inside the current path but OUTSIDE the gate
   loop. Today `r_s = 0` leaves the source-side power path unmodeled, so its T
   dependence aliases into `r_d` (hence γ = 3.3 ≫ package-only). Expected source
   of the "fit looks a bit weird" caveat the user flagged.
6. **CSV header contract mismatch (small, do it early):** the IMCQ digitizer
   export uses `index,V_ds,I_d,V_gs`, which **breaks `load_csv()`** (contract is
   `vgs,vds,ids`). Both scripts carry a local `load_imcq()`. Decide: rename the
   CSVs to the contract, or extend `data.py` to accept the digitizer spelling.
7. Unchanged: third-quadrant channel (sign-symmetric sync-rect conduction),
   HV leakage/BV element, T_j thermal node, CI, diode 5 A re-digitization,
   ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF points).

# Notes / gotchas (carry forward)

- Tool name: `engineering__execute_python` (**not** `files__execute_python`).
  Hard 120 s wall per call — long fits must be launched with
  `subprocess.Popen(..., DETACHED_PROCESS)` writing a log, then poll with
  `tasklist /FI "PID eq <pid>"`.
- `engineering__execute_python` writes scratch files into the workspace; there is
  **no `/tmp`** on this box — write into `outputs/` or the project `out/`.
- `openpyxl` is available; `names.xlsx` emits a harmless "Data Validation
  extension is not supported" warning. Sheet1 header:
  `Old name / LTspice name | Description - FET Channel | Type | Suffix | New name`.
- `out/` is gitignored; the two new `scripts/` are **untracked** — user decides
  whether to commit.
- Test totals drift (8 Tk digitizer tests self-skip headless) — don't pin counts.
- pF·V = pC, not C. Plateau split at V_gp (v_gd = 0), not ΔV_gs.
  `PchipInterpolator` has no `derivatives=` kwarg in the installed scipy.
- LTspice cap-export gotchas (never use `table()` in a B source or in `Q=` —
  its symbolic derivative is garbage ~1e11 A; use the incremental ramp-sum;
  `.meas tran FIND` strictly inside staircase segments) — full record in
  `CAP_SPEC.md`.
- Scratch from session 18: `outputs/imcq_main_run.log`,
  `outputs/IMCQ120R004M2H_datasheet.pdf`, `pdf-pages/*.png` (datasheet figure
  crops used for the corrections in PART 2).
