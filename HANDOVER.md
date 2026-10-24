# Handover — ImprovedVDMOS session 20: `k_beta_trap` gain LANDED + CSV header contract CLOSED
## NEXT SESSION = digitize the IMCQ transfer characteristic + joint multi-temperature fit

Session 19 landed the interface-trap-filling gain `v_ov_eff2`/`k_beta_trap`
(committed `8574ce2`); session 20 closed the CSV header contract (queue item 6,
**uncommitted at snapshot time**: `data.py`, both IMCQ scripts, `tests/test_data.py`,
README/data-README/STATUS/HANDOVER). Durable findings live in `STATUS.md`
(session-19 and session-20 sections), README ("Interface-trap-filling gain",
"Gate-overdrive charge cap", "Package-R solver"), `CAP_SPEC.md` + `reference.py`
(LTspice 26.0.1 cap gotchas — do not re-derive) and
`out/IMCQ120_TEMPERATURE_REPORT.md` (rev 2). This file replaces the session-18
handover; its PART 1 (the `k_beta` spec) was implemented as specified except for
two user-directed changes: the parameter is named **`k_beta_trap`** (names.xlsx
`beta_trap` row already existed; `k_beta` never entered the package) and
**`BOUNDS["k_beta_trap"] = (0.0, 2.0)` with the lower bound exactly 0.0** plus a
hard `ValueError` guard for negative values in `model.py` (clamp alternative
rejected — it would silently evaluate the non-monotone branch).

---

# PART 1 — THE TASK: transfer-curve digitization + joint multi-temperature fit

**Why (the two open fit-quality problems both trace to data the IV curves cannot see):**

1. **High-v_gs overshoot.** With the gain, the model predicts **3650 A vs ~2300 A**
   at v_gs = 18, V_DS = 20 V (datasheet transfer figure, p. 9). The output curves
   cannot see that region — they exit the 1600 A plot frame (plot truncation, NOT
   a current limit; the transfer characteristic is convex to ≈4200 A at V_GS = 20 V,
   no gate-charge saturation anywhere in this device). Bending the transfer over
   is what `k_lambda` does, and `k_lambda = 0` was pinned in every per-T fit.
   **Do not quote anything above ~1500 A until this is fixed.**
2. **Cold-end TC failure.** Two-temperature interpolation predicts Rds(on)-T_j to
   max |err| 0.95 % over 25…195 °C but **−15.3 % at −46 °C** — linear TC
   extrapolation dies cold (`r_d(T)` linear goes NEGATIVE below −26 °C).

**Task:**

1. Digitize the **transfer characteristic** (i_d vs V_GS at V_DS = 20 V;
   25/175/200 °C curves, datasheet p. 9) with `curve_digitizer`. The export
   header `index,V_ds,I_d,V_gs` now loads natively through `load_csv()`
   (session 20) — no loader workaround needed. Datasheet PDF:
   `tooling/workspace/outputs/IMCQ120R004M2H_datasheet.pdf` (Infineon rev 1.10);
   figure crops from session 18 in `outputs/pdf-pages/`.
2. Build a **joint multi-temperature fit**: both IV sets + the transfer curves +
   `data/IMCQ120R004M2H_Rdson_Tj_185.2A.csv` + `data/IMCQ120R004M2H_Vth_Tj.csv`,
   with **TC laws as the only extra unknowns** (per-parameter T-dependence
   instead of independent per-T fits):
   - `r_d(T)`: power law `(T/298 K)^γ` (γ ≈ 3.3 measured) — NOT linear (goes
     negative below −26 °C).
   - `v_th0(T)`: linear, slope ≈ −6.4…−8.6 mV/K (model-form dependent; the
     fitted value carries a constant +2.4 V offset vs the datasheet V_th spec —
     never substitute the datasheet number, the softplus bridges the offset).
   - `v_sub`: keep T-independent (measured +0.3 % over 25→175 °C; it is a shape
     knob here, not a thermal voltage).
   - `k_beta_trap(T)`: 0.0570 → 0.0758 1/V (+33 %, a real TC) — linear or
     Arrhenius; the joint fit decides.
   - **Re-free `k_lambda`** (the transfer data is exactly what senses it).
3. **Pin `r_d(T)` externally (from the Rds(on)-T_j CSV) before trusting `k_p(T)`**:
   the `k_p` TC is not identifiable from independent per-T fits (corr(k_p,r_d)
   +0.82…+0.88; its sign flips with the gain, β = +0.44 vs −0.20) — it is
   model-form dependent, not measurable.
4. A **third temperature** (the 200 °C transfer curve helps; a third IV set would
   be better) pins `r_d(T)` curvature — the cold-end −15 % is a TC-curvature
   problem, not a 25/175 °C problem.

**Acceptance checks (independent specs, never fitted):** `g_fs` = 122 S at
I_D = 185.2 A / V_DS = 20 V (model currently 128.5–128.8 S, +5.4 %); ZTC crossing
≈14.2 V / 1500 A (model 14.40 V); Ron chords 3.70/4.70 mΩ @25 °C (model 3.701/4.729);
V_th spec 58.3 mA (model 57.3/58.2 mA); Rds(on)-T_j over the full −46…195 °C span.
The transfer figure itself becomes a fitted input — after the joint fit, the
high-v_gs corner must come down from 3650 A toward ~2300 A at v_gs = 18.

**Watch out:** the 25 °C v_gs = 7 curve has two identical i_d = 19.9638 A points
(v_ds = 6.90 / 9.36 V — no C¹ model can fit that) plus 2 suspect near-origin
points; removing those 4 moves RMSE 8.53→7.03 % but barely moves parameters —
parameters are robust, the headline number is not. Same doctrine for the joint
fit: quote curves and spec checks, not a single RMSE.

---

# PART 2 — WHAT SESSIONS 19–20 DID (details in STATUS.md; do not re-derive)

**Session 19 (`8574ce2`) — `v_ov_eff2`/`k_beta_trap` in the package model:**

- `v_ov_eff2 = v_ov_eff·(1 + k_beta_trap·v_ov_eff)` applied **after** the cap;
  `v_knee` and the i_d prefactor both use `v_ov_eff2`; `k_theta` multiplies
  `v_ov_eff2` (decision: knee set by the final channel charge). Default 0.0 =
  exact identity; linear space; NOT in `DEFAULT_FREE`/`LOG_PARAMS`.
- **Bracket-invariant amendment:** the package-R proof needs the partials
  NON-NEGATIVE, not ≤ 1. The gain lifts `dv_ov_eff2/dv_ov_eff = 1 + 2·k_beta_trap·v_ov_eff`
  to 2.36 (25 °C) / 3.00 (175 °C) at v_gs = 18 — breaks the ≤ 1 half, keeps
  monotonicity; r_s > 0 + gain tested (brentq parity + stall-guard analogue).
- **Composition gotcha (silicon SJ):** composite asymptote is
  `v_ovmax·(1 + k_beta_trap·v_ovmax)`, NOT `v_ovmax`; synthetic cap+gain round
  trip gives corr(k_beta_trap, v_ovmax) = −0.85. The official ST8L65 record
  keeps k_beta_trap out of the free set and stays bit-identical; if ever freed
  on SJ parts, re-quote the record and report the correlation — never mix records.
- **Anchors (real model, cap off, k_lambda = 0, r_s = 0, 8 free incl. k_beta_trap):**
  RMSE 4.61 / 3.30 % (25/175 °C); k_beta_trap 0.0570 → 0.0758 1/V; k_p 25.74;
  v_th0 6.060; v_sub 0.7644; r_d 1.731 → 7.032 mΩ; Ron chords 3.701/4.729 mΩ;
  Vth spec 57.3/58.2 mA. The motivating diagnostic (implied channel overdrive
  slope rising past 1: 0.48 → 1.29 across v_gs 7→15) is reproduced by synthetic
  sets with k_beta_trap > 0 (`tests/test_overdrive_gain.py`).
- **Environment gotcha found during verification:** `scripts/fit_st8l65_nd.py`'s
  official-record drift metric is **scipy-version sensitive** — the fit lands at
  different points of the flat σ=inf `k_lambda`/`v_dbsat` valley (v_dbsat 339.02
  under scipy 1.18.1 vs 1287.03 under 1.16.0; k_lambda 5.13e-10 vs 1.63e-9) and
  the naive max-relative-drift test flags STALE. A/B-verified (`outputs/_s19_drift_ab*.log`):
  every SENSED parameter identical to ~1e-5 either way. Quote the curve, never
  the valley point; don't panic at a STALE flag — A/B the environment first.

**Session 20 (uncommitted) — CSV header contract:**

- `load_csv()` accepts BOTH spellings (`vgs,vds,ids` legacy / `index,V_ds,I_d,V_gs`
  digitizer), per-column, case-insensitive, order-free, extra columns ignored,
  mixable. `save_csv()` still writes the legacy spelling. The `load_imcq()`
  workarounds are deleted from both IMCQ scripts; the gain probe re-run
  reproduces the session-18 anchors through the new loader
  (`outputs/_s20_gainprobe_loadcsv.log`). Pinned by `tests/test_data.py` (8 tests).

**Durable SiC temperature findings (session 18, still valid):**

- Fitted `v_th0` ≠ datasheet V_th: slope −6.41 mV/K ≈ datasheet −6.67 mV/K but a
  constant +2.4 V offset (58.3 mA spec point vs strong-inversion apparent
  threshold).
- `v_sub` temperature-invariant (shape knob, not thermal voltage).
- `r_d` carries the whole Ron(T) TC: 1.893 → 7.251 mΩ (×3.83, γ = 3.30);
  linear `r_d(T)` goes NEGATIVE below −26 °C.
- Charge cap NOT sensed on SiC (it is the silicon SJ knob): v_ovmax → 67–77 V
  vs max sampled v_ov 13.4 V, σ = inf, `dv_ov_eff/dv_ov = 1.00000` everywhere;
  it also cannot help (model already BELOW data at high v_gs; cap only reduces).
- `k_lambda = 0` was the right call for per-T fits (pinned at the 2e-3 search
  guard with σ_rel 35/64 — a bound, not a measurement); cost +0.16/+0.06 pp RMSE
  and the lost high-v_gs output slope — the joint fit (PART 1) re-frees it.
- Datasheet corrections (report rev 2): the ~1600 A "plateau" is plot-frame
  truncation; I_DM = 1433 A; transfer convex to ≈4200 A at V_GS = 20 V.

---

# PART 3 — Queue (updated)

1. **Digitize the IMCQ transfer characteristic + joint multi-temperature fit** —
   PART 1 above. Highest priority.
2. **`v_gs`-dependent C_gd blend (channel charge, CAP_SPEC limitation 1)** —
   unblocked since session 18. Targets: plateau deficits SJ +61 %, SiC +9.7 %,
   trench +42.8/24.5/22.0 %.
3. **Full-device LTspice subckt** (channel B source + diode + r_d/r_s + the caps
   subckt) then double-pulse Eon/Eoff within ~10 %. The B-source expression now
   needs the `k_beta_trap` factor — trivial to express.
4. **Kelvin-aware subckt:** `r_s` inside the current path but OUTSIDE the gate
   loop. Today `r_s = 0` leaves the source-side power path unmodeled, so its T
   dependence aliases into `r_d` (hence γ = 3.3 ≫ package-only). Expected source
   of the "unphysical r_d and r_s" caveat the user flagged.
5. Unchanged: third-quadrant channel (sign-symmetric sync-rect conduction),
   HV leakage/BV element, T_j thermal node, CI, diode 5 A re-digitization,
   ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF points).

Closed: item 1 (`k_beta_trap`, session 19), item 6 (CSV header contract,
session 20) — see STATUS.md sections.

# Notes / gotchas (carry forward)

- Tool name: `engineering__execute_python` (**not** `files__execute_python`).
  Hard 120 s wall per call — long fits must be launched with
  `subprocess.Popen(..., DETACHED_PROCESS)` writing a log, then poll with
  `tasklist /FI "PID eq <pid>"`. Keep the polling sleep itself under the cap.
- Two interpreters on this box: the session interpreter
  (`tooling/engineering313`, scipy 1.18.1) and system Python
  (`AppData\Local\Programs\Python\Python313`, scipy 1.16.0). Fit results can
  differ in unsensed directions (see the drift-metric gotcha in PART 2);
  there is **no `.venv`** in the project — scripts add `src/` to `sys.path`
  themselves.
- `engineering__execute_python` writes scratch files into the workspace; there is
  **no `/tmp`** on this box — write into `outputs/` or the project `out/`.
- `openpyxl` is available; `names.xlsx` emits a harmless "Data Validation
  extension is not supported" warning. Sheet1 header:
  `Old name / LTspice name | Description - FET Channel | Type | Suffix | New name`.
- `out/` is gitignored. Test totals drift (Tk digitizer tests self-skip headless,
  1–2 skips) — don't pin counts.
- pF·V = pC, not C. Plateau split at V_gp (v_gd = 0), not ΔV_gs.
  `PchipInterpolator` has no `derivatives=` kwarg in the installed scipy.
- LTspice cap-export gotchas (never use `table()` in a B source or in `Q=` —
  its symbolic derivative is garbage ~1e11 A; use the incremental ramp-sum;
  `.meas tran FIND` strictly inside staircase segments) — full record in
  `CAP_SPEC.md`.
- Scratch (`tooling/workspace/outputs/`): `verify_gain_anchors.py` +
  `_s19_gain_anchor.log` (session-19 anchor verification), `_s19_drift_ab.py` +
  `_s19_drift_ab*.log` (scipy A/B), `_s20_gainprobe_loadcsv.log` (session-20
  loader smoke test), `imcq_main_run.log`, `IMCQ120R004M2H_datasheet.pdf`,
  `pdf-pages/*.png` (datasheet figure crops).
