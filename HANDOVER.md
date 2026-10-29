# Handover — ImprovedVDMOS session 21: PART 1 CLOSED + PROMOTED (official IMCQ record)

## NEXT SESSION = Kelvin-aware subckt / r_s split (now has quantitative targets)

Session 21 closed HANDOVER PART 1 (queue item 1): the IMCQ120R004M2H transfer
characteristic (V_DS = 20 V, 25/175/200 °C, datasheet p. 9) was digitized
(`data/IMCQ120R004M2H_Transfer.csv`) and a **joint multi-temperature fit** landed as
`scripts/fit_imcq120_joint.py` — one residual vector over IV 25 °C + IV 175 °C +
transfer + Rds(on)-Tj + Vth-Tj (200 pts, 29 temperatures) with the TC LAWS as the only
extra unknowns. It was then **promoted to the OFFICIAL IMCQ record** (user go-ahead,
same session): `reference.py::IMCQ120_JOINT` (an `ImcqJointRecord` — 25 °C base +
TC-law coefficients, `params_at(t_c)` builds the Params), pinned by
`tests/test_imcq120_official.py` (11 tests, record evaluation against all four CSVs,
no re-fit; the two known open gaps are pinned as tolerances). The session-18
`out/IMCQ120_TEMPERATURE_REPORT.md` was **deleted** (user direction — its durable
findings are carried in `out/IMCQ120_JOINT_FIT_REPORT.md` §8, which is now the IMCQ
document); `scripts/fit_imcq120_temperature.py` / `fit_imcq120_gain_probe.py` stay as
historical baselines with supersession notes.
**Uncommitted at snapshot time:** the joint script, the `reference.py` record +
`__init__` exports, `tests/test_imcq120_official.py`, data-README entry, README,
the STATUS.md session-21 section, and the report deletion (out/ is gitignored). This file replaces the
session-20 handover; its PART 1 was implemented as specified except for three
user-directed changes, all now on record:

1. **`r_d(T)` is QUADRATIC in T(°C)** (Trs1/Trs2 law form), NOT the power law /
   exponential of the temperature report — the user's read of the digitized curve
   was right by a factor of 7: standalone chord fits give 0.067 % / 0.065 % RMSE
   over −46…196 °C (power law 9.5/9.1 %, exponential 5.9/6.3 %). The quadratic
   curvature IS the cold-end flattening; joint-fit Rds(on) check max error 1.08 %
   (was −15.3 % at −46 °C). Never linear, never power law on this part.
2. **`k_theta` freed** — lands at 0.109 1/V and is the transfer-bending mechanism
   (v_knee → 1/k_theta = 9.2 V makes i_sat ~ linear in v_ov_eff2 at high gate
   drive). Pinning it off doubles the overall RMSE (1.93 → 4.14 %).
3. **`k_lambda` re-freed — and the data says NO**: 4.4e−10, σ = inf, ablation
   identical. The PART 1 hypothesis "bending the transfer is what k_lambda does"
   is DISCONFIRMED. Knee saturation + package-r_d divider + the transfer data do it.

Headline: overall RMSE **1.93 %** (IV 1.73 / transfer 3.39 / Rds(on) 0.65 / Vth 2.06);
the session-18 8.53/6.13 % floor is gone. High-v_gs corner v_gs = 18 @ 20 V:
**2118 A** at 175 °C (digitized ~2070; gain fit predicted 3650) — the "don't quote
above 1500 A" gate is LIFTED. ZTC 25/175 at 14.92 V / 1703 A (matches the digitized
output-curve crossing ~14.9 V; the datasheet-figure ~14.2 V read stays +5 % off —
the figure read is itself soft). Ron chords 3.722/4.652 mΩ @25 °C. Vth spec 57.0/59.1 mA.

TC laws (25 °C base, all in `out/imcq120_joint_params.json`): k_p25 9.318 A/V²,
beta_p 0.378 (k_p ~ T^−β); v_th0 5.746 V, slope −7.16 mV/K (σ 0.5 %; datasheet
−6.67; fitted offset ≈ +1.5 V — never substitute the datasheet V_th); v_sub 0.7933
invariant; k_beta_trap 0.4687 1/V with +3.40e−4 1/V/K LINEAR law (Arrhenius ties at
1.92 vs 1.93 %, Ea σ = inf → linear kept, the joint fit decided); r_d 2.3985 mΩ with
a = 6.727e−3 /K, b = 5.367e−5 /K² (vertex ≈ −38 °C — valid over the data span only);
k_kneesharp 1.465; k_eta pinned at the 0.05 cap (σ inf), v_dbsat 320 (σ inf),
k_lambda ≈ 0 — DIBL/CLM block inert (bounds, not measurements).

**Identifiability doctrine CONFIRMED (old PART 1.3):** the `no_transfer` ablation
(IV-only, same TC laws) fits IV *better* (1.59 %) but its r_d(T) TC is degenerate
(a_rd σ_rel 295, corr(a_rd,b_rd) −1.00) and its OUT-OF-SAMPLE Rds(on)-Tj fails at
20.2 % RMSE / 34.1 % cold. The chord data pins r_d(T); only then is beta_p real.
Composition gotchas on record: k_p25~k_bt25 −0.95 (the pair is a model-form split —
the composite k_p·v_ov_eff2 is the record; quote curves, never the pair),
beta_p~s_bkt +0.97 (channel TC shared). Gain derivative dv_ov_eff2/dv_ov_eff =
1 + 2·k_bt·v_ov reaches **12.5 at v_gs = 18** — legal for the package-R bracket only
because r_s = 0 (monotonicity, not ≤ 1, is what the proof needs); re-check before
any r_s > 0 use of this record.

**The two open quantitative gaps BOTH point at the Kelvin source split (queue 2):**
g_fs 129.5 S vs 122 spec (+6.2 %), and the 25 °C transfer corner −5…−7 % at
v_gs 15–19 (r_d25 = 2.40 mΩ drops 7.2 V at 3 kA — the source-side power path
between Kelvin point and channel is unmodeled, today fully aliased into r_d).

---

# Queue (updated)

1. **`v_gs`-dependent C_gd blend (channel charge, CAP_SPEC limitation 1)** —
   targets unchanged: plateau deficits SJ +61 %, SiC +9.7 %, trench +42.8/24.5/22.0 %.
2. **Kelvin-aware subckt: `r_s` inside the current path but OUTSIDE the gate loop**
   — PROMOTED: the session-21 joint fit quantified what it owes (g_fs +6.2 %,
   25 °C corner −7 %, γ = 3.3 aliasing explanation). Today r_s = 0 leaves the
   source-side power path unmodeled. Needs the model to support asymmetric r_s
   (gate loop vs current loop) — a model/terminal question, decide the design first.
3. **Full-device LTspice subckt** (channel B source + diode + r_d/r_s + the caps
   subckt) then double-pulse Eon/Eoff within ~10 %. The B-source expression needs
   the `k_beta_trap` factor — trivial to express. Depends on the Kelvin decision (2).
4. Closed (session 21, user go-ahead): **IMCQ record bookkeeping** — the joint fit
   IS the official record now: `reference.py::IMCQ120_JOINT` +
   `tests/test_imcq120_official.py` + drift check in the script; the
   k_p25~k_bt25 −0.95 composition caveat is carried in the record docstring and
   README ("Interface-trap-filling gain").
5. Unchanged: third-quadrant channel (sign-symmetric sync-rect conduction),
   HV leakage/BV element, T_j thermal node, CI, diode 5 A re-digitization,
   ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF points).

Closed: item "digitize transfer + joint fit" (session 21 — this file's PART 1, now
the official `IMCQ120_JOINT` record), item 4 record bookkeeping (session 21),
`k_beta_trap` (session 19), CSV header contract (session 20) — see STATUS.md.

# Notes / gotchas (carry forward)

- Tool name: `engineering__execute_python` (**not** `files__execute_python`).
  Hard 120 s wall per call — long fits must be launched with
  `subprocess.Popen(..., DETACHED_PROCESS)` writing a log, then poll with
  `tasklist /FI "PID eq <pid>"`. Keep the polling sleep itself under the cap
  (sleep(60) with the default 60 s timeout FAILS — use sleep ≤ 45 s or pass
  `timeout_seconds`). Use `python -u` for a streaming log (session 21 did).
- Two interpreters on this box: the session interpreter
  (`tooling/engineering313`, scipy 1.18.1 — what session 21 used) and system
  Python (`AppData\Local\Programs\Python\Python313`, scipy 1.16.0). Fit results
  can differ in unsensed directions; there is **no `.venv`** — scripts add
  `src/` to `sys.path` themselves.
- `engineering__execute_python` writes scratch into the workspace; **no `/tmp`**
  — write into `outputs/` or the project `out/`.
- The joint script's bounds table is SCRIPT-LOCAL (`FREE_SPEC` in
  `fit_imcq120_joint.py`), mirroring `params.vector_bounds` semantics (log
  params mapped through log — a raw-bounds bug cost one wasted run; the pinned-check
  compares physical values, not vector coords). r_d25 upper bound = datasheet
  25 °C chord (3.7 mΩ) as a search guard against the degenerate "channel ≈ 0"
  branch; k_eta cap 0.05 as in the per-T protocol.
- The OFFICIAL record's TC-law evaluation lives in `reference.py::ImcqJointRecord
  .params_at` (package side) AND in the script's dict-based `params_at` (variant
  side) — deliberately two implementations; they are cross-checked every script
  run (record-RMSE re-eval OK/STALE) and by the pinned test. Keep them in sync.
  The drift check compares SENSED parameters only (session-19 lesson).
- Transfer CSV: header `index,V_gs,I_d,condition` — `load_csv()` does NOT read it
  (no V_ds column); the script-local loader parses `condition` (Qg-style). The
  175/200 °C curves share the crossing point (12.025 V, 999.1 A) once per curve;
  no data below 12 V for them (digitization limit, kept not cleaned).
- `openpyxl` available; `names.xlsx` emits a harmless Data Validation warning.
- `out/` is gitignored. Test totals drift (Tk digitizer tests self-skip headless)
  — don't pin counts. Package `src/` was NOT touched by session 21 — the joint
  fit is script-level (TC laws compose Params per temperature).
- pF·V = pC, not C. Plateau split at V_gp (v_gd = 0), not ΔV_gs.
  `PchipInterpolator` has no `derivatives=` kwarg in the installed scipy.
- LTspice cap-export gotchas (never `table()` in a B source or in `Q=` — garbage
  symbolic derivative ~1e11 A; use the incremental ramp-sum; `.meas tran FIND`
  strictly inside staircase segments) — full record in `CAP_SPEC.md`.
- Scratch (`tooling/workspace/outputs/`): `verify_gain_anchors.py` +
  `_s19_gain_anchor.log`, `_s19_drift_ab*.log` (scipy A/B),
  `_s20_gainprobe_loadcsv.log`, `imcq_main_run.log`,
  `IMCQ120R004M2H_datasheet.pdf`, `pdf-pages/*.png` (figure crops).
- Project `out/` additions (session 21): `IMCQ120_JOINT_FIT_REPORT.md`,
  `imcq120_joint_params.json`, `imcq120_joint_{transfer,iv_25C,iv_175C,ron_Tj}.png`,
  `imcq120_rdson_law_compare.png`, `_s21_joint_run.log`.
