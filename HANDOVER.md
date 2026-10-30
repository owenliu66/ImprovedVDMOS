# Handover — ImprovedVDMOS session 22: r_s SPLIT LANDED + queue 2 CLOSED BY EXPERIMENT

## NEXT SESSION = full-device LTspice subckt (queue 3) — the Kelvin decision is MADE

Session 22 implemented the Kelvin-aware package-R split the user asked for
(`r_sk` row added to names.xlsx by the user; semantics decided and landed by
the assistant, pinned by `tests/test_kelvin_split.py`):

```
v_ds' = v_ds − (r_d + r_s)·i_d      current loop: ALL of i_d
v_gs' = v_gs − (r_s − r_sk)·i_d     gate loop: returns at the KELVIN tap
```

`r_sk` = resistance between the kelvin tap and the power source pin (the
downstream segment of `r_s`, carrying only gate current — zero in the DC
model). `r_sk = 0`: no kelvin contact, legacy full degeneration
(bit-identical). `r_sk = r_s`: ideal tap at the die, gate loop clean.
Guard `0 ≤ r_sk ≤ r_s` (ValueError). Solver signature generalized:
`_solve_package_r(..., r_g=...)` — gate-loop coefficient, default `r_s`
(legacy). The bracket proof was RE-CHECKED with `r_g > 0` AND the
trap-filling gain active (gate derivative ~12.5, the HANDOVER-21 warning):
monotonicity survives, the proof only needs non-negative partials. DONE.

**The experiment** (`scripts/fit_imcq120_kelvin.py`, same 200-pt joint
vector/TC laws/acceptance as the official record; `scripts/
fit_imcq120_kelvin_probe.py` for forced-`r_g` probes; report
`out/IMCQ120_KELVIN_SPLIT_REPORT.md`; log `out/_s22_kelvin_run.log`;
params `out/imcq120_kelvin_params.json`; plot
`out/imcq120_kelvin_transfer.png`):

| variant | landed | RMSE % | g_fs S | corner25 A |
|---|---|---|---|---|
| control (record repro) | drift 4e-6 OK | 1.93 | 129.5 | 2990 |
| kelvin_ideal (r_sk = r_s, gate CLEAN) | r_s25 σ_rel 75, corr(r_d25,r_s25) −1.00, r_tot25 = 2.3985 mΩ invariant | 1.93 (identical) | 129.5 | 2990 |
| degenerated (r_sk = 0, degeneration = r_s) | **r_s25 → 0 pinned LOWER** from 3 starts (0.2/0.5/0.9 mΩ + stronger-channel) | 1.93 | 129.5 | 2990 |
| split_free (r_s25, f_sk = r_sk/r_s free) | r_s25 → 0 pinned; f_sk σ inf | 1.93 | 129.5 | 2990 |

Forced probe: r_g = 0.25/0.5/1.0 mΩ → RMSE 2.11/2.28/2.71 %, corner25
residual −7.8/−8.7/−10.0 %, g_fs 130.3/130.5/131.4 S (it goes UP, not down).

**VERDICTS (on record — do not re-litigate):**
1. The current-loop split is DC-invisible: only `r_d + r_s` (+ total TC) is
   sensed; the r_d/r_s/r_sk split rides a ±1.00-correlated degenerate valley.
2. Gate-loop degeneration ("the datasheet static measurement did NOT use
   the kelvin pin") is REJECTED by the data — boundary optimum r_g = 0 from
   every start; the 25 °C corner needs MORE channel drive and degeneration
   removes exactly that.
3. g_fs +6.2 % is a spec-vs-digitized-slope softness: the transfer curve
   itself pins the terminal slope at 185 A to ~130 S; NO series-R topology
   can move it without breaking the curve fit (the probe re-tunes the
   channel and lands back at ~130 S).
4. The 25 °C corner −6 % is a channel-shape compromise (k_p25~k_bt25~k_theta
   squeezed by the hot curves + chords), not a series-path gap. Candidate
   future mechanism (NOT fitted): the gain form itself — v(1 + k_bt·v) is
   quadratic/unbounded; the corner wants the 25 °C gain to bend less at
   high v_ov than the linear k_bt(T) law allows.
5. The official record `IMCQ120_JOINT` is UNCHANGED (r_s = r_sk = 0); the
   two gap tolerances in `tests/test_imcq120_official.py` stay pinned with
   updated wording. `ImcqJointRecord` now carries `r_sk` (fixed 0) in sync
   with the script's params_at.

**What the subckt (queue 3) must do with this:** the kelvin pin carries only
gate CURRENT — r_sk is an AC element. The 4-terminal subckt: `D G S K`,
`r_d` drain leg, channel B-source between internal D'/S', `r_s` S'→S,
`r_sk` S'→K, gate caps (caps subckt) referenced G→internal, driver return
at K. Then double-pulse Eon/Eoff within ~10 %. The B-source expression
needs the `k_beta_trap` factor — trivial to express (README model block).

# Queue (updated)

1. **`v_gs`-dependent C_gd blend (channel charge, CAP_SPEC limitation 1)** —
   targets unchanged: plateau deficits SJ +61 %, SiC +9.7 %, trench +42.8/24.5/22.0 %.
2. CLOSED (session 22): **Kelvin-aware r_s split** — `Params.r_sk` landed
   (semantics above, tests pinned); the joint-vector refit shows the two
   open IMCQ gaps are NOT series-R topology problems (verdicts above).
3. **NEXT: full-device LTspice subckt** (channel B source + diode +
   r_d/r_s/r_sk wiring + the caps subckt, 4 terminals D G S K) then
   double-pulse Eon/Eoff within ~10 %. UNBLOCKED — the Kelvin decision is
   made (session 22). Depends on nothing else.
4. Unchanged: third-quadrant channel (sign-symmetric sync-rect conduction),
   HV leakage/BV element, T_j thermal node, CI, diode 5 A re-digitization,
   ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF points).

Closed: item 2 Kelvin split (session 22), "digitize transfer + joint fit"
(session 21, now the official `IMCQ120_JOINT` record), record bookkeeping
(session 21), `k_beta_trap` (session 19), CSV header contract (session 20)
— see STATUS.md.

# Notes / gotchas (carry forward)

- Tool name: `engineering__execute_python` (**not** `files__execute_python`).
  Hard 120 s wall per call — long fits must be launched with
  `subprocess.Popen(..., DETACHED_PROCESS)` writing a log, then poll with
  `tasklist /FI "PID eq <pid>"`. Keep the polling sleep itself under the cap
  (sleep(60) with the default 60 s timeout FAILS — use sleep ≤ 45 s or pass
  `timeout_seconds`). Use `python -u` for a streaming log.
- Two interpreters on this box: the session interpreter
  (`tooling/engineering313`, scipy 1.18.1 — what sessions 21–22 used, and
  what `engineering__python_environment` selects) and system Python
  (`AppData\Local\Programs\Python\Python313`, scipy 1.16.0). Fit results can
  differ in unsensed directions; there is **no `.venv`** — scripts add
  `src/` to `sys.path` themselves.
- `engineering__execute_python` writes scratch into the workspace; **no
  `/tmp`** — write into `outputs/` or the project `out/`.
- `fit_imcq120_kelvin.py` IMPORTS `fit_imcq120_joint` and PATCHES its
  module-global `params_at` (kelvin-aware) + extends `FREE_SPEC`/`VARIANTS`
  — importing the kelvin module is what activates the patch (the probe
  script does exactly that). `JF.acceptance` needs the base dict to carry
  `_rd_T/_rd_vg/_rd_chord` (from `JointData.rd_*`) — a KeyError otherwise.
- The joint script's printed `sigma_rel` is the GN marginal sd in VECTOR
  space (log params → relative; linear params → RAW units) and is INFLATED
  1–2 orders for linear knobs on this problem (session-22 calibration:
  `s_vth` printed "1σ" refit explodes RMSE 1.93 → 9.48 %; true σ ≈ 2–5 % of
  the slope; `r_d25` −0.5 mΩ fixed → 2.91 %). For linear-knob uncertainty
  quote FIX-AND-REFIT penalties, not the print. The joint report §2 carries
  the correction note; the record's "σ 0.5 %" for s_vth was a misread.
- The joint script's bounds table is SCRIPT-LOCAL (`FREE_SPEC`), mirroring
  `params.vector_bounds` semantics (log params mapped through log). r_d25
  upper bound = datasheet 25 °C chord (3.7 mΩ) as the degenerate-branch
  guard; the kelvin script caps r_s25 at the same chord and parameterizes
  `r_sk = f_sk·r_s25` with `f_sk ∈ [0,1]` (bounds cannot reference another
  parameter; the model guard 0 ≤ r_sk ≤ r_s is the hard line).
- The OFFICIAL record's TC-law evaluation lives in `reference.py::ImcqJointRecord
  .params_at` AND in the joint script's dict-based `params_at` — deliberately
  two implementations, cross-checked every script run (OK/STALE) and by the
  pinned test. Keep them in sync (session 22 added r_sk to both).
  The drift check compares SENSED parameters only (session-19 lesson).
- Transfer CSV: header `index,V_gs,I_d,condition` — `load_csv()` does NOT
  read it; the script-local loader parses `condition` (Qg-style). The
  175/200 °C curves share the crossing point (12.025 V, 999.1 A) once per
  curve; no data below 12 V for them (digitization limit, kept not cleaned).
- `openpyxl` available; `names.xlsx` emits a harmless Data Validation
  warning. names.xlsx now carries the r_sk row (user-added, uncommitted at
  snapshot time along with the session-22 code).
- pF·V = pC, not C. Plateau split at V_gp (v_gd = 0), not ΔV_gs.
  `PchipInterpolator` has no `derivatives=` kwarg in the installed scipy.
- LTspice cap-export gotchas (never `table()` in a B source or in `Q=` —
  garbage symbolic derivative ~1e11 A; use the incremental ramp-sum;
  `.meas tran FIND` strictly inside staircase segments) — full record in
  `CAP_SPEC.md`.
- Scratch (`tooling/workspace/outputs/`): `_s22_kelvin_smoke.py` (kelvin
  machinery smoke), plus session-19/20/21 items (verify_gain_anchors.py,
  _s19/_s20 logs, IMCQ120R004M2H_datasheet.pdf, pdf-pages/*.png).
- Project `out/` additions (session 22): `IMCQ120_KELVIN_SPLIT_REPORT.md`,
  `imcq120_kelvin_params.json`, `imcq120_kelvin_transfer.png`,
  `_s22_kelvin_run.log`. `out/` is gitignored.
- Test totals drift (Tk digitizer tests self-skip headless) — don't pin
  counts. Session 22 touched the PACKAGE (params.py/model.py/reference.py
  record comments) — the r_sk=0 fast path and all pre-split anchors are
  pinned bit-identical (test_kelvin_split.py + test_package_r.py).
