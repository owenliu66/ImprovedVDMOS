# Handover — ImprovedVDMOS session 24: RECOVERY DATASHEET EXTRACTION PROTOCOL + t_minor KNOB

## NEXT SESSION = IMPROVE THE RECOVERY FIT FURTHER

Session 24 built the datasheet-driven extraction pipeline
(`src/vdmos_fit/extract_recov.py`, `extract()` two-knob + `extract3k()`
with the user's t_minor knob), verified it on ISC040N10NM8, and mapped the
two-knob limitation quantitatively. Durable record:
`out/RECOV_EXTRACTION_PROTOCOL.md` (bench design, gotchas, t_minor scan +
fit tables), `out/RECOV_3k_waveforms.png` / `RECOV_3k_hi_waveforms.png`,
`out/recov_extract3k_{run,eq,hi}/extract3k_report.txt`. Session-23 material
unchanged in STATUS.md.

**Key session-24 findings (do not re-derive):**
- Fast single-event bench (seeded `.ic V(X1:q_j)={i_test*t_minor/q_0}`,
  stiff `Vv N001 0 V_R` + D5 clamp + 1p across D5, C2 Ki A = Cds
  charge-master from `data/ISC040N10NM8_cds_expr.txt`, r_max=0, gear,
  dTmax=t_rr_target/10) reproduces the 500-period bench to <0.5 %;
  full extraction ~15 s.
- `.meas` windows terminate at fractions of MEASURED I_RM (chained
  `i_pk`), not of I_F — the old −frac·I_F targets fail at low di/dt.
- t_minor knob (`ClampD Roff=tau*t_minor/q_0`): stored charge =
  I_F·t_minor/tau (q_0-free), identity with two-knob at t_minor=q_0.
  Qrr_lo/Qrr_hi ratio has a FLOOR ≈ 0.30 (t_minor 15–35 ns); datasheet
  needs 0.146. Single storage state cannot decouple the two di/dt points.
- Best fits: hi-priority 3-knob (35.3n, 3.51, 34.1n) → Qrr_hi −0.5 %;
  equal-weight (46.1n, 9.84, 24.4n) → Qrr_lo +28 %, t_lo +3.6 % but
  Qrr_hi −42 %. Two-knob C (34.8n/5) remains the hard-switching record.

**Next session — investigate how to further improve the fit.** Ideas
(unevaluated, pick from data):
1. Second independent storage state (dual lifetime: fast/slow
   recombination channel) to break the 0.30 Qrr-ratio floor — the only
   structural route to both di/dt points; estimate identifiability first.
2. di/dt-dependent q_0 or tau law (q_0·(1+c·di/dt)) — empirical, check
   against the 3-point datasheet if available.
3. Revisit the D1 DC branch / Is,N fit (V_SD typ 0.82 vs 0.91 @ 25 A) —
   session 23 noted q_0/tau insensitive, but the plateau voltage split
   could shift the snap point slightly.
4. Check whether the datasheet t_rr definition is really t_rr25 (IEC) —
   the +7…+10 % residual might be a definition mismatch, not a model error.
5. If a second part (ST8L65/IMCQ) carries recovery specs, run the
   protocol there as a cross-part check of the model family.
6. Then: `reference.py` `*_RECOV` record + pinned test (same doctrine as
   channel/diode/caps records); ltspice.py Cgs/Cgd-only export flag; the
   ISC040 channel record (queue item 3 unchanged).

**Gotchas added this session** (full list in protocol §7): ASCII raw parse
(index line + one value per line; column k = variable k); single-value
`.step` rejected; `.step` nesting last=outer; failed `.meas` rows print
`failed`; dTmax mandatory; `V(X1:node)` unquoted in `.ic`/`.meas`.

---

# Handover — ImprovedVDMOS session 23: RECOVERY NETWORK INTEGRATED + FITTED (ISC040)

Session 23 integrated the user's behavioral body-diode recovery network
(`DiodeRecov`) into the full-device LTspice subckt and fitted its two knobs
on ISC040N10NM8 — **LTspice-side only, NO package Python touched**. Durable
record: `out/ISC040N10NM8_recov.subckt` (template, Cds ramp-sum inline) +
`out/RECOV_REPORT.md` (semantics, sweep tables, integration decisions, A/Bs);
session-23 section in STATUS.md. Session-22 material is unchanged there
(session-22 section — verdicts on record, do not re-litigate).

**What landed (all simulation-verified, LTspice 26.0.1, `method=gear
reltol=1e-4`):**
- Fit: **q_0 = 35 n, tau = 5** → Qrr 249.6 nC vs datasheet typ 247 (+1.1 %),
  t_rr(25 % I_RM) 27.6 ns vs 25 (+10 %) @ V_R = 50 V, I_F = 12.5 A,
  1000 A/µs. q_0 sets timing (Qrr/t_rr ≈ linear), tau sets shape
  (plateau→tail departure, snap hardness). Variants: 34n/3 datasheet-strict,
  37n/5 waveform-strict.
- Model semantics decoded: ClampD `Roff=tau` closes the storage loop —
  V(q_j) = I_F volts at steady state (tau-independent), q_j time constant =
  Roff·C1 = q_0 (35 ns; the 1 Ω normalization ties the "charge at 1 A" and
  "carrier lifetime" readings), stored charge = q_0·I_F/tau ∝ I_F. The
  −I_RM plateau is carried by B2, not D1 (I(D1) ~ −1 nA); plateau ends at
  V(q_j) ≈ I_RM·tau/V_R (predicted 1.25 V, observed 1.33 V).
- Cds placement: caps delta exports **Cgs/Cgd only**; the measured C_ds
  charge-master lives INSIDE the recovery branch at **Ki–A** (electrically
  D–S through B3+R1/R2, fc ≈ 3 GHz). D–S placement crashes at hard forward
  re-turn-on (timestep 1.25e-19; D1 has no Cjo, B3+R1 isolate the branch);
  Ki–Ai is singular (user's finding, reproduced).
- Half-bridge stability (user fix, now template default): `R2 = 1 µΩ` + B3
  `+1µ` floor — breaks near-singular B3/R1/R2/Ron bridge loops; waveform
  impact < 0.1 % (A/B over r2 × floor, `DiodeTest_fixcheck.net`).

**Known two-knob limit (quantified, do not chase):** datasheet Qrr *rises*
with di/dt (36 nC @100 A/µs → 247 nC @1000 A/µs); no single stored-charge
state follows both — at 100 A/µs the fit gives Qrr 77 nC (≈2× high). The fit
is tuned to the hard-switching condition; the 100 A/µs numbers are "defined
by design" in the datasheet.

**Extraction protocol (next session) — where to start:**
- Datasheet inputs available on modern parts: t_rr + Qrr (rarely I_RM) at
  ≥2 di/dt points (ISC040: 100 and 1000 A/µs), V_SD@I_F/T_j for the D1 DC
  branch (already the `diode.py`/`ST8L65_DIODE` record path), C_oss/Q_oss
  for the Cds branch (already the caps record).
- Order of operations: fit the D1 DC branch FIRST (existing protocol), then
  pin Qrr at the HIGH-di/dt point with q_0 (one knob, ≈linear), then shape
  (t_rr / I_RM snap step) with tau, then CHECK the low-di/dt point and
  REPORT the residual as the model-family limitation. Extra knobs
  (di/dt-dependent tau or q_0, second storage node, soft B2 factor) only if
  the residual demands them — decide from data, not in advance.
- Bench: `tooling/workspace/outputs/vdmos_recov/DiodeTest_sweep.net`
  (settled 500-period; `.meas TD` RELATIVE to the shifted axis — gotcha
  below). k_chg = 5, r_max = 0.1 were held fixed at the user's values —
  probe whether they are identifiable at all (likely weak: k_chg only moves
  the forward-recovery tail, r_max only the turn-on overshoot).
- T_j: datasheets give recovery at 25 °C only; carrier lifetime has a strong
  TC (rule of thumb ≈ +50 %/100 °C on Qrr). Decide whether to add a q_0(T)
  law or stay isothermal — do not invent a TC without data.
- Once the protocol is settled, consider a `reference.py` `*_RECOV` record
  + a pinned test, same doctrine as the channel/diode/caps records.

# Queue (updated)

1. **`v_gs`-dependent C_gd blend (channel charge, CAP_SPEC limitation 1)** —
   targets unchanged: plateau deficits SJ +61 %, SiC +9.7 %, trench +42.8/24.5/22.0 %.
2. CLOSED (session 22): **Kelvin-aware r_s split** — `Params.r_sk` landed;
   the joint-vector refit shows the two open IMCQ gaps are NOT series-R
   topology problems (verdicts in STATUS session-22 section).
3. **Full-device LTspice subckt** (channel B source + diode + r_d/r_s/r_sk
   wiring + caps subckt, 4 terminals D G S K) then double-pulse Eon/Eoff
   within ~10 %. **PARTLY LANDED (session 23):** the diode-recovery half is
   done and fitted (`out/ISC040N10NM8_recov.subckt`). Still needed: the
   channel B-source expression from `model.py` (needs the `k_beta_trap`
   factor — trivial to express, README model block); a `ltspice.py` flag to
   export **Cgs/Cgd only** (Cds now lives in the recovery branch — never
   export both, double counting); the ISC040 channel record (no official
   channel fit exists for that part yet). Kelvin wiring (session 22, still
   the plan): `D G S K`, `r_d` drain leg, channel B-source between internal
   D'/S', `r_s` S'→S, `r_sk` S'→K, gate caps referenced G→internal, driver
   return at K.
4. **NEW: recovery-parameter extraction from datasheet info** — see NEXT
   SESSION above; protocol design + first application (ISC040 re-extraction
   as the round-trip check, then ST8L65/IMCQ if their datasheets carry
   recovery specs).
5. Unchanged: third-quadrant channel (sign-symmetric sync-rect conduction),
   HV leakage/BV element, T_j thermal node, CI, diode 5 A re-digitization,
   ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF points).

Closed: item 2 Kelvin split (session 22), "digitize transfer + joint fit"
(session 21, official `IMCQ120_JOINT`), record bookkeeping (session 21),
`k_beta_trap` (session 19), CSV header contract (session 20) — see STATUS.md.

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
- Recovery-net gotchas (session 23): `.tran` Tstart shifts the raw/meas time
  axis to 0 — `.meas TD` must be RELATIVE (the original bench's TD=0 worked
  by accident); `Roff=tau` unbraced works in 26.0.1 but write `{tau}` for
  portability; **gear required** (trapezoidal rings badly on the snap); B3
  needs the +1 µΩ floor + R2 = 1 µΩ in any bridge topology; Ki–A cap
  placement is mandatory (D–S crashes re-turn-on, Ki–Ai singular); ideal-
  switch benches show a ~1.2 kA Cds dump impulse at midpoint swings — real
  loop inductance (~1–5 nH) domesitcs it; LTspice `.raw` "real" flag =
  float32 records (24 var × 4 B + 2 B flag + 2 B pad in the stepped file —
  parse carefully if hand-parsing).
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
  snapshot time along with the session-22 code). The recovery parameters
  (q_0/tau/k_chg/r_max) are deliberately NOT in names.xlsx — LTspice-side
  netlist values, outside the package naming convention.
- pF·V = pC, not C. Plateau split at V_gp (v_gd = 0), not ΔV_gs.
  `PchipInterpolator` has no `derivatives=` kwarg in the installed scipy.
- LTspice cap-export gotchas (never `table()` in a B source or in `Q=` —
  garbage symbolic derivative ~1e11 A; use the incremental ramp-sum;
  `.meas tran FIND` strictly inside staircase segments) — full record in
  `CAP_SPEC.md`.
- Scratch (`tooling/workspace/outputs/vdmos_recov/`): all recovery benches —
  `DiodeTest_base/500p/sweep/fine/integ/integ2/didt100/fixcheck.net`,
  `HalfBridge_conv.net`, `cds_expr.txt`. Project `out/` additions
  (session 23): `ISC040N10NM8_recov.subckt`, `RECOV_REPORT.md`. `out/` is
  gitignored.
- README slimmed to a pure model reference (session 23, LaTeX display math).
  Three implementation details that lived ONLY in README prose are parked
  here (all pinned by tests — the prose is gone, the behavior is not):
  (1) CLI: `--show`/`--out/-o` are registered on every subcommand with
  `argparse.SUPPRESS` defaults so they parse before *or* after the
  subcommand — without SUPPRESS the subparser resets them unconditionally
  (`test_cli.py`). (2) `make_synthetic` noise is multiplicative and applied
  WITHOUT clipping the values — flooring the current distorts subthreshold
  up to 2×; near-zero currents are handled in the residual denominator
  (`i_d_floor` in `residuals`), not in the data. (3) `k_p` is NOT the
  Level-1 Kp — it must be re-extracted when k_theta/v_vsat are added
  (28 → ~34.7 at DEFAULT); `v_kneemin` is a numerical floor, keep fixed,
  never fitted; `v_ovmax` lower bound (0.05 V) sits deliberately above 0 —
  between 0 and 0.05 V the cap is a hard clip at nearly zero volts.
- Test totals drift (Tk digitizer tests self-skip headless) — don't pin
  counts. Session 23 touched NO package code (docs + `out/` only), so all
  model anchors are untouched by construction.
