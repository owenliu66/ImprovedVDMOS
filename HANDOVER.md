# Handover — ImprovedVDMOS session 13: BODY DIODE COMPLETE + committed; pending list closed

Session 12 (body diode) is COMPLETE and committed to `main`; its handover
content was replaced by this file. Escape hatch for the session-11 rename:
tag `pre-rename-2026-10`.

## What session 13 did (closing session 12's PENDING list)

1. **tests/test_diode.py fixes:**
   - Fixes (a)–(d) from the session-12 list (forward-floor sign, top-level
     `from vdmos_fit import identifiability`, `.meas` parser regex
     `^(m\d+):\s*I\(V1\)\s*=\s*(\S+)` with `m03/m06/m09` tags, `naive_diode`
     ±700 exponent clip) were **already applied in the working tree** when
     session 12 ended — the handover predated the last edits. Verified, not
     re-applied.
   - Fix (e) APPLIED here: `test_st8l65_diode_record_reproduces` RMSE bound
     tightened 0.050 → **0.049** (record is 0.0478).
   - Fix (f) VERIFIED: LTspice batch mode `LTspice.exe -b -w` works headless
     (rc=0, `xcheck.log` written, all three `.meas` parsed). Exe found at
     `C:\Program Files\ADI\LTspice\LTspice.exe` (not on PATH — the
     `_ltspice_exe()` candidate list covers it). The cross-validation test
     RUNS, skipif never fires on this machine.
2. **Full suite:** `pytest` → **120 passed, 0 failed, 0 skipped** (LTspice
   cross-validation ran; the 8 Tk digitizer window tests also ran on this
   desktop — totals still not pinned). `test_diode.py` alone: 17 passed.
   `ruff check` → clean.
3. **Docs:** README "Body diode" section, Package-R bullet ("the body diode
   carries it"), Scope update, STATUS.md session-12 snapshot with diode out
   of Missing, data/README.md `ST8L65N044M9_diode.csv` entry — all were
   **already written** by session 12's end (the handover's "Docs NOT started"
   was stale). Content checked against the handover spec: complete.
   `names.xlsx`: the `v_sd` row EXISTS at **row 43** (not row 28 as the
   handover anticipated — blank row 28 now separates channel rows 2–27 from
   diode rows 29–43; position is cosmetic, every row's E-column formula is
   row-relative and resolves). Verified all 15 diode-block formulas against
   `params.py`: i_di, v_therm, v_di, i_diffuse, i_recomb, i_s, k_n, i_sr,
   k_nr, i_kf, r_diode, t_nom, k_xti, v_eg, v_sd. Prefix sheet covers
   `temperature` (row 7, inside VLOOKUP range $A$2:$B$7) so t_nom resolves.
4. **Committed to main** (session-12 work + this closure in one commit).

## Evidence numbers (session-13 verification)

- LTspice batch probe (Is=1e-6 N=1.5 Ikf=5 Rs=1e-3 Xti=3 Eg=1.11 Tnom=25,
  .temp 25, .dc V1 0→1): m03 = −2.40095541812e-3, m06 = −3.66624522209,
  m09 = −86.0688171387 A — parsed cleanly by the test regex.
- Official diode record reproduces under the tightened 0.049 bound; per-T
  rms −55 °C 3.40 %, 25 °C 3.85 %, 150 °C 6.49 % (unchanged).

## Queue (unchanged, user's decision 3)

- **Third-quadrant channel**: `sign(v_ds)·i_d(|v_ds|, v_gs)` for
  synchronous-rectifier support; channel stays clamped 0 for v_ds ≤ 0.
- **HV leakage / BV/Ibv**: reverse-blocking floor is currently the Shockley
  i_s(T)+i_sr(T) ≈ 3.7 µA @ −1 V/25 °C (datasheet ~1 µA @650 V) — dedicated
  leakage element still needed; the removed (0 V, 650 V, 1 µA) IV point
  waits on it.
- Q(V)/Qrr, `Tj` node (diode `t_c` is the placeholder interface), LTspice
  `.subckt` export, CI.

## Notes / gotchas

- `.vscode/settings.json` left UNTRACKED (editor-local, not part of the
  session). Add to .gitignore or commit deliberately in a future session.
- `out/` is gitignored: `st8l65_diode_{params,probe}.txt` and
  `st8l65_diode_fit.png` are on disk as evidence but not in git; the record
  itself lives in `vdmos_fit/reference.py` and is drift-checked by
  `scripts/fit_st8l65_diode.py` (`OFFICIAL_TOL`, 3.65e-6 [OK] at session-12
  close) and pinned by `tests/test_diode.py`.
- i_kf "off" = sentinel 0.0, never the 1e-3 bound (active knee).
- Fit protocol lessons stand: log-space packing mandatory; r_diode device
  bound (0, 20 mΩ); k_n fixed 2.0 (i_s↔k_n corr +1.000 degenerate);
  r_diode↔i_kf corr +0.948; i_s weak prefactor (σ_rel 0.88) — quote the
  (k_n, i_s) pair.
- Dev env: no venv in repo; pytest/ruff live in the assistant's
  `engineering313` env (Python 3.13.5).
