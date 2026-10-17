# Handover — ImprovedVDMOS session 14: DIODE TEMPERATURE COEFFICIENTS (a_ikf, a_rdiode) — complete, record promoted, committed

Sessions 12 (body diode) and 13 (pending-list closure) are COMPLETE and
committed; their handover content was replaced by this file. Escape hatch for
the session-11 rename: tag `pre-rename-2026-10`.

## What session 14 did

Added the LTspice diode **linear temperature coefficients** for i_kf and
r_diode (user-flagged: a temperature-independent r_diode was bending the
multi-T curves — it was, TC-off misses LTspice's own TC'd diode by up to 25 %).

1. **Names (names.xlsx, user-updated):** new Prefix type `temp. coeff` ->
   prefix `a` (row 8; VLOOKUP ranges widened to $A$2:$B$8); row 44 `Tikf` ->
   **`a_ikf`**, row 45 `Trs1` -> **`a_rdiode`** (suffix `rdiode`; the earlier
   `rs` suffix was corrected by the user -- a_rdiode scales r_diode, NOT the
   channel r_s; docstring carries the near-collision warning).
2. **LTspice form PROBED, not assumed** (help lists Tikf/Trs1 but no formula):
   `.meas` probes on LTspice 26.0.1 (`Tikf=0.006 Trs1=0.008`, .temp
   -55/25/100, 0.6/0.9 V) pin the law to exactly
   `X(T) = X*(1 + TC*(T - Tnom))`, **Tnom-referenced** (a 27 degC reference
   misses 0.2-0.6 %; TC-off misses up to 25 %). Agreement ~1e-5 relative.
   Pinned by `test_diode.py::test_ltspice_tc_laws` (runs live, headless
   batch mode). Trs2 (quadratic) deliberately NOT modelled.
3. **Implementation:** params.py 2 new fields (default 0.0 = off =
   bit-identical to session 12; linear space; BOUNDS (-5e-3, +5e-2) /degC
   keep the factor positive across -55..150 degC, worst 0.375); diode.py
   `_tc_factor` helper + `_TC_FLOOR = 1e-6` clamp (keeps r/i_kf positive --
   solver bracket proof `dg <= -1` survives absurd extrapolation; tested at
   t_c = 1000 degC); i_kf(T) inside the log-space high-injection divisor;
   r_diode(T) per-point in the series-R solver. Channel untouched (user
   decision).
4. **Tests:** 7 new (TC identity-at-t_nom + sentinel bit-identity, solver vs
   brentq at 100 degC, direction tests x2, extrapolation guard, TC round trip
   rank 6/6 exact recovery, LTspice TC laws). **127 tests green** (1 Tk
   digitizer test self-skip), ruff clean.
5. **RECORD PROMOTED (user go-ahead):** `ST8L65_DIODE` = 6-param TC fit,
   **4.15 %** (per-T 2.26/3.41/5.90), rank 6/6, cond 7.12e3, nothing pinned:
   i_s=3.83498e-6, r_diode=1.49732 m, v_eg=1.08815, i_kf=15.6344,
   a_ikf=+0.00544705, a_rdiode=+0.00341287 (k_n 2.0, k_xti 3.0, i_sr 0 fixed).
   TC-off session-12 record (4.78 %) kept as the script's `tc_off_probe`.
   Drift check 2.3e-6 [OK]; pinned-test bound tightened 0.049 -> **0.043**
   (between 4.15 % record and 4.78 % TC-off), per-T 0.08 -> 0.07, TCs
   asserted interior + positive.

## Probe findings worth remembering

- TCs are the **best-sensed diode parameters after r_diode** (sigma_rel
  0.019/0.028) -- three temperatures pin two linear TCs cleanly. Physical
  values: a_rdiode = +0.34 %/degC (silicon-ish series R), a_ikf = +0.54
  %/degC (high-injection knee rises with T).
- free-k_n probe WITH TCs: lands at k_n = 1.31 (3.00 %, corr(i_s,k_n)
  +0.996, reverse floor drifts to 3.9 nA -- TC-off it was 2.07/4.66 %/0.1-0.4
  mA). The valley's direction moves with the TCs; the unconstrained floor
  does not. k_n stays fixed at 2.0 from physics.
- i_kf pinned off now 6.63 % (was 6.48 % TC-off); i_sr/k_nr still unsensed.
- 150 degC remains the worst curve (5.90 %) -- high-current digitization
  noise, not a missing mechanism at this point.

## Queue (unchanged)

- **Third-quadrant channel**: `sign(v_ds)*i_d(|v_ds|, v_gs)` for
  synchronous-rectifier support; channel stays clamped 0 for v_ds <= 0.
- **HV leakage / BV/Ibv**: reverse floor is the Shockley i_s(T)+i_sr(T)
  ~ 3.8 uA @ -1 V/25 degC (datasheet ~1 uA @650 V) -- dedicated element
  still needed; the removed (0 V, 650 V, 1 uA) IV point waits on it.
- Q(V)/Qrr, `Tj` node (diode `t_c` is the placeholder interface), LTspice
  `.subckt` export (diode part is now LTspice-exact incl. TCs:
  `D(Is=… N=… Ikf=… Tikf=… Rs=… Trs1=… Tnom=… Xti=… Eg=…)`), CI.
- Possible future: Trs2 quadratic Rs TC if a device needs it (LTspice has it).

## Notes / gotchas

- TC=0 gives factor exactly 1.0 in IEEE -> bit-identical; TCs never touch
  the Xti/Eg saturation-current law (separate mechanism).
- i_kf "off" = sentinel 0.0, never the 1e-3 bound (active knee).
- `.vscode/settings.json` still untracked (editor-local).
- `out/` gitignored: st8l65_diode_{params,probe}.txt + fit pngs are on disk
  as evidence; the record itself lives in `vdmos_fit/reference.py`,
  drift-checked by `scripts/fit_st8l65_diode.py` (OFFICIAL_TOL 1e-3).
- Dev env: no venv in repo; pytest/ruff live in the assistant's
  `engineering313` env (Python 3.13.5).
