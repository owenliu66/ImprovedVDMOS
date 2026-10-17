# Handover — ImprovedVDMOS session 15: Trs2 QUADRATIC r_diode TC (b_rdiode) — implemented, LTspice-verified, PROMOTED to the record

Sessions 12-14 (body diode, pending-list closure, TC promotion) are COMPLETE
and committed; their handover content was replaced by this file. Escape hatch
for the session-11 rename: tag `pre-rename-2026-10`.

## What session 15 did

User flag: "the improvement isn't as large as expected, especially the 150 C
branch" -> add LTspice Trs2 (quadratic Rs TC). Names per names.xlsx row 46 +
Prefix row 9 (`temp. coeff^2` -> prefix `b`): **`b_rdiode`** (the user's
message said "b_rdrain" but the sheet resolves Trs2 -> b_rdiode, suffix
matching a_rdiode/r_diode; flagged to user).

1. **LTspice law PROBED** (help lists Trs2 but no formula):
   `Rs(T) = Rs*(1 + Trs1*(T-Tnom) + Trs2*(T-Tnom)^2)`, Tnom-referenced --
   Tikf=0.006 Trs1=0.008 Trs2=2e-5 card, .temp -55/25/100/150, 0.6/0.9 V:
   agreement ~1e-5 relative at all 8 points. Linear-only misses up to 10.6 %
   at 150 C; 27 C-referenced misses 0.1-0.8 %. Pinned by
   `test_diode.py::test_ltspice_tc_laws` (extended: card now has Trs2, sweep
   now includes 150 C).
2. **Implementation:** params.py `b_rdiode` (default 0.0 = pure Trs1 law =
   bit-identical to session 14; linear space; BOUNDS (-1e-4, +1e-3)/degC^2);
   diode.py `_tc_factor(a, t_k, t_nom_k, b=0.0)` -- quadratic polynomial,
   same `_TC_FLOOR` clamp (bracket proof dg <= -1 survives absurd
   extrapolation; tested at t_c=1000 C with b=-1e-4). Channel untouched.
3. **Tests:** +1 new (b_rdiode curvature direction + |dT|-growth of the
   effect); TC identity/solver/guard/round-trip tests extended to b_rdiode
   (round trip now 7 free, exact recovery; unit-column SVD rank 7/7 --
   identifiability's raw report says 6/7, see gotchas).
4. **Probe + PROMOTION (user override, 2026-10-17):**
   4.148 -> **3.906 %** (per-T 2.37/3.16/5.49 vs 2.26/3.41/5.90),
   b_rdiode = -1.784e-5 /degC^2, a_rdiode +0.48 %/degC, a_ikf +0.83 %/degC,
   i_kf 21.97, r_diode 1.896 m. **PROMOTED to ST8L65_DIODE** (7-param free
   set) despite the caveats: the five-knob r_diode/i_kf/a_ikf/a_rdiode/
   b_rdiode block is strongly correlated (|corr| up to 0.985, sigma(i_kf)
   5.9 -> 15.4) -- quote the CURVE, never the split. (identifiability prints
   rank 6/7 -- a tolerance ARTIFACT: b_rdiode's column norm ~1.9e4 inflates
   the relative dead-column cutoff; unit-column SVD is full rank 7/7,
   smallest singular value 3.6e-2 of max. The round-trip test asserts
   rank >= 6 with a comment.) Promotion bookkeeping: script MULTI_START now
   carries b_rdiode starts (four starts, all 3.91 %, drift 3.8e-5 [OK]);
   trs2_probe flipped to trs2_off_probe (b_rdiode pinned 0, reproduces the
   session-14 4.15 % baseline); pinned-test bound 0.043 -> 0.040 (between
   3.91 % record and 4.15 % baseline) + b_rdiode interior assertion.
5. **Root cause of the 150 C residual (the real finding):** it is the 5 A
   point (+13..15 %). Its digitized v_sd (0.562 V) sits ~12 mV ABOVE what any
   smooth fit wants; at 5 A the r_diode drop is only ~10 mV of 562 mV, so NO
   Rs TC (Trs1 or Trs2) can move it. Excluding that one point: 150 C rms
   3.85 % (record) / 3.76 % (Trs2 fit). The 5 A point at 25 C is off the
   OTHER way (-6 mV, -9 %): the leftmost digitized point of each curve is
   the data's weak spot. The 150 C curve itself is smooth (dV/dlnI rises
   monotonically 69 -> 256 mV/decade -- no digitization kink).
6. **Recombination branch re-probed** (the physical candidate for
   low-current/high-T bending): 3.50 % overall but sigma(i_sr) = 645
   (unsensed), rank 7/9, 5 A/150 C point stays +11.7 % -- keep i_sr off.

## Queue (unchanged)

- **Third-quadrant channel**: `sign(v_ds)*i_d(|v_ds|, v_gs)` for
  synchronous-rectifier support; channel stays clamped 0 for v_ds <= 0.
- **HV leakage / BV/Ibv**: reverse floor is the Shockley i_s(T)+i_sr(T)
  ~ 3.85 uA @ -1 V/25 C (datasheet ~1 uA @650 V) -- dedicated element still
  needed; the removed (0 V, 650 V, 1 uA) IV point waits on it.
- Q(V)/Qrr, `Tj` node (diode `t_c` is the placeholder interface), LTspice
  `.subckt` export (diode part is LTspice-exact incl. all TCs:
  `D(Is=… N=… Ikf=… Tikf=… Rs=… Trs1=… Trs2=… Tnom=… Xti=… Eg=…)`), CI.
- Data priority: re-digitize (or measure) the 5 A points at 25/150 C --
  they, not missing physics, set the diode RMSE floor (~3.1-3.2 % excl.
  those points).

## Notes / gotchas

- TC=0 gives factor exactly 1.0 in IEEE -> bit-identical; TCs never touch
  the Xti/Eg saturation-current law (separate mechanism).
- b_rdiode sigma_rel = 0.0003 in the 7-param record fit is NOT real
  information (pinv along the near-degenerate correlated direction) -- the
  whole TC/knee block must be quoted as a curve, never as individual splits.
  identifiability's rank on mixed-unit linear-space
  columns is tolerance-sensitive: b_rdiode's dT^2 scale (~1e4 column norm)
  inflates the relative dead cutoff. Unit-column SVD is the honest check
  (7/7 here). A future session could scale columns inside identifiability --
  would change the channel fit's quoted rank 8/10, so do it deliberately.
- i_kf "off" = sentinel 0.0, never the 1e-3 bound (active knee).
- `.vscode/settings.json` still untracked (editor-local).
- `out/` gitignored: st8l65_diode_{params,probe}.txt + fit pngs on disk as
  evidence; the record lives in `vdmos_fit/reference.py`, drift-checked by
  `scripts/fit_st8l65_diode.py` (OFFICIAL_TOL 1e-3).
- Dev env: no venv in repo; pytest/ruff live in the assistant's
  `engineering313` env (Python 3.13.5).
