# Handover — ImprovedVDMOS session 17: capacitance/charge model COMPLETE
(caps.py + c_gdmax extraction + records + LTspice nonlinear-cap export +
parity + names.xlsx + docs — committed)

Session 16 (design settled, data digitized, docs + data.py plumbing) and
session 17 (implementation + extraction + export + parity) are COMPLETE and
committed. Durable findings live in STATUS.md + CAP_SPEC.md + reference.py.
This file replaces the session-16 handover.

## What session 17 did (all verified, 145 tests green, ruff clean)

1. **Confirmed session 16's partial code**: the `load_qg_csv` header bug was
   already fixed (header keys normalized by removing spaces); `data.py`
   loaders/savers tested (CSV round-trip tests).
2. **`caps.py`** (CAP_SPEC Phase 2): `_BranchTable` (PCHIP on ln C, GL-5
   quadrature, clamp-C policy), `origin_anchor` (first-4-point quadratic
   convention, ANCHOR_N=4 pinned), `CapParams` (c_gd0, s_gd0, c_gdmax;
   k_agd derived; flat-branch limit k_agd=inf), `CapModel` (c_gd/q_gd
   composite with origin-quadratic head below the first rss node — exact
   splice slope — and PCHIP above; c_gs/c_ds as union-node difference
   tables). dq/dv ≡ C identity pinned by central-FD tests.
3. **`extract_cmax` + `scripts/fit_cmax.py`** (Phase 3, corrected split):
   plateau = max-dQ/dV segment; plateau charge is a PURE table-branch
   validation (it ENDS at v_gd = 0); c_gdmax from the post-plateau segment
   via bracketed brentq on q_exp(V_end − i_d·r_ds_on; c_gdmax) = remainder,
   bracket [c_gd0, 2·c_gs(0)], LOWER/UPPER clamp diagnostics.
   **RESULTS (records `CAP_*` in reference.py, drift-checked by
   test_caps.py):** all three parts land at c_gdmax ≈ c_gd0 — the exp
   branch is nearly FLAT, the seed c_gs(0) over-predicts 4–8×.
   ST8L65 2750 pF (≈c_gd0), IMCQ120 3823 pF (LOWER clamp), ISC040 663 pF
   (median of 647/798/663). Plateau deficits (limitation-1 quantifier):
   SJ +61.2 %, SiC +9.7 %, trench +42.8/+24.5/+22.0 %.
4. **`ltspice.py` + `scripts/export_caps.py`** (Phase 4): the delta exports
   as `.subckt <name> G D S` of three **charge-master nonlinear capacitors**
   (`Cxxx n1 n2 Q=<expr>`, `x` = device voltage). Q= is the **incremental
   ramp-sum** of the sampled charge table:
   `q0 + m1·(x−x0) + Σ (m_{k+1}−m_k)·max(0, x−x_k)` — LTspice's symbolic
   dQ/dx is exactly the piecewise-constant C staircase (N_SUB=12 samples
   per digitised interval), and the end-segment slopes carry the clamp-C
   policy for free. Exports at the EXTRACTED records, not the seed.
   **Live parity on all three parts** (PWL-ramp drive, 3 delta configs,
   14 `.meas tran FIND` points each): measured vs exported staircase
   ≤0.3 % (most ≤0.05 %); staircase vs PCHIP ≤0.4 % incl. the SJ cliff.
   `test_caps.py::test_ltspice_caps_parity` runs it in the suite (skipped
   only when LTspice is absent).
5. **names.xlsx**: Prefix rows capacitance→`c`, charge→`q`; Sheet1 rows
   c_gd/c_gs/c_ds, c_gdmax, c_gd0, k_agd, q_gd/q_gs/q_ds (VLOOKUP ranges
   widened to $A$2:$B$11). README naming table + layout tree updated
   (session 16's README claim had NOT been saved — it was stale).
6. **Docs**: CAP_SPEC §Code TODO now records Phase 4 landed + the LTspice
   gotchas below; STATUS snapshot → session 17; data/README Qg entries
   (session 16 had added them — verified).

## LTspice 26.0.1 gotchas (measured, on the record in CAP_SPEC)

- `table()` in a B source: rejected ("Function not supported in this
  context").
- `table()` inside `Q=`: compiles, but its **symbolic derivative is
  garbage** (~1e11 A measured). Never use.
- Absolute-slope ramp sums are WRONG (they sum segment slopes instead of
  stepping them — showed up as exact n× parity factors). Use the
  incremental form.
- `.MODEL x C (C=...)` on the instance: "Unknown parameter".
- Parity measurement: `.meas tran FIND` must land strictly INSIDE
  staircase segments (on a knot you measure the solver's blend), and
  maxstep must keep the Gear-2 companion's q-history window ≪ the
  segment width (5 ns at 1e7 V/s ramps worked).
- Sign convention verified: `Cxxx n1 n2 Q=q(x)` sinks dq/dt from n1;
  driving the delta with Vxx gives I(Vxx) = −dq/dt.

## Queue (updated)

- **Full-device LTspice subckt** — channel B source (vdmos_id is all
  expressible: softplus/tanh/min/max) + diode `D(...)` + r_d/r_s + this
  caps subckt; then A/B switching sim vs Level-1. This is the natural
  next session (Phase 5 of the validation plan: double-pulse Eon/Eoff
  within ~10 %; turn-on time arbitrates the flat-branch finding).
- **v_gs-dependent C_gd blend (channel charge, limitation 1)** — the
  plateau deficits above are its quantified targets (SJ +61 % worst).
  Needs the SiC DC IV fit first (IMCQ120 output curves NOT digitized —
  the repo has only ST8L65 IV).
- Clean the ST8L65 `iss` gridline-snap artifact (four identical
  4869.28716 pF points) on the next re-digitization pass.
- Unchanged from session 15: third-quadrant channel, HV leakage/BV,
  Tj node, CI, diode 5 A re-digitization.

## Notes / gotchas (carry forward)

- Tool name: `engineering__execute_python` (NOT files__execute_python).
- `PchipInterpolator` in the installed scipy has no `derivatives=` kwarg —
  directed endpoint slopes come from the quadratic splice, not PCHIP.
- pF·V = pC, not C. The plateau split is at V_gp (v_gd = 0), not ΔV_gs.
- ln vs log10 interpolation: identical (affine rescale).
- Dev env: no venv in repo; pytest/ruff live in the assistant's
  `engineering313` env (Python 3.13.5). openpyxl available.
- Test totals drift: 8 digitizer window tests drive Tk and self-skip
  headless (0–2 skipped per run) — don't pin the count.
- `out/` gitignored: fit_cmax_params.txt, fit_cmax_cgd.png,
  <part>_caps.subckt, parity_*.cir live there. `.vscode/` still untracked.
- Scratch figures from sessions 16/17 in workspace `outputs/`
  (cv_pchip_tables.png, cv_gd_splice_fixed.png, cv_splice_check.png).
