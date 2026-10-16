# Handover — ImprovedVDMOS (naming-convention rename, COMPLETE)

## What happened

Repo-wide variable rename to the new naming convention (`<type prefix>_<suffix>`,
Greek letters get the `k_` prefix). The scheme was reviewed and locked with the
user over two sessions; `names.xlsx` (committed with the rename) is the single
source of truth. Decisions locked:

- `m_triode` **stays** (LTspice model-card mapping; `vdmos_id_level1` is a
  foreign-model adapter).
- **CSV headers stay `vgs,vds,ids`** — external contract with the digitizer.
  `data.py` has `_CSV_COLUMNS`/`_CSV_HEADERS` mapping the new attribute names
  (`CurveSet.v_gs/v_ds/i_d`) to the legacy header spellings.
- Theta joins the unit exception list: `k_` factors are unitless EXCEPT
  `k_p` [A/V²] and the standard 1/V trio `k_lambda`, `k_eta`, `k_theta`.
- `gk` is a **conductance** [S], not transconductance (xlsx Type fixed).
- `db` = drain **bias**, not drain-to-body (no bulk terminal in this model) —
  documented in the README naming section.

## Final mapping (old → new)

| old | new | old | new |
|---|---|---|---|
| kp | `k_p` | lam | `k_lambda` |
| vth0 | `v_th0` | gk | `g_k` |
| ksub | `v_sub` | vthk | `v_thk` |
| kcap | `v_ovmax` | kk | `v_kovsoft` |
| m | `k_ovsharp` | vk | `v_dbsat` |
| theta | `k_theta` | vov_floor | `v_kneemin` |
| vsat | `v_vsat` | rd / rs | `r_d` / `r_s` |
| nd | `k_kneesharp` | gate (local) | `k_db` |
| eta | `k_eta` | vth / vov (locals) | `v_th` / `v_ov` (+ new `v_ov_eff`) |
| vds / vgs / Id | `v_ds` / `v_gs` / `i_d` | vknee / v (locals) | `v_knee` / `v_dseff` |
| kink_gate (local) | `k_kink` | id_floor (fit kwarg) | `i_d_floor` |

Sentinel semantics unchanged: `v_ovmax = 0` = cap off, `k_kneesharp = 0` =
legacy tanh, `v_kneemin` never fitted. `LOG_PARAMS`/`BOUNDS`/`DEFAULT_FREE`
keys renamed; values identical.

## All work items completed (A–J)

A. README model-equation block rewritten: `k_db`, `v_ov_eff`, `v_dseff`,
   `^k_ovsharp`, `softmin(v_ds, v_knee, k_kneesharp)` form.
B. README "Naming convention" section ADDED (prefix table from names.xlsx
   Prefix sheet, `k_` rule + unit exceptions, `db` = drain-bias caveat,
   sharp-vs-soft duality, `m_triode` carve-out, CSV-header contract).
C. Derivative notation fixed repo-wide: `di_d/dv_ds`, `di_d/dv_gs`,
   `dv_ov_eff/dv_ov`, `dv_dseff/dv_ds`, `dg/di_d` (README, plotting.py,
   test_model/test_gate_cap/test_drain_softness).
D. Greek letters referring to OUR model → `k_theta`/`k_eta`/`k_lambda`
   (README prose + table, STATUS fit tables). λ/Kp/Vov/Vto/mtriode in
   Level-1 SPICE context KEPT (item E).
E. Level-1 SPICE-name prose kept; README `Vov_eff` (our model) → `v_ov_eff`.
F. test_model.py kink local `gate` → `k_kink`.
G. README section renamed "Gate-cap numerics" → "Charge-cap numerics"
   (model.py docstring points there); all `v_ovmax/m` → `v_ovmax/k_ovsharp`
   in README/STATUS/HANDOVER/scripts.
H. Column widths bumped: fit.py `{name:<10s}` → `<12s`; scripts `{name:<8s}`
   → `<12s`; roundtrip_demo `<10s` → `<12s` (k_kneesharp is 11 chars).
I. Mask locals `m = np.isclose(...)` → `msk` in plotting.py and both fit
   scripts.
J. `~$names.xlsx` Excel lock file: already gone (nothing to delete).

Extra identifier cleanup beyond the word-boundary audit (old names embedded in
identifiers, no external references — verified by grep):
`test_nd_*` → `test_kneesharp_*`, `test_no_ndr_with_fitted_nd` →
`..._fitted_kneesharp`, `test_round_trip_recovers_nd` → `..._kneesharp`,
`NDTRUTH` → `KNEE_TRUTH`, `test_kcap_zero` → `test_cap_zero`,
`test_cap_monotone_with_kcap_asymptote` → `..._with_v_ovmax_asymptote`,
`test_subthreshold_slope_matches_ksub` → `..._matches_v_sub`,
`nd_bound_probe` → `kneesharp_bound_probe`, `lam_cap_probe` →
`klambda_cap_probe`, `theta_vsat_probe` → `ktheta_vsat_probe`,
`nd_s/nd0/p_nd` → `knee_s/ks0/p_ks`, plot labels `nd_A/nd_B` →
`knee_A/knee_B`, `Id_plateau` → `i_plateau`, `vov_` → `v_ov`, and the
"ksub is anchored…" output strings in both fit scripts → "v_sub is anchored…".

## Verification results (all run post-edit)

1. `pytest` → **101 passed, 2 skipped** (103 collected — identical per-file
   collection to tag `pre-rename-2026-10`; the 1–2 skips are the flaky Tk
   digitizer tests, STATUS says don't pin the count).
2. `ruff check` → clean (17 E501 line-length violations introduced by the
   rename were wrapped).
3. Fingerprint: all 10 arrays in `out/_rename_baseline.npz` reproduced
   **bit-identically** (`np.array_equal`). Recipes: grid
   `v_gs=linspace(-2,16,181)` × `v_ds=linspace(0,60,241)`;
   default_grid / official_grid (= ST8L65_OFFICIAL incl. its r_d/r_s) /
   level1_grid (default kwargs) / default_pkgR (DEFAULT + r_d=2e-3, r_s=1e-3) /
   official_pkgR (OFFICIAL + r_d=0.024, r_s=0.007) / cap_on (DEFAULT +
   v_ovmax=1.0, k_ovsharp=9.0) / nd_on (k_kneesharp=4.0) / nd_tanh (2.8) /
   official_resid (model − meas, amperes, on the CSV) / example (DEFAULT on the
   example_dataset grid). npz deleted after comparison.
4. Bonus A/B: a 12-probe battery (cap/knee/package-R/DIBL/CLM/subthreshold/
   kink/Level-1 variants) run against a `pre-rename-2026-10` worktree →
   **all bit-identical**.
5. `python scripts/roundtrip_demo.py` → recovery unchanged (0.86 % rel RMSE).
6. `python scripts/fit_st8l65_nd.py` → fit A 3.79 %, fit B 1.52 %,
   official-record drift **1.79e-06 [OK]** (no `[STALE]`).
7. Word-boundary grep for every old name in `src/vdmos_fit`, `tests`,
   `scripts` → only allowed hits: `m_triode`/`mtriode`/`Kp`/`Vov`/`Vto`/
   `Lambda` (Level-1 SPICE names), CSV header strings in `data.py`,
   "python -m", `m` = milli in an output format, mtriode loop variable in
   test_model. `src/curve_digitizer` untouched (its `m` is a transform slope).

## Notes / gotchas

- `vdmos_id_level1()` kwargs are `k_p, v_th0, v_sub, m_triode, k_lambda, r_d,
  r_s` — all call sites updated; any notebook or external snippet using
  `kp=`/`lam=` breaks. Hard break accepted: no alias layer; tag
  `pre-rename-2026-10` is the escape hatch.
- `reference.py` docstring carries a "(re-labelled 2026-11; values unchanged)"
  note — keep it until the next release.
- The regex pass renamed `id_` → `i_d` inside `_solve_package_r` and
  `ids` → `i_d` in scripts/tests; audited, no collisions with `i_d_floor`.
- `plotting.py` `_vov` helper and `test_model.py` `v_ov()` helper take
  `v_th0`/`v_sub` kwargs now.
- Dev environment note: the repo has no venv; pytest/ruff were installed into
  the assistant's `engineering313` tool env for the verification runs.

## Standing project state (unchanged by the rename)

- OFFICIAL reference: fit A of `scripts/fit_st8l65_nd.py` (10-param,
  `k_kneesharp` free), 3.79 % rel RMSE, rank 8/10, cond 2.19e3; record in
  `vdmos_fit.reference.ST8L65_OFFICIAL`, pinned by `tests/test_st8l65_official.py`.
  `k_lambda`/`v_dbsat` σ = inf is expected (data stops at Vds = 20 V).
- Queue after the rename: (1) knee-cluster identifiability via dense knee-region
  data; (2) body diode / reverse polarity + HV leakage — biggest model hole,
  package-R path is v_ds ≥ 0 only; (3) architecture: Q(V), Tj node, LTspice
  export, CI.
- If the official set is ever exported to LTspice: pin `k_lambda = 0` and any
  `v_dbsat` rather than shipping the fitted (unsensed) numbers.
- `.vscode/` remains untracked (user's call). `out/` is gitignored.
