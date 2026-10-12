# Handover — ImprovedVDMOS: cap refactor + docs COMPLETE, committed; next: ksub/vk degeneracy + bounded drain element (session 6)

## 0x. SESSION 6 UPDATE (docs + commit DONE — read this first)

- **§1d docs DONE:** `README.md` (Vov_eff in the equation block, new
  "Gate-overdrive charge cap" subsection, stall-guard paragraph in the
  package-R section, test count, comparison-table row);
  `src/vdmos_fit/STATUS.md` (problem 6 closed, problem 9 updated, new
  "Gate-overdrive cap (kcap, m) — sessions 4+5" section with the fit A/B
  table and the five new open problems, verified state updated);
  `data/README.md` (new "Datasets" section: ST8L65 provenance — digitised
  25 °C pulsed, 0–20 V, spec anchors, kept Vgs=5 duplicate, removed 650 V
  off point, Figure 3 not on disk).
- **Committed** (see `git log`): everything except `.vscode/` (left
  untracked — user's call), including the root `STATUS.md` deletion.
- **Discrepancy vs session-5 §8/§9:** the working tree also contained
  `fit.py` + `test_fit.py` + `test_package_r.py` changes that the session-5
  handover did NOT list — the identifiability dead-column fix (relative
  1e-5-of-max tolerance; correlations/cond over live columns only) and its
  test, plus the Newton two-cycle stall test. All included in the commit.
  Consequence: the §3 "ksub↔vk = ±1.000, cond ~2e9" numbers were measured
  with the OLD tolerance; ksub's column is actually dead (7.8e-7 of max,
  pinning leaves the fit identical). Queue item 2 (fix ksub / tie vk, rerun,
  expect cond ~1e4–1e5) stands and now also re-measures identifiability.
- **Re-verified before commit:** 86 passed, 1 skipped (digitizer test needs
  a Tk display), ruff clean. (Session-5's "85 + 1" predates the two extra
  tests.)
- **Next-session queue = §7 items 2–6 below, unchanged** (item 1 done).

---

## 0. READ THIS FIRST — the tree is GREEN and complete except docs

Session 4's half-applied refactor is **finished**. `params.py`, `model.py`,
`scripts/fit_st8l65.py` and the new `tests/test_gate_cap.py` are all done and
verified. Both pending fits (§2 decision 5: with / without the Vgs = 5 curve)
have been **run** — results in §3 below, which *change* the open-problem list.

- **pytest: 85 passed, 1 skipped** (78 old + 7 new `test_gate_cap.py`).
- **`ruff check .` clean.**
- Everything is **still uncommitted** (§9). Commit is safe now — but see §1c:
  the model contains one change *beyond* session 4's spec (solver stall guard).

Still to do: **§1d docs only** (README + STATUS still describe `kappa`), then
commit. Everything else from session 4's §1 is done.

---

## 1. What was implemented this session (exactly)

### 1a. `model.py` — DONE, plus one addition beyond the session-4 spec

All three session-4 edits applied verbatim:
- `vdmos_id` no longer passes `kappa` to `_solve_package_r`.
- `_channel_current`: `vov = _soft_cap(vov, p.kcap, p.m)` inserted after the
  softplus, **in place**, before `vknee` (variant A — charge and knee saturate
  together).
- `_soft_cap` added (log-space form, `kcap <= 0` → identity, bit-identical).
- `_solve_package_r`: `kappa` kwarg, `kappa < 0` guard, bound (b), the residual
  split, the `dudI` branch and the pole comment all deleted. Docstring updated.
  The `moved` progress test kept.

**NEW (not in session 4's plan) — Newton two-cycle stall guard.** The
no-NDR test (Vds to 200 V, fitted params, rd = 26 mΩ) exposed a *real* solver
stall mode that `kappa` was not the only source of: at Vds ≈ 178 V the
bracketed Newton enters a **two-cycle** — `g(Id)` is nearly flat (slope ≈ −1.2)
at both bracket endpoints with a steep gate-degeneration cliff in between
(`rs·Id` dropping the internal gate across the cap knee), so Newton maps each
endpoint into the other's neighbourhood and the bracket shrinks only ~1e-5 per
iteration (fails at 60 iterations with residual ~854 A). Fix: a per-point
`stall` counter — if the bracket width fails to shrink by 10 % three
consecutive times, Newton is vetoed and bisection forced until progress
resumes. ~6 lines in `_solve_package_r`; does not change any converged result
(all 78 pre-existing tests still pass bit-identically). **This is the honest
replacement for the deleted pole-avoidance machinery: the pole is gone, but
flat-then-cliff `g(Id)` shapes survive the cap and need the stall guard.**

### 1b. `scripts/fit_st8l65.py` — DONE (full rewrite, zero `kappa` refs)

- Free set = session-4 §4: `kp, vth0, ksub, lam, gds, rd, rs, vk, kcap, m`.
- `DEVICE_BOUNDS = {"rd": (0, 0.036), "rs": (0, 0.036), "lam": (1e-6, 0.004)}`
  via `bounds_overrides` — global `BOUNDS` untouched.
- Fixed at start: `theta=0, vsat=1e4, eta=0, gk=0, vthk=1.6, kk=0.1`
  (measured dead, session-4 §3e). `kcap`/`m` freed explicitly per device.
- 4 multi-starts around the session-4 basin + different vk/kcap/m basins.
- The "fully free / widened bounds" run stays deleted (session-4 §6).
- Runs **fit A** (full 41 pts) and **fit B** (Vgs = 5 curve dropped, 36 pts)
  — session-4 §2 decision 5, now DONE.
- Outputs: `out/st8l65_capfit.png`, `out/st8l65_capfit_params.txt`.
  (Session-4's `out/gate_cap_vs_kappa.png` etc. are stale — they show the old
  unconstrained-lam fit; regenerate plots from the new script if needed.)

### 1c. `tests/test_gate_cap.py` — DONE (7 tests, all green)

Covers everything session-4 §1c asked: kcap=0 bit-identity (Params-level +
through the package-R solver), `0 <= dVov_eff/dVov <= 1` across
vov∈[1e-6,50]×kcap∈[0.1,20]×m∈[0.5,40], monotonicity + asymptote, subthreshold
slope = 1/ksub for any kcap>0, **no NDR** for Vds∈[0,200] at Vgs=10 with the
fitted params (this is what forced the stall guard, §1a), no inf/NaN at
kcap=0.05/m=40/vov=50 and vov=1e6, synthetic round-trip (kcap=2.0, m=8 →
recovered to 0.1 % / 0.2 %).

**Corrections to the session-4 test spec (verified numerically, do not
"fix" the tests back):**
1. The asymptote bound `kcap*(kcap/vov)**m` in §1c is **wrong for m < 1**.
   Leading term is `(kcap/m)*(kcap/vov)**m` (Bernoulli); m=0.5 violates the
   §1c constant by 1.34× at u = vov/kcap = 10. Test uses `kcap/m` factor.
2. Near the asymptote the true error falls below float64 resolution of kcap
   (m=40: u^-m ~ 1e-30 << eps), and `err = kcap - capped` is computed by
   cancellation. Test adds an `8*eps*kcap` floor to the bound and allows
   ±8 eps on monotonicity/ceiling. `_soft_cap` itself rounds ~1 ulp *over*
   kcap at extreme u (0.050000000000000050) — harmless, pinned in the test.
3. Round-trip start **must set `gk=0`**: `fit()` takes non-free params from
   `start`, and DEFAULT's `gk=0.4` injected ~7.7 A of phantom kink current
   (RMSE 3166 %). Same gotcha class as session 4 §10 ("plumbing silently
   lies").

### 1d. Docs — **NOT DONE. This is the only remaining implementation task.**

- `README.md`: the Model section still omits the cap entirely; add
  `Vov_eff = Vov/(1+(Vov/kcap)**m)**(1/m)` to the equation block + a short
  "Gate-overdrive charge cap" subsection (physical picture is already written
  in `params.py`/`model.py` docstrings — reuse it). Mention the cap preserves
  `dg/dId <= -1` and that the solver has a stall guard.
- `src/vdmos_fit/STATUS.md`:
  - Problem 6 ("no real datasheet data") — **close it**: `data/ST8L65N044M9_IV.csv` exists (41 pts).
  - Add a "Gate-overdrive cap (kcap, m) — session 4+5" section: cap replaces
    `kappa` (6.88 % vs 7.45 % with full rank, session-4 §3a/§5), the
    bit-identity guarantee, the stall guard, and the **new §3 results below**.
  - Note the root `STATUS.md` deletion is still unstaged (session-4 §0).
  - Update "Verified state": 85 passed + 1 skipped, ruff clean.
- `data/README.md` — check it describes the ST8L65 CSV provenance (digitised,
  25 °C pulsed, 0–20 V, the removed 650 V off point).

---

## 2. Decisions (unchanged from session 4 — settled, do not re-litigate)

1. **Variant A** (cap before Vk). 2. **rd, rs ≤ 36 mΩ per device** (implemented).
3. **ksub free; theta/vsat/eta fixed** (implemented). 4. **lam ≤ 0.004 1/V**
   (implemented). 5. Two fits with/without Vgs=5 (**done, §3**).
6. Vgs=5 duplicate Id points are real data; all data 25 °C pulsed; no
   Rds(on)-vs-Vgs curve in the datasheet.

---

## 3. NEW RESULTS — the lam-capped fits (this session, `out/st8l65_capfit_params.txt`)

Both fits: all 4 multi-starts converged to the **same basin** (no basin
hopping — the cap + capped lam made the landscape tame).

| | fit A (full, 41 pts) | fit B (no Vgs=5, 36 pts) |
|---|---|---|
| RMSE_rel | **9.17 %** | **5.54 %** |
| kp | 63.59 | 50.60 |
| vth0 | 4.6656 | 4.4460 |
| ksub | 0.02336 | 0.02645 |
| lam | **0.004 PINNED UPPER** | **0.004 PINNED UPPER** |
| gds | 6.766e-5 S | 6.757e-5 S |
| rd | 31.92 mΩ | 29.85 mΩ |
| rs | 3.8e-9 **PINNED LOWER** | 0.66 mΩ |
| vk | **4.245** | **5.562** |
| kcap | 2.610 V | 2.904 V |
| m | 6.203 | 10.241 |
| rd+rs | 31.92 mΩ | 30.51 mΩ (chord 34.5 mΩ) |
| plateau Vgs=10 @19.7 V | **226.4 A** (meas 240.2, −5.7 %) | **223.9 A** (−6.8 %) |
| rank / cond | 10/10 / **3.1e9** | 10/10 / **2.2e9** |
| per-curve % | 3.7:0.13, 5:21.26, 6:6.03, 7:2.48, 8:6.04, 9:6.10, 10:6.84 | 3.7:0.0, 6:5.83, 7:2.02, 8:6.12, 9:6.10, 10:6.30 |
| Vgs≥6 only | 5.83 % | 5.62 % |

σ_rel (fit A): kp 3.5, vth0 0.52, ksub 5.8e-4, lam 15.6, gds 4.0, rd 0.027,
rs 0.036, **vk 112**, kcap 5.3, m 5.0. Correlations >0.8: kp↔vth0 +0.98,
kp↔rs +0.96, kp↔kcap −0.97, vth0↔rs +0.95, vth0↔kcap −0.97,
**ksub↔vk ±1.000**, rs↔kcap −0.98.

### What these numbers say (the new open problems, replacing session-4 §7)

1. **lam is pinned at its 0.004 bound in both fits** — the data wants more
   output slope than the cap allows. Consequence: the plateau is now **5–7 %
   LOW** (226/224 A vs 240 A), the mirror image of session 4's problem (then
   β was too steep at high Vgs). Session-4 §7 item 2 predicted the
   degradation: 6.88 % → 9.17 % full / 5.6 % on Vgs≥6. The honest conclusion
   is forming that a **bias-dependent drain element is genuinely needed** —
   but *not* `kappa` (pole, divider physics). Candidates: the `rdrift`+`kappa`
   variant was 2.85 % on Vgs≥6 (session-4 §3a) but reintroduces the pole;
   a *bounded* bias-dependent rd (e.g. `rd*(1+kd*Vds)` with the same
   `dV/dI <= 1` discipline, or JFET-style depletion of rd) is the cleaner
   next experiment. **Not yet tried.**
2. **ksub ↔ vk correlation = ±1.000 and cond jumped 3.8e4 → ~2e9.** ksub is
   no longer "dead" (session-4 §3e): with vk free and lam capped, the softplus
   floor and the tanh gate trade off exactly (σ(ksub)=5.8e-4 is meaningless —
   the pair is degenerate, not pinned). ksub drifted to 0.023–0.026 (mid-box,
   not to a bound, so session-4 §7 item 3's "expect inf σ" did not happen —
   the degeneracy is with vk instead). **Recommendation for next session: fix
   ksub (e.g. 0.05 or from a subthreshold spec if one exists) or tie vk to Vk;
   re-measure cond.** The three-gate coupling of vk (DIBL/CLM/kink all use
   `tanh(Vds/vk)`) is still unexamined (session-4 §7 item 1).
3. **vk rose 0.65 → 4.2–5.6 V** — the knee is now much softer, consistent
   with the user's "knee too sharp" observation (session-4 §7 item 1). The
   capped-lam refit largely fixed the knee *shape* complaint at the cost of
   the plateau level (item 1 above). σ(vk) ≈ 112–150: still barely sensed.
4. **Vgs=5 curve costs 21.3 % on fit A** (session-4 §3f stands: unreproducible
   by any C1 model). Fit B (5.54 %) is the fair model-quality number.
5. **rs pinned at 0 in fit A** (fit B: 0.66 mΩ) — rd/rs split still
   unidentifiable; quote only rd+rs ≈ 30–32 mΩ.

---

## 4. Parameter table — unchanged from session-4 §4

Global `params.BOUNDS`; ST8L65 overrides `rd, rs ≤ 0.036`, `lam ≤ 0.004`.
Free (10): `kp, vth0, ksub, lam, gds, rd, rs, vk, kcap, m`. Fixed:
`theta=0, vsat=1e4, eta=0, gk=0, vthk=1.6, kk=0.1, vov_floor=1e-3`.
`kcap` linear space (0 = off, never log); `m` log. `kappa` DELETED.

---

## 5. Why the cap replaced `kappa` — unchanged, keep for docs

See session-4 text preserved in `model.py`/`params.py` docstrings (already
written there — copy into README/STATUS for §1d): `kappa` was a voltage
divider (`Id → Vds/(rd(1+κVds))`), not a charge limit; the fit made rd carry
75 % of terminal voltage at the 0.483 V Ron anchor. The cap is the charge-side
mechanism (finite free charge in the SJ pillars), leaves rd a plain resistor,
keeps `0 <= dVov_eff/dVov <= 1` (solver bracket invariant survives), and
removes the `u(Id)` pole. Cap beats kappa: 6.88 % vs 7.45 % (session 4,
unconstrained lam), full rank, ~1000× better conditioning.

---

## 6. Do NOT re-run the "fully free, widened bounds" fit — still cancelled

Session-4 §6 stands. The crash (`RuntimeError ... 192 A`) was the `u(Id)`
pole; deleting `kappa` removed it. The stall guard (§1a) handles the residual
two-cycle mode; do not "improve" the solver further without a failing test.

---

## 7. Next-session queue (in priority order)

1. **§1d docs** (README + STATUS + check data/README) — the only unfinished
   session-4 task. Then **commit** (§9).
2. **Break the ksub↔vk degeneracy** (§3 item 2): fix ksub or tie vk↔Vk;
   re-run `scripts/fit_st8l65.py`; cond should drop back to ~1e4–1e5.
3. **Bounded bias-dependent drain element** to un-pin lam and fix the −6 %
   plateau (§3 item 1). Design constraint: keep `dV_element/dId` bounded so
   the bracket invariant survives; do NOT resurrect `kappa` as-is.
4. α/β ladders (session-4 §3c/§3d) against the *new* fit B params — the β
   ladder should now be too *shallow* at high Vgs (mirror of before); verify.
5. Datasheet Figure 3 (typical output curves, 25 °C) was **not saved to disk**
   — ask the user for it again when checking knee shape.
6. Body diode / reverse polarity (STATUS problem 7) still the biggest model
   hole; package-R path is Vds ≥ 0 only.

---

## 8. Files touched this session

**Modified:** `src/vdmos_fit/model.py` (cap + `_soft_cap` + kappa removal +
stall guard), `scripts/fit_st8l65.py` (full rewrite).
**Created:** `tests/test_gate_cap.py`.
**Generated (gitignored):** `out/st8l65_capfit.png`, `out/st8l65_capfit_params.txt`.
**Unchanged this session:** `params.py` (was already done), all other tests.
**Still to touch:** `README.md`, `src/vdmos_fit/STATUS.md`, maybe `data/README.md`.

## 9. Git state — NOTHING COMMITTED YET

```
## main...origin/main
 M HANDOVER.md          (session 4's notes; now replaced by this file)
 D STATUS.md            (root; moved to src/vdmos_fit/STATUS.md — still unstaged since session 3)
 M src/vdmos_fit/model.py
 M src/vdmos_fit/params.py
 M scripts/fit_st8l65.py
?? tests/test_gate_cap.py
?? .vscode/             (leave uncommitted or gitignore — user's call)
```
Commit is safe now (85 passed, 1 skipped, ruff clean). Suggested: stage
everything except `.vscode/`, one commit "remove kappa, add gate-overdrive
charge cap (kcap,m) + solver stall guard; ST8L65 cap fits". Add the root
STATUS.md deletion (`git rm STATUS.md` equivalent) to the same commit.

## 10. Gotchas (carry over + new)

- All session-4 §10 gotchas still apply: `sys.executable -m pytest` (base
  Python lacks pytest/ruff); frozen `Params` + `Params(**d)` splat; print one
  model current before trusting an RMSE; `identifiability(res, free)` needs
  the same `free` as `fit`; package-R is Vds ≥ 0 only; re-anchor START per
  device; kcap linear-space, "off" = not in free set.
- **NEW: `fit()` takes non-free parameters from `start`, not from DEFAULT or
  the truth** — a start with the wrong fixed params (e.g. `gk=0.4` when the
  truth has `gk=0`) silently poisons the fit (3166 % RMSE). Always build
  `start` from the truth's fixed values.
- **NEW: the bracketed Newton can two-cycle** on flat-then-cliff `g(Id)`
  (gate degeneration across the cap knee). The stall guard in
  `_solve_package_r` is load-bearing; `test_no_ndr_with_fitted_cap` fails
  without it (Vds ≈ 178 V).
- **NEW: asymptote of the soft cap is `(kcap/m)(kcap/vov)^m`**, not
  `kcap(kcap/vov)^m`, for m < 1.
- Near-asymptote comparisons of the cap hit float64 cancellation floors
  (eps·kcap); don't chase "violations" at the 1e-16 level.
