# Handover — ImprovedVDMOS: gds leak term DELETED (session 7). Code edits DONE, docs + pytest/ruff + commit PENDING.

## 0. READ THIS FIRST — where session 7 stopped

Task: "delete everything that doesn't make physical sense", starting with the
`gds·Vds` drift-leakage term (reverse leakage will be a body-diode/leakage
element later; not enough data to model it now). Hypothesis: deleting gds lets
ksub capture subthreshold/near-threshold conduction. **Confirmed by probes.**

Status when tokens ran out:
- **Probes (steps 1–3): DONE**, results in `out/` (see §3).
- **Step 5 code edits: DONE** (model/params/data/fit + all 4 test files + both
  scripts). See §1.
- **Step 5 docs: NOT DONE.** README.md, src/vdmos_fit/STATUS.md, data/README.md
  still mention gds. Exact edits listed in §2.
- **pytest + ruff: NOT RUN** since the code edits. Run before commit.
- **Commit 1 (gds deletion): NOT DONE.**
- **Step 6 (on_bounds bug fix + test): NOT DONE.** Plan in §4.
- **Fits/plots: SKIPPED** (user: not enough tokens). `out/st8l65_capfit_params.txt`
  is now STALE (10-param, still has gds). The authoritative gds=0 fit numbers are
  in `out/nogds_probe.txt`. Re-running `scripts/fit_st8l65.py` (now 9-param) is a
  next-session task, not done.

Git: HEAD is session-6 commit `addd7ca` (cap + docs + commit). Working tree has
uncommitted gds-deletion edits to: `src/vdmos_fit/{model,params,data,fit}.py`,
`tests/{test_model,test_package_r,test_fit,test_gate_cap}.py`,
`scripts/{fit_st8l65,roundtrip_demo}.py`. `.vscode/` still untracked (user's call).

---

## 1. What was implemented this session (code, DONE)

Deleted the `gds` parameter and its `gds·Vds` term everywhere:

- **`model.py`**: removed `idrift = p.gds * vds` from `_channel_current`
  (`return ich + ikink`); removed `+ gds*Vds` from the module-docstring equation
  block; changed "theta = vsat = eta = gds = gk = 0" → "…= gk = 0"; added a
  "No drift-leakage term (gds deleted, session 7)" docstring section explaining
  the one-point-parameter finding. Also fixed a pre-existing typo
  `(rdp-rs)*kp*Vov` → `(rd+rs)*kp*Vov`.
- **`params.py`**: removed the `gds` field (was default `1e-3`), its docstring,
  its `BOUNDS` entry `(1e-9, 10.0)`, its `LOG_PARAMS` membership, and its
  `DEFAULT_FREE` entry. DEFAULT_FREE now 10 names (kp,vth0,theta,vsat,eta,lam,
  gk,vthk,rd,rs).
- **`data.py`**: removed `gds=2e-3` from `example_dataset()`'s truth.
- **`fit.py`**: removed the "gds only matters in the off state" clause from the
  `sensitivity_report` docstring.
- **`tests/test_model.py`**: stripped `gds=0.0` from the ohmic-limit,
  level1-limit, subthreshold-slope and ksub tests; **renamed**
  `test_off_state_drift_conductance_defined` → `test_off_state_is_gate_dependent_not_a_leak_resistor`
  (asserts 0 < Id(20,0) < 1e-9 and gate-dependence: Id(20,2) > 1e4·Id(20,0));
  kink one-sided test now asserts off < 1e-12 instead of `== gds*1.0`;
  monotone-in-vgs docstring updated (channel-only, resolvable to ~1e-47 A).
- **`tests/test_package_r.py`**: stripped gds from RTRUE and the brentq test;
  **renamed** `test_off_state_exact_series_fix` → `test_off_state_solves_to_channel_floor`
  (off-state Id is now the softplus floor, self-consistent through package-R,
  0 < Id < 1e-9, terminal ≈ internal); stall-guard test params stripped of gds
  (still passes: verified Id(178,10) ≈ 789 A, max rel resid 8.9e-11 — the
  two-cycle cliff is gate-degeneration physics, NOT the deleted leak).
- **`tests/test_fit.py`**: stripped gds from TRUTH, bounds test, uncertainty-set
  assert, legacy_free (now 8 names), and the dead-parameter test (free tuple now
  7 names; the "genuinely sensed but weak stays live" witness switched from gds
  to eta, σ≈12).
- **`tests/test_gate_cap.py`**: stripped gds from FITTED, CAPTRUTH, the
  subthreshold-untouched test, and the round-trip free tuple (now 9 names).
- **`scripts/fit_st8l65.py`**: FREE now 9 names (kp,vth0,ksub,lam,rd,rs,vk,kcap,m);
  removed `GDS_OFF`; docstring rewritten (off point now anchors ksub/vth0 via
  subthreshold channel); `report()` now prints the Vgs(th) point (3.7,3.7)
  current; tail note updated.
- **`scripts/roundtrip_demo.py`**: stripped gds from TRUTH and the print loop.

### Verified numerically BEFORE the edits (so the next session need not re-derive)
With DEFAULT at gds=0 (simulating the deletion), every reference anchor still
passes its existing tolerance — **no anchor table needs changing**:
- Ron @Vgs=10 = 4.271 mΩ (test wants 4.27±0.02) ✓
- Id,sat(20V) = 31.448 / 91.017 / 157.061 / 225.879 / 297.670 A
  (test table 31.5/91.0/157.1/225.9/297.7, rel 0.01) ✓
- α = 1.0781 (test 1.08±0.01) ✓
- Level-1/improved ratio @Vgs 10/11/12 = 4.2518 / 4.7997 / 5.3334
  (refs 4.25/4.80/5.33, rel 5e-3) ✓
- off Id(20,0) = 9.76e-14 A; monotone-in-Vgs still holds (strictly >0 from
  Vgs=1.5 at all Vds) ✓
- legacy-8 theta↔vsat pair still reported (cond 247) ✓
- dead-param rewrite: ksub still inf/dead, eta live (σ 12.5), rank 6/7 ✓

---

## 2. PENDING docs — exact edits to make (then pytest + ruff + commit 1)

### README.md (user-specified, do exactly these)
1. Comparison table row: `| Off-state slope | undefined | gds·Vds |`
   → `| Off-state slope | undefined | softplus subthreshold floor, gate-dependent |`
2. Equation block: delete the `      + gds·Vds` line.
3. "Setting `θ = vsat = η = gds = gk = 0` makes…" → drop `gds` from that list.
4. Fitting strategy 1: "the low-current shape (vth0, ksub, gds) is not swamped"
   → "(vth0, ksub)".
5. Add a short **"Off state"** note after the "Gate-overdrive charge cap"
   subsection: the `gds·Vds` leak term was deleted (session 7); on the ST8L65
   data it was a one-point parameter (100 % of the 250 µA Vgs(th) point, <0.03 %
   of every other point) that flattened the subthreshold curve to a
   gate-independent floor and stole the threshold anchor from ksub. The physical
   off state is now the gate-dependent softplus floor (slope 1/ksub); real HV
   leakage will come with the body diode / a dedicated leakage element, not a
   fitted resistor.

### src/vdmos_fit/STATUS.md
- **Close problem 5** ("gds σ ≈ 455 % — only the off-state constrains it"):
  mark RESOLVED — the gds term was deleted (session 7), so the problem is moot.
- Add a **"gds leak term deleted — session 7"** section: the one-point-parameter
  evidence, the gate-independent-floor flaw, that ksub becomes live (σ 3.5) and
  cond drops 2.2e9 → 3.7e4 on fit B when gds is removed (see out/nogds_probe.txt),
  and that all DEFAULT anchors survive unchanged (list them, §1 above).
- Update "Verified state" test count **after you run pytest** (was 86 passed,
  1 skipped; the gds edits renamed tests but should not change the count —
  confirm).
- Note the sessions-4+5 fit table's gds row and rank/cond are superseded by the
  gds=0 probe (out/nogds_probe.txt); a fresh 9-param `fit_st8l65.py` run is due.

### data/README.md
In the "Datasets / ST8L65N044M9_IV.csv" section:
- "the only off-state point, anchors gds ≈ 68 µS" → "anchors ksub/vth0 through
  the subthreshold channel (the gds leak term was deleted, session 7)".
- The "Removed point" bullet: "the channel model has no HV leakage term, and one
  gds cannot give both 1 µA at 650 V and 250 µA at 3.7 V" → reword without gds:
  no HV leakage element exists yet; real leakage waits for the body diode.

### After docs
`python -m pytest` (expect ~86 passed, 1 skipped) and `python -m ruff check .`
must be clean, then **commit 1**: "delete gds drift-leakage term; off state is
the gate-dependent softplus floor". Stage everything except `.vscode/`.

---

## 3. Probe results (DONE, in out/) — the evidence for the deletion

- `out/baseline_with_gds.txt` — copy of the with-gds fit A/B (session 6).
- `out/nogds_probe.txt` — multi-start refit with gds fixed at 0 (term still in
  code, just not free), 4 starts × {fit A, fit B}.
- `out/nogds_probe_fitB.json` — best gds=0 fit-B params.
- `out/nogds_subthreshold.png` — log Id(Vgs) @Vds=50 mV, gds=0 vs gds-kept.

**Fit B (36 pts, no Vgs=5) — every criterion passes, all 4 starts → one basin:**
RMSE 5.54 % (identical to with-gds); off point hit exactly (250.0 µA via the
softplus channel, not a leak); **ksub now LIVE σ=3.46** (was dead, σ=inf);
Jacobian rank 9/9, **cond 3.75e4** (was 10/10, cond 2.2e9 — gds was pure
degeneracy); plateau 223.9 A (−6.8 %, unchanged); Vgs≥6 5.62 %. ksub=0.1836,
kp 50.58, vth0 4.4459, vk 5.562, kcap 2.904, m 10.24, rd 29.85 mΩ, rs 0.66 mΩ.
lam still pinned UPPER 0.004 (unchanged open problem).

**Fit A (41 pts, with the anomalous Vgs=5 curve) — the known casualty:** RMSE
9.17 → 9.97 %, off point abandoned (0 µA), ksub dead (σ=inf), rank 8/9. The
Vgs=5 curve's steep 5→6 V gate slope fights the 250 µA anchor; gds was the cheat
that decoupled them. Per decision 6 (Vgs=5 curve is real data but unreproducible
by any C¹ model) fit B is the fair number, so this does not block deletion.

**Subthreshold plot:** gds=0 curve is monotone, floor 1.4e-13 A, crosses 247.9 µA
at Vgs=3.7 V. gds-kept curve has a flat gate-independent 3.38 µA floor (slope
∞ mV/dec) — the physically wrong shape the deletion removes.

**Honest caveat:** gds=0 fit B gives ksub=0.184 V ⇒ ~423 mV/decade swing, not
physical (real ~60–100 mV/dec). With only ONE off-state point the data pins the
250 µA anchor, not the slope; ksub's σ=3.5 overstates what we know. Do NOT
hard-fix ksub to 0.05 to force physicality — that pushes the anchor into vth0.
A near-threshold transfer curve (data wish) would pin the real swing.

---

## 4. Step 6 — on_bounds false-positive bug (PENDING, separate commit)

`FitResult.on_bounds(rel_tol=1e-3)` in `src/vdmos_fit/fit.py` flags a param as
pinned if `abs(val - lo) <= rel_tol * max(span, abs(lo), 1e-12)`. With a wide
span this is a huge absolute tolerance: gds span was (1e-9, 10.0) → span 10 →
tol 1e-2 S, so **any gds < 0.01 S read "pinned LOWER"** even at 6.8e-5 (a false
positive seen in the session-6 fit output). After the gds deletion the bug still
mis-scales for other wide-span params.

Fix: make the pinned test relative to the bound's own magnitude / a per-param
scale, not the full span. E.g. compare `abs(val-lo) <= rel_tol * max(abs(lo),
abs(val), <small floor>)` for the lower bound (and symmetric for upper), so a
param at 6.8e-5 with lo=1e-9 is NOT called pinned, while rs=3.8e-9 with lo=0
(span 0.036) IS (genuinely at zero). Add a test: a param that lands well inside
a wide box must not be flagged; a param genuinely at a bound must be. Verify
lam=0.004 (UPPER 0.004, real) and rs=0 (LOWER 0, real) still flag correctly.
Commit separately from the gds deletion.

---

## 5. Next-session queue after §2/§4 (unchanged from session 6 §7 items 2–6)

1. Re-run `scripts/fit_st8l65.py` (now 9-param, gds gone) → regenerate
   `out/st8l65_capfit.png` + `out/st8l65_capfit_params.txt`; expect fit B ≈ 5.5 %,
   cond ~1e4, ksub live. Update the STATUS fit table with the real new numbers.
2. Break the ksub↔vk degeneracy / decide ksub handling now that it's live but
   one-point-anchored (fix from a subthreshold spec if you can digitise a
   near-threshold transfer curve).
3. Bounded bias-dependent drain element to un-pin lam and fix the −6 % plateau
   (keep dV_element/dId bounded so the package-R bracket invariant survives; do
   NOT resurrect kappa as-is).
4. α/β ladders vs the new fit B params (β ladder should now be too shallow at
   high Vgs — verify).
5. Datasheet Figure 3 (typical output curves, 25 °C) still not on disk — ask the
   user before checking knee shape.
6. Body diode / reverse polarity (STATUS problem 7) — biggest model hole; the
   package-R path is Vds ≥ 0 only. This is now also where HV leakage belongs.

## 6. Gotchas (carry over + new this session)
- `fit()` takes non-free params from `start`, not DEFAULT/truth — build start
  from the truth's fixed values (gk=0 gotcha still live).
- New: the gds term is GONE. Any code/test/script that passes `gds=` to Params
  now raises TypeError. All in-repo uses were removed; watch notebooks/ and any
  external scripts.
- New: DEFAULT anchors (Ron 4.27 mΩ, sat table, α=1.08, L1 ratios) are unchanged
  by the deletion — verified, do not "fix" them.
- New: fit B ksub=0.184 V is an anchor-fit artifact of one off-state point, not a
  physical swing; don't over-constrain it.
- `python -m pytest` / `python -m ruff` with the venv interpreter
  (`.../engineering313/Scripts/python.exe`), not base Python.
- out/ is gitignored; the stale `st8l65_capfit_params.txt` still shows gds —
  regenerate before quoting it.
