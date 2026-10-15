# Handover — ImprovedVDMOS (session 10 complete)

## Current status

Session 10 executed session 9's "next tasks" in full.

**Session-9 (nd) work: committed** as `7c771e1 add knee softness parameter`
(params.py nd field + bounds, model.py shared `_soft_min`, README section,
tests/test_drain_softness.py). 98 tests green (1 headless skip), ruff clean. 

**ST8L65 refit with nd free** — new script `scripts/fit_st8l65_nd.py`
(same fit-A/fit-B recipe as fit_st8l65.py; nd added to the free tuple,
started near 2.8, linear space bound **0.1–40** — relaxed from 0.5 during the
session, see below):

| | nd fit A (**OFFICIAL**) | nd fit B | superseded 9-param fit B |
|---|---|---|---|
| RMSE_rel | **3.79 %** (was 9.97) | **1.52 %** (was 5.54) | 5.54 % |
| nd | 0.4359 (interior) | 0.4081 (interior) | 0 |
| lam | 5.1e-10, pinned LOWER 1e-12 | 1.9e-12, pinned LOWER | 0.004 pinned |
| plateau Vgs=10 | 246.9 A (+2.8 %) | 244.7 A (+1.9 %) | 223.9 A (−6.8 %) |
| rank / cond | 8/10 / 2.19e3 | 8/10 / 1.27e4 | 9/9 / 3.75e4 |

Headline findings (details + full correlation tables in
`src/vdmos_fit/STATUS.md`, "ST8L65 fit with nd free"):

- The data wants a knee **softer than the legacy tanh, and softer than the old
  0.5 guard allowed**: nd lands at 0.436 / 0.408, interior to the relaxed
  (0.1, 40) bound. Probes from starts nd = 0.5 and 0.3 both converge to the
  same interior optimum, so the 0.1 floor is not load-bearing. nd → 0
  collapses the current, so lower is not physical; do not lower it further.
- **Bounds were widened mid-session** (user): kp ≤ 1e4, vsat ≤ 1e4,
  lam 1e-12–1.0, vk ≤ 1300. This released the kp / vk / lam bound-pins. The
  old per-device `lam ≤ 0.004` cap was **removed from DEVICE_BOUNDS** — a
  probe shows it is inert: capped vs free lam give identical RMSE and
  identical nd/kp/kcap/m to 4 sig figs, only vk moves (1300 ↔ 0.249).
- The −6 % plateau deficit and the lam bound-pin (old open problem 1) were
  knee-SHAPE problems: nd free un-pins lam and fixes the plateau. Queue 2
  (bias-dependent drain element) is NOT needed on this data.
- theta↔vsat probe: freeing theta/vsat pushes both to their bounds — the
  data rejects that mechanism; corr(theta,vsat) = −0.848 appears in the
  12-param run but nd is now the correlation hub (nd↔ksub +0.988). The
  pair becomes moot rather than loosening.
- Knee cluster (kp, kcap, m, nd, vk, lam) is strongly entangled; rank 8/10
  (lam and vk dead, σ = inf — **expected**: no data beyond Vds = 20 V, so there
  is no deep-saturation signal for CLM or the drain gate to act on). Quote
  shape + rd+rs only, never individual knee parameters. gds ≥ 0 everywhere for
  the official fit (min −3e-12 S = finite-difference noise at the 1e-14 A
  subthreshold floor, not a physical NDR).
- One multi-start caution: the nd = 6 start stays stuck at 13.5 % — always
  multi-start when freeing nd.

**Saved to `out/`** (gitignored): `st8l65_nd_params.txt` (both fits, same
format as the capfit file), `st8l65_nd_probe.txt` (pre-fit sensitivity_report
+ FitResult.summary + identifiability summary + theta/vsat and nd-bound
probes), `st8l65_nd_fit.png` (full + knee zoom, nd vs official),
`st8l65_nd_percurve.png` (Id vs Vds at each measured Vgs — the Figure 3
knee check, matches at every Vgs), `st8l65_nd_overlay.png`.

**Decision (made 2026-10-13): fit A of `fit_st8l65_nd.py` is now the OFFICIAL
reference** — 10-param with nd free, 3.79 % on the full 41-point set. It
supersedes the 9-param fit B (5.54 %). Rationale: fit A drops no curve, so its
RMSE is the honest device-quality number (fit B's 1.52 % stays as a diagnostic
of the unreproducible Vgs = 5 curve), and it is better conditioned and more
physical (cond 2.19e3 vs 1.27e4; kp 579 vs 2354-on-the-bound; rd/rs both
interior). The parameter record now lives in the package —
`vdmos_fit.reference.ST8L65_OFFICIAL` (single source of truth) — re-exported by
the script's `OFFICIAL` dict, and `main()` checks every run reproduces it to
1e-3 relative and prints `[STALE]` if it does not: re-promote deliberately,
never silently. New `tests/test_st8l65_official.py` (5 tests) pins the record
against the measured CSV — RMSE, Ron chord, plateau, Vgs(th), per-curve, nd
interior, no negative gds — so any model/bounds change that moves the official
fit fails the suite. Suite is now **103 tests green**, ruff clean.
`lam`/`vk` σ = inf is accepted as expected behaviour for this dataset.

## Next tasks (do these first, before the Queue)

1. ~~Decide the "official fit" question~~ — **done**: fit A promoted
   (script `OFFICIAL` dict + STATUS.md + README + this file).
2. Anchor re-check at the official set — **done**: min gds over
   Vds 0–200 V × Vgs −2…16 V is −3.1e-12 S (FD noise at the subthreshold
   floor, no NDR), model finite everywhere, solver converges, and the full
   suite is green with the relaxed nd = 0.1 floor. The vsat-on demo-set caveat
   is untouched (it is a DEFAULT-params property, not an official-fit one).
3. If the official set is ever exported to LTspice, remember `lam` and `vk`
   are unsensed there — pin lam = 0 and any vk rather than shipping the fitted
   numbers.

## Queue

1. Knee-cluster identifiability: the shape is pinned but the kp/kcap/m/nd/vk/lam
   split is not (rank 8/10, |corr| up to 0.999). Dense knee-region data
   (Vds 0.5–8 V, fine Vgs steps) is what would break it.
2. Body diode / reverse polarity + HV leakage — biggest model hole; the
   package-R path is Vds ≥ 0 only.
3. Then the architecture: Q(V), Tj node, LTspice export, CI.

`.vscode/` remains untracked (user's call). `out/` is gitignored; new
contents are the nd run (`st8l65_nd_*`) plus `_nd_fit_run.log`.
