# Handover — ImprovedVDMOS (session 10 complete)

## Current status

Session 10 executed session 9's "next tasks" in full.

**Session-9 (nd) work: committed** as `7c771e1 add knee softness parameter`
(params.py nd field + bounds, model.py shared `_soft_min`, README section,
tests/test_drain_softness.py). 98 tests green (1 headless skip), ruff clean.

**ST8L65 refit with nd free** — new script `scripts/fit_st8l65_nd.py`
(same fit-A/fit-B recipe and DEVICE_BOUNDS as fit_st8l65.py; nd added to
the free tuple, started near 2.8, linear space bound 0.5–40):

| | nd fit A | nd fit B | official 9-param fit B |
|---|---|---|---|
| RMSE_rel | **3.89 %** (was 9.97) | **1.62 %** (was 5.54) | 5.54 % |
| nd | 0.5 PINNED LOWER | 0.5 PINNED LOWER | 0 |
| lam | un-pinned (0.0031) | un-pinned (0.0025) | 0.004 pinned |
| plateau Vgs=10 | 246.3 A (+2.5 %) | 243.1 A (+1.2 %) | 223.9 A (−6.8 %) |
| rank / cond | 9/10 / 1.15e4 | 9/10 / 3.23e4 | 9/9 / 3.75e4 |

Headline findings (details + full correlation tables in
`src/vdmos_fit/STATUS.md`, "ST8L65 fit with nd free"):

- The data wants a knee **softer than the legacy tanh**: nd hits the 0.5
  lower bound from every basin-finding start. Relaxed-bound probe (0.1, 40)
  lands at nd = 0.47 for 1.59 % — the 0.5 sentinel guard costs ~0.03 pt;
  keep it. nd -> 0 collapses the current, so lower is not physical anyway.
- The −6 % plateau deficit and the lam bound-pin (old open problem 1) were
  knee-SHAPE problems: nd free un-pins lam and fixes the plateau. Queue 2
  (bias-dependent drain element) is NOT needed on this data.
- theta↔vsat probe: freeing theta/vsat pushes both to their bounds — the
  data rejects that mechanism; corr(theta,vsat) = −0.848 appears in the
  12-param run but nd is now the correlation hub (nd↔ksub +0.988). The
  pair becomes moot rather than loosening.
- Knee cluster (kp, kcap, m, nd, vk, lam) is strongly entangled; rank 9/10
  (vk dead in fit A, σ_rel = 719 in fit B). Quote shape + rd+rs only, never
  individual knee parameters. gds ≥ 0 everywhere for the nd fit.
- One multi-start caution: the nd = 6 start stays stuck at 13.5 % — always
  multi-start when freeing nd.

**Saved to `out/`** (gitignored): `st8l65_nd_params.txt` (both fits, same
format as the capfit file), `st8l65_nd_probe.txt` (pre-fit sensitivity_report
+ FitResult.summary + identifiability summary + theta/vsat and nd-bound
probes), `st8l65_nd_fit.png` (full + knee zoom, nd vs official),
`st8l65_nd_percurve.png` (Id vs Vds at each measured Vgs — the Figure 3
knee check, matches at every Vgs), `st8l65_nd_overlay.png`.

**Decision left open (user's call):** promote nd fit B (1.62 %) to "official"
despite nd pinned at the bound + rank 9/10, or keep the 9-param fit as the
quoted reference. STATUS.md currently presents nd fit B as best fit and keeps
the 9-param table as the official reference.

## Next tasks (do these first, before the Queue)

1. **Commit the session-10 work** — `scripts/fit_st8l65_nd.py` + STATUS.md
   updates are UNCOMMITTED on disk.
2. Decide the "official fit" question above (promote nd-B or annotate).
3. If promoted: re-check the model anchors/tests against the nd regime
   (nd = 0.5 soft-min at vk ~ 2.6: confirm solver bracket invariant and
   no gds dip in the vsat-on demo set too).

## Queue

1. Digitise a near-threshold transfer curve to pin the real subthreshold
   swing — now the TOP priority: the nd fit pushed ksub to 0.258 V
   (~595 mV/dec), an even more obvious one-point-anchor artifact.
2. ~~Bounded bias-dependent drain element~~ — **dropped for this dataset**
   (nd fixed the plateau and un-pinned lam). Keep the `rd·(1 + kd·Vds)`
   idea in reserve for datasets that re-pin lam.
3. α/β ladders vs the nd fit B (β ladder should now be too shallow at high
   Vgs — recheck against the new fit).
4. Knee-shape check vs Figure 3 — visual check done (out/st8l65_nd_percurve.png);
   formal digitisation of Figure 3 would quantify it.
5. Knee-cluster identifiability: the shape is pinned but the kp/kcap/m/nd/vk/lam
   split is not (rank 9/10, |corr| up to 0.99). Dense knee-region data
   (Vds 0.5–8 V, fine Vgs steps) is what would break it.
6. Body diode / reverse polarity + HV leakage — biggest model hole; the
   package-R path is Vds ≥ 0 only.
7. Then the architecture: Q(V), Tj node, LTspice export, CI.

`.vscode/` remains untracked (user's call). `out/` is gitignored; new
contents are the nd run (`st8l65_nd_*`) plus `_nd_fit_run.log`.
