# Handover — ImprovedVDMOS (session 9 complete)

## Current status

Session 9 added the **drain-knee softness exponent `nd`** — the missing
shape-only knob identified in session 8's knee analysis. It replaces the tanh
soft-switch with the same power-law soft-min used for the gate cap:

```
v = Vk·tanh(Vds/Vk)                     nd = 0  (default, legacy)
v = Vds/(1 + (Vds/Vk)^nd)^(1/nd)        nd > 0
```

Shape-only by construction: `dv/dVds(0) = 1` and `v → Vk` hold for every `nd`,
so Ron and the Id,sat plateau are untouched while the knee width shrinks
monotonically with `nd` (nd ≈ 2.8 reproduces the legacy tanh, L2 rms
0.024·Vk). This decouples knee width from knee position, which theta/vsat
could not do (they are ~84 % correlated in fits). `nd = 0` keeps every
pre-nd anchor bit-identical — pinned by `tests/test_drain_softness.py`
(9 new tests: bit-identity, 0 ≤ dv/dVds ≤ 1 solver-bracket invariant,
asymptotes, shape-only sweep, no-NDR, synthetic round trip recovering nd).

Files touched: `params.py` (nd field, bounds (0.5, 40), deliberately NOT in
LOG_PARAMS — 0 is the sentinel — and NOT in DEFAULT_FREE, same rationale as
kcap/m), `model.py` (cap math refactored into shared `_soft_min`, log-space,
sign-preserving for negative Vds), README ("Drain-knee softness (nd)"
section), new test file. **98 tests (96 pass + 2 headless skips), ruff clean.**

**The session-9 work is UNCOMMITTED on disk** — review and commit first.

Caveat found while testing (pre-existing, not a regression): with the DEFAULT
demo set (vsat = 1.5 ≪ Vk ≈ 5 V) the velocity-saturation divisor creates a
mild gds dip near the knee even at nd = 0 (−1.9 S at Vds ≈ 6 V, Vgs = 10;
the −2.0 S tolerance in `test_improved_model_has_no_negative_gds` was already
absorbing it). A sharp `nd` deepens it. Property of those parameters, not of
the switch — documented in README; refit vsat up or nd down if a device lands
there.

Session-8 state still stands: model/fitter/diagnostics in intended final
shape for the static-I-V scope; official fit ST8L65N044M9 fit B, **5.54 %
relative RMSE, rank 9/9, cond 3.75e4** (table in `src/vdmos_fit/STATUS.md`).
Datasheet Figure 3 on disk: `data/ST8L65N044M9_Id_Vds_Vgs.png`.

## Next tasks (do these first, before the Queue)

1. **Commit the session-9 (nd) work** after a quick review.
2. **Re-run the ST8L65N044M9 fit with the updated model**
   (`scripts/fit_st8l65.py`, fit-B recipe and DEVICE_BOUNDS). Free `nd`
   explicitly — it is not in DEFAULT_FREE; add it to the free tuple, start
   near 2.8 (≈ legacy tanh), fitted in linear space (lower bound 0.5).
   Compare relative RMSE against the official 5.54 % and check whether the
   theta↔vsat ~84 % correlation pair loosens now that knee width has its
   own knob. Watch `on_bounds` (nd pinned at 0.5 or 40 means the data does
   not want the new shape).
3. **Save to `out/`**: the fitted parameters (same format as
   `out/fitted_params.txt`) AND the convergence/identifiability probe —
   `FitResult.summary()` plus `identifiability(res).summary()` (rank, cond,
   sigma_rel, pinned bounds, |corr| > 0.8 pairs) — as a text file. Include
   the pre-fit `sensitivity_report` if it is cheap.
4. **Plot Id vs Vds at each measured Vgs point** (the Vgs grid of the
   digitised CSV in `data/`), measured vs fitted overlay (`plotting.py`
   overlay machinery), PNGs into `out/`. This is also the visual check for
   task 4 of the Queue (knee shape vs Figure 3) — the nd fit should now be
   able to match the knee without dragging the plateau.

## Queue

1. Digitise a near-threshold transfer curve to pin the real subthreshold
   swing (fit B ksub = 0.184 V ≈ 420 mV/dec is a one-point anchor artifact).
2. Bounded bias-dependent drain element to un-pin lam (0.004 UPPER) and fix
   the −6 % plateau — NOT kappa-as-is.
3. α/β ladders vs fit B (β ladder should now be too shallow at high Vgs).
4. Knee-shape check of fit B vs Figure 3 (figure now on disk; see Next
   task 4 — nd is the tool for this).
5. Body diode / reverse polarity + HV leakage — biggest model hole; the
   package-R path is Vds ≥ 0 only.
6. Then the architecture: Q(V), Tj node, LTspice export, CI.

`.vscode/` remains untracked (user's call). `out/` is gitignored; current
contents are the fresh 9-param run plus `nd_knee_sweep.png` (session-9 demo).
