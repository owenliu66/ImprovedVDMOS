# Handover — ImprovedVDMOS (session 8 complete)

## Current status

All session-7/8 work is committed to `main` (`f9e6173` gds leak-term deletion,
`f1868cc` on_bounds fix, `001b6e8` docs + official 9-param fits). The model,
fitter and diagnostics are in their intended final shape for the static-I-V
scope: one C¹ channel expression + package rd/rs, gate-charge cap (kcap, m),
gate-dependent softplus off state; `kappa` and `gds` deleted as unsupported by
the data. Official fit: ST8L65N044M9 fit B, **5.54 % relative RMSE, rank 9/9,
cond 3.75e4** (table in `src/vdmos_fit/STATUS.md`). 89 tests green, ruff
clean. Datasheet Figure 3 is now on disk: `data/ST8L65N044M9_Id_Vds_Vgs.png`.
Documentation was consolidated: README.md and STATUS.md each carry one short
current-status description instead of per-session history sections.

## Queue

1. Digitise a near-threshold transfer curve to pin the real subthreshold
   swing (fit B ksub = 0.184 V ≈ 420 mV/dec is a one-point anchor artifact).
2. Bounded bias-dependent drain element to un-pin lam (0.004 UPPER) and fix
   the −6 % plateau — NOT kappa-as-is.
3. α/β ladders vs fit B (β ladder should now be too shallow at high Vgs).
4. Knee-shape check of fit B vs Figure 3 (figure now on disk).
5. Body diode / reverse polarity + HV leakage — biggest model hole; the
   package-R path is Vds ≥ 0 only.
6. Then the architecture: Q(V), Tj node, LTspice export, CI.

`.vscode/` remains untracked (user's call). `out/` is gitignored; current
contents are the fresh 9-param run.
