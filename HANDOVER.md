# Handover — ImprovedVDMOS: `kappa` (drain-bias-dependent drain R) added + ST8L65N044M9 refit (session 3)

State: **78 tests pass (1 skip), `ruff check` clean, `scripts/fit_st8l65.py` runs end-to-end.**
**NOT yet git-committed** (session-2 package-R work is still uncommitted too — see §6).

## 0. Read this first: user correction that invalidates a session-2 conclusion

I reported "rd+rs = 27.5 mΩ is 20 % **below** the measured 34.5 mΩ Ron chord → the fit
sacrifices Ron" as a **caveat/defect**. **The user corrected this: it is normal and expected.**

> "rd+rs lower than R_dson is normal since there's channel resistance. rd and rs come from
> the package (wirebonds and pins)."

Implications to carry forward:
- `Rds(on)` = **rd + rs + R_channel(intrinsic)**. The intrinsic channel resistance is *not*
  a defect to be absorbed into rd/rs. rd/rs should be interpreted as **package only**
  (wire bonds, leadframe, pins), which is their original meaning.
- **RETRACT the `Params.kappa` docstring claim** that rd means "zero-bias drain-side series
  resistance (package + undepleted drift), not package-only". It is wrong for the same
  reason. Fix in `src/vdmos_fit/params.py` (search for "zero-bias drain-side").
- **RETRACT the CAVEATS line in `scripts/fit_st8l65.py`** that says rd+rs being 20 % below
  the Ron chord means "the drift-R fit trades Ron accuracy for the plateau". It does not.
  A *lower* fitted rd+rs with a *higher* model current at the Ron anchor is still a real
  +8.7 % ohmic-region current error, but that is a channel/knee-shape error, not evidence
  that rd/rs are wrong.
- Session-2's "rd+rs ≈ 36–41 mΩ is the trustworthy output" framing is also suspect for the
  same reason — it was effectively fitting rd+rs ≈ Rds(on), i.e. double-counting channel R.

## 1. What was built this session: `Params.kappa` (in the package, not a prototype)

Smallest model change that can produce the measured gate-independent flattening.
`src/vdmos_fit/params.py`:

- `kappa: float = 0.0` (new field, **default 0 → bit-identical to pre-kappa model**; all
  legacy anchors unchanged: Ron 4.27 mΩ, sat table, α = 1.08).
- `BOUNDS["kappa"] = (0.0, 2.0)`, **linear space** (0 is a legitimate answer).
- **`kappa` deliberately NOT in `DEFAULT_FREE`** — adding it would silently change every
  existing fit/anchor. Must be freed explicitly per device.

`src/vdmos_fit/model.py` — `_solve_package_r(..., kappa=0.0)` now handles it; the model
embedding is

```
u     = Id*rd / (1 - Id*rd*kappa)          # voltage across the drain-side resistance
R_d(u)= rd*(1 + kappa*u)
Vds'  = Vds - u - rs*Id                    Vgs' = Vgs - rs*Id
Id    = f_intrinsic(Vds', Vgs')            (implicit, per point)
```

`kappa == 0` takes the **original code path verbatim** (`dg = -rtot*fd - rs*fg - 1`), so the
fast path and all `rd/rs` tests are untouched.

### Solver bug found and fixed (important — this is why the first run crashed)

Bracketing Newton at `hi = Vds/(rd+rs)` alone lets iterates approach the `u(Id)` pole at
`Id = 1/(rd*kappa)`, where `du/dId → ∞`, `dg → -∞`, Newton steps become vanishingly small,
and the solve **stalls** → `RuntimeError: package-R solve failed after 60 iterations:
max |residual| = 296 A`. Two fixes in `_solve_package_r`:

1. **Tighter physical upper bound.** The channel must support a non-negative internal Vds,
   so `u <= Vds` ⇒ `Id <= Vds/(rd*(1 + kappa*Vds))`. Take `min(Vds/rtot, that)`. This sits
   far below the pole.
2. **No-progress safeguard:** accept a Newton step only if it is in-bracket **and** moves
   by `> 1e-12*(hi-lo)`; otherwise bisect.

Also fixed while editing: the `kappa>0` Jacobian must use
`dg = -((dudI + rs)*fd + rs*fg) - 1` — `rs` drops the internal **drain** too, not just the
gate. (My first edit dropped `rs` from the drain term.)

### Verified safe: the linear form is NOT negative-differential-resistance

Checked numerically: `dId/dVds > 0` everywhere up to 200 V; `Id(650 V) = 397.6 A` vs the
asymptote `1/(rd*kappa) = 404 A`. The pole is a **current asymptote**, not an NDR region.
I prototyped a `tanh`-capped variant `R_d = rd*(1+kappa*ucap*tanh(u/ucap))` — it **gains
nothing** (RMSE 7.66–11.2 % vs 7.45 %, and some variants went slightly NDR). **Dropped.**
Superseded prototype kept at `out/prototype_drift_r.py.superseded` (out/ is gitignored).

## 2. Dataset changed (user deleted the 650 V point)

`data/ST8L65N044M9_IV.csv` is now **41 rows** (was 27 in session 2), Vgs grid
[3.7, 5, 6, 7, 8, 9, 10], Vds ≤ 19.7, Id ≤ 240 A. Curves are **denser** now (4–9 pts each).
The `(0 V, 650 V, 1 µA)` off point is **gone** — correct call, the channel model has no HV
leakage term (one `gds` cannot give 1 µA@650 V *and* 250 µA@3.7 V: 1.5 nS vs 68 µS).

- Only remaining off point: `(Vgs=3.7, Vds=3.7, 250 µA)` → `gds = 67.57 µS` anchored directly.
- **Session-2 RMSE numbers are NOT comparable** (27 pts vs 41 pts). Always quote the
  `kappa=0` control on the *same* data.
- Ohmic anchors Vgs=10: (0.483,14), (1.0,29), (1.26,35) → chord 34.5 / 34.5 / 36.0 mΩ.

### Measured exponents (re-derived, on the 41-pt set)

Gate exponent α at Vds ≈ 19 V, for Vth0 = 4.7:
`5→6: 1.77 | 6→7: 1.64 | 7→8: 1.01 | 8→9: 0.157 | 9→10: 0.072`
→ **bimodal, confirmed quantitatively.** α ≈ 1.6–1.8 at low Vov, collapses to ~0.1 at high
Vov (gate-independent = drift/RESURF limited). Output exponent `d ln Id/d ln Vds` along each
curve: 0.12 (Vgs 5), 0.11 (6), 0.42 (7), 0.40 (8), 0.38 (9), 0.41 (10).

## 3. Fit results (the deliverable) — `scripts/fit_st8l65.py`

Relative weighting, 5 multi-starts, `gds` anchored 67.57 µS, `vsat`/`lam` pinned off
(100 V / 1e-6). Free sets:
`BASE_FREE = (kp, vth0, ksub, theta, eta, gds, rd, rs)`; `DRIFT_FREE = BASE_FREE + (kappa,)`.

| | RMSE_rel (41 pts) | per-curve RMSE % (Vgs 3.7/5/6/7/8/9/10) |
|---|---|---|
| **kappa = 0** (control) | **10.35 %** | 0.02 / 18.61 / 4.78 / 9.98 / 7.02 / 8.97 / 10.93 |
| **kappa free** | **7.45 %** | 0.02 / 18.07 / 1.67 / 1.86 / 4.98 / 4.22 / 6.01 |

- **Vgs ≥ 6 only: 8.72 % → 4.30 % (2.0×).** Overall 1.39×.
- **All 5 starts converge to the same basin** for kappa-free → real optimum, not luck.
  (kappa=0 control splits into two basins: 10.35 % ×3 starts, 14.88 % ×2.)

Best (kappa free, 9 free params, converged, rank 9/9, cond 2.87e7):
```
kp = 90.19   vth0 = 4.7306   ksub = 0.005   theta = 1e-4 (bound)
vsat = 100 (pinned off)   lam = 1e-6 (pinned off)   eta = 4.05e-3
gds = 6.758e-5   gk = 0   rd = 24.90 mΩ   rs = 2.638 mΩ   kappa = 0.10273 1/V
plateau 1/(rd*kappa) = 391 A ; R_d(19.7 V) = 75.3 mΩ (x3.02 of zero-bias)
```
sigma_rel: kp 2.94, vth0 0.351, ksub 8.71, **theta 1.05e5 (NOT sensed)**, eta 2.56, gds 4,
rd 0.124, rs 0.123, **kappa 0.733**. Correlations: `rd↔kappa −0.997`, `rd↔rs −0.978`,
`rs↔kappa +0.971`, `theta↔rd +0.987`, `kp↔vth0 +0.967`.
→ Only the **product** `rd*kappa` (the plateau current) is well determined.

Model vs measured plateau at 19.7 V: model 250.2/254.1/254.1 A (Vgs 8/9/10) vs measured
227.0/236.7/240.2 A (spread 5.6 %). Model is ~7–10 % high there.

Ohmic check (Vgs=10): model +8.7 % / +5.0 % / +7.7 % at the three Ron anchors; fitted
rd+rs = 27.5 mΩ vs 34.5 mΩ chord. **Per §0, the rd+rs shortfall is fine — the +8.7 %
ohmic current error is a channel/knee error, not a package-R error.**

### Remaining genuine defects
1. **Vgs = 5 curve: 18 % and barely improved (18.61 → 18.07).** This is knee *position*,
   not flattening — kappa does not address it. It dominates the overall RMSE.
2. **ksub pinned at its 0.005 bound → 12 mV/dec**, far below the ~60 mV/dec physical floor.
   It is acting as a knee-shape knob, not subthreshold physics. Needs a lower bound of
   ~0.026 (or 0.06) or a separate knee parameter.
3. theta/vsat/lam all pinned off (theta σ = 1.05e5). Fine for this device, but confirm.

## 4. NOT DONE — the user's last request, still pending

**User asked: "Please do a free fit with no constraints and present the results."**
I ran out of tokens before running it. Do this next:

- Free fit, **no constraints**: all params free (including `vsat`, `lam`, `theta`, `gds`,
  `gk`, `vthk`, `kk`, `vk`, `ksub`, `kappa`) — i.e. drop the pinned `vsat=100`/`lam=1e-6`,
  drop the `gds` anchor, and per §0 do **not** constrain rd/rs to Rds(on).
- Use `DEFAULT_FREE + ("kappa",)` or literally all of `PARAM_NAMES` minus the numerical
  `vov_floor` (keep `vov_floor` fixed — it is a numerical floor, never fitted).
- Multi-start (≥5, spread kp/theta/vsat/ksub/kappa widely), relative weighting, `max_nfev`
  ≥ 6000. Report RMSE, per-curve RMSE, all params, `identifiability()`, and the Ron-anchor
  residuals. Compare against the 10.35 % / 7.45 % numbers above.
- Expect `ksub` to hit a bound and `gk/vthk/kk` to be unsensed; report σ = inf rather than
  pretending they converged.

## 5. Files touched this session

- `src/vdmos_fit/model.py` — `kappa` in `_solve_package_r` + docstring; bracket fix;
  no-progress bisection; `vdmos_id` passes `kappa=p.kappa`.
- `src/vdmos_fit/params.py` — `kappa` field + docstring (**needs §0 correction**),
  `BOUNDS["kappa"]`, comment on why kappa is not in `DEFAULT_FREE`.
- `scripts/fit_st8l65.py` — **new**, reproducible ST8L65 fit (control vs kappa-free,
  multi-start, identifiability, per-curve RMSE, plots, writes `out/st8l65_fitted_params.txt`).
  **Its CAVEATS block needs the §0 correction.**
- `scripts/prototype_drift_r.py` → moved to `out/prototype_drift_r.py.superseded`.
- Generated: `out/st8l65_driftfit.png`, `out/st8l65_fitted_params.txt`,
  `out/_bestE.json`, `out/_bestE.npy` (scratch — delete or gitignore).
- `README.md` / `src/vdmos_fit/STATUS.md` **not yet updated for kappa.**

## 6. Git state (unchanged from session 2 — still uncommitted!)

`git log`: `d399939 move status.md`, `efce72d add 2-point grid calibration`,
`cdc2c50 add curve digitizer tool`, `e19e799 fix CLi --out syntax`, `c4543ba initial commit`.
`git status --short`: modified `README.md`, `src/vdmos_fit/{STATUS,fit,model,params}.py`,
`tests/test_fit.py`; deleted root `STATUS.md` (moved to `src/vdmos_fit/STATUS.md`);
untracked `.vscode/`, `HANDOVER.md`, `data/ST8L65N044M9_IV.csv`, `tests/test_package_r.py`,
plus this session's `scripts/fit_st8l65.py`. Nothing staged.
**Two full sessions of work are uncommitted — commit before doing anything else.**

## 7. Gotchas

- **Environment:** the selected interpreter (`engineering313`, Python 3.13.5) has numpy,
  scipy, matplotlib, pandas, sympy, **pytest, ruff**. But `subprocess.run(["python", ...])`
  resolves to *base* Python which **lacks pytest/ruff** — always invoke as
  `[sys.executable, "-m", "pytest"]`. There is **no `.venv`** in the project.
- `Params` is `frozen=True` and has **no `kappa` kwarg in `vdmos_id_level1`**; `kappa` lives
  only on `Params`. Prototype-style code that splats a dict containing `kappa` into `Params(**d)`
  will `TypeError` — pop it first.
- Relative weighting + scalar `id_floor` still silently lets fits zero low-current rows
  (session-2 gotcha, still true).
- `identifiability(res, free)` must get the **same** `free` as `fit`.
- Package-R path is **Vds ≥ 0 only**; intrinsic channel model is garbage in reverse
  (divisor `1+v/vsat` singular near `Vds ≈ -vsat`; `f(-20,10) = +1038 A`). Body diode
  (STATUS problem 7) must define reverse operation.
- DEFAULT params are for a ~4 mΩ low-VT FET — **always re-anchor the start** for other
  devices (`kp ≈ 1/(Ron*Vov)`). For ST8L65 the re-anchored start is in
  `scripts/fit_st8l65.py:START`.
- With `kappa > 0`, `least_squares` can probe parameter combos that push the solver near the
  pole; the bracket fix handles it, but if a future form makes `R_d` non-monotone the
  bracket argument (and `dg <= -1`) breaks.

## 8. Recommended next steps (priority)

1. **Run the unconstrained free fit** (§4) — explicitly requested, not yet done.
2. **Apply the §0 correction** to `params.py` docstring and `fit_st8l65.py` caveats.
3. **git commit** sessions 2 + 3 (§6).
4. Attack the **Vgs = 5 knee position** (18 % residual, the dominant remaining error).
   Ideas: free `ksub` with a physical lower bound, or a separate knee/Vk-shape parameter,
   or check the digitised Vgs=5 curve again (user previously confirmed it is correct:
   Vgs=5 is barely above threshold, Vth ≈ 4.5–4.8 V — consistent with fitted 4.73 V).
5. Then STATUS 7–8: body diode (also fixes reverse-polarity garbage), leakage term for HV
   off-state, Q(V), LTspice `.subckt` export.
