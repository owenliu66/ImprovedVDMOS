# Improved VDMOS — continuous drain-current model + parameter fitter

A compact power-MOSFET model that replaces the built-in LTspice Level-1 VDMOS
two-branch (triode/saturation) drain current with **one C¹ continuous
expression**, plus least-squares parameter extraction from static I-V data.

## Why

> **Level-1 is C1 continuous for any mtriode.** LTspice places the
> triode/saturation knee at `Vds = Vov / mtriode` so the triode slope
> `m·Kp·(Vov − m·Vds)` vanishes exactly where saturation begins. (An early
> version of this project reported a −40 S NDR dip and a gds jump at the knee;
> that was a bug in its own reference formula, since corrected and pinned by
> `tests/test_model.py::test_level1_is_c1_continuous`.) Verified against
> LTspice 26.0.1 `.meas` at Vgs = 5 V, Vto = 3 V, Kp = 2, Lambda = 0:

| mtriode | knee = Vov/m | Ron conductance = m·Kp·Vov | Id,sat |
|---|---|---|---|
| 1 | 2.0 V | 4.00 S | 4.000 A |
| 2 | 1.0 V | 8.00 S | 4.000 A |
| 4 | 0.5 V | 16.0 S | 4.000 A |

So **mtriode is a shape knob only**: it moves the knee and Ron together, and can
never change Id,sat = ½·Kp·Vov². Consequences:

- **Ron, knee shape and saturation level are not independently settable.** Fitting
  Rds(on) forces Ron, which forces the knee, which forces Id,sat. At matched Ron the
  square law then over-predicts saturation current by **4.3× at Vgs = 10 V rising to
  5.3× at Vgs = 12 V** (Vds = 20 V) — the regime that matters for buck/bridge inrush.
  Re-measured against the corrected reference: Level-1's defaults (Kp = 28,
  mtriode = 1.2, λ = 0.02) already match the improved model's Ron within 0.15% over
  Vgs = 10–12 V, and the ratio is *unchanged* by the knee correction because at
  Vds = 20 V both the old and corrected formulas evaluate the same saturation
  branch ½·Kp·Vov²(1+λVds). Pinned by
  `tests/test_model.py::test_matched_ron_square_law_overprediction`.
- **The saturation exponent is locked at α = 2.00.** Real power devices sit at
  α ≈ 1.2–1.6 because of drift-velocity saturation.
- **No subthreshold conduction, no mobility degradation, no DIBL, no kink, no
  defined off-state slope.**

| Quantity | Level-1 | Improved |
|---|---|---|
| C1 continuity at knee | yes (by construction) | yes |
| Id,sat ∝ Vov^α | α = 2.00, fixed | α = 1.08 at DEFAULT (Vgs 6–12 V); tunable across 1.2–1.6 with vsat ≈ 3–12 V, θ secondary |
| Ron vs Id,sat coupling | locked together | independent |
| Subthreshold | none (Id = 0 below Vto) | smooth softplus, slope → 1/ksub |
| Off-state slope | undefined | softplus subthreshold floor, gate-dependent |
| Trench kink | none | one-sided sigmoid |
| Gate-charge saturation (α collapse at high Vgs) | none | soft cap `kcap`, `m` (off by default) |

## Model

```
I_D = kp·Vov_eff·(v − v²/(2Vk)) / (1 + v/vsat) · (1 + λ·Vds·tanh(Vds/Vk))
      + gk·Vds·tanh((Vgs − Vth_k)/kk)·tanh(Vds/Vk)

Vth     = Vth0 − η·Vds·tanh(Vds/Vk)              DIBL, bounded
Vov     = ksub·ln(1 + exp((Vgs − Vth)/ksub))     smooth subthreshold → strong inversion
Vov_eff = Vov/(1 + (Vov/kcap)^m)^(1/m)           soft gate-charge cap (kcap = 0 → off)
Vk      = Vov_eff/(1 + θ·Vov_eff)                vertical-field mobility degradation
v       = Vk·tanh(Vds/Vk)                        soft drain variable (replaces the knee switch)
```

Setting `θ = vsat = η = gk = 0` makes the model a smooth approximation of
the Level-1 equations with `m_triode = 1`: the soft drain variable rounds the
hard corner, so it is exact in the ohmic limit and in deep saturation but up to
~6 % low around the knee. Quantified by
`tests/test_model.py::test_level1_limit_agreement`.

### Gate-overdrive charge cap (kcap, m)

On a super-junction part the measured gate exponent α = d ln Id/d ln Vov is
**bimodal**: ~1.6–1.8 at low overdrive, collapsing to ~0.1 above Vgs ≈ 8 V
(ST8L65N044M9, 25 °C pulsed: 1.68 / 1.61 / 1.00 / 0.16 / 0.07 across
Vgs 5→6→7→8→9→10 at the top of each curve). No setting of θ/vsat/λ can drop α
by an order of magnitude. The physical picture: the lightly doped N pillars
that carry the drain current hold only a **finite amount of free charge** —
once the surface channel is strongly enhanced, driving the gate harder cannot
add much more mobile charge in the pillar, so the *effective* overdrive
saturates even though Vgs keeps rising, as if one of two series channels hit
its ceiling.

The cap implements that as a soft power-law saturation of Vov (equation
above): identity for Vov ≪ kcap (subthreshold and the Ron limit untouched),
asymptote kcap for Vov ≫ kcap, and **0 ≤ dVov_eff/dVov ≤ 1 everywhere** —
which is what keeps the package-R solver's bracket invariant dg/dId ≤ −1
intact. `m` sets the sharpness; the measured collapse is abrupt, so
m ≈ 10 rather than the m = 1 of a tanh. The cap is applied **before Vk** 
with a depleted pillar the pinch-off voltage is set by the
depleted column, not by the surface charge. `kcap = 0` (default) is the exact
identity, so every pre-cap anchor stays bit-identical — pinned by
`tests/test_gate_cap.py`.

The previous mechanism, `kappa` (a drain-bias-dependent series resistance), was
removed: `kappa` is a *voltage divider*, not a charge limit, and the cap beat it
on the fit (6.88 % vs 7.45 % RMSE) with full Jacobian rank and ~1000× better
conditioning.

### Off state

Below threshold the softplus channel is the **only** current: Id falls
gate-dependently with slope 1/ksub down to ~1e-13 A at Vgs = 0. Real
high-voltage leakage will come with the body diode / leakage element, not a
fitted resistor.

### Package resistance (rd, rs)

`rd`/`rs` (bond wires, leadframe, DBC) are embedded **outside** the intrinsic
channel, exactly as LTspice's Level-1 `Rd`/`Rs` model-card parameters do —
`vdmos_id_level1()` takes the same kwargs for a fair A/B:

```
Vds' = Vds - (rd + rs)·Id      Vgs' = Vgs - rs·Id
Id   = f_intrinsic(Vds', Vgs')     (implicit; solved per data point)
```

rs is not a pure series drop: it **degenerates the gate** (gm_eff = gm/(1 + rs·gm)).
The implicit equation is solved with bracketed safeguarded Newton, never
fixed-point iteration — in the ohmic region the fixed-point contraction factor
is (rd+rs)·kp·Vov ≈ 1.6 for DEFAULT, which diverges. The solver also carries a
per-point **stall guard**: with the gate cap, `rs·Id` can drag the internal
gate across the cap knee and make g(Id) flat-then-cliffy, putting Newton into
a two-cycle that shrinks the bracket only ~1e-5 per iteration; if the bracket
width fails to shrink by 10 % three consecutive times, Newton is vetoed and
bisection forced until progress resumes (pinned by
`tests/test_package_r.py::test_solver_survives_newton_two_cycle`).
rd = rs = 0 takes the explicit fast path (bit-identical to the pre-package-R
model; pinned by `tests/test_package_r.py`). Forward polarity only (Vds ≥ 0):
the channel-only model is singular in reverse until the body diode exists.

## Project status

Eight sessions, all committed to `main`. The package provides one C¹ continuous
drain-current model (softplus subthreshold, gate-charge cap kcap/m, mobility
degradation, velocity saturation, DIBL, gated CLM, one-sided kink) embedded with
package rd/rs, plus bounded least-squares extraction with identifiability
diagnostics (dead Jacobian columns report σ = inf; `on_bounds` scales with the
bound's magnitude, not the box width). Two mechanisms were added and then
deleted as unsupported by the data: `kappa` (a voltage divider posing as a
package resistance) and a `gds·Vds` leak (a one-point parameter that flattened
the subthreshold curve) — the off state is now the gate-dependent softplus
floor. Current fit on ST8L65N044M9: **5.54 % relative RMSE (fit B), Jacobian
rank 9/9, cond 3.75e4**, 89 tests green. Remaining gaps: Q(V), body diode +
HV leakage, thermal node, LTspice export — details and open problems in
[`src/vdmos_fit/STATUS.md`](src/vdmos_fit/STATUS.md).

## Layout

```
ImprovedVDMOS/
├── pyproject.toml          metadata, deps, pytest + ruff config
├── src/vdmos_fit/          the importable package
│   ├── model.py            vdmos_id(), vdmos_id_level1()
│   ├── params.py           Params dataclass, bounds, log-space packing
│   ├── data.py             CurveSet, CSV I/O, synthetic data
│   ├── fit.py              least_squares extraction + identifiability
│   ├── plotting.py         output/transfer/gds/overlay plots
│   └── cli.py              `vdmos-fit` command
├── tests/                  pytest suite (model regressions + fit round-trip)
├── scripts/                thin entry points (no logic here)
├── data/                   measured/digitised CSVs (see data/README.md)
├── notebooks/              exploratory work
└── out/                    generated plots + fitted_params.txt (gitignored)
```

## Setup

```bash
cd Projects/ImprovedVDMOS
python -m venv .venv
.venv\Scripts\activate            # Windows
pip install -e ".[dev]"
pytest                            # 89 tests green (digitizer Tk window tests self-skip headless)
vdmos-fit check                   # reference numbers + plots
python scripts/roundtrip_demo.py  # synthetic recovery of known parameters
```

Without installing, the `scripts/` entry points add `src/` to `sys.path`
themselves, so `python scripts/roundtrip_demo.py` works in a bare environment.

## Usage as a library

```python
from vdmos_fit import Params, vdmos_id, fit, load_csv

p = Params(kp=30.0, vth0=3.1, theta=0.02, vsat=1.2)
id_amps = vdmos_id(vds=20.0, vgs=10.0, p=p)

res = fit(load_csv("data/IPW065N04N.csv"), weight="relative")
print(res.summary())
```

## Fitting strategy

1. **Relative weighting** first — treats subthreshold and high-current points
   equally, so the low-current shape (vth0, ksub) is not swamped.
2. **Ampere weighting** second — makes the Ron / Rds(on) corner dominate, which
   is what datasheet specs actually pin down.
3. Positive parameters are optimised in **log space** with physical bounds; this
   keeps the trust-region solver out of singular regions (kp ≤ 0, ksub → 0).
4. Check `identifiability(result).summary()` — parameters with order-unity sigma
   are weakly determined by static data alone. On a typical grid `vthk` carries
   ~2000× less sensitivity than `kp`, and `theta`/`vsat` are ~84 % correlated
   (both flatten the knee), so their individual values are meaningless even when
   the fit is excellent. Fix one from independent data, or get them from Q(V)
   and double-pulse switching instead.
5. `rd` and `rs` are fitted by default (datasheet Rds(on) is a 2-terminal spec
   that already contains them). Expect **rd ↔ rs ≈ −0.999 correlated**: static
   data pins their *sum* tightly (synthetic round trip: 0.3 % on rd+rs) but the
   split leans on rs's gm-degeneration signature (individually ±15 %). Treat
   rd+rs as the measured quantity; split from geometry or 4-terminal data if
   you have it.

## Implementation notes

Design rationale consolidated out of the `src/vdmos_fit` code comments.

### Model limits and the Level-1 A/B reference

- Asymptotics of the unified expression: for `Vds ≪ Vk`,
  `I_D → kp·Vov_eff·Vds·(1 + …)` (ohmic, m_triode = 1); for `Vds ≫ Vk`,
  `I_D → 0.5·kp·Vov_eff²·(0.5 + θ·Vov_eff)/(1 + θ·Vov_eff)² / (1 + Vk/vsat)`.
- `vdmos_id_level1()` implements Level-1 **as actually implemented**: knee at
  `Vds = Vov/mtriode`, verified against LTspice 26.0.1 `.meas` (the C1 table
  in "Why"; min Id over 0–6 V = 4.0 A for every mtriode — no negative-gds
  region). `m_triode` sets `Ron = 1/(m·Kp·Vov)` and the knee position but
  cannot change `Id,sat = 0.5·Kp·Vov²`.
- Level-1's rd/rs get the same external series embedding as the improved
  model (implicit solve, shared solver) for a fair A/B; pass `rd = rs = 0`
  for the bare equation reference.

### Trench kink: one-sided sigmoid

The kink gate is `0.5·(1 + tanh((Vgs − Vthk)/kk))`, not bare `tanh()`. Tanh
is odd, so a bare tanh injects **negative** current below the kink onset
(about −40 mA at Vgs = 0, Vds = 1 V for the nominal parameter set) —
non-physical and a convergence hazard in the off state. The one-sided form
is identical to tanh once Vgs is a few kk above Vthk (tanh already saturates
at 1.0 in float) and decays smoothly to zero below.

### Gate-cap numerics

- The cap is evaluated in **log space** (`logaddexp`) so `(Vov/kcap)^m`
  cannot overflow: m ~ 10–40 is a fitted regime, not a hypothetical one.
  Its derivative is `(1 + u^m)^(−1/m − 1)` with `u = Vov/kcap`, hence
  0 ≤ dVov_eff/dVov ≤ 1 everywhere.
- "Cap off" is reached by **leaving kcap out of the free set**, not by
  driving it to the lower bound: kcap → 0⁺ means "hard cap at zero volts",
  the opposite of off. `kcap = 0` exactly is the identity transform.

### Package-R solver

- Bracket: `Id ∈ [0, Vds/(rd + rs)]` — the terminal resistance is never
  below rd + rs. `g(Id) = f(·) − Id` is strictly decreasing
  (`dg/dId = −(rd+rs)·∂f/∂Vds − rs·∂f/∂Vgs − 1 ≤ −1`) and `f(0, ·) = 0`
  exactly (every model term is gated by Vds or tanh(Vds/Vk)), which also
  makes clamping the internal Vds at 0 harmless when a trial point
  overshoots. Partial derivatives are central differences (the model is C1).
- The solver **raises RuntimeError** on non-convergence; it never silently
  returns an unconverged current.
- Forward polarity only: the velocity-saturation divisor `1 + v/vsat` is
  singular near `Vds ≈ −vsat` and the numerator is not odd, so the bracket
  argument fails for Vds < 0. Reverse operation needs the body diode
  (STATUS problem 7), not this solver.

### Bounds and fitting-space decisions

- `kp` is not the Level-1 kp: it must be re-extracted when θ/vsat are added
  (28 → ~34.7 at DEFAULT).
- `kcap` is excluded from `LOG_PARAMS` because 0 is a legitimate value
  ("cap off") and cannot be logged. Its lower bound (0.05 V) is
  deliberately well above 0 — between 0 and 0.05 V the cap is a hard clip
  at nearly zero volts, so a free kcap must never wander there.
- `rd`/`rs` are fitted in linear space (the physical answer can be exactly
  zero, and datasheet package resistance never approaches 0.5 Ω on a power
  FET). The 0.5 Ω bound is a loose *global* ceiling only; per device, cap
  rd and rs at the datasheet Rds(on) (see `scripts/fit_st8l65.py`
  `DEVICE_BOUNDS`).
- `rd + rs` is always **below** the measured Rds(on) chord — the rest of
  Rds(on) is the intrinsic channel resistance, which the model computes
  separately; rd + rs is the lower bound on Rds(on), approached only as
  Vgs → large. With the gate cap in place there is no bias-dependent drain
  element left to distinguish "package" from "undepleted drift": the two
  are exactly degenerate (correlation +1.000 when split), so they are kept
  as one parameter.
- `kcap`/`m` are deliberately **not** in `DEFAULT_FREE`: the cap is 0 (off)
  for the low-VT anchor device and only needed for super-junction /
  long-drift parts, where it must be freed explicitly (see
  `scripts/fit_st8l65.py`). Adding them to the default set would silently
  change every existing fit and anchor.
- `vov_floor` is a numerical floor on Vk — keep fixed, never fitted.

### Fitter diagnostics

- **`on_bounds`**: a value counts as pinned only when it sits within
  `rel_tol` of a bound measured against
  `max(|bound|, |value|, ZERO_FLOOR)` with ZERO_FLOOR = 1e-5 (model units).
  The old span-based tolerance was a false-positive machine on wide boxes:
  the deleted gds bound (1e-9, 10) S had span 10 S, so any gds below 0.01 S
  was called "pinned LOWER" — including the physically meaningful 6.8e-5 S
  (session-6 fit output). The floor keeps the test meaningful for
  zero-bounded linear-space parameters: rs = 3.8e-9 Ω is numerically at
  zero (pinned), while rs = 0.66 mΩ is a real fitted value (not pinned).
- **`identifiability`**: a zero-sensitivity column reports σ = `inf`, never
  0 — a bare `pinv(Jᵀ J)` returns *zero variance* for a column of zeros,
  which reads as "perfectly determined" when the truth is the exact
  opposite (e.g. `vthk` with kk = 0.1 V is fully saturated across a dataset
  whose lowest Vgs is 2 V, so the kink transition is never sampled and
  vthk is unconstrained). The dead-column test is **relative** to the
  strongest Jacobian column (1e-5 of it), not absolute: on the ST8L65
  super-junction fit `ksub`'s column measured 7.8e-7 of max — provably dead
  (pinning it at 0.007/0.02/0.04/0.06/0.09 gives identical RMSE and
  identical other parameters) — yet an absolute cutoff classified it live
  and reported a confident σ_rel = 0.038, while also manufacturing a
  spurious ksub ↔ vk = +1.000 correlation and inflating cond by ~3 orders
  of magnitude. 1e-5 sits in the empty gap between genuinely live columns
  (~3e-4 of max: θ/gds/vthk on the synthetic grid) and dead ones (~1e-6).
  Rank, cond and correlations are computed over live columns only, so a
  dead column can no longer fake a ±1.000 pair.
- **`sensitivity_report`** (pre-fit): parameters sitting at their
  non-negative lower bound (rd, rs at 0) get a one-sided absolute probe of
  size `rel_step`, and the reported number is the mean *relative* Id change
  it produces. The log-space formula d ln Id / d ln p is meaningless at
  p = 0, and a two-sided probe would request negative resistance.

### Synthetic data

`make_synthetic` noise is multiplicative (datasheet curves are read off
log-ish axes) and applied **without clipping** the values: flooring the
current itself distorts the subthreshold points by up to 2× and makes the
generator worse than the fit it is supposed to validate. Near-zero currents
are handled where they belong — in the residual denominator (`id_floor` in
`residuals`), not in the data.

### CLI

`--show` and `--out/-o` are registered on every subcommand with
`argparse.SUPPRESS` defaults so they parse before *or* after the
subcommand. Without SUPPRESS the subparser would reset them to its own
defaults and clobber a value already parsed before the subcommand —
argparse fills subparser defaults into the namespace unconditionally.

## Scope / not yet implemented

This package is the **channel current plus package series resistance**. The
agreed architecture keeps the rest of the subcircuit as-is: charge-based Q(V)
capacitors (symbolic differentiation, guarantees charge conservation), split
series R for channel/drift (intrinsic; the *package* rd/rs are now in the model
itself), body diode with Qrr, and the `Tj` node for Cauer-ladder co-simulation.
Validation target remains double-pulse Eon/Eoff within ~10 %, which static
fitting alone cannot deliver.

## Prior art

Hefner's physics-based MOSFET model, Infineon Level-100, Wolfspeed SiC models,
and the Vayser / "Modelling Farmers" LTspice lineage.
