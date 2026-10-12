# Improved VDMOS — continuous drain-current model + parameter fitter

A compact power-MOSFET model that replaces the built-in LTspice Level-1 VDMOS
two-branch (triode/saturation) drain current with **one C¹ continuous
expression**, plus least-squares parameter extraction from static I-V data.

## Why

> **Correction (verified against LTspice 26.0.1).** An earlier version of this
> project claimed Level-1 VDMOS produces a negative-differential-resistance
> region (gds ≈ −40 S) and a several-hundred-percent gds jump at the knee. That
> was wrong, and the artifact came from a bug in this project's own reference
> formula, which put the knee at `Vds = Vov` instead of LTspice's
> `Vds = Vov / mtriode`. LTspice chooses that knee location *specifically* so the
> triode slope `m·Kp·(Vov − m·Vds)` vanishes exactly where saturation begins, making
> Level-1 **C1 continuous for any mtriode**. Pinned by
> `tests/test_model.py::test_level1_is_c1_continuous`.

The real Level-1 limitations, confirmed by `.meas` at Vgs = 5 V, Vto = 3 V, Kp = 2,
Lambda = 0:

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
| Off-state slope | undefined | gds·Vds |
| Trench kink | none | one-sided sigmoid |
| Gate-charge saturation (α collapse at high Vgs) | none | soft cap `kcap`, `m` (off by default) |

## Model

```
I_D = kp·Vov_eff·(v − v²/(2Vk)) / (1 + v/vsat) · (1 + λ·Vds·tanh(Vds/Vk))
      + gds·Vds
      + gk·Vds·tanh((Vgs − Vth_k)/kk)·tanh(Vds/Vk)

Vth     = Vth0 − η·Vds·tanh(Vds/Vk)              DIBL, bounded
Vov     = ksub·ln(1 + exp((Vgs − Vth)/ksub))     smooth subthreshold → strong inversion
Vov_eff = Vov/(1 + (Vov/kcap)^m)^(1/m)           soft gate-charge cap (kcap = 0 → off)
Vk      = Vov_eff/(1 + θ·Vov_eff)                vertical-field mobility degradation
v       = Vk·tanh(Vds/Vk)                        soft drain variable (replaces the knee switch)
```

Setting `θ = vsat = η = gds = gk = 0` makes the model a smooth approximation of
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
(variant A): with a depleted pillar the pinch-off voltage is set by the
depleted column, not by the surface charge. `kcap = 0` (default) is the exact
identity, so every pre-cap anchor stays bit-identical — pinned by
`tests/test_gate_cap.py`.

*History:* the previous mechanism was `kappa`, a drain-bias-dependent series
resistance. It is removed: `kappa` is a *voltage divider*, not a charge limit,
and the fit made rd carry 75 % of the terminal voltage already at the ohmic
Ron anchor (Vds = 0.483 V, Id = 14 A). The cap leaves rd a plain resistor and
beat kappa on the fit (6.88 % vs 7.45 % RMSE) with full Jacobian rank and
~1000× better conditioning.

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
pytest                            # 86 passed, 1 skipped (Tk digitizer test needs a display)
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
   equally, so the low-current shape (vth0, ksub, gds) is not swamped.
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
