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

## Model

```
I_D = kp·Vov·(v − v²/(2Vk)) / (1 + v/vsat) · (1 + λ·Vds·tanh(Vds/Vk))
      + gds·Vds
      + gk·Vds·tanh((Vgs − Vth_k)/kk)·tanh(Vds/Vk)

Vth = Vth0 − η·Vds·tanh(Vds/Vk)              DIBL, bounded
Vov = ksub·ln(1 + exp((Vgs − Vth)/ksub))     smooth subthreshold → strong inversion
Vk  = Vov/(1 + θ·Vov)                        vertical-field mobility degradation
v   = Vk·tanh(Vds/Vk)                        soft drain variable (replaces the knee switch)
```

Setting `θ = vsat = η = gds = gk = 0` makes the model a smooth approximation of
the Level-1 equations with `m_triode = 1`: the soft drain variable rounds the
hard corner, so it is exact in the ohmic limit and in deep saturation but up to
~6 % low around the knee. Quantified by
`tests/test_model.py::test_level1_limit_agreement`.

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
pytest                            # 27 tests
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

## Scope / not yet implemented

This package is the **channel current only**. The agreed architecture keeps the
rest of the subcircuit as-is: charge-based Q(V) capacitors (symbolic
differentiation, guarantees charge conservation), split series R for
channel/drift, body diode with Qrr, and the `Tj` node for Cauer-ladder
co-simulation. Validation target remains double-pulse Eon/Eoff within ~10 %,
which static fitting alone cannot deliver.

## Prior art

Hefner's physics-based MOSFET model, Infineon Level-100, Wolfspeed SiC models,
and the Vayser / "Modelling Farmers" LTspice lineage.
