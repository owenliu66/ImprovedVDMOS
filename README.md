# Improved VDMOS — continuous drain-current model + parameter fitter

A compact power-MOSFET model that replaces the built-in LTspice Level-1 VDMOS
two-branch (triode/saturation) drain current with **one C¹ continuous
expression**, plus least-squares parameter extraction from static I-V data.

## Why

LTspice's VDMOS model:

| mtriode | knee = v_ov/m | Ron conductance = m·Kp·v_ov | i_d,sat |
|---|---|---|---|
| 1 | 2.0 V | 4.00 S | 4.000 A |
| 2 | 1.0 V | 8.00 S | 4.000 A |
| 4 | 0.5 V | 16.0 S | 4.000 A |

So **mtriode is a shape knob only**: it moves the knee and Ron together, and can
never change i_d,sat = ½·Kp·Vov². Consequences:

- **Ron, knee shape and saturation level are not independently settable.** Fitting
  Rds(on) forces Ron, which forces the knee, which forces i_d,sat. At matched Ron the
  square law then over-predicts saturation current by **4.3× at v_gs = 10 V rising to
  5.3× at v_gs = 12 V** (v_ds = 20 V) — the regime that matters for buck/bridge inrush.
  Re-measured against the corrected reference: Level-1's defaults (Kp = 28,
  mtriode = 1.2, λ = 0.02) already match the improved model's Ron within 0.15% over
  v_gs = 10–12 V, and the ratio is *unchanged* by the knee correction because at
  v_ds = 20 V both the old and corrected formulas evaluate the same saturation
  branch ½·Kp·Vov²(1+λVds). Pinned by
  `tests/test_model.py::test_matched_ron_square_law_overprediction`.
- **The saturation exponent is locked at α = 2.00.** Real power devices sit at
  α ≈ 1.2–1.6 because of drift-velocity saturation.
- **No subthreshold conduction, no mobility degradation, no DIBL, no kink, no
  defined off-state slope.**

| Quantity | Level-1 | Improved |
|---|---|---|
| C1 continuity at knee | yes (by construction) | yes |
| i_d,sat ∝ v_ov^α | α = 2.00, fixed | α = 1.08 at DEFAULT (v_gs 6–12 V); tunable across 1.2–1.6 with v_vsat ≈ 3–12 V, k_theta secondary |
| Ron vs i_d,sat coupling | locked together | independent |
| Subthreshold | none (i_d = 0 below Vto) | smooth softplus, slope → 1/v_sub |
| Off-state slope | undefined | softplus subthreshold floor, gate-dependent |
| Trench kink | none | one-sided sigmoid |
| Gate-charge saturation (α collapse at high v_gs) | none | soft cap `v_ovmax`, `k_ovsharp` (off by default) |
| Knee width vs knee position | locked (mtriode moves both) | independent shape knob `k_kneesharp` (off by default) |

## Model

```
i_d  = k_p·v_ov_eff·(v_dseff − v_dseff²/(2·v_knee))/(1 + v_dseff/v_vsat)·(1 + k_lambda·v_ds·k_db)
       + g_k·v_ds·k_kink·k_db

k_db     = tanh(v_ds/v_dbsat)                            drain bias: 0 ohmic, ~1 above knee
v_th     = v_th0 − k_eta·v_ds·k_db                       DIBL, bounded
v_ov     = v_sub·ln(1 + exp((v_gs − v_th)/v_sub))        smooth subthreshold → strong inversion
v_ov_eff = v_ov/(1 + (v_ov/v_ovmax)^k_ovsharp)^(1/k_ovsharp)   soft charge cap (v_ovmax = 0 → off)
v_knee   = max(v_ov_eff/(1 + k_theta·v_ov_eff), v_kneemin)  vertical-field mobility degradation
v_dseff  = v_knee·tanh(v_ds/v_knee)   or   softmin(v_ds, v_knee, k_kneesharp)
k_kink   = 0.5·(1 + tanh((v_gs − v_thk)/v_kovsoft))      one-sided trench-kink gate
```

Setting `k_theta = v_vsat = k_eta = g_k = 0` makes the model a smooth approximation of
the Level-1 equations with `m_triode = 1`: the soft drain variable rounds the
hard corner, so it is exact in the ohmic limit and in deep saturation but up to
~6 % low around the knee. Quantified by
`tests/test_model.py::test_level1_limit_agreement`.

### Naming convention

Every model quantity is `<type prefix>_<suffix>` (source of truth: `names.xlsx`):

| prefix | type | examples |
|---|---|---|
| `v_` | voltage [V] | `v_th0`, `v_th`, `v_ov`, `v_ov_eff`, `v_ovmax`, `v_sub`, `v_knee`, `v_kneemin`, `v_dseff`, `v_vsat`, `v_dbsat`, `v_thk`, `v_kovsoft` |
| `i_` | current [A] | `i_d` |
| `r_` | resistance [Ω] | `r_d`, `r_s` |
| `g_` | conductance [S] | `g_k` — trench-kink conductance, *not* transconductance |
| `k_` | factor | `k_db`, `k_kink`, `k_ovsharp`, `k_kneesharp`, `k_p`, `k_theta`, `k_eta`, `k_lambda`, `k_agd` |
| `c_` | capacitance [F] | `c_gd`, `c_gs`, `c_ds`, `c_gd0`, `c_gdmax` (CAP_SPEC) |
| `q_` | charge [C] | `q_gd`, `q_gs`, `q_ds` — the charge-master state variables |

- Greek-letter parameters take the `k_` prefix so no name collides with a
  Python keyword: θ → `k_theta`, η → `k_eta`, λ → `k_lambda` (bare `lambda` is
  a keyword). The `k_` factors are unitless EXCEPT the industry-standard set:
  `k_p` [A/V²] and the 1/V trio `k_lambda`, `k_eta`, `k_theta`.
- `db` = drain **bias**, not drain-to-body — this model has no bulk terminal.
  `k_db = tanh(v_ds/v_dbsat)` is the drain-bias factor that gates DIBL, CLM and
  the kink term (0 in the ohmic region, ~1 above the knee).
- Sharp/soft duality: `*_sharp` names are power-law **exponents** — bigger =
  sharper transition (`k_ovsharp`, `k_kneesharp`); `*_soft` names are voltage
  **scales** — bigger = softer transition (`v_kovsoft`).
- `m_triode` is the one carve-out: `vdmos_id_level1()` is a foreign-model
  adapter for LTspice's Level-1, and keeping the model-card spelling `mtriode`
  *is* the mapping. Everything else in this repo follows the convention.
- CSV on-disk headers stay `vgs,vds,ids` — an external contract with the curve
  digitizer and the already-digitised datasheet files. `data.py` maps the
  attribute names (`CurveSet.v_gs/v_ds/i_d`) to those header spellings.

### Gate-overdrive charge cap (v_ovmax, k_ovsharp)

On a super-junction part the measured gate exponent α = d ln i_d/d ln v_ov is
**bimodal**: ~1.6–1.8 at low overdrive, collapsing to ~0.1 above v_gs ≈ 8 V
(ST8L65N044M9, 25 °C pulsed: 1.68 / 1.61 / 1.00 / 0.16 / 0.07 across
v_gs 5→6→7→8→9→10 at the top of each curve). No setting of k_theta/v_vsat/k_lambda can drop α
by an order of magnitude. The physical picture: the lightly doped N pillars
that carry the drain current hold only a **finite amount of free charge** —
once the surface channel is strongly enhanced, driving the gate harder cannot
add much more mobile charge in the pillar, so the *effective* overdrive
saturates even though v_gs keeps rising, as if one of two series channels hit
its ceiling.

The cap implements that as a soft power-law saturation of v_ov (equation
above): identity for v_ov ≪ v_ovmax (subthreshold and the Ron limit untouched),
asymptote v_ovmax for v_ov ≫ v_ovmax, and **0 ≤ dv_ov_eff/dv_ov ≤ 1 everywhere** —
which is what keeps the package-R solver's bracket invariant dg/di_d ≤ −1
intact. `k_ovsharp` sets the sharpness; the measured collapse is abrupt, so
k_ovsharp ≈ 10 rather than the k_ovsharp = 1 of a tanh. The cap is applied **before v_knee** 
with a depleted pillar the pinch-off voltage is set by the
depleted column, not by the surface charge. `v_ovmax = 0` (default) is the exact
identity, so every pre-cap anchor stays bit-identical — pinned by
`tests/test_gate_cap.py`.

The previous mechanism, `kappa` (a drain-bias-dependent series resistance), was
removed: `kappa` is a *voltage divider*, not a charge limit, and the cap beat it
on the fit (6.88 % vs 7.45 % RMSE) with full Jacobian rank and ~1000× better
conditioning.

### Drain-knee softness (k_kneesharp)

k_theta and v_vsat move the knee position **and** width together — they are ~84 %
correlated in fits because both flatten the knee while setting it, and the
tanh soft-switch ties the transition width to `v_knee` itself. `k_kneesharp` breaks that
coupling: it replaces the tanh with the **same power-law soft-min** used for
the gate cap, `v_dseff = v_ds/(1 + (v_ds/v_knee)^k_kneesharp)^(1/k_kneesharp)`. Shape-only by construction:
`dv_dseff/dv_ds(0) = 1` and `v_dseff → v_knee` hold for every `k_kneesharp`, so **Ron and the i_d,sat
plateau are untouched** while the knee width shrinks monotonically with `k_kneesharp`
(`k_kneesharp → ∞` is a hard min, `k_kneesharp = 1` the harmonic soft-min; `k_kneesharp ≈ 2.8` reproduces
the legacy tanh knee, L2 rms 0.024·v_knee). Same log-space evaluation and the same
`0 ≤ dv_dseff/dv_ds ≤ 1` derivative bound as the cap, so the package-R solver bracket
invariant survives. `k_kneesharp = 0` (default) is the exact legacy tanh — every
pre-k_kneesharp anchor stays bit-identical (`tests/test_drain_softness.py`). Caveat:
with the DEFAULT demo set (`v_vsat = 1.5 ≪ v_knee ≈ 5 V`) the velocity-saturation
divisor creates a mild gds dip near the knee *even at k_kneesharp = 0* (−1.9 S at
v_ds ≈ 6 V, absorbed by the −2.0 S tolerance in
`test_improved_model_has_no_negative_gds`); a sharper `k_kneesharp` deepens it. That is
a property of those parameters, not of the switch — refit `v_vsat` up or `k_kneesharp`
down if a device lands in that regime.

### Off state

Below threshold the softplus channel is the **only** channel current: i_d
falls gate-dependently with slope 1/v_sub down to ~1e-13 A at v_gs = 0. With
the body diode ON (i_s > 0), the diode's reverse saturation i_s(T) + i_sr(T)
becomes the off-state floor at v_ds > 0 (µA scale on ST8L65 —
Shockley-consistent, exactly what LTspice does); the softplus floor only
dominates while the diode is off. Real high-voltage leakage still needs a
dedicated element (BV/Ibv), not a fitted resistor.

### Package resistance (r_d, r_s)

`r_d`/`r_s` (bond wires, leadframe, DBC) are embedded **outside** the intrinsic
channel, exactly as LTspice's Level-1 `Rd`/`Rs` model-card parameters do —
`vdmos_id_level1()` takes the same kwargs for a fair A/B:

```
v_ds' = v_ds - (r_d + r_s)·i_d      v_gs' = v_gs - r_s·i_d
i_d   = f_intrinsic(v_ds', v_gs')     (implicit; solved per data point)
```

r_s is not a pure series drop: it **degenerates the gate** (gm_eff = gm/(1 + r_s·gm)).
The implicit equation is solved with bracketed safeguarded Newton, never
fixed-point iteration — in the ohmic region the fixed-point contraction factor
is (r_d+r_s)·k_p·v_ov ≈ 1.6 for DEFAULT, which diverges. The solver also carries a
per-point **stall guard**: with the gate cap, `r_s·i_d` can drag the internal
gate across the cap knee and make g(i_d) flat-then-cliffy, putting Newton into
a two-cycle that shrinks the bracket only ~1e-5 per iteration; if the bracket
width fails to shrink by 10 % three consecutive times, Newton is vetoed and
bisection forced until progress resumes (pinned by
`tests/test_package_r.py::test_solver_survives_newton_two_cycle`).
r_d = r_s = 0 takes the explicit fast path (bit-identical to the pre-package-R
model; pinned by `tests/test_package_r.py`). Forward polarity only (v_ds ≥ 0):
the channel-only model is singular in reverse — `vdmos_id` never routes
negative v_ds into this solver; the body diode carries it (next section).

### Body diode

The body diode is the **LTspice diode DC branch** (module `diode.py`), wired
at the device TERMINALS so it **bypasses r_d/r_s** — only channel current
flows through the package resistances, and the two implicit solves
(package-R, diode r_diode) stay independent one-dimensional problems:

```
   drain o──[r_d]──(D')══════ channel ══════(S')──[r_s]──┬─o source
         │                                               │
         └───────── diode + r_diode ─────────────────────┘   anode = source

i_d = i_ch(v_ds > 0) − i_di(v_sd = −v_ds)

i_di = i_diffuse + i_recomb
v_di = v_sd − r_diode·i_di                                (implicit)
i_diffuse = i_s(T)·[exp(v_di/(k_n·v_therm)) − 1] / sqrt(1 + i_s(T)/i_kf·[exp(…) − 1])
i_recomb  = i_sr(T)·[exp(v_di/(k_nr·v_therm)) − 1]
i_s(T)  = i_s ·(T/T_nom)^(k_xti/k_n) ·exp[(v_eg·q/k_B/k_n)·(1/T_nom − 1/T)]
          (same law for i_sr(T), divisor k_nr)
i_kf(T)    = i_kf·(1 + a_ikf·(T − T_nom))          (LTspice Tikf law)
r_diode(T) = r_diode·(1 + a_rdiode·(T − T_nom)
                     + b_rdiode·(T − T_nom)²)       (LTspice Trs1 + Trs2 laws)
```

Sign convention: `v_sd = −v_ds`, anode = source (N-FET); positive i_di flows
source → drain. The diode CSV keeps the datasheet source-referenced spelling
(`I_sd,v_sd,T`); the conversion lives in the model/fitter.

LTspice mapping: `Is→i_s  N→k_n  Isr→i_sr  Nr→k_nr  Ikf→i_kf  Rs→r_diode
Tnom→t_nom  Xti→k_xti  Eg→v_eg  Tikf→a_ikf  Trs1→a_rdiode
Trs2→b_rdiode`. Capacitance
(Cjo/Vj/M), transit time, breakdown (BV/Ibv) and noise are out of scope;
instance multiplicity (`m=`)
is a wiring-level parallel-device count, not a model parameter (the
model-card `M` is the junction grading coefficient — capacitance only), so
`M = 1` here. The Trs1+Trs2 pair is the complete LTspice Rs temperature law.

Verified against LTspice 26.0.1 (`.dc`/`.step temp`/`.meas` probes, pinned by
`tests/test_diode.py::test_ltspice_cross_validation`):

- the high-injection divisor matches LTspice to ~1e-5 relative across
  3.4 A … 375 kA;
- the saturation-current temperature law uses **N** (not Berkeley's NTT) as
  the divisor of k_xti and v_eg; the recombination branch uses **Nr**;
- LTspice's own recombination branch deviates from its documented
  `Isr·(exp(V/NrVt)−1)` by ~1 % at 0.02 V growing to ~17 % at 0.3 V
  (V-dependent, exactly linear in Isr — an implementation quirk). We
  implement the *documented* form; the cross-validation test covers the
  diffusion + Ikf + temperature terms only, and the ST8L65 data lives at
  0.56–1.05 V where i_sr is unsensed anyway;
- the TC laws are exactly `X(T) = X·(1 + TC·(T − Tnom))` (Tikf) and
  `Rs(T) = Rs·(1 + Trs1·(T − Tnom) + Trs2·(T − Tnom)²)`, referenced to
  **Tnom** (not 27 °C): a `Tikf=0.006 Trs1=0.008 Trs2=2e-5` probe matches to
  ~1e-5 relative at −55/25/100/150 °C × 0.6/0.9 V; a 27 °C-referenced
  variant misses by 0.1–0.8 %, TC-off by up to 25 %, and linear-only (no
  Trs2) by up to 10.6 % at 150 °C
  (`test_diode.py::test_ltspice_tc_laws`).

Channel for `v_ds ≤ 0` is **0**: the channel formulation is forward-polarity
only (odd k_db, singular velocity divisor). All reverse current is diode
current. Sign-symmetric channel conduction for synchronous rectification
(`sign(v_ds)·i_d(|v_ds|, v_gs)`) is a Queue item, not this model.

Sentinels: `i_s = i_sr = 0` (the default) → diode current exactly 0 →
`vdmos_id` is bit-identical to the pre-diode model (pinned by
`tests/test_diode.py`). `i_kf = 0` → high-injection divisor off (LTspice
`Ikf = inf` equivalent); the divisor applies to the forward branch only —
reverse saturation is physically unaffected by high-level injection. "Off"
is reached by leaving the parameter out of the free set, never by driving it
to the lower bound (bounds start at 1e-22 A / 1e-3 A for that reason).
`a_ikf = a_rdiode = b_rdiode = 0` (the defaults, LTspice's own defaults) →
no temperature dependence, bit-identical to the pre-TC model; the TCs act on
i_kf/r_diode only, never on the Xti/Eg law. TC bounds ((−0.005, +0.05)/°C;
b_rdiode (−1e-4, +1e-3)/°C²) keep the factor sane across −55…150 °C; where
the polynomial would extrapolate negative the factor is clamped at a tiny
floor so the solver bracket proof survives absurd extrapolation temperatures.

Extraction: `fit_diode()` on `DiodeCurveSet` (multi-temperature, per-point
temperature). The free set IS the protocol — see `scripts/fit_st8l65_diode.py`
and STATUS, "Current ST8L65 diode fit": k_n fixed at 2.0 from physics,
i_kf/v_eg/r_diode + all three TCs free, i_sr off (unsensed). **Session 15
promotion (user override):** the official record is the 7-param TC+Trs2 fit —
**3.91 %** (−55 °C 2.37, 25 °C 3.16, 150 °C 5.49): a_ikf = +0.83 %/°C,
a_rdiode = +0.48 %/°C, b_rdiode = −1.78e-5/°C². The five-knob
r_diode/i_kf/a_ikf/a_rdiode/b_rdiode block is strongly correlated
(|corr| ≤ 0.985, σ(i_kf) ≈ 15) — quote the CURVE, never the split. Baselines
live on as script probes: `tc_off_probe` 4.78 % (session 12),
`trs2_off_probe` 4.15 % (session 14).
**Trs2 diagnosis (why the gain is modest):** the stubborn 150 °C
residual is the 5 A point (+13 %): its digitized v_sd is ~12 mV off, in the
exponential region where r_diode carries only ~10 mV of the 562 mV drop — no
series-R TC (linear or quadratic) can fix it; excluding that point the
150 °C rms is ~3.8 % for both the 6- and 7-param fits. See the Trs2
promotion section in STATUS.

## Project status

Fifteen sessions, all committed to `main` (session 13 closed the body-diode
pending list; session 14 added the diode temperature coefficients
`a_ikf`/`a_rdiode` and `b_rdiode` (Trs2), and promoted the TC+Trs2 fit to the
official record). The
package provides one C¹ continuous
drain-current model (softplus subthreshold, gate-charge cap v_ovmax/k_ovsharp, mobility
degradation, velocity saturation, DIBL, gated CLM, one-sided kink) embedded with
package r_d/r_s, **plus the body diode** (LTspice diode DC branch with
high-injection compression and the Xti/Eg temperature law, terminal-wired,
bypassing r_d/r_s), plus bounded least-squares extraction with identifiability
diagnostics (dead Jacobian columns report σ = inf; `on_bounds` scales with the
bound's magnitude, not the box width). Two mechanisms were added and then
deleted as unsupported by the data: `kappa` (a voltage divider posing as a
package resistance) and a `gds·v_ds` leak (a one-point parameter that flattened
the subthreshold curve) — the off state is now the gate-dependent softplus
floor (plus the diode's reverse saturation when the diode is on). Current
**official** fits on ST8L65N044M9: channel **3.79 % relative RMSE
(fit A, 10-param with `k_kneesharp` free, full 41-point set), Jacobian rank 8/10,
cond 2.19e3**; diode **3.91 % on 30 points (−55/25/150 °C), v_eg =
1.082 eV ≈ silicon bandgap, TCs a_ikf = +0.83 %/°C / a_rdiode = +0.48 %/°C /
b_rdiode = −1.78e-5/°C²**,
145 tests green (baselines: linear-TC-only 4.15 %, TC-off 4.78 %). `k_lambda`/`v_dbsat` report σ = inf because the data
stops at v_ds = 20 V and carries no deep-saturation information — expected,
not a fit failure. The Q(V) charge model **landed** (charge-master LUTs +
c_gd splice + c_gdmax extraction from Qg + LTspice nonlinear-cap export with
live parity — see [`src/vdmos_fit/CAP_SPEC.md`](src/vdmos_fit/CAP_SPEC.md)).
Remaining gaps: HV
leakage, thermal node, full-device LTspice subckt — details and open problems in
[`src/vdmos_fit/STATUS.md`](src/vdmos_fit/STATUS.md).

## Layout

```
ImprovedVDMOS/
├── pyproject.toml          metadata, deps, pytest + ruff config
├── src/vdmos_fit/          the importable package
│   ├── model.py            vdmos_id(), vdmos_id_level1(), package-R solver
│   ├── diode.py            vdmos_diode(), LTspice diode DC branch + r_diode solver
│   ├── caps.py             charge-master C–V LUTs + c_gd splice + c_gdmax extraction
│   ├── params.py           Params dataclass, bounds, log-space packing
│   ├── data.py             CurveSet / DiodeCurveSet / CapCurveSet / QgCurveSet, CSV I/O
│   ├── fit.py              least_squares extraction + identifiability (+ fit_diode)
│   ├── reference.py        official fitted device records (channel + diode + caps)
│   ├── ltspice.py          LTspice export: charge-master nonlinear-cap subckt
│   ├── plotting.py         output/transfer/gds/overlay plots
│   ├── cli.py              `vdmos-fit` command
│   ├── STATUS.md           session log, open problems, missing architecture
│   └── CAP_SPEC.md         settled capacitance/charge-model design spec
├── tests/                  pytest suite (model regressions + fit round-trips,
│                           test_diode.py incl. LTspice cross-validation)
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
pytest                            # 145 tests green (digitizer Tk window tests self-skip headless)
vdmos-fit check                   # reference numbers + plots
python scripts/roundtrip_demo.py  # synthetic recovery of known parameters
```

Without installing, the `scripts/` entry points add `src/` to `sys.path`
themselves, so `python scripts/roundtrip_demo.py` works in a bare environment.

## Usage as a library

```python
from vdmos_fit import Params, vdmos_id, fit, load_csv

p = Params(k_p=30.0, v_th0=3.1, k_theta=0.02, v_vsat=1.2)
id_amps = vdmos_id(v_ds=20.0, v_gs=10.0, p=p)

res = fit(load_csv("data/IPW065N04N.csv"), weight="relative")
print(res.summary())
```

## Fitting strategy

1. **Relative weighting** first — treats subthreshold and high-current points
   equally, so the low-current shape (v_th0, v_sub) is not swamped.
2. **Ampere weighting** second — makes the Ron / Rds(on) corner dominate, which
   is what datasheet specs actually pin down.
3. Positive parameters are optimised in **log space** with physical bounds; this
   keeps the trust-region solver out of singular regions (k_p ≤ 0, v_sub → 0).
4. Check `identifiability(result).summary()` — parameters with order-unity sigma
   are weakly determined by static data alone. On a typical grid `v_thk` carries
   ~2000× less sensitivity than `k_p`, and `k_theta`/`v_vsat` are ~84 % correlated
   (both flatten the knee), so their individual values are meaningless even when
   the fit is excellent. Fix one from independent data, or get them from Q(V)
   and double-pulse switching instead.
5. `r_d` and `r_s` are fitted by default (datasheet Rds(on) is a 2-terminal spec
   that already contains them). Expect **r_d ↔ r_s ≈ −0.999 correlated**: static
   data pins their *sum* tightly (synthetic round trip: 0.3 % on r_d+r_s) but the
   split leans on r_s's gm-degeneration signature (individually ±15 %). Treat
   r_d+r_s as the measured quantity; split from geometry or 4-terminal data if
   you have it.

## Implementation notes

Design rationale consolidated out of the `src/vdmos_fit` code comments.

### Model limits and the Level-1 A/B reference

- Asymptotics of the unified expression: for `v_ds ≪ v_knee`,
  `i_d → k_p·v_ov_eff·v_ds·(1 + …)` (ohmic, m_triode = 1); for `v_ds ≫ v_knee`,
  `i_d → 0.5·k_p·v_ov_eff²·(0.5 + k_theta·v_ov_eff)/(1 + k_theta·v_ov_eff)² / (1 + v_knee/v_vsat)`.
- `vdmos_id_level1()` implements Level-1 **as actually implemented**: knee at
  `v_ds = v_ov/mtriode`, verified against LTspice 26.0.1 `.meas` (the C1 table
  in "Why"; min i_d over 0–6 V = 4.0 A for every mtriode — no negative-gds
  region). `m_triode` sets `Ron = 1/(m·Kp·v_ov)` and the knee position but
  cannot change `i_d,sat = 0.5·Kp·Vov²`.
- Level-1's r_d/r_s get the same external series embedding as the improved
  model (implicit solve, shared solver) for a fair A/B; pass `r_d = r_s = 0`
  for the bare equation reference.

### Trench kink: one-sided sigmoid

The kink gate is `0.5·(1 + tanh((v_gs − v_thk)/v_kovsoft))`, not bare `tanh()`. Tanh
is odd, so a bare tanh injects **negative** current below the kink onset
(about −40 mA at v_gs = 0, v_ds = 1 V for the nominal parameter set) —
non-physical and a convergence hazard in the off state. The one-sided form
is identical to tanh once v_gs is a few v_kovsoft above v_thk (tanh already saturates
at 1.0 in float) and decays smoothly to zero below.

### Charge-cap numerics

- The cap is evaluated in **log space** (`logaddexp`) so `(v_ov/v_ovmax)^k_ovsharp`
  cannot overflow: k_ovsharp ~ 10–40 is a fitted regime, not a hypothetical one.
  Its derivative is `(1 + u^k_ovsharp)^(−1/k_ovsharp − 1)` with `u = v_ov/v_ovmax`, hence
  0 ≤ dv_ov_eff/dv_ov ≤ 1 everywhere.
- "Cap off" is reached by **leaving v_ovmax out of the free set**, not by
  driving it to the lower bound: v_ovmax → 0⁺ means "hard cap at zero volts",
  the opposite of off. `v_ovmax = 0` exactly is the identity transform.

### Package-R solver

- Bracket: `i_d ∈ [0, v_ds/(r_d + r_s)]` — the terminal resistance is never
  below r_d + r_s. `g(i_d) = f(·) − i_d` is strictly decreasing
  (`dg/di_d = −(r_d+r_s)·∂f/∂v_ds − r_s·∂f/∂v_gs − 1 ≤ −1`) and `f(0, ·) = 0`
  exactly (every model term is gated by v_ds or tanh(v_ds/v_knee)), which also
  makes clamping the internal v_ds at 0 harmless when a trial point
  overshoots. Partial derivatives are central differences (the model is C1).
- The solver **raises RuntimeError** on non-convergence; it never silently
  returns an unconverged current.
- Forward polarity only: the velocity-saturation divisor `1 + v/v_vsat` is
  singular near `v_ds ≈ −v_vsat` and the numerator is not odd, so the bracket
  argument fails for v_ds < 0. `vdmos_id` never routes negative v_ds here —
  the body diode carries reverse polarity — but the guard stays for direct
  callers (and `vdmos_id_level1`, which has no diode, still fails loudly).

### Bounds and fitting-space decisions

- `k_p` is not the Level-1 k_p: it must be re-extracted when k_theta/v_vsat are added
  (28 → ~34.7 at DEFAULT).
- `v_ovmax` is excluded from `LOG_PARAMS` because 0 is a legitimate value
  ("cap off") and cannot be logged. Its lower bound (0.05 V) is
  deliberately well above 0 — between 0 and 0.05 V the cap is a hard clip
  at nearly zero volts, so a free v_ovmax must never wander there.
- `r_d`/`r_s` are fitted in linear space (the physical answer can be exactly
  zero, and datasheet package resistance never approaches 0.5 Ω on a power
  FET). The 0.5 Ω bound is a loose *global* ceiling only; per device, cap
  r_d and r_s at the datasheet Rds(on) (see `scripts/fit_st8l65.py`
  `DEVICE_BOUNDS`).
- `r_d + r_s` is always **below** the measured Rds(on) chord — the rest of
  Rds(on) is the intrinsic channel resistance, which the model computes
  separately; r_d + r_s is the lower bound on Rds(on), approached only as
  v_gs → large. With the gate cap in place there is no bias-dependent drain
  element left to distinguish "package" from "undepleted drift": the two
  are exactly degenerate (correlation +1.000 when split), so they are kept
  as one parameter.
- `v_ovmax`/`k_ovsharp` are deliberately **not** in `DEFAULT_FREE`: the cap is 0 (off)
  for the low-VT anchor device and only needed for super-junction /
  long-drift parts, where it must be freed explicitly (see
  `scripts/fit_st8l65.py`). Adding them to the default set would silently
  change every existing fit and anchor. Same rationale for `k_kneesharp`: 0 is the
  legacy tanh, it is fitted in linear space (0 cannot be logged), and its
  lower bound 0.1 keeps a free `k_kneesharp` away from the sentinel. The floor was
  relaxed 0.5 → 0.1 when the ST8L65 optimum landed at k_kneesharp ≈ 0.47, below the old
  guard: the data wants a knee softer than 0.5 allows. Do not relax it further
  — k_kneesharp → 0 collapses the drain current.
- `v_kneemin` is a numerical floor on v_knee — keep fixed, never fitted.

### Fitter diagnostics

- **`on_bounds`**: a value counts as pinned only when it sits within
  `rel_tol` of a bound measured against
  `max(|bound|, |value|, ZERO_FLOOR)` with ZERO_FLOOR = 1e-5 (model units).
  The old span-based tolerance was a false-positive machine on wide boxes:
  the deleted gds bound (1e-9, 10) S had span 10 S, so any gds below 0.01 S
  was called "pinned LOWER" — including the physically meaningful 6.8e-5 S
  (session-6 fit output). The floor keeps the test meaningful for
  zero-bounded linear-space parameters: r_s = 3.8e-9 Ω is numerically at
  zero (pinned), while r_s = 0.66 mΩ is a real fitted value (not pinned).
- **`identifiability`**: a zero-sensitivity column reports σ = `inf`, never
  0 — a bare `pinv(Jᵀ J)` returns *zero variance* for a column of zeros,
  which reads as "perfectly determined" when the truth is the exact
  opposite (e.g. `v_thk` with v_kovsoft = 0.1 V is fully saturated across a dataset
  whose lowest v_gs is 2 V, so the kink transition is never sampled and
  v_thk is unconstrained). The dead-column test is **relative** to the
  strongest Jacobian column (1e-5 of it), not absolute: on the ST8L65
  super-junction fit `v_sub`'s column measured 7.8e-7 of max — provably dead
  (pinning it at 0.007/0.02/0.04/0.06/0.09 gives identical RMSE and
  identical other parameters) — yet an absolute cutoff classified it live
  and reported a confident σ_rel = 0.038, while also manufacturing a
  spurious v_sub ↔ v_dbsat = +1.000 correlation and inflating cond by ~3 orders
  of magnitude. 1e-5 sits in the empty gap between genuinely live columns
  (~3e-4 of max: k_theta/gds/v_thk on the synthetic grid) and dead ones (~1e-6).
  Rank, cond and correlations are computed over live columns only, so a
  dead column can no longer fake a ±1.000 pair.
- **`sensitivity_report`** (pre-fit): parameters sitting at their
  non-negative lower bound (r_d, r_s at 0) get a one-sided absolute probe of
  size `rel_step`, and the reported number is the mean *relative* i_d change
  it produces. The log-space formula d ln i_d / d ln p is meaningless at
  p = 0, and a two-sided probe would request negative resistance.

### Synthetic data

`make_synthetic` noise is multiplicative (datasheet curves are read off
log-ish axes) and applied **without clipping** the values: flooring the
current itself distorts the subthreshold points by up to 2× and makes the
generator worse than the fit it is supposed to validate. Near-zero currents
are handled where they belong — in the residual denominator (`i_d_floor` in
`residuals`), not in the data.

### CLI

`--show` and `--out/-o` are registered on every subcommand with
`argparse.SUPPRESS` defaults so they parse before *or* after the
subcommand. Without SUPPRESS the subparser would reset them to its own
defaults and clobber a value already parsed before the subcommand —
argparse fills subparser defaults into the namespace unconditionally.

## Scope / not yet implemented

This package is the **channel current plus package series resistance plus the
body-diode DC branch**. The
agreed architecture keeps the rest of the subcircuit as-is: charge-based Q(V)
capacitors (symbolic differentiation, guarantees charge conservation), split
series R for channel/drift (intrinsic; the *package* r_d/r_s are now in the model
itself), body-diode recovery charge Qrr (the static diode I–V is in; the
reverse-recovery and HV-leakage elements are not), and the `Tj` node for
Cauer-ladder co-simulation.
Validation target remains double-pulse Eon/Eoff within ~10 %, which static
fitting alone cannot deliver.

## Prior art

Hefner's physics-based MOSFET model, Infineon Level-100, Wolfspeed SiC models,
and the Vayser / "Modelling Farmers" LTspice lineage.
