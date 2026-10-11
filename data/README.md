# data/

Measured or digitised I-V data, one CSV per device/condition.

## Format

```csv
vgs,vds,ids
4.0,0.05,0.812
4.0,1.0,15.3
...
```

One row per measured point (flat, not a grid) — digitised datasheet curves are
ragged and different curves have different lengths.

## Suggested filenames

```
<part>_<condition>.csv
IPW065N04N_g175_static.csv
IMZA120R030H_t25_static.csv
```

## Digitising datasheet curves

Plotly/Origin/WebPlotDigitizer export → keep only the three columns above.
Record the source (datasheet rev, figure number, temperature) in a sidecar
`<name>.meta.txt` so fits stay traceable.

## What to capture (in priority order)

1. **Transfer curves** Id(Vgs) at several Vds, including the low-Vds ohmic line
   → fixes kp, vth0, theta.
2. **Rds(on) vs Vgs** spec table → the strongest constraint on Ron; fit with
   `--weight abs` on these points.
3. **Output curves** Id(Vds) at several Vgs → fixes vsat, lam, eta.
4. **Off-state / subthreshold** Id(Vds) at Vgs = 0 → fixes gds.
5. **Kink region** (trench parts): dense Vgs spacing near the kink onset → gk, vthk, kk.

Static data cannot identify capacitance or switching parameters — those come
from Q(V) and double-pulse Eon/Eoff, handled separately.
