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
4. **Off-state / subthreshold** Id(Vds) at Vgs = 0 → fixes the subthreshold
   floor (ksub, vth0).
5. **Kink region** (trench parts): dense Vgs spacing near the kink onset → gk, vthk, kk.

Static data cannot identify capacitance or switching parameters — those come
from Q(V) and double-pulse Eon/Eoff, handled separately.

## Datasets

### ST8L65N044M9_IV.csv

- **Device:** ST8L65N044M9 — 650 V super-junction power MOSFET, 36–44 mΩ
  Rds(on) class.
- **Source:** digitised from the datasheet's typical output characteristics
  (Figure 3, 25 °C, pulsed measurement) plus spec-table anchor points. No
  `.meta.txt` sidecar — this entry is the provenance record.
- **Figure 3 on disk:** `ST8L65N044M9_Id_Vds_Vgs.png` — the datasheet's full
  typical output curves (Vgs = 5, 6, 7, 8, 9, 10 V; Vds 0–19 V; Id up to
  ~240 A; print code GPDP160620251503OCH). Use it for the knee-shape check
  of the fitted curves.
- **Contents:** 41 points, header `index,vds,ids,vgs` (`load_csv` matches
  columns by header, so order/extras don't matter). Seven gate curves:
  Vgs = 3.7, 5, 6, 7, 8, 9, 10 V; Vds range 0–20 V.
- **Spec anchors included:** the Vgs(th) point (Vgs = Vds = 3.7 V,
  Id = 250 µA) — the only off-state point, anchoring ksub/vth0 through the
  subthreshold channel — and the
  Rds(on) points at Vgs = 10 V (0.483 V/14 A, 1 V/29 A, 1.26 V/35 A →
  34.5–36 mΩ chord).
- **Real-data oddity kept, not cleaned:** the Vgs = 5 curve has two points
  with identical Id ≈ 3.937 A at Vds = 1.127 and 5.153 V. Genuine datasheet
  data, but no C¹ model can reproduce it (β = 0 over a 4.6× voltage interval,
  then β ≈ 0.12); `scripts/fit_st8l65.py` reports fits with and without that
  curve (fit A / fit B).
- **Removed point:** (Vgs = 0, Vds = 650 V, Id = 1 µA) was digitised but
  removed — the channel model has no high-voltage leakage element yet, and no
  single conductance could give both 1 µA at 650 V (1.5 nS) and 250 µA at
  3.7 V (68 µS) anyway. Real HV leakage waits for the body diode / a
  dedicated leakage element.
- **No Rds(on)-vs-Vgs curve exists in this datasheet** (checked); the
  rd/rs-vs-channel split stays weakly identified without it.
