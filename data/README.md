# data/

Measured or digitised I-V data, one CSV per device/condition.

## Format

```csv
v_gs,v_ds,i_d
4.0,0.05,0.812
4.0,1.0,15.3
...
```

One row per measured point (flat, not a grid) — digitised datasheet curves are
ragged and different curves have different lengths.

Body-diode files use their own header contract, `index,I_sd,v_sd,T` (T in
degC) — the datasheet diode spec is current-referenced and multi-temperature;
`load_diode_csv()` reads it into a `DiodeCurveSet`.

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

1. **Transfer curves** i_d(v_gs) at several v_ds, including the low-v_ds ohmic line
   → fixes k_p, v_th0, k_theta.
2. **Rds(on) vs v_gs** spec table → the strongest constraint on Ron; fit with
   `--weight abs` on these points.
3. **Output curves** i_d(v_ds) at several v_gs → fixes v_vsat, k_lambda, k_eta.
4. **Off-state / subthreshold** i_d(v_ds) at v_gs = 0 → fixes the subthreshold
   floor (v_sub, v_th0).
5. **Kink region** (trench parts): dense v_gs spacing near the kink onset → g_k, v_thk, v_kovsoft.

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
  typical output curves (v_gs = 5, 6, 7, 8, 9, 10 V; v_ds 0–19 V; i_d up to
  ~240 A; print code GPDP160620251503OCH). Use it for the knee-shape check
  of the fitted curves.
- **Contents:** 41 points, header `index,v_ds,i_d,v_gs` (`load_csv` matches
  columns by header, so order/extras don't matter). Seven gate curves:
  v_gs = 3.7, 5, 6, 7, 8, 9, 10 V; v_ds range 0–20 V.
- **Spec anchors included:** the v_gs(th) point (v_gs = v_ds = 3.7 V,
  i_d = 250 µA) — the only off-state point, anchoring v_sub/v_th0 through the
  subthreshold channel — and the
  Rds(on) points at v_gs = 10 V (0.483 V/14 A, 1 V/29 A, 1.26 V/35 A →
  34.5–36 mΩ chord).
- **Real-data oddity kept, not cleaned:** the v_gs = 5 curve has two points
  with identical i_d ≈ 3.937 A at v_ds = 1.127 and 5.153 V. Genuine datasheet
  data, but no C¹ model can reproduce it (β = 0 over a 4.6× voltage interval,
  then β ≈ 0.12); `scripts/fit_st8l65.py` reports fits with and without that
  curve (fit A / fit B).
- **Removed point:** (v_gs = 0, v_ds = 650 V, i_d = 1 µA) was digitised but
  removed — the channel model has no high-voltage leakage element yet, and no
  single conductance could give both 1 µA at 650 V (1.5 nS) and 250 µA at
  3.7 V (68 µS) anyway. Real HV leakage waits for a dedicated leakage
  element (the body diode landed in session 12 and sets a µA-scale
  Shockley floor, but not the 650 V blocking characteristic).
- **No Rds(on)-vs-v_gs curve exists in this datasheet** (checked); the
  r_d/r_s-vs-channel split stays weakly identified without it.

### ST8L65N044M9_diode.csv

- **Device:** ST8L65N044M9 — body-diode forward characteristics.
- **Source:** digitised from the datasheet body-diode forward voltage curves
  at three temperatures. No `.meta.txt` sidecar — this entry is the
  provenance record.
- **Contents:** 30 points, header `index,I_sd,v_sd,T` (`load_diode_csv`
  matches columns by header, case-insensitive; `index` ignored). Ten current
  points 5–55 A at each of T = 25, 150, −55 °C; v_sd 0.56–1.05 V.
- **Sign convention:** source-referenced (anode = source): `i_sd` flows
  source → drain, `v_sd = −v_ds` in the model's drain-referenced convention.
- **Identifiability limit to remember:** the data starts at 5 A — there is
  **no low-current tail**, so the saturation prefactors (i_s, and i_sr even
  more) are weakly determined; only the (k_n, i_s) pair, r_diode, i_kf,
  v_eg and the Tikf/Trs1/Trs2 temperature coefficients (a_ikf, a_rdiode,
  b_rdiode) reproduce
  the curves. The temperature law (Xti/Eg) IS pinned by the
  three temperatures (v_eg ≈ 1.08 eV ≈ silicon). Fit protocol:
  `scripts/fit_st8l65_diode.py`.
