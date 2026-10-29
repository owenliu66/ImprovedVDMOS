# data/

Measured or digitised I-V data, one CSV per device/condition.

## Format

I-V files: `load_csv()` accepts **both** external header spellings, matched
per column, case-insensitively, order-free (extra columns like `index` are
ignored):

```csv
vgs,vds,ids            <- legacy SPICE-style (what save_csv() writes)
4.0,0.05,0.812
...

index,V_ds,I_d,V_gs    <- digitizer/model naming (what curve_digitizer exports)
1,0.05,0.812,4.0
...
```

One row per measured point (flat, not a grid) — digitised datasheet curves are
ragged and different curves have different lengths.

Body-diode files use their own header contract, `index,I_sd,v_sd,T` (T in
degC) — the datasheet diode spec is current-referenced and multi-temperature;
`load_diode_csv()` reads it into a `DiodeCurveSet`.

Capacitance files use the header contract `index,V_ds,Cap (pF),Type` with
`Type ∈ {iss,oss,rss}` — one row per digitized point, all three branches in
one file (the delta decomposition c_gd = Crss, c_gs = Ciss − Crss,
c_ds = Coss − Crss is the model contract, see `src/vdmos_fit/CAP_SPEC.md`).
Gate-charge files use `index,Q_g(nC),V_gs,condition` (condition =
`I_d=..A, V_ds=..V`, parsed per row so one file can hold several test
curves; rows are polyline corner points).

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
- **Contents:** 41 points, header `index,Vds,Ids,Vgs` (`load_csv` matches
  columns by header, so order/extras/case don't matter). Seven gate curves:
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

### IMCQ120R004M2H_IV_25C.csv / IMCQ120R004M2H_IV_175C.csv

- **Device / source:** Infineon IMCQ120R004M2H CoolSiC 1200 V (Kelvin
  source), typical output curves digitized from the datasheet (rev 1.10)
  at 25 °C and 175 °C. Official model: `reference.py::IMCQ120_JOINT`
  (session-21 joint multi-temperature fit — protocol
  `scripts/fit_imcq120_joint.py`, report
  `out/IMCQ120_JOINT_FIT_REPORT.md`, pinned by
  `tests/test_imcq120_official.py`); the session-18 per-temperature report
  `out/IMCQ120_TEMPERATURE_REPORT.md` was deleted in session 21 (durable
  findings carried in the joint report §8) and
  `scripts/fit_imcq120_temperature.py` / `fit_imcq120_gain_probe.py` stay as
  historical baselines.
- **Contents:** 78 / 73 points, header `index,V_ds,I_d,V_gs` — the
  curve_digitizer export spelling, read natively by `load_csv()` since
  session 20 (the old per-script `load_imcq()` workaround is gone).
  Gate curves v_gs = 4.2…20 V; v_ds up to ~21 V (curves exit the datasheet
  1600 A plot frame — that is plot truncation, not a current limit).
- **Related IMCQ spec CSVs** (`Rdson_Tj_185.2A`, `Vth_Tj`) use the same
  digitizer export style and are consumed directly by the joint-fit script
  and the pinned record test (not through `load_csv`); they pin `r_d(T)`
  (quadratic TC law) and `v_th0(T)` in the official record.

### IMCQ120R004M2H_Transfer.csv

- **Device / source:** Infineon IMCQ120R004M2H, transfer characteristic
  i_d vs V_GS at V_DS = 20 V, T_j = 25/175/200 °C — digitized from the
  datasheet (rev 1.10, p. 9) with `curve_digitizer` for the session-21
  joint multi-temperature fit (`scripts/fit_imcq120_joint.py`,
  `out/IMCQ120_JOINT_FIT_REPORT.md`). No `.meta.txt` sidecar — this entry
  is the provenance record.
- **Contents:** 23 rows, header `index,V_gs,I_d,condition` (condition =
  `"V_ds=20V,T_j=25C"` style, parsed by the script-local loader — same
  pattern as the Qg `condition` column; `load_csv()` does NOT read this
  file, there is no `V_ds` column).
- **Digitization limit kept, not cleaned:** the 175 °C and 200 °C curves
  are indistinguishable below V_gs = 12 V on the datasheet figure, so no
  points were digitized there; the two curves cross at 12.025 V / 999.1 A
  and that crossing point is carried once per curve (rows 22–23 are the
  same point for T_j = 175 and 200 °C). Residuals exist only where data
  exists.
- **Role in the fit:** this is the data that senses the high-v_gs shape
  (it fixed the 3650 A gain-fit overshoot at v_gs = 18 / V_DS = 20) and
  the third temperature that pins the `r_d(T)` law together with the
  Rds(on)-Tj chords. The joint fit that consumed it is now the **official
  IMCQ record** (`reference.py::IMCQ120_JOINT`, pinned by
  `tests/test_imcq120_official.py`). The figure itself is a fitted input —
  after the joint fit quote the model-vs-digitized residuals (±3.7 %), not
  the figure read.

### ST8L65N044M9_Cap_Vds.csv / IMCQ120R004M2H_Cap_Vds.csv / ISC040N10NM8_Cap_Vds.csv

- **Devices / source:** digitized typical-capacitance curves (C = f(V_DS),
  V_GS = 0) from the three datasheets — ST8L65N044M9 super-junction 650 V
  (Figure 7, f = 250 kHz, print code GPDP160620251505CVR); IMCQ120R004M2H
  CoolSiC 1200 V (f = 100 kHz); ISC040N10NM8 trench 100 V (Diagram 11,
  f = 1 MHz). No `.meta.txt` sidecars — these entries are the provenance
  record. All 25 °C typical; no temperature dependence available.
- **Contents:** header `index,V_ds,Cap (pF),Type`, Type ∈ {iss,oss,rss};
  46 / 36 / 28 points. Digitized with `curve_digitizer` (log-log axis
  calibration).
- **Model contract:** see `src/vdmos_fit/CAP_SPEC.md` (settled design:
  charge-master LUTs, PCHIP on ln C vs linear v, c_gd(v_gd) splice at
  v_gd = 0 with exponential-saturation branch).
- **Data oddities kept, not cleaned:** (a) ST8L65 `iss` has four identical
  points (4869.28716 pF — gridline snap; harmless for PCHIP, clean on next
  re-digitization); (b) the SiC plot's x-axis starts at 1 V, so the origin
  anchor (C0, slope) is an extrapolation (±10–15 %); (c) the SJ `rss` branch
  is genuinely non-monotonic (rises again above ~200 V) — real SJ physics,
  not a digitization error.
- **Cross-checks on record:** trench `rss` fits the graded-junction law
  C₀(1+v/V_bi)^−M with M = 0.721, V_bi = 6.41 V (max residual 10.6 %);
  the same law fails on SiC/SJ (table approach confirmed). Near-origin
  quadratic anchors (first-4-point convention) are tabulated in CAP_SPEC.md.

### ST8L65N044M9_Vgs_Qg.csv / IMCQ120R004M2H_Vgs_Qg.csv / ISC040N10NM8_Vgs_Qg.csv

- **Devices / source:** digitized gate-charge curves (V_GS vs Q_g at fixed
  I_d, V_DS) from the same three datasheets. Polyline corner points only
  (the datasheet specs are piecewise-linear): 4 / 4 / 10 rows (ISC040
  carries three test conditions: V_DS = 20/50/80 V at I_d = 12.5 A).
  ST8L65: I_d = 29 A, V_DS = 400 V. IMCQ120: I_d = 185.2 A, V_DS = 800 V
  (labeled "turn-on pulse" in the datasheet — a dynamic test, so the
  quasi-static decomposition is approximate; noted in CAP_SPEC).
- **Model contract:** c_gdmax extraction + plateau validation via
  `extract_cmax()` / `scripts/fit_cmax.py`; records `CAP_*` in
  reference.py. Plateau detection = max-dQ/dV segment; the (0,0) origin
  corner is prepended automatically when the digitized curve starts at the
  plateau.
