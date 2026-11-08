# Handover — ImprovedVDMOS session 28: ALL SESSION-27 PENDING DONE (diode record LANDED + subckt D1 = DC set + pinned tests)

## ⚠️ START HERE — session 28 (this session, all verified)

DONE:
1. **`ISC040_DIODE` record LANDED** in `reference.py` (+ `ISC040_DIODE_FREE` /
   `_RMSE_REL` / `_GROUP_RMSE_REL`), exported from `__init__.py`. Reproducible
   protocol `scripts/fit_isc040_diode.py` (4 multi-starts, local `T_j (C)`
   loader). Probes reproduce the session-27 findings: k_n=2.0 pinned → 35.2 %,
   TC-off → 3.31 %, freeing b_rdiode/a_ikf or i_sr buys nothing (0.514 %).
   ALSO FIXED: session 27 had left a DUPLICATE `ISC040_CHANNEL`/
   `IscChannelRecord` block in reference.py (the second silently overrode the
   first); deduped keeping the full-precision values from
   out/isc040_channel_params.json. Constants renamed
   `ISC040_CHANNEL_RMSE25_REL` / `ISC040_CHANNEL_RMSE_TC_REL`.
2. **Full subckt D1 branch UPDATED to the DC-fit set.** `gen_full_subckt.py`
   now emits the record values on the Xrec instance line and the model card
   `.model MyDiode D(Is=Is,Rs=0,N=N,Ikf=Ikf,Eg=Eg,Xti=Xti,Tnom=Tnom)` plus
   `R1 N001 Ki {Rs*(1+Trs1*(temp-Tnom))}`. RECONCILIATION NOTE: the
   session-27 wording ("D1 Rs=0.96m ... R1 = r_diode is fine") would
   DOUBLE-COUNT the series R (~1.92 mΩ extra → ~0.7 V error at 367 A); the
   implemented split carries r_diode(T) in R1 (Trs1 law, Tnom-referenced,
   verified to ~1e-5 by test_diode::test_ltspice_tc_laws) with Rs=0 on D1 —
   exactly ONCE — and the branch is now temperature-complete (the `temp`
   reserved variable works in R value expressions; verified).
   `out/ISC040N10NM8_full.subckt` regenerated (35.1 KB). Recovery knobs still
   from `ISC040_RECOV` (that record UNCHANGED — it documents the extraction-
   run bench conditions).
3. **BENCHES** (all in out/, run via `"C:\Program Files\ADI\LTspice\LTspice.exe"
   -b -w file.cir`):
   - `parity_isc040_full_meas.cir` (NEW .meas variant of the parity bench,
     reltol 1e-8, S/K grounded): worst 333 ppm vs the RAW record — fully
     explained by the documented 1 µΩ Rs floor (it carries the full channel
     current: I·1µΩ gate degeneration ~1 mV at 1 kA) plus .6g netlist
     formatting; vs the NETLIST-CONSTANT Python model (floor + formatted
     values) worst **0.05 ppm**. The session-27 "<0.001 %" claim therefore
     holds only against a floor-aware comparison.
   - `diode_fwd_25.cir` / `diode_fwd_175.cir` (NEW): .dc LIST at the exact
     CSV V_sd points (no FIND interpolation error), gate 0 V. LTspice vs
     netlist-constant Python: **≤27 ppm at BOTH 25 and 175 °C** (residual =
     LTspice's internal physical constants, same level as the test_diode
     cross-validation); LTspice vs CSV RMSE 0.450 %/0.578 % = the fit's own
     residual (0.445/0.575) — the subckt branch adds nothing on top. The
     175 °C run confirms the Eg/Xti/Tnom card TC + R1 Trs1 law end-to-end.
   - `recov_d1_check.cir` (NEW, clean-ASCII rewrite of recov_r1_sens.cir
     which carries stray 0xA7/0xB5 bytes): `.step param d1 list 0 1` switches
     the D1 card between the OLD bench set (Is=1u N=1.4, Ikf=1e9) and the NEW
     DC set; R1=0.961328m both. didt=1000 A/µs:
     q_rri 248.171 → 248.059 nC, t_rri 24.684 → 24.702 ns,
     i_fi_min 18.436 A — recovery INSENSITIVE to the D1 Is/N/Ikf swap
     (≤0.075 %), and the record's achieved 248.20/24.687 is reproduced at
     the new R1. Is/N were untested in recovery before; now verified.
4. **PINNED TESTS LANDED**: `tests/test_isc040_official.py` — 11 tests, all
   pass: channel 25 C 132-pt RMSE reproduction + TC (chord + Vth-VOLTAGE
   residuals exactly as fit_tc), acceptance pins (Ron chords 3.35/3.56/3.89
   mΩ @Vgs 15/10/8, Ron@Vgs15.7-T 2.12/3.34/5.20/7.76 mΩ @−55/25/100/175,
   Vth@520µA 3.52/3.06/2.50/2.07 V — data comparisons at interpolated reads,
   since the Rdson_Vgs CSV Vgs points are digitized ~15.7/10.6/8.8 not round
   numbers), TC-law invariants (r_d ×2.293 @175 C, s_vth −7.486 mV/K, v_sub
   T-invariant, k_p falls with T), transfer monotonicity, record sanity;
   diode RMSE reproduction + interior/fixed-structure pins (k_n ~1.05 NOT
   2.0, i_sr=0, a_ikf=b_rdiode=0) + the i_s↔k_n pair anchor
   vdmos_diode(1 V, 25 C) = 152.6 A; plus a structure pin that the GENERATED
   subckt's D1 branch matches the record (skips when out/ absent).
5. **FULL SUITE + RUFF**: pytest 199 tests → **1 failure, PRE-EXISTING and
   unrelated**: `test_caps.py::test_ltspice_caps_parity` c1 (got
   −0.00734197 vs pin −0.00736535, 0.32 % against the 0.2 % pin). It fails
   IDENTICALLY at HEAD 6004afa with this session's changes stashed. Root
   cause: the harness samples FIND currents at digitized-interval midpoints,
   which coincide with one of the 12 fine stair breakpoints (knife-edge);
   LTspice 26.1.1's adaptive step sequence lands the sample on the other
   stair segment than the version the pin was written on. A fixed-1 ns-step
   experiment moved the failure from c1 to a2 — harness fragility, NOT a
   record regression (the caps records are untouched; git diff shows zero
   cap-related lines). OWNER DECISION NEEDED: re-baseline the pins on
   26.1.1 or redesign the sample points off the fine breakpoints — do NOT
   just loosen rel. RUFF: all checks pass (fixed: `__init__` import sort,
   quoted annotations in the two fit scripts, `zip(..., strict=True)`, the
   long Xrec line split in the generator).

PENDING next session:
- The caps-parity pin decision above (the only red test).

## Session 27 (history) — START HERE as written then

DONE (all verified):
1. **ISC040N10NM8 channel record LANDED** — closes the queue item that was
   BLOCKED ON DATA (data/ now has ISC040N10NM8_IV_25C.csv + Rdson/Vth CSVs).
   `reference.py::ISC040_CHANNEL` (`IscChannelRecord`, `params_at(T)`),
   exported in `__init__.py` (imports verified clean). Protocol
   `scripts/fit_isc040_channel.py` (reproducible, 5 multi-starts → same basin).
2. **FULL 4-terminal subckt GENERATED + LTspice-parity VERIFIED** —
   `scripts/gen_full_subckt.py` → `out/ISC040N10NM8_full.subckt` (D G S K,
   ~33.8 KB, self-contained: channel B-source + caps + DiodeRecov + package R).
   Channel B-source reproduces Python `vdmos_id` to **<0.001 %** at every point
   (bench `out/parity_isc040_full.cir`, LTspice 26.1.1, clean Newton solve).
   (Session-28 note: vs the RAW record the worst is 333 ppm — the documented
   1 µΩ Rs floor + .6g formatting; vs the netlist-constant model 0.05 ppm.)
3. **Diode DC branch FITTED** from `data/ISC040N10NM8_Diode.csv` (25/175 °C,
   26 pts) — best set **RMSE 0.51 %** (25 °C 0.45 %, 175 °C 0.58 %).
4. **Recovery R1-insensitivity MEASURED** (see below) — recovery q_rri/t_rri
   unchanged across R1 = 10 µ/1 m/10 m, so the DC-fit diode params can go into
   D1 without disturbing the pinned recovery record.

PENDING (do these first next session, in order): → **ALL THREE DONE in
session 28** — see the session-28 START HERE block at the top for the results;
the wording below is kept as the task definition (note: item 2's literal
"Rs=0.96m on D1 + R1 = r_diode" would double-count the series R — session 28
implemented D1 Rs=0 + R1 = r_diode(T) with the Trs1 law, single count).
1. **Add the diode record to `reference.py`** (`ISC040_DIODE`, `Params` with
   the fitted DC set below) + export in `__init__.py`. NOT done yet.
2. **Update the full subckt's D1 branch** to the DC-fit diode params (currently
   `gen_full_subckt.py` still emits the OLD bench diode Is=1u N=1.4 Rs=10m via
   `ISC040_RECOV.is_diode/rs_diode/n_diode`). Reconcile: D1 = DC-fit set
   (Is=3.77e-12, N=1.053, Rs=0.96m, Ikf=357, Eg=1.11, Trs1=0.00102, Tnom=25,
   Xti=3); R1 = r_diode (0.96 m) is fine (recovery insensitive, measured).
   Then re-run `out/parity_isc040_full.cir` AND a forward-V_SD check vs the
   diode CSV, and re-run the recovery bench to confirm q_rri≈248.2/t_rri≈24.69
   still hold with the new D1 (only Rs was stepped so far — Is/N not yet tested
   in the recovery bench; expected insensitive because I(D1)≈−1 nA carries no
   plateau, but VERIFY).
3. **Pinned tests**: `tests/test_isc040_official.py` (channel record vs CSVs,
   like test_imcq120_official/test_st8l65_official) + a diode record test.
   Run full pytest + ruff.

## ISC040_CHANNEL record (OFFICIAL, promoted this session)

Model family: **SJ charge-cap** (v_ovmax/k_ovsharp free, k_theta=0, v_vsat=1e4
off, k_lambda=0, g_k=0) — same family as ST8L65N044M9. Beats the velocity-sat
form (1.33 % vs 1.73 % 25 C RMSE) and is genuinely sensed (ceilings the Vgs=15
conductance; without it the Ron chord is ~13 % low).

Protocol (user's 3 notes): (1) NO transfer curve (datasheet Vds conditions
vague) → fit from output curves + Rds(on) chords only; (2) NOT Kelvin → r_sk=0,
and the data REJECTS gate-loop degeneration (r_s→0 from every start, monotonically
worse as r_s grows — the IMCQ session-22 finding), so the whole current-loop
series R is in r_d, r_s=r_sk=0; (3) fit 25 C FIRST and LOCK, then TC (no high-T
IV data — TC comes from Rds(on)-vs-Vgs at −55/25/100/175 °C + Vth-vs-Tj).

Stage 1 (25 C, locked): free = k_p, v_th0, v_sub, v_ovmax, k_ovsharp,
k_kneesharp, k_eta, v_dbsat, r_d. Data = IV_25C (69) + Rds(on)-vs-Id_25C (55,
converted to V_ds=R·I chords) + Rds(on)-vs-Vgs 25 C subset (8, I_test=25 A).
Relative weight. r_d ≤ 3.4 mΩ (datasheet 25 C chord). **RMSE 1.33 %** (132 pts).
Stage 2 (TC, 25 C locked): free = s_vth, beta_p, a_rd, b_rd. TC laws:
  k_p(T)=k_p25·(T_K/298.15)^−beta_p; v_th0(T)=v_th0_25+s_vth·(T−25);
  r_d(T)=r_d25·(1+a_rd·ΔT+b_rd·ΔT²) [quadratic Trs1/Trs2, NEVER linear];
  v_sub/v_ovmax/k_ovsharp/k_kneesharp/k_eta/v_dbsat T-independent.
  Data = Rds(on)-vs-Vgs all 4 T (chord residual as CURRENT at R_data·25 A) +
  Vth-vs-Tj both currents (residual on the Vth VOLTAGE: solve model gate V at
  V_ds=V_gs giving I_d, vs datasheet Vth — fits the Vth(T) curve/slope, avoids
  the exponential-current sensitivity the single-softplus subthreshold can't
  reproduce). **TC overall 2.97 %** (Rds(on) 2.35 %, Vth-voltage 4.04 %).

Record values (ISC040_CHANNEL):
  k_p25=222.59, v_th0_25=4.28667, v_sub=0.244088, v_ovmax=5.77459,
  k_ovsharp=2.19649, k_kneesharp=1.78701, k_eta=0.0313065, v_dbsat=0.52529,
  r_d25=0.00248295, s_vth=−0.00748634, beta_p=2.2004, a_rd=0.0056206,
  b_rd=2.00079e−05. Fixed: k_theta=0, v_vsat=1e4, k_lambda=0, g_k=0,
  r_s=r_sk=0.
Acceptance (independent specs, never fitted): Ron@25A 3.35/3.56/3.89 mΩ @Vgs
15/10/8 (data 3.33/3.52/3.87); Ron@Vgs15 vs T 2.12/3.34/5.20/7.76 mΩ @−55/25/
100/175 (data 2.05/3.33/5.23/8.00); Vth@520µA vs T 3.52/3.06/2.50/2.07 V
@−52/22/107/173 (data 3.37/3.05/2.58/2.07). s_vth −7.49 mV/K matches datasheet
Vth_Tj slope (~−6.3 mV/K, within model-form). β_p essential (without it Rdson
RMSE→13 %, channel-dominated low-Vgs breaks) and identifiable ONLY because
Rds(on)-Vgs(T) pins r_d(T) at high Vgs (IMCQ doctrine).
Identifiability caveats (quote curves, never split): r_d↔r_s corr −1.0 (only
sum sensed); k_p↔v_th0↔v_ovmax↔r_d cluster |corr|≤0.86 (composite channel
conductance is the record); β_p weak (σ≈4); v_dbsat σ≈77 (inert, k_db(5V)~0.99
but data stops at Vds≈5 V).

## Diode DC fit (data/ISC040N10NM8_Diode.csv, 25/175 °C, 26 pts) — record LANDED session 28 (`ISC040_DIODE`; protocol `scripts/fit_isc040_diode.py`)

**KEY FINDING: k_n=2.0 (the ST8L65 fixed value) is WRONG here.** With k_n free
RMSE drops 35 %→3.31 %; k_n lands ≈1.05 (near-ideal diffusion, unlike the SJ
part). i_s↔k_n corr +0.999 (degenerate prefactor pair — quote the pair/curve,
not i_s alone; i_s σ≈12). i_kf (high-injection knee) and v_eg (silicon ~1.11 eV)
well-sensed. Recombination branch NOT needed (i_sr free is worse, 4.5 %).
**Linear r_diode TC (a_rdiode) is the decisive TC: 3.31 %→0.51 %**, strongly
identified (σ 0.013). Quadratic b_rdiode and a_ikf buy nothing and go degenerate
with only 2 temperatures (a_rdiode↔b_rdiode corr −1.0) → DO NOT free them.

Best set (free = i_s, k_n, r_diode, i_kf, v_eg, a_rdiode; fixed k_xti=3.0,
i_sr=0, k_nr=2.0, t_nom=25):
  i_s=3.7679e−12, k_n=1.0529, r_diode=0.00096133, i_kf=357.18, v_eg=1.1097,
  a_rdiode=0.0010248.  RMSE 0.51 % (25 C 0.45 %, 175 C 0.58 %).
  sigma: i_s 12.3, k_n 0.45, r_diode 0.005, i_kf 147, v_eg 0.055, a_rdiode 0.013.
  corr: i_s↔k_n +0.999, r_diode↔i_kf +0.994.
NOTE: r_diode ≈ 0.96 mΩ is ~10× SMALLER than the recovery branch's R1=10 mΩ
(the old bench set). Reconcile in the subckt (see PENDING 2).

## Recovery R1-insensitivity (MEASURED this session — enables the D1 update)

`out/recov_r1_sens.cir` (record recovery knobs tau=3.78928e−4, t_minor=22.3635n,
t_pole=0.832389n, D1 unchanged Is=1u N=1.4, didt=1000 A/µs, `.step param Rs list
10m 1m 10u` → Rs drives R1 only, D1 has Rs=0):
  q_rri = 248.12 / 248.12 / 248.16 nC;  t_rri = 24.684 / 24.684 / 24.687 ns;
  i_fi_min ≈ 18.435 A  (for R1 = 10µ / 1m / 10m).
→ **Recovery is INSENSITIVE to R1** (record 248.2/24.69 reproduced at all three).
Consistent with session 23 (q_0/tau insensitive to V_SD; plateau carried by B2,
I(D1)≈−1 nA). So R1 can be set to the DC r_diode (0.96 m) and D1 to the DC-fit
set without moving recovery. (Is/N not yet stepped in the recovery bench —
expected insensitive, but VERIFY when updating D1.) → **VERIFIED session 28**:
`out/recov_d1_check.cir` A/B steps old vs DC-fit D1 at R1 = 0.961 m: q_rri
248.171 → 248.059 nC, t_rri 24.684 → 24.702 ns (≤0.075 %), record achieved
248.20/24.687 reproduced.

## Files this session
- `scripts/fit_isc040_channel.py` — channel 2-stage fit (NEW, reproducible).
- `scripts/gen_full_subckt.py` — full 4-term subckt generator (NEW).
- `reference.py` — added `ISC040_CHANNEL`/`IscChannelRecord` + `ISC040_CHANNEL_*`
  constants; `__all__` + `__init__.py` imports/exports updated (verified clean).
- `out/ISC040N10NM8_full.subckt` — the deliverable (D G S K).
- `out/parity_isc040_full.cir` — channel DC parity bench (verified <0.001 %).
- `out/recov_r1_sens.cir` — recovery R1-sensitivity bench (measured above).
- `out/isc040_iv_25C_fit.png`, `out/isc040_channel_tc.png`,
  `out/isc040_channel_params.json` — fit artifacts.

## Gotchas learned this session (do not re-derive)
- **Diode CSV temp column is `T_j (C)`, NOT `T`** → `load_diode_csv()` (expects
  `T`) FAILS on ISC040N10NM8_Diode.csv; use a local loader (see
  fit_isc040_channel.py / exploration). Same for the Rdson/Vth CSVs (script-local
  loaders, like the IMCQ spec CSVs).
- **Rdson_Vgs figure is a 25 A chord** (verified: Vgs=15.7→3.33 mΩ matches IV_25C
  row Vgs=15,Vds=0.08325,I=25). Convert to V_ds=R·I for the fit.
- **LTspice rejects R=0** ("Resistance must not be zero") and a **0 V short in
  series with the channel current source makes the DC matrix singular** (node sx
  no conductance to ground). Fix in gen_full_subckt.py: for r_s=r_sk=0 emit
  `Rs sx S 1e-6` / `Rsk S K 1e-6` (1 µΩ floor, negligible vs r_d=2.48 mΩ) and
  tie K→S by node in the bench (not a 3rd 0 V source → avoids a V-source loop).
  Added `Rleak dx sx 100Meg` to silence the benign "node dx floating on a current
  source" initial-check warning (<1e−6 relative effect).
- **`.include` in a run dir needs an ABSOLUTE path** (relative → "File not found").
- **Xrec instantiation must be ONE line** (backslash continuations in the
  generated subckt → "sub-circuit cannot be instantiated").
- Channel B-source expression (overflow-safe softplus `log(1+e^z)=max(z,0)+
  log(1+e^−|z|)`, forward clamp `max(V(dx,sx),0)`, k_theta/k_beta_trap/g_k/
  k_lambda=0 collapsed) is in `gen_full_subckt.py::channel_expr` — verified exact.

## Unchanged queue (from session 26)
- `test_ltspice_caps_parity` c1 under LTspice 26.1.1 (0.32 % vs 2e−3 tol) —
  pre-existing, re-baseline with evidence.
- Optional physics follow-ups (two-section transit ladder, 300 A/µs DPT,
  T_j recovery TC — no data, don't invent), v_gs-dependent Cgd blend (plateau
  deficit), HV leakage/BV element, T_j thermal node, third-quadrant channel.

---

# Handover — ImprovedVDMOS session 26: extract_recov.py REWRITE COMPLETE — POLE-GAUGE PROTOCOL LANDED END-TO-END

**The session-25c "⚠️ START HERE: extract_recov.py IS BROKEN" warning is
RESOLVED.** The rewrite is finished, verified against the bench, and the
whole queue (items 1–5) is done except the ISC040 channel record (blocked
on data — see below).

## What landed this session (all verified)

1. **`src/vdmos_fit/extract_recov.py` complete** (pole gauge): RecovSpec,
   NETLIST_TMPL (brace bugs fixed; C3 ground node restored), single
   `.step param k` mixed-radix index sweep, log-log multilinear
   surrogates, multi-start surrogate solve, TRF polish (diff_step=0.03,
   tau free ≥1e-6, no dTmax, no .ic), fix-and-refit identifiability at
   PERTURBED values, terminal-path reference numbers, waveform plot.
   VERIFIED FIRST as the handover demanded: rendered netlist is
   structurally identical to `DiodeTest_pole.net` (only the intended
   parameterization + t_rri25 addition), and a smoke run reproduces
   q_rri 248.542 / t_rri 24.643 / I_RM 18.514 at (3e-4, 22.283n, 0.8201n).
2. **Full ISC040 extraction** → optimum **tau = 3.78928e-4,
   t_minor = 22.3635 n, t_pole = 0.832389 n, RMS 1.113 %**; three polish
   starts (surrogate / best-grid / seed) all land in the same valley
   (tau 3.79–3.94e-4). Report
   `out/recov_extract_isc040_run/extract_report.txt` (identifiability
   rerun appended), `extract_results.json`, figure
   `out/RECOV_pole_extract_waveforms.png`.
   Identifiability (fix knob at optimum ×0.5/×2, refit the other two):
   tau +0.13/+0.08 pp (FLAT valley — quote the coupled set, never tau
   alone); t_minor +9.6/+86 pp (dominant); t_pole +8.1/+5.4 pp (sharp).
   NOTE: the 5×5×5 grid (125-step index sweep) takes ~10 min per di/dt
   point, not the hoped 30 s — extreme tau corners (1e-5, B2 gm ~MS)
   cost iterations. Surrogate max resid 8–18 % over the wide box; the
   polish fixes it, but a tighter box or quadratic terms would speed a
   cross-part reuse.
3. **`reference.py`: `RecovRecord` + `ISC040_RECOV`** (tau
   3.7892757739719363e-4, t_minor 2.2363459877611823e-8, t_pole
   8.32388859768649e-10, rms_rel 0.011132; DC branch Is=1u Rs=10m
   N=1.4; k_chg=5 r_max=0), exported from `vdmos_fit`, pinned by
   `tests/test_recov_record.py` (4 tests; pins record values, pole-gauge
   netlist structure — q_0 must never return, no dTmax, no .ic — and the
   index-step decomposition; needs no LTspice run).
4. **`out/ISC040N10NM8_recov.subckt` regenerated in the pole gauge** via
   the new `scripts/gen_recov_subckt.py` (record-driven — never hand-edit
   the subckt). VALIDATED: a standalone `.include` bench reproduces the
   record exactly (q_rri 248.200 nC, t_rri 24.687 ns, I_RM 18.438 A).
5. **`ltspice.py` `cgs_cgd_only=True` flag** + `scripts/export_caps.py`
   `CGS_CGD_ONLY = {"ISC040N10NM8"}`; `out/ISC040N10NM8_caps.subckt`
   regenerated Cgs/Cgd-only (Cds lives in the recovery branch — never
   both). Test added (`test_export_cgs_cgd_only_flag`).
6. **Docs**: README recovery section rewritten to the pole gauge (+
   naming bullet `t_minor/t_pole/tau/k_chg/r_max`, + extract_recov.py in
   the layout tree); STATUS.md snapshot + new "sessions 25a–25c" section
   + Missing bullet + session-23 heading marked superseded (integration
   decisions stand); stale gotchas corrected in `RECOV_INTRINSIC_FIT.md`
   §4, `RECOV_POLE_FIT.md` §8, `RECOV_EXTRACTION_PROTOCOL.md` §2/§7, and
   the 25a/25b/24 HANDOVER notes (chained `{i_fi_min}` stepped-table
   claim RETRACTED; "dTmax mandatory" = pole-free era only).

## Suite state

Full pytest: 185 tests, ONE failure —
`tests/test_caps.py::test_ltspice_caps_parity` point c1 at 0.32 % vs the
2e-3 tolerance. **PRE-EXISTING**: reproduces identically with the
untouched HEAD `ltspice.py` (verified by stash-test); the box upgraded
LTspice 26.0.1 → 26.1.1 since the parity was pinned (session 22 log
passes). Left untouched (out of scope) — investigate or re-baseline with
evidence next session. Ruff clean on every touched file.

## NEXT SESSION — remaining queue

1. **ISC040 channel record — BLOCKED ON DATA**: `data/` has no ISC040
   transfer/IV curves (only Cap_Vds, Vgs_Qg, cds_expr). Digitize a
   transfer curve (curve_digitizer) first, then the standard
   fit_datasheet protocol.
2. `test_ltspice_caps_parity` c1 under LTspice 26.1.1: diagnose
   (gear-integration vs `.meas FIND` timing artifact?) or re-baseline
   the tolerance with evidence. STATUS parity record says ≤0.3 %; the
   test pins 0.2 %.
3. Unchanged optional physics follow-ups: two-section transit-time ladder
   if the last 1–3 % of t_rr matters; out-of-sample DPT at 300 A/us
   (model: Qrr ≈ 104 nC, t_rr ≈ 29 ns, I_RM ≈ 6.7 A vs ~90 nC power-law
   interpolation); T_j dependence (no data — do not invent a TC law);
   cross-part check if ST8L65/IMCQ datasheets carry recovery specs;
   B2 `1/tau` → explicit stiff gate.
4. k_chg/r_max stay fixed (not identifiable from recovery data). The
   ISC040 D1 DC branch is the bench set (Is=1u, Rs=10m, N=1.4), NOT the
   diode.py record path — refine separately against V_SD specs if that
   matters (session 23 noted q_0/tau-insensitivity of V_SD).
5. Extraction cost: consider a tighter grid box (tau lo 1e-5 → ~3e-5
   still brackets the 3.8e-4 optimum with margin) or part-to-part
   surrogate seeding if running other parts.

---

# Handover — ImprovedVDMOS session 25c: dTmax DROPPED, FITTER FIXED, DRIVER REWRITE IN PROGRESS

## ⚠️ START HERE: `src/vdmos_fit/extract_recov.py` IS BROKEN / INCOMPLETE

Session 25c was mid-rewrite when the context window ran out. The file now
contains a valid docstring + `RecovSpec` dataclass + `NETLIST_TMPL` +
`find_ltspice()` + `_cds_clause()` + `_spice_num()`, then a **broken stub**
at the end:

```python
def make_netlist(spec: Spec := None, *, didt: float, steps: str = "",
                 label: str = "bench", **kw) -> str:
    raise TypeError("use make_netlist(spec, didt=...)")
```

`Spec` does not exist (the dataclass is `RecovSpec`) and `spec: Spec := None`
is a syntax error — **the module will not import.** Everything after
`_spice_num()` (netlist builder, index-step emitter, log parser, surrogates,
inverse solve, polish, identifiability, plotting, `extract()` entry point)
was never written. The old session-24 implementation is recoverable from
git (`git show HEAD:src/vdmos_fit/extract_recov.py`) — it is the two-knob/
t_minor-gauge version with the obsolete `Vv`/`D5`/`.ic` bench and `I(R4)`
plot references, so use it for structure only, not for topology.

**The `NETLIST_TMPL` in the new file has NEVER been executed.** It was
written from the verified scratch bench (`outputs/vdmos_recov/
DiodeTest_pole.net`) but the `{{`/`{{{{` brace-doubling for `.format()` is
UNTESTED — e.g. `PULSE({{-i_test}} {{{{i_test*9}}}} ...)` renders as
`{{i_test*9}}` (double braces) which is wrong for a single-pass format.
First action next session: render the template, diff it against
`DiodeTest_pole.net`, and run it before trusting anything.

## What session 25c established (all measured — do not re-derive)

User asked three things: (a) why didn't the fitter reach tau=3e-4, (b) tau
must stay free because 3e-4 is device-specific, (c) dTmax is no longer
needed since reltol=1e-4 forces small steps. Then asked to update
`extract_recov.py`, `README.md`, `STATUS.md`. **README.md and STATUS.md are
NOT yet updated — still pending.**

### (c) dTmax is genuinely unnecessary now — user is right

At the optimum (tau=3e-4, t_minor=22.283 n, t_pole=0.8201 n), 1000 and
100 A/us, vs a dTmax=0.2 n reference:

| dTmax | q_rri err | t_rri err | iters (100 A/us) |
|---|---|---|---|
| none | +0.01 % | +0.01 % | 4088 |
| 2.5 n | +0.01 % | +0.00 % | 2788 |
| 1.0 n | +0.01 % | +0.00 % | 4825 |
| 0.5 n | +0.01 % | −0.01 % | 8368 |

All agree to ≤0.01 %, and **dropping dTmax is FASTER** (4088 iters vs 8368
at 0.5 n). The session-24 "dTmax mandatory" finding belonged to the
**pole-free** model, where the snap was a hard discontinuity; the pole
rounds it and the adaptive controller resolves it. New benches: no dTmax.

### (a)+(b) the fitter CAN find tau — my session-25b "fix tau" advice was wrong

With the tau bound extended to 1e-5 and `max_nfev=200`, TRF from seed
(tau0=0.1, tm0=17 n, tp0=2 n) reaches **tau=3.982e-4, t_minor=22.361 n,
t_pole=0.8372 n, RMS 1.119 %** (nfev=51, status=3). So tau does not need
fixing — the earlier stall was a too-small step budget + bound, exactly as
the user's Newton intuition said. **tau stays a free knob** (device-specific
optimum, user is right).

Refined optimum: **tau ≈ 3.98e-4, t_minor ≈ 22.36 n, t_pole ≈ 0.837 n,
RMS 1.119 %** (profile-fit reference: 3e-4 / 22.283 n / 0.8201 n, 1.122 %).

### Finite-difference step: bigger is better here (counterintuitive, measured)

`probe_fitter.py` configs from seed s1:

| config | result |
|---|---|
| A trf, diff_step=0.03, x_scale=jac, max_nfev=200 | **tau 3.98e-4, RMS 1.119 %** ✓ |
| B trf, diff_step=1e-4, ftol=xtol=gtol=1e-12 | tau stuck 1.185e-1, **RMS 5.566 %** ✗ |
| C trf, x_scale="log" | **ValueError** — scipy `least_squares` accepts only `'jac'` or array_like; configs C/D/E never ran |

Reason: with diff_step=1e-4 the cost differences (~1e-6) fall **below the
simulator's own measurement noise** (~1e-4 relative on q_rri/t_rri), so the
Jacobian is noise-dominated. diff_step=0.03 puts the cost difference ~1e-2,
far above noise. **Tighter tolerances do not help for the same reason.**
Driver default: `diff_step=0.03`, generous `max_nfev`, bounds tau ≥1e-6.

### My dTmax-noise hypothesis was WRONG (recorded so it isn't re-tried)

Cost-surface roughness along a 13-point log t_minor sweep (RMS of 2nd
differences, normalized): dTmax=2.5 n → 0.1181 %, dTmax=none → 0.1174 %
(1000 A/us); 0.3363 % vs 0.3355 % (100 A/us). **Identical** — dTmax does not
roughen the cost surface. The session-25b stall was path-dependence in a
flat valley + small step budget, nothing else.

### CORRECTION to a session-25a gotcha (it was my bookkeeping error)

Chained `{i_fi_min}` measurements **DO work in stepped runs.** A fresh
stepped test emitted both `t_rri` (23.915/23.525/24.901/24.509 ns) and
`t_rri25` (23.123/22.559/24.118/23.560 ns) as full 4-row blocks. The
session-25a "dependent block dropped from stepped tables" claim came from a
netlist that **never contained a `t_rri25` line** (verified:
`grid_hi.net`/`grid_lo.net` → "netlist has t_rri25 line: False"). Delete
that gotcha from `RECOV_INTRINSIC_FIT.md` §4 and the 25a/25b handover
notes. `t_rri` vs `t_rri25` differ by ~0.8 ns at these settings.

### Verified driver building blocks (use these, they work)

- **Single index sweep with mixed-radix decomposition** — removes the
  `.step` nesting gotcha entirely and is deterministic:
  ```
  .step param k list 0 1 2 3 4 5 6 7 8
  .param nA=mod(k,3)
  .param nB=int(k/3)
  .param tau_p=exp(ln(1e-4)+nA*ln(3))
  .param tm_p=exp(ln(20n)+nB*ln(1.1))
  .param tp_p=0.8201n
  ```
  Global step params DO drive values **inside** the subckt including the
  `C3 q_j {tp_p}` value (rc 0, values respond). GOTCHA: every referenced
  param must be defined or LTspice dies at parse with
  `No such parameter defined.` pointing at `C3 q_j 0 {tp_p}` (rc 1, no
  measurements at all).
- **`plotwinsize=0` changes nothing** on these `.meas` results (identical to
  6 digits: q_rri 248.5420, t_rri 24.64342, I_RM 18.51401) — optional, not
  needed for measurement accuracy here.
- Bench runs from the DC operating point with **no `.ic`** (no DC current in
  R3 ⇒ V(q_j)=V(q_j0)); `solver=alt` + `method=gear reltol=1e-4`; tstop
  `{tr*5+50n}`; run time 0.13–0.20 s per point.

## Files this session

- `outputs/vdmos_recov/probe_tau_bound.{py,log}` — tau valley profile
  (fix tau, refit the other two): tau 2e-2/3e-3/1e-3/**3e-4**/1e-4/3e-5/1e-5
  → RMS 3.12/1.78/1.27/**1.12**/1.57/2.74/3.41 %; interior minimum at
  ~3e-4, degrades below. Analytic law verified to ~1 %:
  `Δt(clamp→I_RM) = t_pole·ln(V(q_j)|clamp / V_th)`, `V_th = I_RM·tau/ΔV`,
  so `d(Δt)/d(ln tau) = −t_pole` (measured 1.15/3.22/5.67 ns vs analytic
  1.15/3.25/5.76 ns at tau 2e-2/3e-4/1e-5). Not numerical: 1534/1451/1740
  iterations, B2 gm up to 1.25 MS stays solvable.
- `outputs/vdmos_recov/probe_fitter.{py,log}` — fitter config study (A/B ran;
  C crashed on invalid `x_scale="log"`; D lm and E profile+refine never ran —
  worth finishing, E is the robust route if TRF proves seed-sensitive).
- `outputs/vdmos_recov/fit_pole.py` — bounds fixed (tau ≥1e-5, diff_step
  0.01, max_nfev 200, prints `status`/`active_mask`). NOTE: its
  `diff_step=0.01` is now known to be sub-optimal vs 0.03 (see above).
- `out/RECOV_POLE_FIT.md` §4b added (tau interior minimum + correction of
  §4's "saturating toward tau→0" wording + fitter guidance). §4 and the
  HANDOVER 25b bullet still say "fix tau" — **now superseded by this
  section; needs correcting.**

## NEXT SESSION — ordered

1. **Fix `extract_recov.py`** (it does not import). Finish the rewrite in the
   pole gauge: `RecovSpec` (per-part datasheet points + grid + bounds),
   index-decomposed single `.step`, no dTmax, no `.ic`, intrinsic
   `q_rri`/`t_rri` doctrine, log-log surrogates, multi-start surrogate solve,
   polish with `diff_step=0.03` + generous budget + tau free (≥1e-6),
   fix-and-refit identifiability, terminal-path reference numbers, waveform
   plot. Verify the template renders to the working bench before anything else.
2. **README.md** — the recovery section is stale (still describes the
   two-knob `C1=q_0/tau`, `Roff=tau`, `q_0`/`tau` model, "Two knobs", the
   35 n/5 fit, and "Datasheet-driven extraction protocol: open"). Must be
   rewritten to the pole gauge: `C1={t_minor}`, `Roff=1`, `E1`+`R3`+`C3`
   pole, three free knobs (t_minor, t_pole, tau), q_0 deleted as a gauge,
   snap condition `V(q_j)=I_RM·tau/ΔV`, the Δt law, and the ISC040 result
   (tau 3.98e-4, t_minor 22.36 n, t_pole 0.837 n, RMS 1.12 % on all four
   datasheet numbers). Also the naming-convention bullet lists recovery
   params as (`q_0`, `tau`, `k_chg`, `r_max`) — must become
   (`t_minor`, `t_pole`, `tau`, `k_chg`, `r_max`).
3. **STATUS.md** — add a session-25c section after the session-23 recovery
   section (which is now stale: it documents q_0/tau semantics, the
   35 n/5 fit, and the "Qrr rises with di/dt is unreproducible" verdict —
   that verdict is DEAD, the pole reproduces both di/dt points). Update the
   "Missing" bullet on body-diode recovery (extraction protocol now exists).
4. Correct the stale gotchas: delete the "chained `{i_fi_min}` dropped in
   stepped tables" claim (25a §4 + 25b notes); mark "dTmax mandatory"
   (24 §2/§7) as pole-free-era only.
5. Then the unchanged queue: regenerate `out/ISC040N10NM8_recov.subckt` in
   the pole gauge; `reference.py` `*_RECOV` record + pinned test; `ltspice.py`
   Cgs/Cgd-only export flag; ISC040 channel record.
6. Optional physics follow-ups (unevaluated): two-section transit-time ladder
   if the last 1–3 % of t_rr matters; out-of-sample check at 300 A/us
   (model predicts Qrr ≈104 nC, t_rr ≈29 ns, I_RM ≈6.7 A — a real DPT
   measurement would discriminate against the ~90 nC power-law interpolation).

---

# Handover — ImprovedVDMOS session 25b: TRANSIT-TIME POLE (idea B) — FIT SOLVED

## NEXT SESSION = ADOPT THE 2-KNOB PROTOCOL + REGENERATE THE TEMPLATE

User's idea B landed and WORKS: an RC pole between the storage node and the
conduction-driving node models the finite carrier transit/sweep-out time.
Bench adopted as `outputs/vdmos_recov/DiodeTest_pole.net` (user's netlist,
`solver=alt`). Durable record: `out/RECOV_POLE_FIT.md` + `out/RECOV_pole_fit.png`.

**Headline: RMS relative error over all four datasheet numbers**
(no pole: 17.45 %) -> (pole, 3-knob: 3.20 %) -> (pole, tau FIXED at the
stiff limit, only 2 fitted knobs: **1.12 %**).
The Qrr-ratio problem is DEAD: the pole-free model's 0.175 floor (0.30 in
the old terminal gauge) vs datasheet 0.146 is now tunable through — the
hi-priority weighting hits BOTH Qrr rows to 0.3 %.

**Recommended parameterization / extraction protocol (supersedes 25a):**
- `C1={t_minor}`, `ClampD Roff=1` -> storage tau = t_minor, V_q0,ss = I_F,
  stored charge = I_F·t_minor. **q_0 deleted (it was a pure gauge).**
- `E1` unity buffer + `R3`+`C3` -> EXACT pole tau_p = R3·C3 (no C1 loading).
- **tau: interior minimum at ~3e-4, NOT a monotone approach to the stiff
  limit.** Earlier wording ("saturating toward tau->0") was WRONG and is
  corrected in `RECOV_POLE_FIT.md` sec. 4b. Valley profile (tau fixed,
  other two refit): tau 2e-2/3e-3/1e-3/**3e-4**/1e-4/3e-5/1e-5 ->
  RMS 3.12/1.78/1.27/**1.12**/1.57/2.74/3.41 %. Degrades below 3e-4, so
  3e-4 is a genuine optimum, not an arbitrary stopping point.
  Reason (analytic, verified to 1 %): the clamp->I_RM delay is
  dt = tau_p*ln(V(q_j)|clamp / V_th) with V_th = I_RM*tau/(V_Ki-V_Ai), so
  d(dt)/d(ln tau) = -tau_p — shrinking tau ADDS tail time that must be paid
  by shortening tau_p (sharply identified, wall below ~0.6 n) and t_minor
  (off ~22 n costs the Qrr ratio). Minimum = where those costs balance.
  Not numerical: 1534/1451/1740 iterations at tau 2e-2/3e-4/1e-5, no
  convergence trouble, B2 gm up to 1.25 MS stays solvable.
- **FITTED BOUNDS GOTCHA:** `fit_pole.py` originally had `tau >= 0.02` and
  seedA converged to `tau=0.020000009` — pinned at the bound, costing ~2 %
  RMS. Bounds now extended to 1e-5 (diff_step 0.01, max_nfev 200). With
  bounds fixed, the 3-knob fit descends but terminates `status=3` (ftol) at
  seed-dependent tau (1.28e-3/4.35e-3/1.74e-3, RMS 1.38/2.00/1.51 %) — the
  tau valley is flat and coupled to BOTH other knobs, so coordinate-wise
  finite differences stop short of the bottom. **Treat status==3 as not
  converged; prefer fixing tau=3e-4 and fitting two knobs.**
- Fit exactly TWO physical knobs: **t_minor = 22.283 n** (carrier lifetime)
  and **t_pole = 0.8201 n** (sweep-out memory) at tau=3e-4 ->
  Qrr 248.5 (+0.6 %), t_rr 24.64 (−1.4 %) @1000 A/us;
  Qrr 35.7 (−0.8 %), t_rr 29.42 (+1.4 %) @100 A/us.
- t_pole is SHARPLY identified (fix-and-refit: 0.3 n->15.6 %, 0.6 n->10.8 %,
  optimum->1-3 %, 2 n->11.6 %, 4 n->27.5 %, 8 n->52.8 %, ->0 = 17.45 %).

**Mechanism (verified on raw):** q_j0 clamps to 0 at 54.6 ns (hi) / 399.4 ns
(lo); at that instant V(q_j) still holds +0.29 V / +0.12 V and then decays
with 1.17 ns = R3·C3 exactly. So **I_RM occurs ~1-2 ns AFTER charge
exhaustion** — the transit-memory signature. Snap condition verified exact:
V(q_j)·(V_Ki−V_Ai)/tau = 17.79 A = I_RM (17.76 A); B2 carries the plateau
(I(D1) ~ −1 uA).
**Physical anchor:** t_pole ~0.8-1.2 n is diffusive — L = sqrt(2·D_a·t_pole)
= 1.4-2.7 um (D_a 12-36 cm^2/s), the right scale for the depletion-adjacent
region of a 100 V SJ device. Drift at v_sat would cover 82-98 um, so it is
NOT a drift-transit time. Do not over-interpret.
**Out-of-sample prediction (falsifiable):** at 300 A/us the model predicts
Qrr ~104 nC, t_rr ~29 ns, I_RM ~6.7 A (t_rr nearly constant 29.9/28.9/25.0 ns
at 100/300/1000). A power-law interpolation of the two datasheet rows gives
~90 nC -> model ~15 % high; discriminate with a real DPT at 300 A/us.

**Next session:**
1. Regenerate `out/ISC040N10NM8_recov.subckt` in the pole gauge with the
   recommended set (tau fixed, t_minor/t_pole fitted). Do NOT mix gauges.
2. If the last 1-3 % of t_rr matters: upgrade the single pole to a
   **two-section transit-time ladder** (a distributed drift region gives a
   pole chain, not one pole) — adds phase lag without changing DC gain,
   which is what the residual t_lo bias wants.
3. Consider replacing B2's `1/tau` with an explicit stiff gate: the model is
   converging to an ideal switch controlled by the filtered charge node.
4. Do NOT add di/dt-dependent q_0/tau laws — the pole already produces the
   Qrr-rises-with-di/dt behavior from one physical state.
5. Then the session-24 queue (reference.py `*_RECOV` record, ltspice.py
   Cgs/Cgd-only export flag, ISC040 channel record) unchanged.

**Gotchas added this session** (full list in `RECOV_POLE_FIT.md` §8): the
pole does not break `.ic`-free DC-op-point seeding (no DC current in R3 ->
V(q_j)=V(q_j0)); ~~forward-referenced `i_fi_min` resolves in single runs but
chained dependent `.meas` blocks are still dropped in `.step` tables~~
(CORRECTED in 25c: chained blocks DO work in stepped runs — bookkeeping
error); user extended tstop to `{tr*5+50n}` for the longer tail.

---

# Handover — ImprovedVDMOS session 25: INTRINSIC (C_oss-free) RECOVERY MEASUREMENT + REFIT

## (Superseded by session 25b — the pole model replaces the 3-knob no-pole fit; kept for history)

Session 25 adopted the user's updated DiodeTest netlist
(`outputs/vdmos_recov/DiodeTest_user.net`: Tau removed where unnecessary,
R4 removed, `ClampD Roff=t_minor/q_0`, `C1={q_0}`, `B1=-I(R2)`, B3 off,
C2 charge-master) and refit recovery per the C_oss-contamination decision
(datasheet Qrr does NOT include Q_oss → fit `q_rri`/`t_rri` measured on the
intrinsic branch `I(x1:R2)`; C2 displacement bypasses R2 at node A).
Durable record: `out/RECOV_INTRINSIC_FIT.md` + `out/RECOV_intrinsic_waveforms.png`.

**Key session-25 findings (do not re-derive):**
- New parameterization semantics (sim-verified): storage τ = t_minor
  (Roff·C1, tau-free), V_q,ss = I_F·t_minor/q_0, stored charge = I_F·t_minor;
  tau is a pure B2 gain/snap-shape knob. Model-line instance-param scoping
  (`Roff=t_minor/q_0`) works in 26.1.1.
- Bench runs from DC op point, no `.ic` seed, 0.13 s; loop di/dt = didt exactly.
- q_0 NOT identifiable (<2 % over 30–55 n) — fix 41.1 n. t_minor dominates
  Qrr and t_rr; tau stretches t_rr and softens the snap.
- Qrr_lo/Qrr_hi ratio floor now ≈0.175 at t_minor=20 n (was 0.30 terminal-path);
  datasheet 0.146 still unreachable at the tm≈36 n that Qrr_hi=247 requires
  (ratio ≈0.24 there).
- Fits (targets 247 nC/25 ns + 36 nC/29 ns, equal relative weights unless
  noted): seed (0.1, 30 n) → hi −21 %/−21 %, lo +21 %/+1 %; equal (0.584,
  29.4 n) → hi −24 %/−17 %, lo +12 %/+10 %; hi-prio (0.939, 36.3 n) →
  hi +0.0 %/−4 %, lo +70 %/+35 %. **Recommend hi-prio set for the
  hard-switching record.**
- C_oss contamination at 100 A/µs ≈28 nC / ≈15 ns on the terminal path —
  the intrinsic measurement is what makes the lo point comparable at all.

**Next session ideas:**
1. Regenerate `out/ISC040N10NM8_recov.subckt` from the user's new subckt
   semantics + chosen fit set (template still carries the old two-knob
   `Roff={tau}`/`C1={q_0/tau}` form — do NOT mix the two parameterizations).
2. t_rr definition check: datasheet t_rr likely IEC t_rr25; user's t_rri is
   10 % of I_RM. At tau≈0.9 the snap is soft enough that the 10 %/25 %
   split is ~1 ns — decide the canonical definition and refit if needed.
3. Structural route to the residual Qrr_lo tension remains a second storage
   state (dual lifetime); the ratio-floor numbers above bound what a
   single-state model can do.
4. Then the session-24 queue items (reference.py `*_RECOV` record, ltspice.py
   Cgs/Cgd-only export flag, ISC040 channel record) unchanged.

**Gotchas added this session** (full list in `RECOV_INTRINSIC_FIT.md` §4):
`.tran` dTmax is the FOURTH positional arg (2nd position silently rebinds
Tstop → `No. Points: 0` + "Unexpected file format"); `.meas` INTEG log
values already carry the `*1e9` (nC) — don't rescale; ~~chained `{i_f_min}`
measurements are dropped from stepped-run log tables~~ (CORRECTED in 25c:
they DO work — the claim came from a netlist that never had the line); log
value = token before `FROM` then `split("=")[-1]`; event windows anchor on
the FIRST falling zero-crossing of `Ix(x1:a)`.

---

# Handover — ImprovedVDMOS session 24: RECOVERY DATASHEET EXTRACTION PROTOCOL + t_minor KNOB

## (Superseded by session 25 — intrinsic-path doctrine; kept for history)

Session 24 built the datasheet-driven extraction pipeline
(`src/vdmos_fit/extract_recov.py`, `extract()` two-knob + `extract3k()`
with the user's t_minor knob), verified it on ISC040N10NM8, and mapped the
two-knob limitation quantitatively. Durable record:
`out/RECOV_EXTRACTION_PROTOCOL.md` (bench design, gotchas, t_minor scan +
fit tables), `out/RECOV_3k_waveforms.png` / `RECOV_3k_hi_waveforms.png`,
`out/recov_extract3k_{run,eq,hi}/extract3k_report.txt`. Session-23 material
unchanged in STATUS.md.

**Key session-24 findings (do not re-derive):**
- Fast single-event bench (seeded `.ic V(X1:q_j)={i_test*t_minor/q_0}`,
  stiff `Vv N001 0 V_R` + D5 clamp + 1p across D5, C2 Ki A = Cds
  charge-master from `data/ISC040N10NM8_cds_expr.txt`, r_max=0, gear,
  dTmax=t_rr_target/10) reproduces the 500-period bench to <0.5 %;
  full extraction ~15 s.
- `.meas` windows terminate at fractions of MEASURED I_RM (chained
  `i_pk`), not of I_F — the old −frac·I_F targets fail at low di/dt.
- t_minor knob (`ClampD Roff=tau*t_minor/q_0`): stored charge =
  I_F·t_minor/tau (q_0-free), identity with two-knob at t_minor=q_0.
  Qrr_lo/Qrr_hi ratio has a FLOOR ≈ 0.30 (t_minor 15–35 ns); datasheet
  needs 0.146. Single storage state cannot decouple the two di/dt points.
- Best fits: hi-priority 3-knob (35.3n, 3.51, 34.1n) → Qrr_hi −0.5 %;
  equal-weight (46.1n, 9.84, 24.4n) → Qrr_lo +28 %, t_lo +3.6 % but
  Qrr_hi −42 %. Two-knob C (34.8n/5) remains the hard-switching record.

**Next session — investigate how to further improve the fit.** Ideas
(unevaluated, pick from data):
1. Second independent storage state (dual lifetime: fast/slow
   recombination channel) to break the 0.30 Qrr-ratio floor — the only
   structural route to both di/dt points; estimate identifiability first.
2. di/dt-dependent q_0 or tau law (q_0·(1+c·di/dt)) — empirical, check
   against the 3-point datasheet if available.
3. Revisit the D1 DC branch / Is,N fit (V_SD typ 0.82 vs 0.91 @ 25 A) —
   session 23 noted q_0/tau insensitive, but the plateau voltage split
   could shift the snap point slightly.
4. Check whether the datasheet t_rr definition is really t_rr25 (IEC) —
   the +7…+10 % residual might be a definition mismatch, not a model error.
5. If a second part (ST8L65/IMCQ) carries recovery specs, run the
   protocol there as a cross-part check of the model family.
6. Then: `reference.py` `*_RECOV` record + pinned test (same doctrine as
   channel/diode/caps records); ltspice.py Cgs/Cgd-only export flag; the
   ISC040 channel record (queue item 3 unchanged).

**Gotchas added this session** (full list in protocol §7): ASCII raw parse
(index line + one value per line; column k = variable k); single-value
`.step` rejected; `.step` nesting last=outer; failed `.meas` rows print
`failed`; ~~dTmax mandatory~~ (superseded for the pole model, 25c: no dTmax
needed — pole-free era only); `V(X1:node)` unquoted in `.ic`/`.meas`.

---

# Handover — ImprovedVDMOS session 23: RECOVERY NETWORK INTEGRATED + FITTED (ISC040)

Session 23 integrated the user's behavioral body-diode recovery network
(`DiodeRecov`) into the full-device LTspice subckt and fitted its two knobs
on ISC040N10NM8 — **LTspice-side only, NO package Python touched**. Durable
record: `out/ISC040N10NM8_recov.subckt` (template, Cds ramp-sum inline) +
`out/RECOV_REPORT.md` (semantics, sweep tables, integration decisions, A/Bs);
session-23 section in STATUS.md. Session-22 material is unchanged there
(session-22 section — verdicts on record, do not re-litigate).

**What landed (all simulation-verified, LTspice 26.0.1, `method=gear
reltol=1e-4`):**
- Fit: **q_0 = 35 n, tau = 5** → Qrr 249.6 nC vs datasheet typ 247 (+1.1 %),
  t_rr(25 % I_RM) 27.6 ns vs 25 (+10 %) @ V_R = 50 V, I_F = 12.5 A,
  1000 A/µs. q_0 sets timing (Qrr/t_rr ≈ linear), tau sets shape
  (plateau→tail departure, snap hardness). Variants: 34n/3 datasheet-strict,
  37n/5 waveform-strict.
- Model semantics decoded: ClampD `Roff=tau` closes the storage loop —
  V(q_j) = I_F volts at steady state (tau-independent), q_j time constant =
  Roff·C1 = q_0 (35 ns; the 1 Ω normalization ties the "charge at 1 A" and
  "carrier lifetime" readings), stored charge = q_0·I_F/tau ∝ I_F. The
  −I_RM plateau is carried by B2, not D1 (I(D1) ~ −1 nA); plateau ends at
  V(q_j) ≈ I_RM·tau/V_R (predicted 1.25 V, observed 1.33 V).
- Cds placement: caps delta exports **Cgs/Cgd only**; the measured C_ds
  charge-master lives INSIDE the recovery branch at **Ki–A** (electrically
  D–S through B3+R1/R2, fc ≈ 3 GHz). D–S placement crashes at hard forward
  re-turn-on (timestep 1.25e-19; D1 has no Cjo, B3+R1 isolate the branch);
  Ki–Ai is singular (user's finding, reproduced).
- Half-bridge stability (user fix, now template default): `R2 = 1 µΩ` + B3
  `+1µ` floor — breaks near-singular B3/R1/R2/Ron bridge loops; waveform
  impact < 0.1 % (A/B over r2 × floor, `DiodeTest_fixcheck.net`).

**Known two-knob limit (quantified, do not chase):** datasheet Qrr *rises*
with di/dt (36 nC @100 A/µs → 247 nC @1000 A/µs); no single stored-charge
state follows both — at 100 A/µs the fit gives Qrr 77 nC (≈2× high). The fit
is tuned to the hard-switching condition; the 100 A/µs numbers are "defined
by design" in the datasheet.

**Extraction protocol (next session) — where to start:**
- Datasheet inputs available on modern parts: t_rr + Qrr (rarely I_RM) at
  ≥2 di/dt points (ISC040: 100 and 1000 A/µs), V_SD@I_F/T_j for the D1 DC
  branch (already the `diode.py`/`ST8L65_DIODE` record path), C_oss/Q_oss
  for the Cds branch (already the caps record).
- Order of operations: fit the D1 DC branch FIRST (existing protocol), then
  pin Qrr at the HIGH-di/dt point with q_0 (one knob, ≈linear), then shape
  (t_rr / I_RM snap step) with tau, then CHECK the low-di/dt point and
  REPORT the residual as the model-family limitation. Extra knobs
  (di/dt-dependent tau or q_0, second storage node, soft B2 factor) only if
  the residual demands them — decide from data, not in advance.
- Bench: `tooling/workspace/outputs/vdmos_recov/DiodeTest_sweep.net`
  (settled 500-period; `.meas TD` RELATIVE to the shifted axis — gotcha
  below). k_chg = 5, r_max = 0.1 were held fixed at the user's values —
  probe whether they are identifiable at all (likely weak: k_chg only moves
  the forward-recovery tail, r_max only the turn-on overshoot).
- T_j: datasheets give recovery at 25 °C only; carrier lifetime has a strong
  TC (rule of thumb ≈ +50 %/100 °C on Qrr). Decide whether to add a q_0(T)
  law or stay isothermal — do not invent a TC without data.
- Once the protocol is settled, consider a `reference.py` `*_RECOV` record
  + a pinned test, same doctrine as the channel/diode/caps records.

# Queue (updated)

1. **`v_gs`-dependent C_gd blend (channel charge, CAP_SPEC limitation 1)** —
   targets unchanged: plateau deficits SJ +61 %, SiC +9.7 %, trench +42.8/24.5/22.0 %.
2. CLOSED (session 22): **Kelvin-aware r_s split** — `Params.r_sk` landed;
   the joint-vector refit shows the two open IMCQ gaps are NOT series-R
   topology problems (verdicts in STATUS session-22 section).
3. **Full-device LTspice subckt** (channel B source + diode + r_d/r_s/r_sk
   wiring + caps subckt, 4 terminals D G S K) then double-pulse Eon/Eoff
   within ~10 %. **PARTLY LANDED (session 23):** the diode-recovery half is
   done and fitted (`out/ISC040N10NM8_recov.subckt`). Still needed: the
   channel B-source expression from `model.py` (needs the `k_beta_trap`
   factor — trivial to express, README model block); a `ltspice.py` flag to
   export **Cgs/Cgd only** (Cds now lives in the recovery branch — never
   export both, double counting); the ISC040 channel record (no official
   channel fit exists for that part yet). Kelvin wiring (session 22, still
   the plan): `D G S K`, `r_d` drain leg, channel B-source between internal
   D'/S', `r_s` S'→S, `r_sk` S'→K, gate caps referenced G→internal, driver
   return at K.
4. **NEW: recovery-parameter extraction from datasheet info** — see NEXT
   SESSION above; protocol design + first application (ISC040 re-extraction
   as the round-trip check, then ST8L65/IMCQ if their datasheets carry
   recovery specs).
5. Unchanged: third-quadrant channel (sign-symmetric sync-rect conduction),
   HV leakage/BV element, T_j thermal node, CI, diode 5 A re-digitization,
   ST8L65 `iss` gridline-snap artifact (four identical 4869.28716 pF points).

Closed: item 2 Kelvin split (session 22), "digitize transfer + joint fit"
(session 21, official `IMCQ120_JOINT`), record bookkeeping (session 21),
`k_beta_trap` (session 19), CSV header contract (session 20) — see STATUS.md.

# Notes / gotchas (carry forward)

- Tool name: `engineering__execute_python` (**not** `files__execute_python`).
  Hard 120 s wall per call — long fits must be launched with
  `subprocess.Popen(..., DETACHED_PROCESS)` writing a log, then poll with
  `tasklist /FI "PID eq <pid>"`. Keep the polling sleep itself under the cap
  (sleep(60) with the default 60 s timeout FAILS — use sleep ≤ 45 s or pass
  `timeout_seconds`). Use `python -u` for a streaming log.
- Two interpreters on this box: the session interpreter
  (`tooling/engineering313`, scipy 1.18.1 — what sessions 21–22 used, and
  what `engineering__python_environment` selects) and system Python
  (`AppData\Local\Programs\Python\Python313`, scipy 1.16.0). Fit results can
  differ in unsensed directions; there is **no `.venv`** — scripts add
  `src/` to `sys.path` themselves.
- `engineering__execute_python` writes scratch into the workspace; **no
  `/tmp`** — write into `outputs/` or the project `out/`.
- Recovery-net gotchas (session 23): `.tran` Tstart shifts the raw/meas time
  axis to 0 — `.meas TD` must be RELATIVE (the original bench's TD=0 worked
  by accident); `Roff=tau` unbraced works in 26.0.1 but write `{tau}` for
  portability; **gear required** (trapezoidal rings badly on the snap); B3
  needs the +1 µΩ floor + R2 = 1 µΩ in any bridge topology; Ki–A cap
  placement is mandatory (D–S crashes re-turn-on, Ki–Ai singular); ideal-
  switch benches show a ~1.2 kA Cds dump impulse at midpoint swings — real
  loop inductance (~1–5 nH) domesitcs it; LTspice `.raw` "real" flag =
  float32 records (24 var × 4 B + 2 B flag + 2 B pad in the stepped file —
  parse carefully if hand-parsing).
- `fit_imcq120_kelvin.py` IMPORTS `fit_imcq120_joint` and PATCHES its
  module-global `params_at` (kelvin-aware) + extends `FREE_SPEC`/`VARIANTS`
  — importing the kelvin module is what activates the patch (the probe
  script does exactly that). `JF.acceptance` needs the base dict to carry
  `_rd_T/_rd_vg/_rd_chord` (from `JointData.rd_*`) — a KeyError otherwise.
- The joint script's printed `sigma_rel` is the GN marginal sd in VECTOR
  space (log params → relative; linear params → RAW units) and is INFLATED
  1–2 orders for linear knobs on this problem (session-22 calibration:
  `s_vth` printed "1σ" refit explodes RMSE 1.93 → 9.48 %; true σ ≈ 2–5 % of
  the slope; `r_d25` −0.5 mΩ fixed → 2.91 %). For linear-knob uncertainty
  quote FIX-AND-REFIT penalties, not the print. The joint report §2 carries
  the correction note; the record's "σ 0.5 %" for s_vth was a misread.
- The joint script's bounds table is SCRIPT-LOCAL (`FREE_SPEC`), mirroring
  `params.vector_bounds` semantics (log params mapped through log). r_d25
  upper bound = datasheet 25 °C chord (3.7 mΩ) as the degenerate-branch
  guard; the kelvin script caps r_s25 at the same chord and parameterizes
  `r_sk = f_sk·r_s25` with `f_sk ∈ [0,1]` (bounds cannot reference another
  parameter; the model guard 0 ≤ r_sk ≤ r_s is the hard line).
- The OFFICIAL record's TC-law evaluation lives in `reference.py::ImcqJointRecord
  .params_at` AND in the joint script's dict-based `params_at` — deliberately
  two implementations, cross-checked every script run (OK/STALE) and by the
  pinned test. Keep them in sync (session 22 added r_sk to both).
  The drift check compares SENSED parameters only (session-19 lesson).
- Transfer CSV: header `index,V_gs,I_d,condition` — `load_csv()` does NOT
  read it; the script-local loader parses `condition` (Qg-style). The
  175/200 °C curves share the crossing point (12.025 V, 999.1 A) once per
  curve; no data below 12 V for them (digitization limit, kept not cleaned).
- `openpyxl` available; `names.xlsx` emits a harmless Data Validation
  warning. names.xlsx now carries the r_sk row (user-added, uncommitted at
  snapshot time along with the session-22 code). The recovery parameters
  (q_0/tau/k_chg/r_max) are deliberately NOT in names.xlsx — LTspice-side
  netlist values, outside the package naming convention.
- pF·V = pC, not C. Plateau split at V_gp (v_gd = 0), not ΔV_gs.
  `PchipInterpolator` has no `derivatives=` kwarg in the installed scipy.
- LTspice cap-export gotchas (never `table()` in a B source or in `Q=` —
  garbage symbolic derivative ~1e11 A; use the incremental ramp-sum;
  `.meas tran FIND` strictly inside staircase segments) — full record in
  `CAP_SPEC.md`.
- Scratch (`tooling/workspace/outputs/vdmos_recov/`): all recovery benches —
  `DiodeTest_base/500p/sweep/fine/integ/integ2/didt100/fixcheck.net`,
  `HalfBridge_conv.net`, `cds_expr.txt`. Project `out/` additions
  (session 23): `ISC040N10NM8_recov.subckt`, `RECOV_REPORT.md`. `out/` is
  gitignored.
- README slimmed to a pure model reference (session 23, LaTeX display math).
  Three implementation details that lived ONLY in README prose are parked
  here (all pinned by tests — the prose is gone, the behavior is not):
  (1) CLI: `--show`/`--out/-o` are registered on every subcommand with
  `argparse.SUPPRESS` defaults so they parse before *or* after the
  subcommand — without SUPPRESS the subparser resets them unconditionally
  (`test_cli.py`). (2) `make_synthetic` noise is multiplicative and applied
  WITHOUT clipping the values — flooring the current distorts subthreshold
  up to 2×; near-zero currents are handled in the residual denominator
  (`i_d_floor` in `residuals`), not in the data. (3) `k_p` is NOT the
  Level-1 Kp — it must be re-extracted when k_theta/v_vsat are added
  (28 → ~34.7 at DEFAULT); `v_kneemin` is a numerical floor, keep fixed,
  never fitted; `v_ovmax` lower bound (0.05 V) sits deliberately above 0 —
  between 0 and 0.05 V the cap is a hard clip at nearly zero volts.
- Test totals drift (Tk digitizer tests self-skip headless) — don't pin
  counts. Session 23 touched NO package code (docs + `out/` only), so all
  model anchors are untouched by construction.
