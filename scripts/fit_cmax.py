"""Extract c_gdmax from the digitised gate-charge curves (CAP_SPEC protocol).

Data: data/<part>_Cap_Vds.csv + data/<part>_Vgs_Qg.csv for the three parts.

PROTOCOL (CAP_SPEC.md "Cmax protocol", corrected split -- the plateau ENDS at
v_gd = 0, so the plateau charge is a pure table-branch VALIDATION and carries
no c_gdmax information; c_gdmax is extracted from the post-plateau segment):

  * segment 1 (0 -> plateau start): Qgs check against the c_gs = iss - rss
    table (limitation 4: table shape used against v_gs).
  * segment 2 (Miller plateau): Qgd check against int C_rss dv + c_gs tilt
    term. A POSITIVE deficit = channel charge the off-state table cannot see
    (limitation 1) -- expected to be worst on the super-junction part, whose
    plateau is nearly flat (huge effective Miller capacitance).
  * segment 3 (post-plateau, v_gd > 0): remainder charge after the c_gs term
    is q_exp(V_end - i_d*r_ds_on; c_gdmax); solve 1-D, bracket
    [c_gd0, 2*c_gs(0)]. Non-bracketed -> clamp + diagnostic (spec step 4).

r_ds_on (ohmic tail voltage i_d*r_ds_on) is a datasheet-class assumption,
documented per part; a V_on = 0 probe reports the sensitivity.

Quote the CURVE, not the split: c_gdmax is correlated with the origin slope
s_gd0 (the noisiest number in the model) -- see CAP_SPEC caveats.

Outputs (out/):
  fit_cmax_params.txt    per-part/per-condition extraction records
  fit_cmax_cgd.png       composite c_gd(v_gd) per part at the extracted cmax

    python scripts/fit_cmax.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import (  # noqa: E402
    CapModel,
    build_cap_model,
    extract_cmax,
    load_cap_csv,
    load_qg_csv,
)

# Datasheet-class Rds(on) for the ohmic tail voltage (documented assumptions):
#   ST8L65N044M9: 36 mOhm chord from the digitised IV anchors (data/README).
#   IMCQ120R004M2H: 4 mOhm typ (part-code class "004").
#   ISC040N10NM8: 40 mOhm (part-code class "040").
PARTS = {
    "ST8L65N044M9": {"r_ds_on": 0.036},
    "IMCQ120R004M2H": {"r_ds_on": 0.004},
    "ISC040N10NM8": {"r_ds_on": 0.040},
}


def run_part(name: str, r_ds_on: float) -> tuple[CapModel, list]:
    caps = load_cap_csv(ROOT / "data" / f"{name}_Cap_Vds.csv")
    qg = load_qg_csv(ROOT / "data" / f"{name}_Vgs_Qg.csv")
    model = build_cap_model(caps)
    fits = list(extract_cmax(model, qg, r_ds_on=r_ds_on))
    probe = extract_cmax(model, qg, r_ds_on=0.0)  # V_on sensitivity probe
    return model, fits, probe


def main() -> None:
    out_dir = ROOT / "out"
    out_dir.mkdir(exist_ok=True)
    lines: list[str] = []
    fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
    extracted: dict[str, float] = {}

    for ax, (name, cfg) in zip(axs, PARTS.items(), strict=True):
        model, fits, probe = run_part(name, cfg["r_ds_on"])
        a = model.anchor
        c_gs0 = model.c_gs0()
        lines.append(f"=== {name} ===")
        lines.append(
            f"anchor: c_gd0={a.c_gd0*1e12:.1f} pF  s_gd0={a.s_gd0*1e12:.1f} pF/V  "
            f"b2={a.b2*1e12:.1f} pF/V^2   seed c_gs(0)={c_gs0*1e12:.1f} pF   "
            f"bracket [{a.c_gd0*1e12:.0f}, {2*c_gs0*1e12:.0f}] pF"
        )
        for fit, pr in zip(fits, probe, strict=True):
            raw = "clamped" if fit.c_gdmax_raw is None else f"{fit.c_gdmax_raw*1e12:.1f} pF"
            lines.append(
                f"  I_d={fit.i_d:g} A  V_DS={fit.v_ds:g} V  "
                f"V_gp=[{fit.v_gp_lo:.2f},{fit.v_gp_hi:.2f}] V"
            )
            lines.append(
                f"    Qgs  meas {fit.qgs_meas*1e9:7.2f}  model {fit.qgs_model*1e9:7.2f} nC"
                f"   ({(fit.qgs_model-fit.qgs_meas)/fit.qgs_meas*100:+.1f} %)"
            )
            lines.append(
                f"    Qgd  meas {fit.qgd_meas*1e9:7.2f}  model {fit.qgd_model*1e9:7.2f} nC"
                f"   deficit {fit.plateau_deficit*1e9:+7.2f} nC"
                f"   ({fit.plateau_deficit/fit.qgd_meas*100:+.1f} %)"
            )
            lines.append(
                f"    seg3 meas {fit.q_seg3_meas*1e9:7.2f}  cgs model {fit.cgs_seg3_model*1e9:7.2f}"
                f"  remainder {fit.remainder*1e9:7.2f} nC"
            )
            lines.append(
                f"    c_gdmax raw {raw}  -> {fit.c_gdmax*1e12:.1f} pF"
                f"   {fit.clamp or ''}   (V_on=0 probe: {pr.c_gdmax*1e12:.1f} pF {pr.clamp or ''})"
            )
        good = [f for f in fits if f.c_gdmax_raw is not None]
        chosen = float(np.median([f.c_gdmax for f in good])) if good else fits[0].c_gdmax
        extracted[name] = chosen
        lines.append(f"  record c_gdmax = {chosen*1e12:.1f} pF"
                     f"  ({'extracted' if good else 'clamped ' + str(fits[0].clamp)})")
        lines.append("")

        # composite c_gd(v_gd) at the extracted value
        m = CapModel(model.tables, model.anchor,
                     type(model.params)(a.c_gd0, a.s_gd0, chosen))
        v_lo = np.linspace(-min(model.tables.rss.v[-1], 600), 0, 400)
        v_hi = np.linspace(0, 25, 200)
        ax.plot(v_lo, m.c_gd(v_lo) * 1e12, lw=1.8)
        ax.plot(v_hi, m.c_gd(v_hi) * 1e12, lw=1.8)
        ax.plot(0, a.c_gd0 * 1e12, "k*", ms=10)
        ax.set_xlim(-60, 25)
        ax.set_title(f"{name}: c_gd(v_gd), c_gdmax={chosen*1e12:.0f} pF")
        ax.set_xlabel("v_gd [V]")
        ax.set_ylabel("c_gd [pF]")
        ax.grid(alpha=0.3)

    txt = "\n".join(lines)
    print(txt)
    (out_dir / "fit_cmax_params.txt").write_text(txt, encoding="utf-8")
    fig.tight_layout()
    fig.savefig(out_dir / "fit_cmax_cgd.png", dpi=130)
    print(f"saved: {out_dir / 'fit_cmax_params.txt'}\n       {out_dir / 'fit_cmax_cgd.png'}")
    print("\nextracted c_gdmax [F]:", {k: f"{v:.4e}" for k, v in extracted.items()})


if __name__ == "__main__":
    main()
