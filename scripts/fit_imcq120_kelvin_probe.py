"""Constrained probe: force gate-loop degeneration r_g = r_s (r_sk = 0) at
fixed values and refit everything else. Quantifies what the 'datasheet static
measurement did NOT use the kelvin pin' hypothesis costs in fit quality.

Session 22 companion to fit_imcq120_kelvin.py (run after the main experiment).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src"))

import fit_imcq120_joint as JF  # noqa: E402  (official joint machinery)
import fit_imcq120_kelvin as K  # noqa: E402  (import patches JF.params_at + registry)
from vdmos_fit.reference import IMCQ120_JOINT  # noqa: E402

start = {n: float(getattr(IMCQ120_JOINT, n)) for n in K.REC_FIELDS}
data = JF.JointData(include=K.ALL_GROUPS)
start.update({"ea_bkt": 250.0, "gamma_rd": 3.3,
              "_rd_T": data.rd_T, "_rd_vg": data.rd_vg,
              "_rd_chord": data.rd_chord})

for r_g in (0.25e-3, 0.5e-3, 1.0e-3):
    name = f"probe_rg_{r_g * 1e3:.2f}"
    JF.VARIANTS[name] = {"free": JF.PRIMARY_FREE, "rd_law": "quad",
                         "bkt_law": "lin", "groups": K.ALL_GROUPS}
    st = {**start, "r_s25": r_g, "f_sk": 0.0}  # r_s25 NOT free -> fixed r_g
    r = JF.fit_variant(name, st, seeds=None)
    acc = JF.acceptance(r["base"], "quad", "lin")
    gap = K.gap_report(r["base"], r)
    print(f"\n=== {name}: r_g = {r_g * 1e3:.2f} mOhm FIXED ===")
    print(f"   RMSE {r['rmse_overall']:.2f}% (iv {r['rmse_group']['iv']:.2f} "
          f"trf {r['rmse_group']['transfer']:.2f} rdson {r['rmse_group']['rdson']:.2f} "
          f"vth {r['rmse_group']['vth']:.2f})")
    print(f"   g_fs {acc['g_fs_S']:.1f} S ({(acc['g_fs_S'] - 122) / 122 * 100:+.1f} % vs spec)"
          f"   corner25 {gap['corner25_i18_A']:.0f} A "
          f"({gap['corner25_i18_A'] / 3213 * 100 - 100:+.1f} % vs 3213)"
          f"   corner25 res mean {gap['corner25_res_mean_pct']:+.2f} %")
    print(f"   corner175 {acc['corner_18V']['175C']:.0f} A   rdson RMSE {acc['rdson_rmse']:.2f}%"
          f"   chords25 {acc['ron_25_vg18'] * 1e3:.3f}/{acc['ron_25_vg15'] * 1e3:.3f} mO")
    print(f"   r_d25 {r['base']['r_d25'] * 1e3:.4f} mO  k_p25 {r['base']['k_p25']:.3f}"
          f"  k_bt25 {r['base']['k_bt25']:.4f}  k_theta {r['base']['k_theta']:.4f}")
