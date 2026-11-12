"""ISC040 channel, GUI protocol (load_tables + fit_channel_ref), velocity-
saturation family, knee-form head-to-head via the new start_override path.

start_override carries k_kneeform (a knob the preset never knew) through the
fit verbatim as the first start -- exactly what the GUI now does -- so the
scale form no longer gets silently reset to the preset default (mobility).
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import vdmos_id  # noqa: E402
from vdmos_fit.ingest import IVTable, RdsSpec  # noqa: E402
from vdmos_fit.presets import PRESETS  # noqa: E402
from vdmos_fit.release import load_tables, merge_ivs  # noqa: E402
from vdmos_fit.stages import fit_channel_ref  # noqa: E402

DATA = ROOT / "data"
PATHS = [DATA / "ISC040N10NM8_IV_25C.csv",
         DATA / "ISC040N10NM8_Rdson_Id_25C.csv",
         DATA / "ISC040N10NM8_Rdson_Vgs.csv",
         DATA / "ISC040N10NM8_Vth_Tj.csv"]

# velocity-saturation family WITH the DIBL cap free (the GUI family the user
# fits): k_theta, v_vsat, k_kneesharp and v_diblmax/k_diblsharp free; k_lambda
# left at 0 to match the experiment comparison.
FREE = ("k_p", "v_th0", "v_sub", "v_ovmax", "k_ovsharp", "k_theta", "v_vsat",
        "k_kneesharp", "k_eta", "v_diblmax", "k_diblsharp", "v_dbsat", "r_d")


def main():
    tables, _report = load_tables(PATHS, 25.0)
    iv = merge_ivs([t for t in tables if isinstance(t, IVTable)])
    rds = next(t for t in tables if isinstance(t, RdsSpec))
    p = PRESETS["trench"]

    def base_start(kneeform, ktheta, vsat):
        s = {**p.channel_fixed, **p.channel_start}
        s.update(k_kneeform=kneeform, k_theta=ktheta, v_vsat=vsat,
                 v_diblmax=7.0, k_diblsharp=0.5)
        return s

    def run(tag, kneeform, ktheta, vsat, kt_lo, kt_hi):
        so = base_start(kneeform, ktheta, vsat)
        res = fit_channel_ref(iv, rds, p, ref_temp=25.0, free_override=FREE,
                              bounds_override={"k_theta": (kt_lo, kt_hi)},
                              start_override=so)
        q = res.params
        print(f"=== {tag} ===  RMSE {res.rmse_rel*100:.3f}%  conv {res.success}")
        print(f"  k_theta {q.k_theta:.5g}  v_vsat {q.v_vsat:.4g}  "
              f"v_ovmax {q.v_ovmax:.4g}  k_kneesharp {q.k_kneesharp:.4g}  "
              f"v_dbsat {q.v_dbsat:.4g}")
        print(f"  off i_d(V_gs=0, 80V) = {float(vdmos_id(80.0, 0.0, q)):.3g} A"
              f"   (data 1e-7)")
        print(f"  on_bounds: {res.on_bounds()}")
        return res

    vom = base_start(0.0, 1e-6, 5.0)["v_ovmax"]
    run("MOBILITY form (k_kneeform=0, k_theta free to -1/v_ovmax)",
        0.0, 1e-6, 5.0, -1.0 / vom, 5.0)
    run("SCALE form (k_kneeform=1, v_knee = k_theta*v_ov_eff2)",
        1.0, 1.0, 5.0, 0.0, 5.0)


if __name__ == "__main__":
    main()
