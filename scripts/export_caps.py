"""Export the capacitance delta of the three parts as LTspice subckts.

Each subckt (pins G D S) implements the three-capacitor delta as
charge-master B sources -- I = d(q(v))/dt with q(v) an LTspice table()
sampled from the Python CapModel (see vdmos_fit/ltspice.py for the
rationale: piecewise-linear q -> piecewise-constant C staircase, clamp
policy carried in the table samples).

The c_gdmax used is the EXTRACTED record from reference.py (CAP_*), not
the c_gs(0) seed.

Outputs (out/):
  <part>_caps.subckt    .subckt <part>_caps G D S -- three B caps

    python scripts/export_caps.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit import (  # noqa: E402
    CAP_IMCQ120,
    CAP_ISC040,
    CAP_ST8L65,
    CapModel,
    build_cap_model,
    cap_charge_tables,
    load_cap_csv,
    write_cap_subckt,
)

PARTS = {
    "ST8L65N044M9": CAP_ST8L65,
    "IMCQ120R004M2H": CAP_IMCQ120,
    "ISC040N10NM8": CAP_ISC040,
}


def main() -> None:
    out_dir = ROOT / "out"
    for name, params in PARTS.items():
        model = build_cap_model(load_cap_csv(ROOT / "data" / f"{name}_Cap_Vds.csv"))
        # re-splice at the extracted c_gdmax record (not the seed)
        model = CapModel(model.tables, model.anchor, params)
        path = write_cap_subckt(model, out_dir / f"{name}_caps.subckt",
                                name=f"{name.lower()}_caps")
        qs = cap_charge_tables(model)
        sizes = {k: len(v) for k, (v, _) in qs.items()}
        spans = {k: (f"{v[0]:.1f}", f"{v[-1]:.1f}") for k, (v, _) in qs.items()}
        print(f"{name}: {path.name}  table points {sizes}  v-spans {spans}")
        print(f"  c_gd0={params.c_gd0*1e12:.1f} pF  c_gdmax={params.c_gdmax*1e12:.1f} pF"
              f"  k_agd={params.k_agd:.4g} 1/V")


if __name__ == "__main__":
    main()
