"""Release CLI (session 30): CSVs -> fitted records -> subckt in one shot.

    python scripts/build_part.py MYFET100 data/MYFET100_*.csv ^
        --preset trench --ref 25 [--rdson-itest 25] [--recovery]

Presets: trench | sj | sic. --recovery runs the LTspice-grid extraction
(needs a REC CSV: Q_rr, t_rr, I_RM, I_sd, di/dt columns, and LTspice
installed); without it the FAST subckt exports with a plain diode and no
tt. Outputs: out/<PART>_full|fast.subckt, out/<PART>_caps.subckt,
out/<PART>_record.json, data/<PART>_cds_expr.txt.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from vdmos_fit.release import build_part  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("part", help="part name (file prefix + subckt name)")
    ap.add_argument("csvs", nargs="+", help="input CSVs (IV/CAP/QG/RDS/VTH/REC)")
    ap.add_argument("--preset", default="trench",
                    choices=("trench", "sj", "sic"))
    ap.add_argument("--ref", type=float, default=25.0,
                    help="reference temperature [degC]")
    ap.add_argument("--rdson-itest", type=float, default=None,
                    help="default Rds(on) chord test current [A]")
    ap.add_argument("--recovery", action="store_true",
                    help="extract recovery knobs (needs REC CSV + LTspice)")
    ap.add_argument("--v_r", type=float, default=50.0,
                    help="recovery test reverse voltage [V]")
    ap.add_argument("--out", default=str(ROOT / "out"))
    args = ap.parse_args()

    res = build_part(args.csvs, args.part, preset=args.preset,
                     ref_temp=args.ref, rdson_i_test=args.rdson_itest,
                     out_dir=args.out, with_recovery=args.recovery,
                     v_r=args.v_r)
    print(f"== {args.part} ({args.preset}, ref {args.ref:g} C) ==")
    for k, v in res.rmse.items():
        print(f"  rmse {k}: {v:.4g}")
    for w in res.warnings:
        print(f"  note: {w}")
    for f in res.files:
        print(f"  wrote {f}")


if __name__ == "__main__":
    main()
