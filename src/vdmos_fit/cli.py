"""Command-line entry point: `vdmos-fit`.

Usage
-----
    vdmos-fit check                       # model sanity checks + plots
    vdmos-fit fit data.csv -o out/        # extract parameters from a CSV
    vdmos-fit demo                        # synthetic round-trip demo

The shared flags (--show, --out/-o) may appear before OR after the subcommand:
`vdmos-fit --out x check` and `vdmos-fit check --out x` are equivalent.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib
import numpy as np

from .data import example_dataset, load_csv
from .fit import fit, identifiability, sensitivity_report
from .model import vdmos_id, vdmos_id_level1
from .params import DEFAULT
from .plotting import gds_continuity, mtriode_shape, output_curves, overlay_data, transfer_curves


def _save(fig, outdir: Path, name: str, show: bool) -> None:
    if outdir != Path("."):
        outdir.mkdir(parents=True, exist_ok=True)
        path = outdir / name
        fig.savefig(path, dpi=150, bbox_inches="tight")
        print(f"wrote {path}")
    if show:
        matplotlib.pyplot.show()


def cmd_check(args) -> int:
    print(
        DEFAULT.summary()
        if hasattr(DEFAULT, "summary")
        else "params: " + ", ".join(f"{k}={v:g}" for k, v in DEFAULT.as_dict().items())
    )
    ron = 0.05 / float(vdmos_id(0.05, 10.0, DEFAULT))
    print(f"Ron(Vgs=10V) = {ron * 1e3:.2f} mOhm   (LTspice A/B reference: 4.27 mOhm)")
    for vg in (4, 6, 8, 10, 12):
        a = float(vdmos_id(20.0, vg, DEFAULT))
        b = float(vdmos_id_level1(20.0, vg))
        print(f"  Vgs={vg:2d}  improved={a:7.1f} A   Level-1={b:8.1f} A   ratio={a / b:5.2f}")
    vg = np.linspace(6, 12, 7)
    vov = 0.05 * np.log1p(np.exp((vg - 3.0) / 0.05))
    alpha = np.polyfit(np.log(vov), np.log(vdmos_id(np.full_like(vg, 20.0), vg, DEFAULT)), 1)[0]
    print(f"alpha = {alpha:.2f}   (Level-1 gives 2.00)")

    out = Path(args.out)
    for fn, name in (
        (output_curves, "output_curves.png"),
        (transfer_curves, "transfer_curves.png"),
        (gds_continuity, "gds_continuity.png"),
        (mtriode_shape, "mtriode_shape.png"),
    ):
        fig = fn().figure
        _save(fig, out, name, args.show)
    return 0


def cmd_fit(args) -> int:
    data = load_csv(args.csv)
    print(
        f"loaded {args.csv}: {len(data)} points, "
        f"Vgs {data.vgs.min():.2f}..{data.vgs.max():.2f} V, "
        f"Vds {data.vds.min():.2f}..{data.vds.max():.2f} V"
    )
    res = fit(data, start=DEFAULT, weight=args.weight, verbose=args.verbose)
    print(res.summary())
    print()
    print(identifiability(res).summary())
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "fitted_params.txt").write_text(res.summary(), encoding="utf-8")
    _save(overlay_data(data, res.params).figure, out, "fit_vs_data.png", args.show)
    return 0 if res.success else 1


def cmd_demo(args) -> int:
    data = example_dataset()
    print("pre-fit sensitivity (mean |d ln Id / d ln p|):")
    for k, v in sensitivity_report(DEFAULT, data).items():
        print(f"  {k:<8s} {v:10.3g}")
    print()
    res = fit(data, start=DEFAULT)
    print(res.summary())
    print()
    print(identifiability(res).summary())
    _save(overlay_data(data, res.params).figure, Path(args.out), "demo_fit.png", args.show)
    return 0


def _build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(prog="vdmos-fit", description=__doc__)
    ap.add_argument("--show", action="store_true", help="open interactive plots")
    ap.add_argument("-o", "--out", default="out", help="output directory (default: out/)")
    sub = ap.add_subparsers(dest="cmd", required=True)

    # Shared flags repeated on every subcommand so they parse in either position.
    # SUPPRESS defaults are essential: without them the subparser would reset
    # --out/--show to their defaults and clobber a value already parsed before
    # the subcommand (argparse fills subparser defaults into the namespace
    # unconditionally).
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("--show", action="store_true", default=argparse.SUPPRESS,
                        help="open interactive plots")
    common.add_argument("-o", "--out", default=argparse.SUPPRESS,
                        help="output directory (default: out/)")

    sub.add_parser("check", parents=[common], help="model sanity checks + reference plots")

    fp = sub.add_parser("fit", parents=[common], help="extract parameters from a measured CSV")
    fp.add_argument("csv", help="CSV with columns vgs,vds,ids")
    fp.add_argument("--weight", choices=("relative", "abs"), default="relative")
    fp.add_argument("-v", "--verbose", action="store_true")

    sub.add_parser("demo", parents=[common], help="synthetic round-trip fit")
    return ap


def main(argv: list[str] | None = None) -> int:
    args = _build_parser().parse_args(argv)
    return {"check": cmd_check, "fit": cmd_fit, "demo": cmd_demo}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
