"""Measured / digitised I-V datasets and synthetic test data."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .model import vdmos_id
from .params import Params

__all__ = ["CurveSet", "load_csv", "save_csv", "make_synthetic", "example_dataset"]


@dataclass(frozen=True)
class CurveSet:
    """Flat arrays of (Vgs, Vds, Id) triplets -- one row per measured point.

    Flat form (rather than 2-D grids) is deliberate: datasheet curves digitised
    from graphs are ragged, and different curves can have different lengths.
    """

    vgs: np.ndarray
    vds: np.ndarray
    ids: np.ndarray
    label: str = ""

    def __post_init__(self) -> None:
        n = len(self.vgs)
        if not (len(self.vds) == len(self.ids) == n):
            raise ValueError("vgs, vds and ids must have the same length")

    def __len__(self) -> int:
        return len(self.vgs)

    def subset(self, mask: np.ndarray) -> CurveSet:
        return CurveSet(self.vgs[mask], self.vds[mask], self.ids[mask], self.label)

    def vgs_grid(self, tol: float = 1e-9) -> list[float]:
        return sorted({round(float(v), 9) for v in self.vgs})

    def vds_grid(self, tol: float = 1e-9) -> list[float]:
        return sorted({round(float(v), 9) for v in self.vds})


def load_csv(path: str | Path) -> CurveSet:
    """Read a CSV with header `vgs,vds,ids` (column order flexible, case-insensitive)."""
    path = Path(path)
    with path.open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError(f"{path} has no data rows")
    cols = {k.strip().lower(): k for k in rows[0] if k}
    for need in ("vgs", "vds", "ids"):
        if need not in cols:
            raise ValueError(f"{path} missing column '{need}'; header is {list(rows[0])}")
    get = lambda r, c: float(r[cols[c]])  # noqa: E731
    return CurveSet(
        np.array([get(r, "vgs") for r in rows]),
        np.array([get(r, "vds") for r in rows]),
        np.array([get(r, "ids") for r in rows]),
        label=path.stem,
    )


def save_csv(data: CurveSet, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["vgs", "vds", "ids"])
        w.writerows(zip(data.vgs, data.vds, data.ids, strict=True))


def make_synthetic(
    truth: Params,
    vgs_vals: np.ndarray,
    vds_vals: np.ndarray,
    *,
    noise: float = 0.01,
    seed: int = 0,
) -> CurveSet:
    """Full grid of (Vgs, Vds) points from `truth`, with multiplicative noise.

    Noise is multiplicative and applied *without* clipping: near-zero currents
    are handled in the residual denominator (`id_floor` in
    :func:`vdmos_fit.fit.residuals`), not in the data -- rationale in README,
    "Synthetic data". Used to verify the fitter recovers known parameters.
    """
    vg, vd = np.meshgrid(np.asarray(vgs_vals, float), np.asarray(vds_vals, float), indexing="ij")
    vg, vd = vg.ravel(), vd.ravel()
    ids = vdmos_id(vd, vg, truth)
    rng = np.random.default_rng(seed)
    ids = ids * (1.0 + noise * rng.standard_normal(ids.size))
    return CurveSet(vg, vd, ids, label="synthetic")


def example_dataset(noise: float = 0.01, seed: int = 1) -> CurveSet:
    """Typical power-FET static grid: transfer curves at several Vds plus output curves."""
    vgs = np.linspace(2.0, 14.0, 25)
    vds = np.array([0.05, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 40.0])
    truth = Params(
        kp=30.0,
        vth0=3.1,
        theta=0.02,
        vsat=1.2,
        eta=0.03,
        lam=0.015,
        gk=0.3,
        vthk=1.5,
    )
    return make_synthetic(truth, vgs, vds, noise=noise, seed=seed)
