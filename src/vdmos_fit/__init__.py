"""Improved VDMOS compact model + static parameter extraction."""

from .data import CurveSet, example_dataset, load_csv, make_synthetic, save_csv
from .fit import (
    FitResult,
    Identifiability,
    fit,
    identifiability,
    parameter_uncertainty,
    residuals,
    sensitivity_report,
)
from .model import vdmos_id, vdmos_id_level1
from .params import (
    BOUNDS,
    DEFAULT,
    DEFAULT_FREE,
    LOG_PARAMS,
    Params,
    pack,
    unpack,
    vector_bounds,
)
from .reference import (
    ST8L65_FIT_A_FREE,
    ST8L65_OFFICIAL,
    ST8L65_OFFICIAL_RMSE_REL,
)

__version__ = "0.1.0"

__all__ = [
    "BOUNDS",
    "CurveSet",
    "DEFAULT",
    "DEFAULT_FREE",
    "FitResult",
    "Identifiability",
    "LOG_PARAMS",
    "Params",
    "ST8L65_FIT_A_FREE",
    "ST8L65_OFFICIAL",
    "ST8L65_OFFICIAL_RMSE_REL",
    "example_dataset",
    "fit",
    "identifiability",
    "load_csv",
    "make_synthetic",
    "parameter_uncertainty",
    "pack",
    "residuals",
    "save_csv",
    "sensitivity_report",
    "unpack",
    "vdmos_id",
    "vdmos_id_level1",
    "vector_bounds",
    "__version__",
]
