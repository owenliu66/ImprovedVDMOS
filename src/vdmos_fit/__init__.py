"""Improved VDMOS compact model + static parameter extraction."""

from .data import (
    CurveSet,
    DiodeCurveSet,
    example_dataset,
    load_csv,
    load_diode_csv,
    make_synthetic,
    make_synthetic_diode,
    save_csv,
    save_diode_csv,
)
from .diode import vdmos_diode
from .fit import (
    FitResult,
    Identifiability,
    diode_residuals,
    diode_sensitivity_report,
    fit,
    fit_diode,
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
    ST8L65_DIODE,
    ST8L65_DIODE_FREE,
    ST8L65_DIODE_RMSE_REL,
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
    "DiodeCurveSet",
    "FitResult",
    "Identifiability",
    "LOG_PARAMS",
    "Params",
    "ST8L65_DIODE",
    "ST8L65_DIODE_FREE",
    "ST8L65_DIODE_RMSE_REL",
    "ST8L65_FIT_A_FREE",
    "ST8L65_OFFICIAL",
    "ST8L65_OFFICIAL_RMSE_REL",
    "diode_residuals",
    "diode_sensitivity_report",
    "example_dataset",
    "fit",
    "fit_diode",
    "identifiability",
    "load_csv",
    "load_diode_csv",
    "make_synthetic",
    "make_synthetic_diode",
    "parameter_uncertainty",
    "pack",
    "residuals",
    "save_csv",
    "save_diode_csv",
    "sensitivity_report",
    "unpack",
    "vdmos_diode",
    "vdmos_id",
    "vdmos_id_level1",
    "vector_bounds",
    "__version__",
]
