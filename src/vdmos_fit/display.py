"""GUI number formatting (session 30, Phase E).

The user's displayNum helper, FIXED: the original took log10 of the raw
value, so negative and zero entries crashed (math domain error). Sign is
handled separately and abs() feeds the magnitude; zero short-circuits.
"""

from __future__ import annotations

import math

_PREFIXES = [(1e12, "T"), (1e9, "G"), (1e6, "M"), (1e3, "k"), (1.0, ""),
             (1e-3, "m"), (1e-6, "\u00b5"), (1e-9, "n"), (1e-12, "p"),
             (1e-15, "f")]


def displayNum(x: float, sig: int = 4) -> str:
    """Engineering-suffix display: 0.00248295 -> '2.483m', -123456 ->
    '-123.5k', 0 -> '0'. Negatives/zero are safe (the crash fix)."""
    x = float(x)
    if math.isnan(x):
        return "nan"
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if x == 0.0:
        return "0"
    sign = "-" if x < 0.0 else ""
    a = abs(x)
    for factor, suffix in _PREFIXES:
        if a >= factor:
            return f"{sign}{a / factor:.{sig}g}{suffix}"
    return f"{sign}{a:.{sig - 1}e}"      # below femto: scientific


_SUFFIX = {suf: fac for fac, suf in _PREFIXES if suf}


def parseEng(text: str) -> float:
    """Inverse of displayNum: '2.483m' -> 0.002483, '-123.5k' -> -123500,
    '1e-3' -> 0.001, '0' -> 0. Accepts u as micro; raises ValueError on
    garbage (the GUI manual-edit path)."""
    t = text.strip().replace("\u03bc", "\u00b5").replace("u", "\u00b5")
    if not t:
        raise ValueError("empty number")
    sign = 1.0
    if t[0] in "+-":
        sign = -1.0 if t[0] == "-" else 1.0
        t = t[1:]
    if t and t[-1] in _SUFFIX and "e" not in t and "E" not in t:
        num, fac = t[:-1], _SUFFIX[t[-1]]
    else:
        num, fac = t, 1.0
    return sign * float(num) * fac
