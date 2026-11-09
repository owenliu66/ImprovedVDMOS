"""displayNum (session 30 Phase E): the negative/zero crash fix pins, and
parseEng -- the manual-edit inverse (double-click a value cell in the GUI).
"""

import pytest

from vdmos_fit.display import displayNum, parseEng


@pytest.mark.parametrize(("x", "want"), [
    (0.0, "0"),
    (-0.0, "0"),
    (2.48295e-3, "2.483m"),
    (-2.48295e-3, "-2.483m"),      # the log10 domain-error crash fix
    (-123456.0, "-123.5k"),
    (123456.0, "123.5k"),
    (1.0, "1"),
    (-1.0, "-1"),
    (1e-9, "1n"),
    (-7.4863e-3, "-7.486m"),
    (0.025, "25m"),
    (5.0, "5"),
    (2.0e-18, "2.000e-18"),        # below femto: scientific, no crash
])
def test_displayNum(x, want):
    assert displayNum(x) == want


def test_specials():
    assert displayNum(float("inf")) == "inf"
    assert displayNum(float("-inf")) == "-inf"
    assert displayNum(float("nan")) == "nan"


@pytest.mark.parametrize(("txt", "want"), [
    ("0", 0.0),
    ("2.483m", 2.483e-3),
    ("-2.483m", -2.483e-3),
    ("-123.5k", -123500.0),
    ("1k", 1000.0),
    ("1e-3", 1e-3),          # scientific: no suffix interpretation
    ("-1.5", -1.5),
    ("+3", 3.0),
    ("2u", 2e-6),            # u accepted as micro
    ("1n", 1e-9),
    (" 5M ", 5e6),
    ("1T", 1e12),
])
def test_parseEng(txt, want):
    assert parseEng(txt) == pytest.approx(want)


def test_parseEng_roundtrip():
    for x in (0.0, 2.48295e-3, -123456.0, 1e-9, 5.0, 2e-18, -7.4863e-3):
        assert parseEng(displayNum(x)) == pytest.approx(x, rel=1e-3)


def test_parseEng_garbage():
    for bad in ("", "abc", "m", "1.2.3", "k5"):
        with pytest.raises(ValueError):
            parseEng(bad)
