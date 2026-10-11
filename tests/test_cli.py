"""CLI argument parsing: shared flags must work before OR after the subcommand.

Parser-level tests only -- no plots, no fitting -- so they stay fast.
"""

from __future__ import annotations

import pytest

from vdmos_fit.cli import _build_parser


def test_out_flag_after_subcommand():
    args = _build_parser().parse_args(["check", "--out", "x"])
    assert args.cmd == "check" and args.out == "x" and args.show is False


def test_out_flag_before_subcommand():
    args = _build_parser().parse_args(["--out", "x", "check"])
    assert args.out == "x"


def test_out_default():
    args = _build_parser().parse_args(["check"])
    assert args.out == "out"


def test_subparser_does_not_clobber_pre_subcommand_value():
    # The SUPPRESS-default pattern: a value set before the subcommand must
    # survive the subparser's (suppressed) default.
    args = _build_parser().parse_args(["--show", "--out", "pre", "demo"])
    assert args.show is True and args.out == "pre"


def test_short_o_alias_on_both_sides():
    assert _build_parser().parse_args(["-o", "a", "demo"]).out == "a"
    assert _build_parser().parse_args(["fit", "d.csv", "-o", "b"]).out == "b"


def test_show_flag_after_subcommand():
    args = _build_parser().parse_args(["demo", "--show"])
    assert args.show is True


def test_fit_specific_options_still_parse():
    args = _build_parser().parse_args(
        ["fit", "data.csv", "--weight", "abs", "-v", "-o", "y"]
    )
    assert args.csv == "data.csv" and args.weight == "abs"
    assert args.verbose is True and args.out == "y"


def test_missing_subcommand_errors():
    with pytest.raises(SystemExit):
        _build_parser().parse_args(["--out", "x"])
