"""Test :mod:`~triumvirate.winconv`.

"""
import doctest
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from sympy import Expr, sympify

import triumvirate.winconv as winconv
from triumvirate.winconv import (
    # BispecWinConv,
    # PowspecWinConv,
    # ThreePCFWinConv,
    # ThreePointWindow,
    # TwoPCFWinConv,
    WinConvBase,
    # calc_threept_ic,
    calc_threept_winconv_coeff,
)


def test_doctest():
    # Doctest is included for `Multipole`, `WinConvTerm` and `WinConvFormulae`.
    assert doctest.testmod(winconv, raise_on_error=True).failed == 0


@pytest.mark.parametrize(
    "multipole, multipole_Q, multipole_Z, expected_coeff",
    [
        ('000', '110', '110', sympify('1/3')),
        ('110', '000', '110', 1),
        ('110', '000', '000', 0.),
    ]
)
def test_calc_threept_winconv_coeff(multipole, multipole_Q, multipole_Z,
                                    expected_coeff):

    coeff = calc_threept_winconv_coeff(
        multipole, multipole_Q, multipole_Z,
        symb=isinstance(expected_coeff, Expr)
    )

    if isinstance(expected_coeff, Expr):
        validity = (coeff == expected_coeff)
    else:
        validity = np.isclose(coeff, expected_coeff)
    assert validity


@pytest.mark.parametrize(
    "formulae, window_sampts, window_multipoles, exception_context",
    [
        (
            "wilson+16",
            np.array([0.1, 0.2, 0.3]),
            {
                '0': np.array([1.0, 2.0, 3.0]),
                '2': np.array([0.5, 0.6, 0.7]),
            },
            pytest.raises(ValueError),  # missing multipole
        ),
        (
            "wilson+16",
            np.array([0.1, 0.2, 0.3]),
            {
                '0': np.array([1.0, 2.0, 3.0]),
                '2': np.array([0.5, 0.6, 0.7]),
                '4': np.array([0.5, 0.6, 0.7]),
                '6': np.array([0.5, 0.6, 0.7]),
                '8': np.array([0.5, 0.6, 0.7]),
            },
            does_not_raise(),
        ),
        (
            "sugiyama+18",
            np.array([0.1, 0.2, 0.3]),
            {
                '000': np.array(
                    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
                ),
                '202': np.array(
                    [[0.5, 0.6, 0.7], [0.8, 0.9, 1.0], [1.1, 1.2, 1.3]]
                ),
                '110': np.array(
                    [[0.5, 0.6, 0.7], [0.8, 0.9, 1.0], [1.1, 1.2, 1.3]]
                ),
                '112': np.array(
                    [[0.5, 0.6, 0.7], [0.8, 0.9, 1.0], [1.1, 1.2, 1.3]]
                ),
                '022': np.array(
                    [[0.5, 0.6, 0.7], [0.8, 0.9, 1.0], [1.1, 1.2, 1.3]]
                ),
            },
            does_not_raise(),
        ),
        (
            "sugiyama+18",
            np.array([0.1, 0.2, 0.3]),
            {
                '000': np.array(
                    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
                ),
                '202': np.array(
                    [[0.5, 0.6, 0.7], [0.8, 0.9, 1.0], [1.1, 1.2, 1.3]]
                ),
            },
            pytest.raises(ValueError),  # missing multipole
        ),
    ]
)
def test_winconvbase(formulae, window_sampts, window_multipoles,
                     exception_context):

    with exception_context:
        _ = WinConvBase(formulae, window_sampts, window_multipoles)
