"""Test :mod:`~triumvirate.winconv`.

"""
from fractions import Fraction

import numpy as np
import pytest
from sympy.core.numbers import Number

from triumvirate.winconv import (
    Multipole,
    calc_threept_winconv_coeff,
    WinConvFormulae,
    # ThreePointWindow,
    # calc_threept_ic,
    TwoPointWinConvBase,
    # TwoPCFWinConv,
    # PowspecWinConv,
    ThreePointWinConvBase,
    # ThreePCFWinConv,
    # BispecWinConv,
)


@pytest.fixture
def multipole():
    return Multipole('110')


def test_multipole___init__(multipole):
    assert multipole._multipole == (1, 1, 0)
    assert multipole.multipole == (1, 1, 0)


def test_multipole___str__(multipole):
    assert str(multipole) == '110'


def test_multipole___neq__(multipole):

    multipole_other = Multipole((1, 1, 0))
    assert multipole_other == multipole

    multipole_another = Multipole('0')
    assert multipole_another != multipole


@pytest.mark.parametrize(
    "multipole, multipole_Q, multipole_Z, expected_coeff",
    [
        (Multipole('000'), Multipole('110'), Multipole('110'), 1./3.),
        (Multipole('110'), Multipole('000'), Multipole('110'), Number('1')),
        (Multipole('110'), Multipole('000'), Multipole('000'), 0.),
    ]
)
def test_calc_threept_winconv_coeff(multipole, multipole_Q, multipole_Z,
                                    expected_coeff):

    coeff = calc_threept_winconv_coeff(
        multipole, multipole_Q, multipole_Z,
        symb=isinstance(expected_coeff, Number)
    )

    if isinstance(expected_coeff, Number):
        validity = (coeff == expected_coeff)
    else:
        validity = np.isclose(coeff, expected_coeff)
    assert validity


@pytest.fixture
def winconvformulae():
    formulae_dict = {
        '000': [
            ('000', '000', 1),
            ('000', 'ic', -1),
            ('110', '110', Fraction(1, 3)),
        ],
    }
    return WinConvFormulae(formulae_dict)


# Implicit tests of :class:`~triumvirate.winconv.WinConvTerm`
def test_winconvformulae_getitem(winconvformulae):
    term = winconvformulae['000'][0]
    assert term.ind_Q == '000'
    assert term.ind_Z == '000'
    assert term.coeff == 1


# Implicit tests of all methods in one
@pytest.mark.parametrize(
    "formulae_dict,"
    "expected_multipoles, expected_multipoles_Q, expected_multipoles_Z,"
    "expected_latex_expr",
    [
        (
            {
                '000': [
                    ('000', '000', 1),
                    ('000', 'ic', -1),
                    ('110', '110', Fraction(1, 3)),
                ],
            },
            ['000',],
            ['000', '110',],
            ['000', 'ic', '110',],
            r"Q_\mathrm{000} \zeta_\mathrm{000} "
            r"- Q_\mathrm{000} \zeta_\mathrm{ic} "
            r"+ \frac{1}{3} Q_\mathrm{110} \zeta_\mathrm{110}"
        ),
    ]
)
def test_winconvformulae(formulae_dict, expected_multipoles,
                         expected_multipoles_Q, expected_multipoles_Z,
                         expected_latex_expr):

    formulae = WinConvFormulae(formulae_dict)

    assert set(formulae.multipoles) == set(expected_multipoles)
    assert set(formulae.multipoles_Q) == set(expected_multipoles_Q)
    assert set(formulae.multipoles_Z) == set(expected_multipoles_Z)

    latex_expr = formulae.get_latex_expr('000')
    assert latex_expr == expected_latex_expr


@pytest.mark.parametrize(
    "formulae, window_sampts, window_multipoles, exception",
    [
        (
            "wilson+16", np.array([0.1, 0.2, 0.3]),
            {
                '0': np.array([1.0, 2.0, 3.0]),
                '2': np.array([0.5, 0.6, 0.7]),
            },
            ValueError,
        ),
        (
            "wilson+16", np.array([0.1, 0.2, 0.3]),
            {
                '0': np.array([1.0, 2.0, 3.0]),
                '2': np.array([0.5, 0.6, 0.7]),
                '4': np.array([0.5, 0.6, 0.7]),
                '6': np.array([0.5, 0.6, 0.7]),
                '8': np.array([0.5, 0.6, 0.7]),
            },
            ValueError,
        ),
    ]
)
def test_twopointwinconvbase(formulae, window_sampts, window_multipoles,
                             exception):
    # Create instance.
    if exception is not None:
        with pytest.raises(exception):
            winconv = ThreePointWinConvBase(
                formulae, window_sampts, window_multipoles
            )
        return

    winconv = TwoPointWinConvBase(formulae, window_sampts, window_multipoles)

    # Check formulae attribute.
    assert isinstance(winconv._formulae, WinConvFormulae)

    # Check representative multipole.
    assert isinstance(winconv._multipole_rep, (str, int))

    # Check window multipoles.
    assert np.allclose(winconv._rQ_in, window_sampts)
    assert np.allclose(winconv._Q_in['0'], window_multipoles['0'])
    assert np.allclose(winconv._Q_in['2'], window_multipoles['2'])


@pytest.mark.parametrize(
    "formulae, window_sampts, window_multipoles, exception",
    [
        (
            "sugiyama+18", np.array([0.1, 0.2, 0.3]),
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
            None,
        ),
        (
            "sugiyama+18", np.array([0.1, 0.2, 0.3]),
            {
                '000': np.array(
                    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
                ),
                '202': np.array(
                    [[0.5, 0.6, 0.7], [0.8, 0.9, 1.0], [1.1, 1.2, 1.3]]
                ),
            },
            ValueError,
        ),
    ]
)
def test_threepointwinconvbase(formulae, window_sampts, window_multipoles,
                               exception):
    # Create instance.
    if exception is not None:
        with pytest.raises(exception):
            winconv = ThreePointWinConvBase(
                formulae, window_sampts, window_multipoles
            )
        return

    winconv = ThreePointWinConvBase(formulae, window_sampts, window_multipoles)

    # Check formulae attribute.
    assert isinstance(winconv._formulae, WinConvFormulae)

    # Check representative multipole.
    assert isinstance(winconv._multipole_rep, (str, int))

    # Check window multipoles.
    assert np.allclose(winconv._rQ_in, window_sampts)
    for multipole in window_multipoles:
        assert np.allclose(
            winconv._Q_in[multipole], window_multipoles[multipole]
        )
