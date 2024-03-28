"""
Window Convolution (:mod:`~triumvirate.winconv`)
==========================================================================

.. versionadded:: 0.4.0


Perform window convolution of two- and three-point statistics.

.. autosummary::
    Multipole
    calc_threept_winconv_coeff
    WinConvTerm
    WinConvFormulae
    ThreePointWindow
    calc_threept_ic
    TwoPointWinConvBase
    TwoPCFWinConv
    PowspecWinConv
    ThreePointWinConvBase
    ThreePCFWinConv
    BispecWinConv

"""
import warnings
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from sympy.physics.wigner import wigner_3j, wigner_9j

from triumvirate._arrayops import MixedSignError, SpacingError, _check_1d_array
from triumvirate.transforms import (
    DoubleSphericalBesselTransform,
    SphericalBesselTransform,
    resample_lglin,
    resample_lin,
)


NAMED_FORMULAE = {
    # Wilson et al. (2016) [arXiv:1511.07799]
    'wilson+16': {
        0: [
            (0, 0, 1),
            (2, 2, Fraction(1, 5)),
            (4, 4, Fraction(1, 9)),
        ],
        2: [
            (2, 0, 1),
            (0, 2, 1),
            (2, 2, Fraction(2, 7)),
            (4, 2, Fraction(2, 7)),
            (2, 4, Fraction(2, 7)),
            (4, 4, Fraction(100, 693)),
            (6, 4, Fraction(25, 143))
        ],
        4: [
            (4, 0, 1),
            (2, 2, Fraction(18, 35)),
            (4, 2, Fraction(20, 77)),
            (6, 2, Fraction(45, 143)),
            (0, 4, 1),
            (2, 4, Fraction(20, 77)),
            (4, 4, Fraction(162, 1001)),
            (6, 4, Fraction(20, 143)),
            (8, 4, Fraction(490, 2431)),
        ]
    },
    # Sugiyama et al. (2018) [arXiv:1803.02132]
    'sugiyama+18': {
        '000': [
            ('000', '000', 1),
            ('000', 'ic', -1),
            ('110', '110', Fraction(1, 3)),
        ],
        '110': [
            ('000', '110', 1),
            ('110', '000', 1),
            ('110', 'ic', -1),
        ],
        '202': [
            ('000', '202', 1),
            ('202', '000', 1),
            ('202', 'ic', -1),
            ('110', '112', Fraction(1, 3)),
            ('112', '110', Fraction(1, 3)),
        ],
        '112': [
            ('000', '112', 1),
            ('112', '000', 1),
            ('112', 'ic', -1),
            ('022', '110', Fraction(2, 5)),
            ('202', '110', Fraction(2, 5)),
            ('110', '022', Fraction(2, 5)),
            ('110', '202', Fraction(2, 5)),
        ],
    },
}
"""Named formulae for window convolution.

This is a `dict` where each key corresponds to a named formula below:

- ``'wilson+16'`` (for plane-parallel two-point statistics):
  Wilson et al., 2016. MNRAS 464(3), 3121 [arXiv:1511.07799];
- ``'sugiyama+18'`` (for plane-parallel three-point statistics):
  Sugiyama et al., 2018. MNRAS 484(1), 364 [arXiv:1803.02132].

The value for each key is itself a `dict`; see
:class:`~triumvirate.winconv.WinConvFormulae` for more details.

"""


class Multipole:
    """Representation of a multipole.

    Parameters
    ----------
    multipole : str or int or tuple of int of length 3
        Multipole.

    """

    def __init__(self, multipole):

        if isinstance(multipole, str):
            if len(multipole) == 1:
                self._multipole = int(multipole)
            elif len(multipole) == 3:
                self._multipole = tuple(map(int, multipole))
            else:
                raise ValueError("Multipole string must be of length 1 or 3.")
        elif isinstance(multipole, (int, tuple)):
            self._multipole = multipole
        else:
            raise ValueError("Unknown multipole type.")

    def __str__(self):
        if isinstance(self._multipole, int):
            return str(self._multipole)
        return ''.join(map(str, self._multipole))

    def __eq__(self, other):
        return self._multipole == other._multipole

    def __ne__(self, other):
        return not self.__eq__(other)

    @property
    def multipole(self):
        """Identifiable Multipole representation.

        """
        return self._multipole


def _N_prefactor(multipole=None, ell1=None, ell2=None, L=None):
    """Calculate the Wigner-3j--related pre-factor.

    Parameters
    ----------
    multipole : :class:`~winconv.Multipole`, optional
        Multipole.  If `None` (default), use the individual
        multipole degrees.
    ell1, ell2, L : int, optional
        Multipole degrees (default is `None`).  Only used if
        `multipole` is `None`.

    Returns
    -------
    int
        Multiplicative factor.

    """
    if multipole is not None:
        ell1, ell2, L = multipole.multipole

    return (2 * ell1 + 1) * (2 * ell2 + 1) * (2 * L + 1)


def _H_factor(ell1, ell2, ell3):
    """Calculate the Wigner-3j factor with zero spherical orders.

    Parameters
    ----------
    ell1, ell2, ell3 : int
        Multipole degrees.
    symb : bool, optional
        If `True` (default is `False`), return the coefficient as a
        symbolic expression.

    Returns
    -------
    :class:`sympy.core.numbers.Number`
        Wigner-3j factor.

    """
    return wigner_3j(ell1, ell2, ell3, 0, 0, 0)


def calc_threept_winconv_coeff(multipole, multipole_Q, multipole_Z,
                               symb=False):
    """Calculate the window convolution coefficient for a specific
    window convolution term.

    Parameters
    ----------
    multipole : str or tuple of int or :class:`~winconv.Multipole`
        Convolved multipole.
    multipole_Q : str or tuple of int or :class:`~winconv.Multipole`
        Window function multipole.
    multipole_Z : str or tuple of int or :class:`~winconv.Multipole`
        Unwindowed correlation function multipole.
    symb : bool, optional
        If `True` (default is `False`), return the coefficient as a
        symbolic expression.

    Returns
    -------
    float or :class:`sympy.core.numbers.Number`
        Window convolution coefficient.

    """
    if not isinstance(multipole, Multipole):
        multipole = Multipole(multipole)
    if not isinstance(multipole_Q, Multipole):
        multipole_Q = Multipole(multipole_Q)
    if not isinstance(multipole_Z, Multipole):
        multipole_Z = Multipole(multipole_Z)

    ell1, ell2, L = multipole.multipole
    ell1_, ell2_, L_ = multipole_Z.multipole
    ell1__, ell2__, L__ = multipole_Q.multipole

    H_frac = (
        _H_factor(ell1, ell2, L)
        / _H_factor(ell1_, ell2_, L_)
        / _H_factor(ell1__, ell2__, L__)
        * _H_factor(ell1, ell1_, ell1__)
        * _H_factor(ell2, ell2_, ell2__)
        * _H_factor(L, L_, L__)
    )

    coeff = (
        _N_prefactor(ell1=ell1, ell2=ell2, L=L)
        * wigner_9j(ell1__, ell2__, L__, ell1_, ell2_, L_, ell1, ell2, L)
        * H_frac
    )

    if symb:
        return coeff
    return float(coeff)


class WinConvTerm(namedtuple('WinConvTerm', ['ind_Q', 'ind_Z', 'coeff'])):
    r"""Window convolution term as a tuple of three factors.

    Attributes
    ----------
    ind_Q : str or int or tuple of int
        Window function multipole index/indices.
    ind_Z : str or int or tuple of int
        Unwindowed correlation funciton multipole index/indices.
    coeff : float
        Numerical coefficient for the convolution term.

    Examples
    --------
    The following `term` corresponds to the term
    :math:`- 1 \cdot Q_{000} \zeta_{\mathrm{ic}}`:

    >>> term = WinConvTerm('000', 'ic', -1)
    >>> print(term.coeff, term.ind_Q, term.ind_Z)
    -1 000 ic

    Alternative, it can be constructed as:

    >>> term = WinConvTerm((0, 0, 0), 'ic', -1)

    """

    def __init__(self, ind_Q, ind_Z, coeff):
        if coeff == 0:
            warnings.warn(
                "Coefficient of window convolution term is zero.",
                RuntimeWarning
            )


class WinConvFormulae:
    r"""Window convolution formulae.

    The full set of formulae are encoded as a dictionary, where each key
    corresponds to a windowed correlation function (CF) multipole, and
    each value is a sequence of (named) tuples each corresponding to
    a window convolution term.

    Parameters
    ----------
    formulae : |formulae_dicttype|
        Window convolution formulae.

    Attributes
    ----------
    formulae : |formulae_type|
        Window convolution formulae.
    multipoles : list of str or list of (tuple of) int
        Windowed CF multipoles.
    multipoles_Q : list of str or list of (tuple of) int
        Required window function multipoles.
    multipoles_Z : list of str or list of (tuple of) int
        Required unwindowed CF multipoles.

    Examples
    --------
    The following `formulae` corresponds to a single formula for the
    three-point CF monopole (0, 0, 0) consisting of three terms, namely
    :math:`Q_{000} \zeta_{000} - Q_{000} \zeta_{\mathrm{ic}} + \frac{1}{3}
    Q_{110} \zeta_{110}`:

    >>> from fractions import Fraction
    >>> formulae_dict = {
    ...     '000': [
    ...         ('000', '000', 1),
    ...         ('000', 'ic', -1),
    ...         ('110', '110', Fraction(1, 3)),
    ...     ],
    ... }
    >>> formulae = WinConvFormulae(formulae_dict)

    It encodes the formula for the windowed (0, 0, 0) monopole only:

    >>> print(formulae.multipoles)
    ['000']

    The window function multipoles required are:

    >>> print(formulae.multipoles_Q)
    ['000', '110']

    The unwindowed three-point CF multipoles required are:

    >>> print(formulae.multipoles_Z)
    ['000', 'ic', '110']

    The formular can be printed as a LaTeX string (without the maths-mode
    delimiters):

    >>> print(formulae.get_latex_expr('000'))
    Q_\\mathrm{000} \\zeta_\\mathrm{000} - Q_\\mathrm{000} \\zeta_\\mathrm{ic} + \\frac{1}{3} Q_\\mathrm{110} \\zeta_\\mathrm{110}


    .. |formulae_dicttype| replace:: \
        dict of \
        {|pole_type|: sequence of tuples of \
        (|pole_type|, |pole_type|, |coeff_type|)}

    .. |formulae_type| replace:: \
        dict of \
        {|pole_type|: list of :class:`~triumvirate.winconv.WinConvTerm`}

    .. |pole_type| replace:: str or int or tuple of int

    .. |coeff_type| replace:: int, float or :class:`fractions.Fraction`

    """  # noqa: E501

    def __init__(self, formulae):

        self.formulae = {
            multipole: [WinConvTerm(*term) for term in formula]
            for multipole, formula in formulae.items()
        }

        self.multipoles = list(formulae.keys())
        self.multipoles_Q = list(set([
            term.ind_Q
            for formula in self.formulae.values()
            for term in formula
        ]))
        self.multipoles_Z = list(set([
            term.ind_Z
            for formula in self.formulae.values()
            for term in formula
        ]))

    def __getitem__(self, multipole):
        """Get the formula for a specific windowed CF multipole.

        Parameters
        ----------
        multipole : str or (tuple of) int
            Windowed CF multipole.

        Returns
        -------
        list of :class:`~triumvirate.winconv.WinConvTerm`
            List of window convolution terms for `multipole`.

        """
        return self.formulae[multipole]

    def get_latex_expr(self, multipole, symbol=r"\zeta"):
        """Get the LaTeX string representing the window convolution
        formula for a specific windowed CF multipole.

        Parameters
        ----------
        multipole : str or (tuple of) int
            Windowed CF multipole.
        symbol : str, optional
            Symbol for the unwindowed CF (default is ``r"\\zeta"``).

        Returns
        -------
        str
            LaTeX string representing the window convolution formula
            (without the maths-mode delimiters).

        """
        terms_str = []
        for idx, term in enumerate(self[multipole]):
            if term.coeff == 0:
                continue
            elif term.coeff < 0:
                coeff_sgn = '-'
            elif idx > 0:
                coeff_sgn = '+'
            else:
                coeff_sgn = ''

            if abs(term.coeff) == 1:
                coeff_num = ''
            elif isinstance(term.coeff, Fraction):
                coeff_num = r"\frac{{{}}}{{{}}}".format(
                    abs(term.coeff.numerator), abs(term.coeff.denominator)
                )
            else:
                coeff_num = f"{abs(term.coeff)}"

            coeff_str = (coeff_sgn + ' ' + coeff_num).strip()

            if isinstance(term.ind_Q, str):
                Qind_str = term.ind_Q
            elif isinstance(term.ind_Q, int):
                Qind_str = str(term.ind_Q)
            else:
                Qind_str = ''.join(map(str, term.ind_Q))

            if isinstance(term.ind_Z, str):
                Zind_str = term.ind_Z
            elif isinstance(term.ind_Z, int):
                Zind_str = str(term.ind_Z)
            else:
                Zind_str = ''.join(map(str, term.ind_Z))

            term_str_ = "{} Q_\\mathrm{{{}}} {}_\\mathrm{{{}}}".format(
                coeff_str, Qind_str, symbol, Zind_str
            ).strip()

            terms_str.append(term_str_)

        return ' '.join(terms_str)


class ThreePointWindow:
    """Three-point window function in configuration space.

    By construction, the window function is discretely sampled as a
    symmetric matrix where the upper triangular part is stored as a
    flattened array in row-major order, where each value corresponds
    to a pair of separation sample points.

    Parameters
    ----------
    paired_sampts : 2-d array of float
        Paired separation window function sample points.
    flat_multipoles : dict of {str or tuple of int: 1-d array of float}
        Flattened window function multipole samples (each key is
        a multipole).
    alpha_contrast : float, optional
        Alpha contrast factor (default is 1.) for normalising the higher
        number density of the random catalogue, e.g. 0.1 for a
        random catalogue with 10 times the number density of the
        data catalogue.  The square of this factor is applied to the
        amplitude of the window function.
    make_loglin : bool or int, optional
        Whether to resample the window function multipole samples
        logarithmically if they are not already (default is `True`).
        If an integer, it is the number of points to resample to.

    Attributes
    ----------
    multipoles : list of str or list of tuple of int
        Window function multipoles.
    r : 1-d array of float
        window function separation sample points.
    Q : dict of {str or tuple of int: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    Q_diag : dict of {str or tuple of int: 1-d array of float}
        Diagonal of the window function multipole samples (each key
        is a multipole).
    sources : list of str
        Sources of the window function multipole samples as either 'array'
        or the file path.

    Raises
    ------
    ValueError
        If `paired_sampts` is not a 2-d array.
    ValueError
        If `paired_sampts` and `flat_multipoles` have different lengths.
    ValueError
        If `paired_sampts` and `flat_multipoles` do not correspond to an
        upper triangular matrix.

    """

    def __init__(self, paired_sampts, flat_multipoles, alpha_contrast=1.,
                 make_loglin=True):

        # Check dimensions.
        if paired_sampts.ndim != 2:
            raise ValueError(
                "Paired separation sample points must be a 2-d array."
            )
        for flat_multipole_ in flat_multipoles.values():
            if len(paired_sampts) != len(flat_multipole_):
                raise ValueError(
                    "Paired separation sample points and flattened multipole "
                    "samples must have the same length."
                )

        nbins = (np.sqrt(8*len(paired_sampts) + 1) - 1) / 2.
        if not np.isclose(nbins, int(nbins)):
            raise ValueError(
                "Paired separation sample points do not correspond to "
                "an upper triangular matrix."
            )
        nbins = int(nbins)

        # Reshape arrays.
        triu_indices = np.triu_indices(nbins)

        rQ_in = paired_sampts[:nbins, 1]

        Q_in = {}
        for multipole, Q_flat in sorted(flat_multipoles.items()):
            Q_mat = np.zeros((nbins, nbins))
            Q_mat[triu_indices] = Q_flat
            Q_mat.T[triu_indices] = Q_flat
            Q_in[multipole] = alpha_contrast**2 * Q_mat

        self.r = rQ_in
        self.Q = Q_in

        self.Q_diag = {}
        self.multipoles = []
        for multipole, Q_in_ in self.Q.items():
            self.Q_diag[multipole] = np.diag(Q_in_)
            self.multipoles.append(multipole)

        self.sources = ['array'] * len(self.multipoles)

        # Check sample spacing.
        try:
            _check_1d_array(self.r, check_loglin=True)
        except (MixedSignError, SpacingError,):
            if make_loglin:
                self.resample(
                    size=make_loglin if type(make_loglin) is int else None
                )
            else:
                warnings.warn(
                    "Window function separation sample points are not "
                    "logarithmically spaced.",
                    RuntimeWarning
                )

    @classmethod
    def load_from_textfiles(cls, filepaths, subtract_shotnoise=True,
                            usecols=(1, 4, 6, 8), alpha_contrast=1.,
                            make_loglin=True):
        """Load window function from plain-text files as saved by
        :func:`~triumvirate.threept.compute_3pcf_window` with
        ``form='full'``.

        Parameters
        ----------
        filepaths : dict of {str or tuple of int: str}
            Window function file paths (values) for each multipole (key).
        subtract_shotnoise : bool, optional
            Whether to subtract the shot-noise contribution from the
            window function (default is `True`).
        usecols : tuple of int, optional
            Zero-indexed columns to read from the file (default is
            ``(1, 4, 6, 8)``).   If `subtract_shotnoise` is `True`,
            must be of length 4 and correspond to the effective separation
            sample points (two columns), and the raw window function and
            shot-noise (two columns); otherwise, must be of length 3 and
            correspond to the effective separation sample points (two
            columns) and the raw window function (one column).
        alpha_contrast : float, optional
            Alpha contrast factor (default is 1.) for normalising the
            higher number density of the random catalogue, e.g. 0.1 for a
            random catalogue with 10 times the number density of the
            data catalogue.
        make_loglin : bool, optional
            Whether to resample the window function multipole samples
            logarithmically if they are not already (default is `True`).

        Returns
        -------
        :class:`~triumvirate.winconv.ThreePointWindow`
            Three-point window function.


        .. attention::

            The effective separation sample points are assumed to be
            the same for all window function files.

        """
        r1_eff, r2_eff = None, None

        flat_multipoles = {}
        sources = []
        for multipole, filepath in sorted(filepaths.items()):
            if subtract_shotnoise:
                if len(usecols) != 4:
                    raise ValueError(
                        "Must provide four columns to read from the file "
                        "when subtracting shot noise."
                    )
                r1_eff_, r2_eff_, Qflat_raw, Qflat_shot = np.loadtxt(
                    filepath, unpack=True, usecols=usecols
                )
                flat_multipoles[multipole] = Qflat_raw - Qflat_shot
            else:
                if len(usecols) != 3:
                    raise ValueError(
                        "Must provide three columns to read from the file "
                        "when not subtracting shot-noise."
                    )
                r1_eff_, r2_eff_, Qflat_raw = np.loadtxt(
                    filepath, unpack=True, usecols=usecols
                )
                flat_multipoles[multipole] = Qflat_raw

            if r1_eff is not None:
                if not np.allclose(r1_eff_, r1_eff):
                    warnings.warn(
                        "Effective separation sample points "
                        "for the first dimension "
                        "are not the same for all window function files."
                    )
            else:
                r1_eff = r1_eff_
            if r2_eff is not None:
                if not np.allclose(r2_eff_, r2_eff):
                    warnings.warn(
                        "Effective separation sample points "
                        "for the second dimension "
                        "are not the same for all window function files."
                    )
            else:
                r2_eff = r2_eff_

            sources.append(filepath)

        paired_sampts = np.column_stack((r1_eff, r2_eff))

        self = cls(
            paired_sampts, flat_multipoles,
            alpha_contrast=alpha_contrast, make_loglin=make_loglin
        )

        self.sources = sources

        return self

    @classmethod
    def load_from_dict(cls, dict_obj):
        """Load window function from a dictionary.

        Parameters
        ----------
        dict_obj : dict
            Dictionary of 2-d window function multipoles
            (key 'Q') and sample points (key 'r'), where the multipoles
            are themselves a single dictionary with each key corresponding
            to a multipole.

        Returns
        -------
        :class:`~triumvirate.winconv.ThreePointWindow`
            Three-point window function.

        """
        self = cls.__new__(cls)

        self.r = dict_obj['r']
        self.Q = dict_obj['Q']

        self.Q_diag = {
            multipole: np.diag(Q_in_)
            for multipole, Q_in_ in self.Q.items()
        }

        self.multipoles = sorted(list(self.Q.keys()))

        try:
            self.sources = dict_obj['sources']
        except KeyError:
            self.sources = ['array'] * len(self.multipoles)

        return self

    @classmethod
    def load_from_file(cls, filepath):
        """Load window function multipoles from a compressed NumPy file.

        Parameters
        ----------
        filepath : str
            File path to load the window function multipoles.

        Returns
        -------
        :class:`~triumvirate.winconv.ThreePointWindow`
            Three-point window function.

        """
        self = cls.__new__(cls)

        win_obj = np.load(filepath)

        self.r = win_obj['r']
        self.Q = win_obj['Q']
        self.Q_diag = win_obj['Q_diag']

        self.multipoles = win_obj['multipoles']
        self.sources = win_obj['sources']

        return self

    def save_to_file(self, filepath):
        """Save window function to a compressed NumPy file.

        Parameters
        ----------
        filepath : str
            File path to save the window function multipoles.

        """
        win_obj = {
            'r': self.r,
            'Q': self.Q,
            'Q_diag': self.Q_diag,
            'multipoles': self.multipoles,
            'sources': self.sources,
        }

        np.savez(filepath, **win_obj)

    def subselect(self, rrange=None, idxrange=None, inplace=False):
        """Select a range of separation sample points.

        Parameters
        ----------
        rrange : tuple of float or None, optional
            Range of separation sample points (default is `None`), e.g.
            ``(None, 300.)``.
        idxrange : tuple of int or None, optional
            Range of separation sample point indices (default is `None`),
            e.g. ``(None, -1)``, equivalent to ``slice(None, -1)``.
        inplace : bool, optional
            If `False` (default), return a new window function object;
            otherwise, modify the window function in place.

        Returns
        -------
        :class:`~triumvirate.winconv.ThreePointWindow` or None
            If returned, three-point window function for the selected
            range of separation sample points.

        """
        if rrange is None:
            sel_range = np.full_like(self.r, True, dtype=bool)
        else:
            rmin, rmax = rrange
            rmin = rmin or 0.
            rmax = rmax or np.inf
            sel_range = np.logical_and(rmin <= self.r, self.r <= rmax)

        if idxrange is None:
            sel_idx = np.full_like(self.r, True, dtype=bool)
        else:
            sel_idx = np.full_like(self.r, False, dtype=bool)
            sel_idx[slice(*idxrange)] = True

        selector = sel_range & sel_idx

        if not inplace:
            return self.__class__.load_from_dict({
                'r': self.r[selector],
                'Q': {
                    multipole: Q_in_[selector, :][:, selector]
                    for multipole, Q_in_ in self.Q.items()
                }
            })

        self.r = self.r[selector]
        self.Q = {
            multipole: Q_in_[selector, :][:, selector]
            for multipole, Q_in_ in self.Q.items()
        }
        self.Q_diag = {
            multipole: np.diag(Q_in_)
            for multipole, Q_in_ in self.Q.items()
        }

    def resample(self, spacing='lglin', size=None):
        """Resample window function multipole samples.

        Parameters
        ----------
        spacing : {'lglin', 'lin'}, optional
            Spacing type for the resampled separation sample points
            either 'lglin' (default) (logarithmic-linear) or 'lin'
            (linear).
        size : (tuple of) int, optional
            Number(s) of points to resample to (default is `None`).

        """
        _rQ, _Q = None, {}
        if spacing == 'lglin':
            for multipole, Q_in_ in self.Q.items():
                _rQ, _Q[multipole] = resample_lglin(self.r, Q_in_, size=size)
        elif spacing == 'lin':
            for multipole, Q_in_ in self.Q.items():
                _rQ, _Q[multipole] = resample_lin(self.r, Q_in_, size=size)
        else:
            raise ValueError("Unknown spacing type: {}".format(spacing))

        self.r, self.Q = _rQ[0], _Q

        self.Q_diag = {
            multipole: np.diag(Q_in_)
            for multipole, Q_in_ in self.Q.items()
        }


def _check_conv_range(r_conv, r):
    """Check if the convolution range is fully covered by the
    sample points.

    Parameters
    ----------
    r_conv : 1-d array of float
        Convolution range sample points.
    r : 1-d array of float
        Separation sample points.

    Returns
    -------
    bool
        `True` if within range; `False` otherwise.

    """
    RTOL = 1.e-5

    def _check(rc, rs):
        return (1 - RTOL) * rs.min() <= rc.min() \
            and (rc.max() <= (1 + RTOL) * rs.max())

    if isinstance(r_conv, dict):
        if isinstance(r, dict):
            multipoles_common = set(r_conv.keys()) & set(r.keys())
            return all(
                _check(r_conv[multipole], r[multipole])
                for multipole in multipoles_common
            )
        else:  # STYLE: unpythonic but logically clearer
            return all(_check(r_conv[multipole], r) for multipole in r_conv)
    else:  # STYLE: unpythonic but logically clearer
        return _check(r_conv, r)


def _integrate_2d_samples(x, y, z):
    r"""Integrate 2-d samples of function :math:`z(x, y)`,

    .. math:: \int \int z(x, y) x^2 y^2 \, \mathrm{d}x \mathrm{d}y \,,

    using Simpson's rule.

    Parameters
    ----------
    x, y : 1-d array of float
        Sample points for both dimensions.
    z : 2-d array of float
        Sample values.

    Returns
    -------
    float
        Integrated value.

    """
    return simpson(simpson(z * y * y, y, axis=-1) * x * x, x, axis=-1)


def calc_threept_ic(window_sampts, window_multipoles, r_in, zeta_in,
                    r_common=None, approx=False):
    """Compute three-point clustering statistic integral constraint.

    Parameters
    ----------
    window_sampts :  1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str or tuple of int: 2-d array of float}
        Window function multipole samples (each key is a multipole).
        Must contain the ``'000'`` or ``(0, 0, 0)`` multipole.
    r_in : dict of {str or tuple of int: 1-d array of float}
        Separation sample points for the input 3PCF multipoles
        (each key is a multipole).
    zeta_in : dict of {str or tuple of int: 2-d array of float}
        Input 3PCF multipole samples (each key is a multipole)
        at sample points `r_in`.
    r_common : 1-d array of float, optional
        Common separation sample points.  If `None` (default), it is the
        same as `r_in`.
    approx : bool, optional
        If `True` (default is `False`), include only the leading-order
        term with ``'000'`` or ``(0, 0, 0)`` multipole as
        an approximation.

    Returns
    -------
    float
        Integral constraint.


    .. attention::

        If the window convolution formulae use `str`, `int` or
        `tuple` of `int` to represent a multipole, the corresponding
        `dict` of window function multipole samples must use exactly
        the same key type.  In other words, the type representation of
        a multipole must be consistent across all inputs.

    """
    # Infer multipole key type.
    keytype = type(next(iter(window_multipoles.keys())))

    monopole = '000' if keytype is str else (0, 0, 0)
    if approx:
        multipoles_common = [monopole]
    else:
        multipoles_common = set(window_multipoles.keys()) & set(zeta_in.keys())
        if monopole not in multipoles_common:
            raise ValueError(
                "The '000' or (0, 0, 0) multipole is required for calculating "
                "the integral constraint."
            )

    if r_common is None:
        r_common = r_in
    elif not isinstance(r_common, dict):
        r_common = {multipole: r_common for multipole in multipoles_common}

    _Q_in = {
        _multipole: RectBivariateSpline(
            window_sampts, window_sampts, window_multipoles[_multipole]
        )(r_common[_multipole], r_common[_multipole])
        for _multipole in multipoles_common
    }
    _zeta_in = {
        _multipole: RectBivariateSpline(
            r_in[_multipole], r_in[_multipole], zeta_in[_multipole].real
        )(r_common[_multipole], r_common[_multipole])
        for _multipole in multipoles_common
    }

    ic_nom = 0.
    for _multipole in multipoles_common:
        ell1, ell2, ELL = map(int, _multipole)

        # `float` conversion necessary as `sympy` returns `sympy.Float`
        # and slows down the computation.
        _H2 = float(_H_factor(ell1, ell2, ELL)) ** 2
        _N = _N_prefactor(ell1=ell1, ell2=ell2, L=ELL)

        ic_nom += _N * _H2 * _integrate_2d_samples(
            r_common[_multipole], r_common[_multipole],
            _Q_in[_multipole] * _zeta_in[_multipole]
        )

    ic_denom = _integrate_2d_samples(
        r_common[_multipole], r_common[_multipole], _Q_in[monopole]
    )

    return ic_nom / ic_denom


class TwoPointWinConvBase:
    """Generic window convolution of two-point statistics.

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'wilson+16').
    window_sampts : 1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str or int: 1-d array of float}
        Window function multipole samples (each key is a multipole).


    .. attention::

        If the window convolution formulae use `str` or `int` to represent
        a multipole, the corresponding `dict` of window function multipole
        samples must use exactly the same key type.  In other words, the
        type representation of a multipole must be consistent across
        all inputs.

    """

    def __init__(self, formulae, window_sampts, window_multipoles):

        # Constuct formulae object.
        if isinstance(formulae, str):
            try:
                formulae = NAMED_FORMULAE[formulae.lower()]
            except KeyError:
                raise ValueError(f'Unknown named formulae: {formulae}.')
        if not isinstance(formulae, WinConvFormulae):
            formulae = WinConvFormulae(formulae)
        self._formulae = formulae

        if formulae.multipoles_Z[0] != 'ic':
            self._multipole_rep = formulae.multipoles_Z[0]
        else:
            self._multipole_rep = formulae.multipoles_Z[-1]

        # Store window multipoles.
        self._rQ_in = window_sampts
        self._Q_in = window_multipoles

        if not set(self._formulae.multipoles_Q).issubset(self._Q_in.keys()):
            raise ValueError(
                "The window function multipole samples must contain all "
                "the multipole keys required by the convolution formulae."
            )


class TwoPCFWinConv(TwoPointWinConvBase):
    """Window convolution of two-point correlation function (2PCF).

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'wilson+16').
    window_sampts : 1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str or int: 1-d array of float}
        Window function multipole samples (each key is a multipole).
    r_in : 1-d array of float or dict of {str or int: 1-d array-like}
        Separation sample points for the input 2PCF multipoles.
    r_out : 1-d array of float, optional
        Separation sample points for the output 2PCF multipoles.
        If `None` (default), it is set to `r_in`.

    Attributes
    ----------
    r_in : dict of {str or int: 1-d array of float}
        Separation sample points for the input 2PCF multipoles.
    r_out : dict of {str or int: 1-d array of float}
        Separation sample points for the output 2PCF multipoles.


    .. attention::

        Integral constraint is ignored.

    .. attention::

        If the window convolution formulae use `str` or `int` to represent
        a multipole, the corresponding `dict` of window function multipole
        samples and any 2PCF multipole samples to be convolved must use
        exactly the same key type.  In other words, the type
        representation of a multipole must be consistent across
        all inputs.

    .. attention::

        The window function multipoles are assumed to be sampled at
        the same separation sample points.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 r_in, r_out=None):

        super().__init__(formulae, window_sampts, window_multipoles)

        if isinstance(r_in, dict):
            self.r_in = r_in
        else:
            self.r_in = {
                multipole: r_in
                for multipole in set(
                    formulae.multipoles_Z + formulae.multipoles_Q
                )
            }
        if r_out is None:
            self.r_out = r_in
        else:
            self.r_out = {
                multipole: r_out
                for multipole in set(
                    formulae.multipoles_Z +
                    formulae.multipoles_Q +
                    formulae.multipoles
                )
            }
        if not _check_conv_range(self.r_out, self._rQ_in):
            warnings.warn(
                "The convolution range `r_out` is not fully covered "
                "by the window function separation sample points. "
                "The window function samples may be extrapolated."
            )
        if not _check_conv_range(self.r_out, self.r_in):
            warnings.warn(
                "The convolution range `r_out` is not fully covered "
                "by the input 3PCF separation sample points. "
                "The input 3PCF samples may be extrapolated."
            )

        self._Q_out = {
            _multipole: InterpolatedUnivariateSpline(
                self._rQ_in, _Qpole_in
            )(self.r_out[_multipole])
            for _multipole, _Qpole_in in self._Q_in.items()
        }

    def convolve(self, xi_in):
        """Convolve 2PCF multipoles.

        Parameters
        ----------
        xi_in : dict of {str or int: 1-d array of float}
            Input 3PCF multipole samples (each key is a multipole)
            at sample points :attr:`r_in`.

        Returns
        -------
        xi_conv_out : dict of {str or int: 1-d :class:`numpy.ndarray`}
            Output windowed 3PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_out`.

        """
        # Interpolate input multipoles at output sample points.
        xi_out = {
            multipole: InterpolatedUnivariateSpline(
                self.r_in[multipole], Zpole_in
            )(self.r_out[multipole])
            for multipole, Zpole_in in xi_in.items()
        }

        # Perform convolution as multiplication.
        xi_conv_out = {}
        for multipole in self._formulae.multipoles:
            xi_conv_out[multipole] = np.sum(
                [
                    term.coeff * self._Q_out[term.ind_Q] * xi_out[term.ind_Z]
                    for term in self._formulae[multipole]
                ],
                axis=0
            )

        return xi_conv_out


class PowspecWinConv(TwoPointWinConvBase):
    """Window convolution of power spectrum.

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'wilson+16').
    window_sampts : 1-d array of float
        Window function separation sample points.  If not logarithmically
        spaced, these are respaced with the same range and number
        of points; `window_multipoles` are resampled accordingly.
    window_multipoles : dict of {str or int: 1-d array of float}
        Window function multipole samples (each key is a multipole).
        If `window_sampts` needs to be logarithmically respaced, these are
        resampled accordingly.
    k_in : 1-d array of float
        Wavenumber sample points for the input power spectrum.  Must be
        logarithmically spaced.
    k_out : 1-d array of float, optional
        Wavenumber sample points for the output power spectrum.  If `None`
        (default), it is automatically derived.
    r_common : 1-d array of float, optional
        Separation sample points for multipoles of both the
        window function and the intermediate 2PCF (transformed from
        `pk_in`).  If `None` (default), they are automatically derived;
        otherwise, `r_common` must be logarithmically spaced.
    transform_kwargs : dict, optional
        FFTLog transform keyword arguments to be passed to
        :class:`~triumvirate.transforms.SphericalBesselTransform`
        (default is `None`).

    Attributes
    ----------
    k_in : 1-d array of float
        Logarithmically spaced wavenumber sample points for the
        input bispectrum.
    k_out : 1-d array of float
        Wavenumber sample points for the output bispectrum.
    r_common : 1-d array of float
        Logarithmically spaced separation sample points for multipoles
        of both the window function and the intermediate 2PCF
        (transformed from `pk_in`).

    Raises
    ------
    :exc:`~triumvirate._arrayops.SpacingError`
        When `k_in` or `r_common` is not logarithmically spaced.


    .. attention::

        If the window convolution formulae use `str` or `int` to represent
        a multipole, the corresponding `dict` of window function multipole
        samples and any power spectrum multipole samples to be convolved
        must use exactly the same key type.  In other words, the type
        representation of a multipole must be consistent across
        all inputs.

    .. attention::

        Integral constraint is ignored.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 k_in, k_out=None, r_common=None, transform_kwargs=None):

        super().__init__(formulae, window_sampts, window_multipoles)

        # Instantiate default transform parameters.
        self._transform_kwargs = transform_kwargs or {}
        # if 'pivot' not in self._transform_kwargs:
        #     self._transform_kwargs['pivot'] = 1.
        # if 'lowring' not in self._transform_kwargs:
        #     self._transform_kwargs['lowring'] = True
        if 'extrap' not in self._transform_kwargs:
            self._transform_kwargs['extrap'] = 4
        # if 'extrap_exp' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap_exp'] = 2.
        # if 'extrap_layer' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap_layer'] = 'outer'
        if 'threaded' not in self._transform_kwargs:
            self._transform_kwargs['threaded'] = False

        # Set input power spectrum wavenumber sample points.
        try:
            _check_1d_array(k_in, check_loglin=True)
        except SpacingError:
            raise SpacingError(
                "Input power spectrum wavenumber sample points are not "
                "logarithmically spaced."
            )

        self.k_in = k_in

        # Instantiate forward transforms for each unwindowed
        # power spectrum multipole.
        self._bsjt = {
            multipole_Z: SphericalBesselTransform(
                degree=int(multipole_Z),
                bias=0,
                sample_pts=self.k_in,
                **self._transform_kwargs
            )
            for multipole_Z in self._formulae.multipoles_Z
            if multipole_Z != 'ic'
        }

        # Set intermediate 2PCF separation sample points.
        if r_common is None:
            r_common_min = max(
                self._bsjt[multipole_Z]._post_sampts.min()
                for multipole_Z in self._formulae.multipoles_Z
                if multipole_Z != 'ic'
            )
            r_common_max = min(
                self._bsjt[multipole_Z]._post_sampts.max()
                for multipole_Z in self._formulae.multipoles_Z
                if multipole_Z != 'ic'
            )
            size_r_common = min(
                self._bsjt[multipole_Z]._post_sampts.size
                for multipole_Z in self._formulae.multipoles_Z
                if multipole_Z != 'ic'
            )
            self.r_common = np.logspace(
                np.log10(r_common_min), np.log10(r_common_max), size_r_common,
                base=10.
            )
        else:
            try:
                _check_1d_array(r_common, check_loglin=True)
            except SpacingError:
                raise SpacingError(
                    "Common separation sample points are not "
                    "logarithmically spaced."
                )
            self.r_common = r_common

        # Instantiate backward transforms for each windowed 2PCF multipole.
        self._fsjt = {
            multipole: SphericalBesselTransform(
                degree=int(multipole),
                bias=0,
                sample_pts=self.r_common,
                **self._transform_kwargs
            )
            for multipole in self._formulae.multipoles
        }

        # Set output power spectrum wavenumber sample points.
        if k_out is None:
            k_out_min = max(
                self._fsjt[multipole]._post_sampts.min()
                for multipole in self._formulae.multipoles
            )
            k_out_max = min(
                self._fsjt[multipole]._post_sampts.max()
                for multipole in self._formulae.multipoles
            )
            size_k_out = min(
                self._fsjt[multipole]._post_sampts.size
                for multipole in self._formulae.multipoles
            )
            self.k_out = np.logspace(
                np.log10(k_out_min), np.log10(k_out_max), size_k_out,
                base=10.
            )
        else:
            self.k_out = k_out

        # Resample window function multipoles at common separation sample
        # points.  This overrides the attributes set by the base class.
        self._rQ_in = self.r_common
        self._Q_in = {}
        for multipole_Q, Qpole_in in window_multipoles.items():
            self._Q_in[multipole_Q] = RectBivariateSpline(
                window_sampts, window_sampts, Qpole_in
            )(self.r_common, self.r_common)

    def convolve(self, pk_in):
        """Convolve power spectrum multipoles.

        Parameters
        ----------
        pk_in : dict of {str or int: 1-d array of float}
            Input power spectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_in`.

        Returns
        -------
        pk_conv_out : dict of {str or int: 1-d :class:`numpy.ndarray`}
            Output windowed power spectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_out`.

        """
        # Backward transform power spectrum multipoles to 2PCF multipoles.
        xi_in = {}
        for multipole_Z in self._formulae.multipoles_Z:
            if multipole_Z == 'ic':
                continue
            r_in, xi_in[multipole_Z] = self._bsjt[multipole_Z].\
                transform_cosmo_multipoles(-1, pk_in[multipole_Z])

        # Perform 2PCF convolution.
        _winconv = TwoPCFWinConv(
            self._formulae, self._rQ_in, self._Q_in,
            r_in,
            r_out=self.r_common
        )
        xi_conv = _winconv.convolve(xi_in)

        # Forward transform 2PCF multipoles to power spectrum multipoles.
        pk_conv_out = {}
        for multipole in self._formulae.multipoles:
            k_conv, pk_conv_out[multipole] = self._fsjt[multipole].\
                transform_cosmo_multipoles(1, xi_conv[multipole])

        # Resample at output wavenumber sample points.
        pk_conv_out = {
            multipole: InterpolatedUnivariateSpline(
                k_conv, Ppole_conv_out
            )(self.k_out)
            for multipole, Ppole_conv_out in pk_conv_out.items()
        }

        return pk_conv_out


class ThreePointWinConvBase:
    """Generic window convolution of three-point statistics.

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'sugiyama+18').
    window_sampts : 1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str or tuple of int: 2-d array of float}
        Window function multipole samples (each key is a multipole).


    .. attention::

        All convolution assumes that the 2-d window function multipole
        samples are square matrices, i.e. they are sampled
        at the same sample points in both dimensions.  This also applies
        to any integral constraint components.

    .. attention::

        If the window convolution formulae use `str`, `int` or
        `tuple` of `int` to represent a multipole, the corresponding
        `dict` of window function multipole samples must use exactly
        the same key type.  In other words, the type representation of
        a multipole must be consistent across all inputs.

    """

    def __init__(self, formulae, window_sampts, window_multipoles):

        # Constuct formulae object.
        if isinstance(formulae, str):
            try:
                formulae = NAMED_FORMULAE[formulae.lower()]
            except KeyError:
                raise ValueError(f'Unknown named formulae: {formulae}.')
        if not isinstance(formulae, WinConvFormulae):
            formulae = WinConvFormulae(formulae)
        self._formulae = formulae

        if formulae.multipoles_Z[0] != 'ic':
            self._multipole_rep = formulae.multipoles_Z[0]
        else:
            self._multipole_rep = formulae.multipoles_Z[-1]

        # Store window multipoles.
        self._rQ_in = window_sampts
        self._Q_in = window_multipoles

        if not set(self._formulae.multipoles_Q).issubset(self._Q_in.keys()):
            raise ValueError(
                "The window function multipole samples must contain all "
                "the multipole keys required by the convolution formulae."
            )


# STYLE: Standard naming convention is not always followed below.
class ThreePCFWinConv(ThreePointWinConvBase):
    """Window convolution of three-point correlation function (3PCF).

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'sugiyama+18').
    window_sampts : 1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str or tuple of int: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    r_in : 1-d array of float or dict of {str or tuple: 1-d array-like}
        Separation sample points for the input 3PCF, either the same
        for all multipoles or a dictionary of sample points for each
        multipole.
    r_out : 1-d array of float, optional
        Separation sample points for all output 3PCF multipoles.
        If `None` (default), it is set to `r_in`.

    Attributes
    ----------
    r_in : dict of {str or tuple: 1-d array of float}
        Separation sample points for the input 3PCF.
    r_out : dict of {str or tuple: 1-d array of float}
        Separation sample points for the output 3PCF.


    .. attention::

        All convolution assumes that the 2-d window function and
        3PCF multipole samples are square matrices, i.e. they are sampled
        at the same sample points in both dimensions.  This also applies
        to any integral constraint components.

    .. attention::

        If the window convolution formulae use `str`, `int` or
        `tuple` of `int` to represent a multipole, the corresponding
        `dict` of window function multipole samples and any 3PCF
        multipole samples to be convolved must use exactly the same
        key type.  In other words, the type representation of a multipole
        must be consistent across all inputs.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 r_in, r_out=None):

        super().__init__(formulae, window_sampts, window_multipoles)

        if isinstance(r_in, dict):
            self.r_in = r_in
        else:
            self.r_in = {
                multipole: r_in
                for multipole in set(
                    formulae.multipoles_Z + formulae.multipoles_Q
                )
            }
        if r_out is None:
            self.r_out = r_in
        else:
            self.r_out = {
                multipole: r_out
                for multipole in set(
                    formulae.multipoles_Z +
                    formulae.multipoles_Q +
                    formulae.multipoles
                )
            }
        if not _check_conv_range(self.r_out, self._rQ_in):
            warnings.warn(
                "The convolution range `r_out` is not fully covered "
                "by the window function separation sample points. "
                "The window function samples may be extrapolated."
            )
        if not _check_conv_range(self.r_out, self.r_in):
            warnings.warn(
                "The convolution range `r_out` is not fully covered "
                "by the input 3PCF separation sample points. "
                "The input 3PCF samples may be extrapolated."
            )

        self._Q_out = {
            _multipole: RectBivariateSpline(
                self._rQ_in, self._rQ_in, _Qpole_in
            )(self.r_out[_multipole], self.r_out[_multipole])
            for _multipole, _Qpole_in in self._Q_in.items()
        }
        self._Qdiag_out = {
            _multipole: InterpolatedUnivariateSpline(
                self._rQ_in, np.diag(_Qpole_in)
            )(self.r_out[_multipole])
            for _multipole, _Qpole_in in self._Q_in.items()
        }

    def convolve(self, zeta_in, ic=None):
        r"""Convolve 3PCF multipoles.

        Parameters
        ----------
        zeta_in : dict of {str or tuple of int: 2-d array of float}
            Input 3PCF multipole samples (each key is a multipole)
            at sample points :attr:`r_in`.
        ic : float, optional
            Integral constraint :math:`\bar{\zeta}`.  If `None` (default)
            and required by the convolution formula, it is calculated
            from the input 3PCF and the window function.

        Returns
        -------
        zeta_conv_out : dict of {str or tuple of int: 2-d array of float}
            Output windowed 3PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_out`.

        Raises
        ------
        ValueError
            When the dimensions of the input 3PCF multipole samples do not
            match the separation sample points of the input 3PCF.

        """
        # Interpolate input multipoles at output sample points.
        try:
            zeta_out = {
                multipole: RectBivariateSpline(
                    self.r_in[multipole], self.r_in[multipole], Zpole_in.real
                )(self.r_out[multipole], self.r_out[multipole])
                for multipole, Zpole_in in zeta_in.items()
            }
        except IndexError:  # raised by :mod:`scipy.interpolate._fitpack2`
            raise ValueError(
                "The input 3PCF multipole samples must be square matrices "
                "matching the dimension of the input 3PCF "
                "separation sample points."
            )

        # Enforce integral constraint.
        if 'ic' in self._formulae.multipoles_Z:
            if ic is None:
                ic = calc_threept_ic(
                    self._rQ_in, self._Q_in, self.r_in, zeta_in,
                    r_common=self.r_out, approx=False
                )
            zeta_out['ic'] = ic

        # Perform convolution as multiplication.
        zeta_conv_out = {}
        for multipole in self._formulae.multipoles:
            zeta_conv_out[multipole] = np.sum(
                [
                    term.coeff * self._Q_out[term.ind_Q] * zeta_out[term.ind_Z]
                    for term in self._formulae[multipole]
                ],
                axis=0
            )

        return zeta_conv_out

    def convolve_diag(self, zeta_diag_in, ic=None):
        r"""Convolve diagonal 3PCF multipoles.

        Parameters
        ----------
        zeta_diag_in : dict of {str or tuple of int: 1-d array of float}
            Input diagonal 3PCF multipole samples (each key
            is a multipole) at sample points :attr:`r_in`.
        ic : float, optional
            Integral constraint :math:`\bar{\zeta}`.  Cannot be `None`
            (default) if required by the convolution formula.

        Returns
        -------
        zeta_diag_conv_out : dict of {str or tuple of int: 1-d array}
            Output windowed diagonal 3PCF multipole samples (each key
            is a multipole) at sample points :attr:`r_out`.

        Raises
        ------
        ValueError
            When the integral constraint is required by the convolution
            formulae but `ic` is `None`.

        """
        # Interpolate input multipoles at output sample points.
        try:
            zeta_diag_out = {
                multipole: InterpolatedUnivariateSpline(
                    self.r_in[multipole], Zpole_diag_in
                )(self.r_out[multipole])
                for multipole, Zpole_diag_in in zeta_diag_in.items()
            }
        except ValueError:  # raised by :mod:`scipy.interpolate._fitpack2`
            raise ValueError(
                "The input 3PCF multipole samples, including any integral "
                "constraint, must be vectors matching the dimension "
                "of the input 3PCF separation sample points."
            )

        # Enforce integral constraint.
        if 'ic' in self._formulae.multipoles_Z:
            if ic is None:
                raise ValueError(
                    "The integral constraint is required by the convolution "
                    "formulae but not provided."
                )
            zeta_diag_out['ic'] = ic

        # Perform convolution as multiplication.
        zeta_diag_conv_out = {}
        for multipole in self._formulae.multipoles:
            zeta_diag_conv_out[multipole] = np.sum(
                [
                    term.coeff
                    * self._Qdiag_out[term.ind_Q] * zeta_diag_out[term.ind_Z]
                    for term in self._formulae[multipole]
                ],
                axis=0
            )

        return zeta_diag_conv_out


class BispecWinConv(ThreePointWinConvBase):
    """Window convolution of bispectrum.

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'sugiyama+18').
    window_sampts : 1-d array of float
        Window function separation sample points.  If not logarithmically
        spaced, these are respaced with the same range and number
        of points; `window_multipoles` are resampled accordingly.
    window_multipoles : dict of {str or tuple of int: 2-d array of float}
        Window function multipole samples (each key is a multipole).
        If `window_sampts` needs to be logarithmically respaced, these are
        resampled accordingly.
    k_in : 1-d array of float
        Wavenumber sample points for the input bispectrum.  Must be
        logarithmically spaced.
    k_out : 1-d array of float, optional
        Wavenumber sample points for the output bispectrum.  If `None`
        (default), they are chosen to correspond to the separation
        sample points `r_common`.
    r_common : 1-d array of float, optional
        Separation sample points for multipoles of both the
        window function and the intermediate 3PCF (transformed from
        `B_in`).  If `None` (default), they are automatically derived;
        otherwise, `r_common` must be logarithmically spaced.
    transform_kwargs : dict, optional
        FFTLog transform keyword arguments to be passed to
        :class:`~triumvirate.transforms.DoubleSphericalBesselTransform`
        (default is `None`).

    Attributes
    ----------
    k_in : 1-d array of float
        Logarithmically spaced wavenumber sample points for the
        input bispectrum.
    k_out : 1-d array of float
        Wavenumber sample points for the output bispectrum.
    r_common : 1-d array of float
        Logarithmically spaced separation sample points for multipoles
        of both the window function and the intermediate 3PCF
        (transformed from `B_in`).

    Raises
    ------
    :exc:`~triumvirate._arrayops.SpacingError`
        When `k_in` or `r_common` is not logarithmically spaced.


    .. attention::

        All convolution assumes that the 2-d window function and
        bispectrum multipole samples are square matrices, i.e. they are
        sampled at the same sample points in both dimensions.  This also
        applies to any integral constraint components.

    .. attention::

        If the window convolution formulae use `str`, `int` or
        `tuple` of `int` to represent a multipole, the corresponding
        `dict` of window function multipole samples and any bispectrum
        multipole samples to be convolved must use exactly the same
        key type.  In other words, the type representation of a multipole
        must be consistent across all inputs.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 k_in, k_out=None, r_common=None, transform_kwargs=None):

        super().__init__(formulae, window_sampts, window_multipoles)

        # Initialise default transform parameters.
        self._transform_kwargs = transform_kwargs or {}
        # if 'pivot' not in self._transform_kwargs:
        #     self._transform_kwargs['pivot'] = 1.
        # if 'lowring' not in self._transform_kwargs:
        #     self._transform_kwargs['lowring'] = True
        if 'extrap' not in self._transform_kwargs:
            self._transform_kwargs['extrap'] = 5
        # if 'extrap_exp' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap_exp'] = 2.
        # if 'extrap2d' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap2d'] = False
        if 'threaded' not in self._transform_kwargs:
            self._transform_kwargs['threaded'] = False

        # Set input bispectrum wavenumber sample points.
        try:
            _check_1d_array(k_in, check_loglin=True)
        except SpacingError:
            raise SpacingError(
                "Input bispectrum wavenumber sample points are not "
                "logarithmically spaced."
            )
        self.k_in = k_in

        # Instantiate forward transforms for each unwindowed
        # bispectrum multipole.
        self._bdsjt = {
            multipole_Z: DoubleSphericalBesselTransform(
                degrees=(int(multipole_Z[0]), int(multipole_Z[1])),
                biases=(0, 0),
                sample_pts=self.k_in,
                **self._transform_kwargs
            )
            for multipole_Z in self._formulae.multipoles_Z
            if multipole_Z != 'ic'
        }

        # Set intermediate 3PCF separation sample points.
        if r_common is None:
            r_common_min = max(
                self._bdsjt[multipole_Z]._post_sampts.min()
                for multipole_Z in self._formulae.multipoles_Z
                if multipole_Z != 'ic'
            )
            r_common_max = min(
                self._bdsjt[multipole_Z]._post_sampts.max()
                for multipole_Z in self._formulae.multipoles_Z
                if multipole_Z != 'ic'
            )
            size_r_common = min(
                self._bdsjt[multipole_Z]._post_sampts.size
                for multipole_Z in self._formulae.multipoles_Z
                if multipole_Z != 'ic'
            )
            self.r_common = np.logspace(
                np.log10(r_common_min), np.log10(r_common_max), size_r_common,
                base=10.
            )
        else:
            try:
                _check_1d_array(r_common, check_loglin=True)
            except SpacingError:
                raise SpacingError(
                    "Common separation sample points are not "
                    "logarithmically spaced."
                )
            self.r_common = r_common

        # Instantiate backward transforms for each windowed 3PCF multipole.
        self._fdsjt = {
            multipole: DoubleSphericalBesselTransform(
                degrees=(int(multipole[0]), int(multipole[1])),
                biases=(0, 0),
                sample_pts=self.r_common,
                **self._transform_kwargs
            )
            for multipole in self._formulae.multipoles
        }

        # Set output bispectrum wavenumber sample points.
        if k_out is None:
            k_out_min = max(
                self._fdsjt[multipole]._post_sampts.min()
                for multipole in self._formulae.multipoles
            )
            k_out_max = min(
                self._fdsjt[multipole]._post_sampts.max()
                for multipole in self._formulae.multipoles
            )
            size_k_out = min(
                self._fdsjt[multipole]._post_sampts.size
                for multipole in self._formulae.multipoles
            )
            self.k_out = np.logspace(
                np.log10(k_out_min), np.log10(k_out_max), size_k_out,
                base=10.
            )
        else:
            self.k_out = k_out

        # Resample window function multipoles at common separation sample
        # points.  This overrides the attributes set by the base class.
        self._rQ_in = self.r_common
        self._Q_in = {}
        for multipole_Q, Qpole_in in window_multipoles.items():
            self._Q_in[multipole_Q] = RectBivariateSpline(
                window_sampts, window_sampts, Qpole_in
            )(self.r_common, self.r_common)

    def convolve(self, B_in, ic=None, ret_zeta=False):
        r"""Convolve bispectrum multipoles.

        Parameters
        ----------
        B_in : dict of {str or tuple of int: 2-d array of float}
            Input bispectrum multipole samples (each key is a multipole)
            at sample points :attr:`k_in`.
        ic : float, optional
            Integral constraint :math:`\bar{\zeta}`.  If `None` (default)
            and required by the convolution formula, it is calculated
            from the input bispectrum and the window function.
        ret_zeta : bool, optional
            Whether to return the intermediate 3PCF multipoles
            (default is `False`).

        Returns
        -------
        B_conv_out : dict of {str or tuple of int: 2-d array of float}
            Output windowed bispectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_out`.

        """
        # Backward transform bispectrum multipoles to 3PCF multipoles.
        zeta_in = {}
        for multipole_Z in self._formulae.multipoles_Z:
            if multipole_Z == 'ic':
                continue
            r_in, zeta_in[multipole_Z] = self._bdsjt[multipole_Z].\
                transform_cosmo_multipoles(-1, B_in[multipole_Z])

        # Perform 3PCF convolution.
        _winconv = ThreePCFWinConv(
            self._formulae, self._rQ_in, self._Q_in,
            r_in[0][:, 0],
            r_out=self.r_common
        )

        zeta_conv = _winconv.convolve(zeta_in, ic=ic)

        # Forward transform 3PCF multipoles to bispectrum multipoles.
        B_conv_out = {}
        for multipole in self._formulae.multipoles:
            k_conv, B_conv_out[multipole] = self._fdsjt[multipole].\
                transform_cosmo_multipoles(1, zeta_conv[multipole])

        # Resample at output wavenumber sample points.
        B_conv_out = {
            multipole:
                RectBivariateSpline(
                    k_conv[0][:, 0], k_conv[-1][0, :], Bpole_conv_out.real
                )(self.k_out, self.k_out) +
                1j * RectBivariateSpline(
                    k_conv[0][:, 0], k_conv[-1][0, :], Bpole_conv_out.imag
                )(self.k_out, self.k_out)
            for multipole, Bpole_conv_out in B_conv_out.items()
        }

        if ret_zeta:
            return zeta_conv, B_conv_out
        return B_conv_out
