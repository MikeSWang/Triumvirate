"""
Window Convolution (:mod:`~triumvirate.winconv`)
==========================================================================

.. versionadded:: 0.4.0


Perform window convolution of two- and three-point statistics.

.. autosummary::
    WinConvTerm
    WinConvFormulae
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
from sympy.physics.wigner import wigner_3j

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
    multipoles : list of str
        Windowed CF multipoles.
    multipoles_Q : array of str
        Required window function multipoles.
    multipoles_Z : array of str
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
        self.multipoles_Q = np.unique([
            term.ind_Q
            for formula in self.formulae.values()
            for term in formula
        ])
        self.multipoles_Z = np.unique([
            term.ind_Z
            for formula in self.formulae.values()
            for term in formula
        ])

    def __getitem__(self, multipole):
        """Get the formula for a specific windowed CF multipole.

        Parameters
        ----------
        multipole : str
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
        multipole : str
            Windowed CF multipole.
        symbol : str, optional
            Symbol for the unwindowed CF multipole (default is
            ``r"\\zeta"``).

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
    flat_multipoles : dict of {str: 1-d array of float}
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
    multipoles : list of str
        Window function multipoles.
    r : 1-d array of float
        window function separation sample points.
    Q : dict of {str: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    Q_diag : dict of {str: 1-d array of float}
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
        for degrees, Q_flat in sorted(flat_multipoles.items()):
            Q_mat = np.zeros((nbins, nbins))
            Q_mat[triu_indices] = Q_flat
            Q_mat.T[triu_indices] = Q_flat
            Q_in[degrees] = alpha_contrast**2 * Q_mat

        self.r = rQ_in
        self.Q = Q_in

        self.Q_diag = {}
        self.multipoles = []
        for degrees, Q_in_ in self.Q.items():
            self.Q_diag[degrees] = np.diag(Q_in_)
            self.multipoles.append(degrees)

        self.sources = ['array'] * len(self.multipoles)

        # Check sample spacing.
        try:
            _check_1d_array(self.r, check_loglin=True)
        except (MixedSignError, SpacingError,):
            if make_loglin:
                self.resample(
                    size=make_loglin if isinstance(make_loglin, int) else None
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
        filepaths : dict of {str: str}
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
                r1_eff, r2_eff, Qflat_raw, Qflat_shot = np.loadtxt(
                    filepath, unpack=True, usecols=usecols
                )
                flat_multipoles[multipole] = Qflat_raw - Qflat_shot
            else:
                if len(usecols) != 3:
                    raise ValueError(
                        "Must provide three columns to read from the file "
                        "when not subtracting shot-noise."
                    )
                r1_eff, r2_eff, Qflat_raw = np.loadtxt(
                    filepath, unpack=True, usecols=usecols
                )
                flat_multipoles[multipole] = Qflat_raw
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
            (key 'Q') and sample points (key 'r').

        Returns
        -------
        :class:`~triumvirate.winconv.ThreePointWindow`
            Three-point window function.

        """
        self = cls.__new__(cls)

        self.r = dict_obj['r']
        self.Q = dict_obj['Q']

        self.Q_diag = {
            degrees: np.diag(Q_in_)
            for degrees, Q_in_ in self.Q.items()
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
                    degrees: Q_in_[selector, :][:, selector]
                    for degrees, Q_in_ in self.Q.items()
                }
            })

        self.r = self.r[selector]
        self.Q = {
            degrees: Q_in_[selector, :][:, selector]
            for degrees, Q_in_ in self.Q.items()
        }
        self.Q_diag = {
            degrees: np.diag(Q_in_)
            for degrees, Q_in_ in self.Q.items()
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
            for degrees, Q_in_ in self.Q.items():
                _rQ, _Q[degrees] = resample_lglin([self.r]*2, Q_in_, size=size)
        elif spacing == 'lin':
            for degrees, Q_in_ in self.Q.items():
                _rQ, _Q[degrees] = resample_lin([self.r]*2, Q_in_, size=size)
        else:
            raise ValueError("Unknown spacing type: {}".format(spacing))

        self.r, self.Q = _rQ[0], _Q

        self.Q_diag = {
            degrees: np.diag(Q_in_)
            for degrees, Q_in_ in self.Q.items()
        }


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
    window_multipoles : dict of {str: 2-d array of float}
        Window function multipole samples (each key is a multipole).
        Must contain the ``'000'`` multipole.
    r_in : 1-d array of float
        Separation sample points for the input 3PCF.
    zeta_in : dict of {str: 2-d array of float}
        Input 3PCF multipole samples (each key is a multipole)
        at sample points `r_in`.
    r_common : 1-d array of float, optional
        Common separation sample points.  If `None` (default), it is the
        same as `r_in`.
    approx : bool, optional
        If `True` (default is `False`), include only the leading-order
        term with ``'000'`` multipoles as an approximation.

    Returns
    -------
    float
        Integral constraint.

    """
    if r_common is None:
        r_common = r_in

    if approx:
        multipoles_common = ['000']
    else:
        multipoles_common = set(window_multipoles.keys()) & set(zeta_in.keys())
        if '000' not in multipoles_common:
            raise ValueError(
                "The '000' multipole is required for calculating "
                "the integral constraint."
            )

    _Q_in = {
        _multipole: RectBivariateSpline(
            window_sampts, window_sampts, window_multipoles[_multipole]
        )(r_common, r_common)
        for _multipole in multipoles_common
    }
    _zeta_in = {
        _multipole: RectBivariateSpline(
            r_in, r_in, zeta_in[_multipole]
        )(r_common, r_common)
        for _multipole in multipoles_common
    }

    ic_nom = 0.
    for _multipole in multipoles_common:
        ell1, ell2, ELL = map(int, _multipole)

        _N = (2 * ell1 + 1) * (2 * ell2 + 1) * (2 * ELL + 1)
        _H2 = wigner_3j(ell1, ell2, ELL, 0, 0, 0) ** 2

        ic_nom += _N * _H2 * _integrate_2d_samples(
            r_common, r_common, _Q_in[_multipole] * _zeta_in[_multipole]
        )

    ic_denom = _integrate_2d_samples(r_common, r_common, _Q_in['000'])

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
    window_multipoles : dict of {str: 1-d array of float}
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

        # Store window multipoles.
        self._rQ_in = window_sampts
        self._Q_in = window_multipoles


class TwoPCFWinConv(TwoPointWinConvBase):
    """Window convolution of two-point correlation function (2PCF).

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'wilson+16').
    window_sampts : 1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str: 1-d array of float}
        Window function multipole samples (each key is a multipole).
    r_in : 1-d array of float
        Separation sample points for the input 2PCF.
    r_out : 1-d array of float, optional
        Separation sample points for the output 2PCF.  If `None`
        (default), it is the same as `r_in`.

    Attributes
    ----------
    r_in : 1-d array of float
        Separation sample points for the input 2PCF.
    r_out : 1-d array of float
        Separation sample points for the output 2PCF.


    .. attention::

        If the window convolution formulae use `str` or `int` to represent
        a multipole, the corresponding `dict` of window function multipole
        samples and any 2PCF multipole samples to be convolved must use
        exactly the same key type.  In other words, the type
        representation of a multipole must be consistent across
        all inputs.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 r_in, r_out=None):

        super().__init__(formulae, window_sampts, window_multipoles)

        self.r_in = r_in
        if r_out is None:
            self.r_out = r_in
        else:
            self.r_out = r_out

        self._Q_out = {
            _multipole: InterpolatedUnivariateSpline(
                self._rQ_in, _Qpole_in
            )(self.r_out)
            for _multipole, _Qpole_in in self._Q_in.items()
        }

    def convolve(self, xi_in):
        """Convolve 2PCF multipoles.

        Parameters
        ----------
        xi_in : dict of {str: 2-d array of float}
            Input 3PCF multipole samples (each key is a multipole)
            at sample points :attr:`r_in`.

        Returns
        -------
        xi_conv_out : dict of {str: 2-d :class:`numpy.ndarray`}
            Output windowed 3PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_out`.

        """
        # Interpolate input multipoles at output sample points.
        xi_out = {
            multipole: InterpolatedUnivariateSpline(
                self.r_in, Zpole_in
            )(self.r_out)
            for multipole, Zpole_in in xi_in.items()
        }

        # Perform convolution as multiplication.
        xi_conv_out = {}
        for multipole in self._formulae.multipoles:
            xi_conv_out[multipole] = np.add.reduce([
                term.coeff * self._Q_out[term.ind_Q] * xi_out[term.ind_Z]
                for term in self._formulae[multipole]
            ])

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
    window_multipoles : dict of {str: 1-d array of float}
        Window function multipole samples (each key is a multipole).
        If `window_sampts` needs to be logarithmically respaced, these are
        resampled accordingly.
    k_in : 1-d array of float
        Wavenumber sample points for the input power spectrum.  Must be
        logarithmically spaced.
    k_out : 1-d array of float, optional
        Wavenumber sample points for the output power spectrum.  If `None`
        (default), it is the same as `k_in`.
    r_common : 1-d array of float, optional
        Separation sample points for multipoles of both the
        window function and the intermediate 2PCF (transformed from
        `pk_in`).  If `None` (default), the window function separation
        sample points are used (after logarithmic respacing if necessary);
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

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 k_in, k_out=None, r_common=None, transform_kwargs=None):

        # Different to other child classes of `TwoPointWinConvBase`,
        # Hankel-like transforms are involved and the input arrays
        # need to be remoulded.
        if r_common is None:
            try:
                _check_1d_array(window_sampts, check_loglin=True)
                _rQ_in, _Q_in = window_sampts, window_multipoles
            except SpacingError:
                _Q_in = {}
                for multipole_Q, Qpole_in in window_multipoles.items():
                    _rQ_in, _Q_in[multipole_Q] = resample_lglin(
                        window_sampts, Qpole_in
                    )
                warnings.warn(
                    "Window function separation sample points are not "
                    "logarithmically spaced.  Window function multipoles "
                    "have been resampled accordingly.",
                    RuntimeWarning
                )
            finally:
                self.r_common = _rQ_in
        else:
            try:
                _check_1d_array(r_common, check_loglin=True)
            except SpacingError:
                raise SpacingError(
                    "Common separation sample points are not "
                    "logarithmically spaced."
                )

            _rQ_in = r_common
            _Q_in = {}
            for multipole_Q, Qpole_in in window_multipoles.items():
                _Q_in[multipole_Q] = InterpolatedUnivariateSpline(
                    window_sampts, Qpole_in
                )(r_common)

            self.r_common = r_common

        try:
            _check_1d_array(k_in, check_loglin=True)
        except SpacingError:
            raise SpacingError(
                "Input power spectrum wavenumber sample points are not "
                "logarithmically spaced."
            )

        self.k_in = k_in
        if k_out is None:
            self.k_out = k_in
        else:
            self.k_out = k_out

        super().__init__(formulae, _rQ_in, _Q_in)

        # Instantiate transforms for each multipole for repeated use.
        self._transform_kwargs = transform_kwargs or {}
        # if 'lowring' not in self._transform_kwargs:
        #     self._transform_kwargs['lowring'] = True
        if 'extrap' not in self._transform_kwargs:
            self._transform_kwargs['extrap'] = 3
        # if 'extrap_exp' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap_exp'] = 2.
        # if 'extrap_layer' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap_layer'] = 'outer'

        self._bsjt = {
            multipole_Z: SphericalBesselTransform(
                degree=int(multipole_Z),
                bias=0,
                sample_pts=self.k_in,
                **self._transform_kwargs
            )
            for multipole_Z in self._formulae.multipoles_Z
        }  # backward sj-transform for unwindowed power spectrum multipoles
        self._fsjt = {
            multipole: SphericalBesselTransform(
                degree=int(multipole),
                bias=0,
                sample_pts=self.r_common,
                **self._transform_kwargs
            )
            for multipole in self._formulae.multipoles
        }  # forward sj-transform for windowed 2PCF multipoles

    def convolve(self, pk_in):
        """Convolve power spectrum multipoles.

        Parameters
        ----------
        pk_in : dict of {str: 1-d array of float}
            Input power spectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_in`.

        Returns
        -------
        pk_conv_out : dict of {str: 1-d :class:`numpy.ndarray`}
            Output windowed power spectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_out`.

        """
        # Backward transform bispectrum multipoles to 3PCF multipoles.
        xi_in = {}
        for multipole_Z, Ppole_in in pk_in.items():
            r_in, xi_in[multipole_Z] = \
                self._bsjt[multipole_Z].transform(Ppole_in)

        # Resample at common separation sample points.
        xi_in = {
            multipole_Z: InterpolatedUnivariateSpline(
                r_in, Zpole_in
            )(self.r_common)
            for multipole_Z, Zpole_in in xi_in.items()
        }

        # Convolve 2PCF multipoles.
        xi_conv = {}
        for multipole in self._formulae.multipoles:
            xi_conv[multipole] = np.add.reduce([
                term.coeff * self._Q_in[term.ind_Q] * xi_in[term.ind_Z]
                for term in self._formulae[multipole]
            ])

        # Forward transform 2PCF multipoles to power spectrum multipoles.
        pk_conv_out = {}
        for multipole, Zpole_conv in xi_conv.items():
            k_conv, pk_conv_out[multipole] = \
                self._fsjt[multipole].transform(Zpole_conv)

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
    window_multipoles : dict of {str: 2-d array of float}
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

        # Store window multipoles.
        self._rQ_in = window_sampts
        self._Q_in = window_multipoles


class ThreePCFWinConv(ThreePointWinConvBase):
    """Window convolution of three-point correlation function (3PCF).

    Parameters
    ----------
    formulae : str, |formulae_dicttype| or |formulae_type|
        Window convolution formulae.  If a string, it is assumed to be
        a named formula (e.g. 'sugiyama+18').
    window_sampts : 1-d array of float
        Window function separation sample points.
    window_multipoles : dict of {str: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    r_in : 1-d array of float
        Separation sample points for the input 3PCF.
    r_out : 1-d array of float, optional
        Separation sample points for the output 3PCF.  If `None`
        (default), it is the same as `r_in`.

    Attributes
    ----------
    r_in : 1-d array of float
        Separation sample points for the input 3PCF.
    r_out : 1-d array of float
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

        self.r_in = r_in
        if r_out is None:
            self.r_out = r_in
        else:
            self.r_out = r_out

        self._Q_out = {
            _multipole: RectBivariateSpline(
                self._rQ_in, self._rQ_in, _Qpole_in
            )(self.r_out, self.r_out)
            for _multipole, _Qpole_in in self._Q_in.items()
        }
        self._Qdiag_out = {
            _multipole: InterpolatedUnivariateSpline(
                self._rQ_in, np.diag(_Qpole_in)
            )(self.r_out)
            for _multipole, _Qpole_in in self._Q_in.items()
        }

    def convolve(self, zeta_in):
        """Convolve 3PCF multipoles.

        Parameters
        ----------
        zeta_in : dict of {str: 2-d array of float}
            Input 3PCF multipole samples (each key is a multipole)
            at sample points :attr:`r_in`.

        Returns
        -------
        zeta_conv_out : dict of {str: 2-d :class:`numpy.ndarray`}
            Output windowed 3PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_out`.

        """
        # Interpolate input multipoles at output sample points.
        zeta_out = {
            multipole: RectBivariateSpline(
                self.r_in, self.r_in, Zpole_in
            )(self.r_out, self.r_out)
            for multipole, Zpole_in in zeta_in.items()
        }

        # Perform convolution as multiplication.
        zeta_conv_out = {}
        for multipole in self._formulae.multipoles:
            zeta_conv_out[multipole] = np.add.reduce([
                term.coeff * self._Q_out[term.ind_Q] * zeta_out[term.ind_Z]
                for term in self._formulae[multipole]
            ])

        return zeta_conv_out

    def convolve_diag(self, zeta_diag_in):
        """Convolve diagonal 3PCF multipoles.

        Parameters
        ----------
        zeta_diag_in : dict of {str: 1-d array of float}
            Input diagonal 3PCF multipole samples (each key
            is a multipole) at sample points :attr:`r_in`.

        Returns
        -------
        zeta_diag_conv_out : dict of {str: 1-d :class:`numpy.ndarray`}
            Output windowed diagonal 3PCF multipole samples (each key
            is a multipole) at sample points :attr:`r_out`.

        """
        # Interpolate input multipoles at output sample points.
        zeta_diag_out = {
            multipole: InterpolatedUnivariateSpline(
                self.r_in, Zpole_diag_in
            )(self.r_out)
            for multipole, Zpole_diag_in in zeta_diag_in.items()
        }

        # Perform convolution as multiplication.
        zeta_diag_conv_out = {}
        for multipole in self._formulae.multipoles:
            zeta_diag_conv_out[multipole] = np.add.reduce([
                term.coeff
                * self._Qdiag_out[term.ind_Q] * zeta_diag_out[term.ind_Z]
                for term in self._formulae[multipole]
            ])

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
    window_multipoles : dict of {str: 2-d array of float}
        Window function multipole samples (each key is a multipole).
        If `window_sampts` needs to be logarithmically respaced, these are
        resampled accordingly.
    k_in : 1-d array of float
        Wavenumber sample points for the input bispectrum.  Must be
        logarithmically spaced.
    k_out : 1-d array of float, optional
        Wavenumber sample points for the output bispectrum.  If `None`
        (default), it is the same as `k_in`.
    r_common : 1-d array of float, optional
        Separation sample points for multipoles of both the
        window function and the intermediate 3PCF (transformed from
        `bk_in`).  If `None` (default), the window function separation
        sample points are used (after logarithmic respacing if necessary);
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
        (transformed from `bk_in`).

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

        # Different to other child classes of `ThreePointWinConvBase`,
        # Hankel-like transforms are involved and the input arrays
        # need to be remoulded.
        if r_common is None:
            try:
                _check_1d_array(window_sampts, check_loglin=True)
                _rQ_in, _Q_in = window_sampts, window_multipoles
            except SpacingError:
                _Q_in = {}
                for multipole_Q, Qpole_in in window_multipoles.items():
                    _rQ_in, _Q_in[multipole_Q] = resample_lglin(
                        window_sampts, Qpole_in
                    )
                warnings.warn(
                    "Window function separation sample points are not "
                    "logarithmically spaced.  Window function multipoles "
                    "have been resampled accordingly.",
                    RuntimeWarning
                )
            finally:
                self.r_common = _rQ_in
        else:
            try:
                _check_1d_array(r_common, check_loglin=True)
            except SpacingError:
                raise SpacingError(
                    "Common separation sample points are not "
                    "logarithmically spaced."
                )

            _rQ_in = r_common
            _Q_in = {}
            for multipole_Q, Qpole_in in window_multipoles.items():
                _Q_in[multipole_Q] = RectBivariateSpline(
                    window_sampts, window_sampts, Qpole_in
                )(r_common, r_common)

            self.r_common = r_common

        try:
            _check_1d_array(k_in, check_loglin=True)
        except SpacingError:
            raise SpacingError(
                "Input bispectrum wavenumber sample points are not "
                "logarithmically spaced."
            )

        self.k_in = k_in
        if k_out is None:
            self.k_out = k_in
        else:
            self.k_out = k_out

        super().__init__(formulae, _rQ_in, _Q_in)

        # Instantiate transforms for each multipole for repeated use.
        self._transform_kwargs = transform_kwargs or {}
        # if 'lowring' not in self._transform_kwargs:
        #     self._transform_kwargs['lowring'] = True
        if 'extrap' not in self._transform_kwargs:
            self._transform_kwargs['extrap'] = 3
        # if 'extrap_exp' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap_exp'] = 2.
        # if 'extrap2d' not in self._transform_kwargs:
        #     self._transform_kwargs['extrap2d'] = False

        self._bdsjt = {
            multipole_Z: DoubleSphericalBesselTransform(
                degrees=(int(multipole_Z[0]), int(multipole_Z[1])),
                biases=(0, 0),
                sample_pts=self.k_in,
                **self._transform_kwargs
            )
            for multipole_Z in self._formulae.multipoles_Z
        }  # backward double sj-transform for unwindowed bispectrum multipoles
        self._fdsjt = {
            multipole: DoubleSphericalBesselTransform(
                degrees=(int(multipole[0]), int(multipole[1])),
                biases=(0, 0),
                sample_pts=self.r_common,
                **self._transform_kwargs
            )
            for multipole in self._formulae.multipoles
        }  # forward double sj-transform for windowed 3PCF multipoles

    def convolve(self, bk_in):
        """Convolve bispectrum multipoles.

        Parameters
        ----------
        bk_in : dict of {str: 2-d array of float}
            Input bispectrum multipole samples (each key is a multipole)
            at sample points :attr:`k_in`.

        Returns
        -------
        bk_conv_out : dict of {str: 2-d :class:`numpy.ndarray`}
            Output windowed bispectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_out`.

        """
        # Backward transform bispectrum multipoles to 3PCF multipoles.
        zeta_in = {}
        for multipole_Z, Bpole_in in bk_in.items():
            r_in, zeta_in[multipole_Z] = \
                self._bdsjt[multipole_Z].transform(Bpole_in)

        # Resample at common separation sample points.
        zeta_in = {
            multipole_Z: RectBivariateSpline(
                r_in[0][:, 0], r_in[-1][0, :], Zpole_in
            )(self.r_common, self.r_common)
            for multipole_Z, Zpole_in in zeta_in.items()
        }

        # Convolve 3PCF multipoles.
        zeta_conv = {}
        for multipole in self._formulae.multipoles:
            zeta_conv[multipole] = np.add.reduce([
                term.coeff * self._Q_in[term.ind_Q] * zeta_in[term.ind_Z]
                for term in self._formulae[multipole]
            ])

        # Forward transform 3PCF multipoles to bispectrum multipoles.
        bk_conv_out = {}
        for multipole, Zpole_conv in zeta_conv.items():
            k_conv, bk_conv_out[multipole] = \
                self._fdsjt[multipole].transform(Zpole_conv)

        # Resample at output wavenumber sample points.
        bk_conv_out = {
            multipole: RectBivariateSpline(
                k_conv[0][:, 0], k_conv[-1][0, :], Bpole_conv_out
            )(self.k_out, self.k_out)
            for multipole, Bpole_conv_out in bk_conv_out.items()
        }

        return bk_conv_out
