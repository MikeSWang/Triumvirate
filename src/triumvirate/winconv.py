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


.. |coeff_type| replace:: \
    float or int or :class:`fractions.Fraction` or \
    str or :class:`sympy.core.expr.Expr`

.. |multipole_type| replace:: \
    str or int or tuple or :class:`~triumvirate.winconv.Multipole`

.. |multipole_3pt_type| replace:: \
    str or tuple or :class:`~triumvirate.winconv.Multipole`

.. |wcterm_type| replace:: \
    :class:`~triumvirate.winconv.WinConvTerm` or \
    tuples of (|multipole_type|, |multipole_type|, |coeff_type|)

.. |formulae_dict_type| replace:: \
    dict of {|multipole_type|: sequence of |wcterm_type|}

.. |formulae_type| replace:: \
    dict of \
    {:class:`~triumvirate.winconv.Multipole`: \
    list of :class:`~triumvirate.winconv.WinConvTerm`}

.. |sampts_type| replace:: \
    1-d array of float or dict of {|multipole_type|: 1-d array of float}

.. |sampts_3pt_type| replace:: \
    1-d array of float or \
    dict of {|multipole_3pt_type|: 1-d array of float}

"""
# STYLE: Standard naming convention is not always followed in this module
# for readability and clarity.
import warnings
from dataclasses import dataclass
from fractions import Fraction

import numpy as np
try:
    from scipy.integrate import simpson           # scipy>=1.6
except ImportError:
    from scipy.integrate import simps as simpson  # scipy<1.6
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline
from sympy import latex, sympify
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
  Wilson et al., 2016.
  `MNRAS 464(3), 3121 <https://doi.org/10.1093/mnras/stw2576>`_
  [`arXiv:1511.07799 <https://arxiv.org/abs/1511.07799>`_];
- ``'sugiyama+18'`` (for plane-parallel three-point statistics):
  Sugiyama et al., 2018.
  `MNRAS 484(1), 364 <https://doi.org/10.1093/mnras/sty3249>`_
  [`arXiv:1803.02132 <https://arxiv.org/abs/1803.02132>`_].

The value for each key is itself a `dict` that is used to construct an
instance of :class:`~triumvirate.winconv.WinConvFormulae`.

"""

SPLINE = 3
"""Default spline order for interpolation.

.. seealso:: :mod:`scipy.interpolate`

"""

_EXT = 3
"""Default extrapolation order for diagonal window interpolation.

.. seealso:: :mod:`scipy.interpolate`

"""


class ConvolutionRangeWarning(UserWarning):
    """Warning issued when the convolution range is not fully covered by
    the sample points.

    """
    pass


@dataclass
class Multipole:
    """Multipole indexing representation.

    Parameters
    ----------
    multipole : |multipole_type|
        Multipole representation.

    Attributes
    ----------
    index : tuple of int
        Multipole index or indices (composite index).
    abstr : str
        Abbreviated/undelimited string representation, e.g.
        '000' for (0, 0, 0).

    Examples
    --------
    The following variables all correspond to the multipole (0, 0, 0):

    >>> multipole_a = Multipole('000')
    >>> multipole_b = Multipole((0, 0, 0))
    >>> multipole_c = Multipole([0, 0, '0'])
    >>> assert multipole_a == multipole_b == multipole_c
    >>> print(multipole_a)
    0,0,0
    >>> print(multipole_a.abstr)
    000


    .. caution::

        When using a string to represent the multipole, it is assumed
        that each multipole degree is a single digit unless the string
        is delimited by commas.  If any of the multipole degrees is more
        than one digit, then the undelimited string representation will
        not be unique.

    """

    def __init__(self, multipole):

        if isinstance(multipole, self.__class__):
            self.index = multipole.index
            self.abstr = multipole.abstr
        else:
            try:
                if isinstance(multipole, str):
                    multipole = multipole.split(',')
                    if len(multipole) == 1:
                        multipole = multipole.pop()
                self.index = tuple(map(int, multipole))
            except TypeError:
                try:
                    self.index = (int(multipole),)
                except TypeError:
                    raise ValueError("Invalid `multipole` type or value.")
            self.abstr = ''.join(map(str, self.index))

    def __repr__(self):
        return f"{self.__class__.__name__}({self.index})"

    def __str__(self):
        return ','.join(map(str, self.index))

    def __hash__(self):
        return hash(self.index)

    def __lt__(self, other):
        self._iscomp(other)

        if self.index[-1] < other.index[-1]:
            return True
        if self.index[-1] > other.index[-1]:
            return False

        for idx, idx_other in zip(self.index[:-2], other.index[:-2]):
            if idx < idx_other:
                return True
            if idx > idx_other:
                return False

        return False

    def __gt__(self, other):
        self._iscomp(other)

        if self.index[-1] > other.index[-1]:
            return True
        if self.index[-1] < other.index[-1]:
            return False

        for idx, idx_other in zip(self.index[:-2], other.index[:-2]):
            if idx > idx_other:
                return True
            if idx < idx_other:
                return False

        return False

    def __le__(self, other):
        return self < other or self == other

    def __ge__(self, other):
        return self > other or self == other

    def __getitem__(self, idx):
        return self.index[idx]

    def _iscomp(self, other):
        if len(self.index) != len(other.index):
            raise ValueError(
                "Cannot compare multipole indices of different lengths."
            )


def _N_prefactor(multipole=None, ell1=None, ell2=None, L=None):
    """Calculate the Wigner-3j--related pre-factor.

    Parameters
    ----------
    multipole : :class:`~triumvirate.winconv.Multipole`, optional
        Multipole indices.  If `None` (default), use the individual
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
        ell1, ell2, L = multipole.index

    return (2 * ell1 + 1) * (2 * ell2 + 1) * (2 * L + 1)


def _H_factor(ell1, ell2, ell3):
    """Calculate the Wigner-3j factor with zero spherical orders.

    Parameters
    ----------
    ell1, ell2, ell3 : int
        Multipole degrees.

    Returns
    -------
    :class:`sympy.core.expr.Expr`
        Wigner-3j factor.

    """
    return wigner_3j(ell1, ell2, ell3, 0, 0, 0)


def calc_threept_winconv_coeff(multipole, multipole_Q, multipole_Z,
                               symb=False):
    """Calculate the window convolution coefficient for a specific
    window convolution term.

    Parameters
    ----------
    multipole : |multipole_3pt_type|
        Windowed multipole indices.
    multipole_Q : |multipole_3pt_type|
        Window function multipole indices.
    multipole_Z : |multipole_3pt_type|
        Unwindowed correlation function multipole indices.
    symb : bool, optional
        If `True` (default is `False`), return the coefficient as a
        symbolic expression.

    Returns
    -------
    float or :class:`sympy.core.expr.Expr`
        Window convolution coefficient.  When `multipole_Z` is 'ic',
        return -1 if `multipole` and `multipole_Q` are the same,
        otherwise 0.

    """
    multipole = Multipole(multipole)
    multipole_Q = Multipole(multipole_Q)

    if multipole_Z == 'ic':
        return -1 if multipole_Q == multipole else 0
    multipole_Z = Multipole(multipole_Z)

    ell1, ell2, L = multipole.index
    ell1_, ell2_, L_ = multipole_Z.index
    ell1__, ell2__, L__ = multipole_Q.index

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

    return coeff if symb else float(coeff)


@dataclass
class WinConvTerm:
    r"""Window convolution term as a tuple of three multiplicative
    factors including the coefficient, window function multipole
    and unwindowed correlation function (CF) multipole.

    Parameters
    ----------
    multipole_Q : |multipole_type|
        Window function multipole index/indices.
    multipole_Z : |multipole_type|
        Unwindowed-CF multipole index/indices.
    coeff : |coeff_type|
        Numerical coefficient for the convolution term.

    Attributes
    ----------
    coeff : |coeff_type|
        Numerical coefficient for the convolution term.
    multipole_Q : :class:`~triumvirate.winconv.Multipole`
        Window function multipole index/indices.
    multipole_Z : :class:`~triumvirate.winconv.Multipole`
        Unwindowed CF multipole index/indices.
    index_Q : tuple of int
        Window function multipole index/indices as a tuple.
    index_Z : tuple of int
        Unwindowed CF multipole index/indices as a tuple.

    Examples
    --------
    The following variables all correspond to the term
    :math:`- 1 \cdot Q_{000} \zeta_{\mathrm{ic}}`:

    >>> term_a = WinConvTerm((0, 0, 0), 'ic', -1)
    >>> term_b = WinConvTerm('000', 'ic', Fraction(-1))
    >>> assert term_a == term_b
    >>> print(term_a)
    -1;0,0,0;ic

    """
    # Special pesudo-indices.
    _SPEC_INDICES = ('ic',)

    def __init__(self, multipole_Q, multipole_Z, coeff):

        self.multipole_Q = Multipole(multipole_Q)
        self.index_Q = self.multipole_Q.index

        if multipole_Z not in self._SPEC_INDICES:
            self.multipole_Z = Multipole(multipole_Z)
            self.index_Z = self.multipole_Z.index
        else:
            self.multipole_Z = multipole_Z
            self.index_Z = multipole_Z

        if isinstance(coeff, str):
            coeff = sympify(coeff)
        if coeff == 0:
            warnings.warn(
                "Coefficient of window convolution term is zero.",
                RuntimeWarning
            )
        self.coeff = coeff

    def __repr__(self):
        return f"{self.__class__.__name__}({self})"

    def __str__(self):
        return f"{self.coeff};{self.multipole_Q};{self.multipole_Z}"

    def get_latex_str(self, symbol=r"\zeta"):
        """Get the LaTeX string representing the window convolution
        formula for a specific windowed CF multipole.

        Parameters
        ----------
        symbol : str, optional
            Symbol for the unwindowed CF (default is ``r"\\zeta"``).

        Returns
        -------
        str
            LaTeX string representing the window convolution formula
            (without the maths-mode delimiters).

        """
        if self.coeff == 0:
            return ''

        if self.coeff < 0:
            sgn_str = '-'
        else:
            sgn_str = ''

        if (coeff_abs := abs(self.coeff)) == 1:
            coeff_str = ''
        else:
            coeff_str = latex(sympify(coeff_abs))

        coeff_str = (sgn_str + ' ' + coeff_str).strip()

        index_Q_str = ''.join(map(str, self.index_Q))
        index_Z_str = ''.join(map(str, self.index_Z))

        term_str = "{} Q_\\mathrm{{{}}} {}_\\mathrm{{{}}}".format(
            coeff_str, index_Q_str, symbol, index_Z_str
        ).strip()

        return term_str


class WinConvFormulae:
    """Window convolution formulae.

    The full set of formulae are encoded as a dictionary, where each key
    corresponds to a windowed correlation function (CF) multipole, and
    each value is a sequence of tuples each corresponding to
    a window convolution term.

    Parameters
    ----------
    formulae : |formulae_dict_type|
        Window convolution formulae.

    Attributes
    ----------
    formulae : |formulae_type|
        Window convolution formulae.
    multipoles : list of :class:`~triumvirate.winconv.Multipole`
        Windowed CF multipole index/indices.
    multipoles_Q : list of :class:`~triumvirate.winconv.Multipole`
        Required window function multipole index/indices.
    multipoles_Z : list of :class:`~triumvirate.winconv.Multipole` or str
        Required unwindowed CF multipole index/indices including
        pseudo-index such as 'ic' for the integral constraint.

    Examples
    --------
    The following `formulae` corresponds to a single formula for the
    three-point CF monopole (0, 0, 0) consisting of three terms, namely
    :math:`Q_{000} \\zeta_{000} - Q_{000} \\zeta_{\\mathrm{ic}}
    + \\frac{1}{3} Q_{110} \\zeta_{110}`:

    >>> formulae_dict = {
    ...     '000': [
    ...         ('000', '000', 1),
    ...         ('000', 'ic', -1),
    ...         ('110', '110', '1/3'),
    ...     ],
    ... }
    >>> formulae = WinConvFormulae(formulae_dict)

    It encodes the formula for the windowed (0, 0, 0) monopole only:

    >>> print(formulae.multipoles)
    [Multipole((0, 0, 0))]

    The window function multipoles required are:

    >>> print(formulae.multipoles_Q)
    [Multipole((0, 0, 0)), Multipole((1, 1, 0))]

    The unwindowed three-point CF multipoles required are:

    >>> print(formulae.multipoles_Z)
    [Multipole((0, 0, 0)), Multipole((1, 1, 0)), 'ic']

    The formular can be printed as a LaTeX string (without the maths-mode
    delimiters):

    >>> print(formulae.get_latex_expr('000'))
    Q_\\mathrm{000} \\zeta_\\mathrm{000} \
- Q_\\mathrm{000} \\zeta_\\mathrm{ic} \
+ \\frac{1}{3} Q_\\mathrm{110} \\zeta_\\mathrm{110}


    """

    def __init__(self, formulae):

        self._formulae = {}
        for multipole, formula in formulae.items():
            if not isinstance(multipole, Multipole):
                multipole = Multipole(multipole)

            terms = []
            for term in formula:
                if not isinstance(term, WinConvTerm):
                    term = WinConvTerm(*term)
                terms.append(term)

            self._formulae[multipole] = terms

        self.multipoles = sorted(self._formulae.keys())
        self.multipoles_Q = sorted(set([
            term.multipole_Q
            for formula in self._formulae.values()
            for term in formula
        ]))

        _multipoles_Z = []
        _multipoles_extras = []
        for formula in self._formulae.values():
            for term in formula:
                if term.multipole_Z not in WinConvTerm._SPEC_INDICES:
                    _multipoles_Z.append(term.multipole_Z)
                else:
                    _multipoles_extras.append(term.multipole_Z)
        self._multipoles_Z_true = sorted(set(_multipoles_Z))
        self._multipoles_Z_peusdo = sorted(set(_multipoles_extras))
        self.multipoles_Z = self._multipoles_Z_true + self._multipoles_Z_peusdo

    def __getitem__(self, multipole):
        """Get the formula for a specific windowed CF multipole.

        Parameters
        ----------
        multipole : str or int or tuple or \
                    :class:`~triumvirate.winconv.Multipole`
            Windowed CF multipole.

        Returns
        -------
        list of :class:`~triumvirate.winconv.WinConvTerm`
            List of window convolution terms for `multipole`.

        """
        if not isinstance(multipole, Multipole):
            multipole = Multipole(multipole)

        return self._formulae[multipole]

    def get_latex_expr(self, multipole, symbol=r"\zeta"):
        """Get the LaTeX string representing the window convolution
        formula for a specific windowed CF multipole.

        Parameters
        ----------
        multipole : str or int or tuple or \
                    :class:`~triumvirate.winconv.Multipole`
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
            term_str = term.get_latex_str(symbol=symbol)
            if not term_str:
                continue
            if idx > 0 and not term_str.startswith('-'):
                term_str = '+' + ' ' + term_str
            terms_str.append(term_str.strip())

        return ' '.join(terms_str)


class ThreePointWindow:
    """Three-point window function in configuration space.

    By construction, the window function is discretely sampled as a
    symmetric matrix where the upper triangular part is stored as a
    flattened array in row-major order, and each matrix element
    corresponds to a pair of separation sample points.

    Parameters
    ----------
    paired_sampts : 2-d array of float
        Paired separation window function sample points, i.e. each row
        is :math:`(r_1, r_2)` where :math:`r_1` and :math:`r_2` are
        the separation sample points for the first and second dimensions.
    flat_multipoles : dict of {|multipole_3pt_type|: 1-d array of float}
        Flattened window function multipole samples (each key is
        a multipole) evaluated at `paired_sampts`.
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
    multipoles : list of :class:`~triumvirate.winconv.Multipole`
        Window function multipoles.
    r : 1-d array of float
        window function separation sample points.
    Q : dict of {:class:`~triumvirate.winconv.Multipole`: 2-d array}
        Window function multipole samples (each key is a multipole).
    Q_diag : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
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


    .. attention::

        All window function multipole samples are assumed to be
        evaluated at the same separation sample points.

    """

    def __init__(self, paired_sampts, flat_multipoles, alpha_contrast=1.,
                 make_loglin=True):

        # Check dimensions.
        paired_sampts = np.squeeze(paired_sampts)
        if paired_sampts.ndim != 2:
            raise ValueError(
                "Paired separation sample points must be a 2-d array."
            )

        for flat_multipole_ in flat_multipoles.values():
            flat_multipole_ = np.squeeze(flat_multipole_)
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
        multipoles_in = []
        for multipole, Qpole_flat in flat_multipoles.items():
            Qpole_mat = np.zeros((nbins, nbins))
            Qpole_mat[triu_indices] = Qpole_flat
            Qpole_mat.T[triu_indices] = Qpole_flat

            multipole_in = Multipole(multipole)
            Q_in[multipole_in] = alpha_contrast**2 * Qpole_mat
            multipoles_in.append(multipole_in)

        self.r = rQ_in
        self.Q = Q_in
        self.multipoles = sorted(multipoles_in)

        self.Q_diag = {
            multipole: np.diag(Qpole)
            for multipole, Qpole in self.Q.items()
        }

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
                    "logarithmically spaced."
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
        filepaths : dict of {|multipole_3pt_type|: str}
            Window function file paths (values) for each multipole (key).
        subtract_shotnoise : bool, optional
            Whether to subtract the shot-noise contribution from the
            window function (default is `True`).
        usecols : sequence of int, optional
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
                r1_eff_, r2_eff_, Qpole_flat_raw, Qpole_flat_shot = np.loadtxt(
                    filepath, usecols=usecols, unpack=True
                )
                flat_multipoles[multipole] = Qpole_flat_raw - Qpole_flat_shot
            else:
                if len(usecols) != 3:
                    raise ValueError(
                        "Must provide three columns to read from the file "
                        "when not subtracting shot-noise."
                    )
                r1_eff_, r2_eff_, Qpole_flat_raw = np.loadtxt(
                    filepath, usecols=usecols, unpack=True
                )
                flat_multipoles[multipole] = Qpole_flat_raw

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
        dict_obj : dict of {str: dict or 1-d array of float}
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

        self.Q = {}
        self.Q_diag = {}
        multipoles = []
        for multipole_in, Qpole_in in dict_obj['Q'].items():
            multipole = Multipole(multipole_in)
            self.Q[multipole] = Qpole_in
            self.Q_diag[multipole] = np.diag(Qpole_in)
            multipoles.append(multipole)

        self.multipoles = sorted(multipoles)

        try:
            self.sources = dict_obj['sources']
        except KeyError:
            self.sources = ['array'] * len(self.multipoles)

        return self

    @classmethod
    def load_from_file(cls, filepath):
        """Load window function multipoles from a compressed NumPy file
        as saved by
        :meth:`~triumvirate.winconv.ThreePointWindow.save_to_file`.

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

        win_obj = np.load(filepath, allow_pickle=True)

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
        """Subselect samples based on a range of separation sample points.

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
            If returned (default is `None`), three-point window function
            for the subselected range of separation sample points.

        """
        if rrange is None:
            sel_range = np.full_like(self.r, True, dtype=bool)
        else:
            rmin, rmax = rrange
            rmin = rmin or 0.
            rmax = rmax or np.inf
            sel_range = np.logical_and(rmin <= self.r, self.r <= rmax)

        sel_idx = np.full_like(self.r, idxrange is None, dtype=bool)
        if idxrange is not None:
            sel_idx[slice(*idxrange)] = True

        selector = np.logical_and(sel_range, sel_idx)

        if not inplace:
            instance = self.__class__.load_from_dict({
                'r': self.r[selector],
                'Q': {
                    multipole: Q_in_[selector, :][:, selector]
                    for multipole, Q_in_ in self.Q.items()
                }
            })
            instance.sources = self.sources
            return instance

        self.r = self.r[selector]
        self.Q = {
            multipole: Qpole[selector, :][:, selector]
            for multipole, Qpole in self.Q.items()
        }
        self.Q_diag = {
            multipole: np.diag(Qpole)
            for multipole, Qpole in self.Q.items()
        }

    def resample(self, spacing='lglin', size=None, spline=SPLINE):
        """Resample window function multipole samples.

        Parameters
        ----------
        spacing : {'lglin', 'lin'}, optional
            Spacing type for the resampled separation sample points
            either 'lglin' (default) (logarithmic-linear) or 'lin'
            (linear).
        size : (tuple of) int, optional
            Number(s) of points to resample to (default is `None`).
        spline : int, optional
            Spline order for interpolation
            (default is :py:const:`SPLINE`).
            See :mod:`scipy.interpolate` for more details.

        """
        _rQ, _Q = None, {}
        if spacing == 'lglin':
            for multipole, Qpole in self.Q.items():
                _rQ, _Q[multipole] = resample_lglin(
                    self.r, Qpole, size=size, spline=spline
                )
        elif spacing == 'lin':
            for multipole, Qpole in self.Q.items():
                _rQ, _Q[multipole] = resample_lin(
                    self.r, Qpole, size=size, spline=spline
                )
        else:
            raise ValueError("Unknown spacing type: {}".format(spacing))

        self.r, self.Q = _rQ[0], _Q

        self.Q_diag = {
            multipole: np.diag(Qpole)
            for multipole, Qpole in self.Q.items()
        }


def _check_conv_range(r_conv, r_samp):
    """Check if the convolution range is consistent with the
    sample points.

    Parameters
    ----------
    r_conv : dict of {:class:`~triumvrate.winconv.Multipole`: 1-d array} \
             or 1-d array of float
        Convolution range sample points.
    r_samp : dict of {:class:`~triumvrate.winconv.Multipole`: 1-d array} \
             or 1-d array of float
        Separation sample points.

    Returns
    -------
    bool
        `True` if `r_conv` within `r_samp` range; `False` otherwise.


    .. note:: There is a small relative tolerance for the range check.

    """
    RTOL = 1.e-5

    def _check_coverage(rc, rs):
        # Check if the convolution range `rc` is fully covered by the
        # sample points `rs`.
        try:
            return (1 - RTOL) * rs.min() <= rc.min() \
                and (rc.max() <= (1 + RTOL) * rs.max())
        except AttributeError as err:
            raise ValueError(
                "{}. `rc` is of type {} and `rs` is of type {}."
                .format(err, type(rc), type(rs))
            )

    # STYLE: Un-Pythonic if-block and return statements for clarity.
    if isinstance(r_conv, dict):
        if isinstance(r_samp, dict):
            multipoles_common = set(r_conv.keys()) & set(r_samp.keys())
            return all(
                _check_coverage(r_conv[multipole], r_samp[multipole])
                for multipole in multipoles_common
            )
        else:
            return all(
                _check_coverage(r_conv[multipole], r_samp)
                for multipole in r_conv.keys()
            )
    else:
        return _check_coverage(r_conv, r_samp)


def _find_conv_range(*coords, spacing='lglin'):
    """Find a convolution range consistent with two sets of sample points.

    Parameters
    ----------
    *coords : 1-d array of float
        Sample points.
    spacing : {'lglin', 'lin'}, optional
        Spacing for the convolution range separation sample points,
        either 'lglin' (default) (logarithmic-linear) or 'lin' (linear).

    Returns
    -------
    coord_conv : 1-d array of float
        Convolution range sample points.

    """
    coord_min = max(np.min(coord) for coord in coords)
    coord_max = min(np.max(coord) for coord in coords)
    nsamp = min(len(coord) for coord in coords)

    if spacing == 'lglin':
        coord_conv = np.logspace(
            np.log10(coord_min), np.log10(coord_max), nsamp, base=10
        )
    elif spacing == 'lin':
        coord_conv = np.linspace(coord_min, coord_max, nsamp)
    else:
        raise ValueError("Unknown spacing type: {}".format(spacing))

    return coord_conv


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


def calc_threept_ic(window_sampts, window_multipoles, r, zeta,
                    r_common=None, approx=False, spline=SPLINE):
    """Compute three-point clustering statistic integral constraint.

    Parameters
    ----------
    window_sampts : |sampts_3pt_type|
        Window function multipole separation sample points, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    window_multipoles : dict of {|multipole_3pt_type|: 2-d array of float}
        Window function multipole samples (each key is a multipole).
        Must contain the (0, 0, 0) monopole.
    r : |sampts_3pt_type|
        Separation sample points for the input 3PCF multipoles, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    zeta : dict of {|multipole_3pt_type|: 2-d array of float}
        Input 3PCF multipole samples (each key is a multipole)
        at sample points `r_in`.
    r_common : 1-d array of float, optional
        Common separation sample points.  If `None` (default), it is the
        same as `r_in`.
    approx : bool, optional
        If `True` (default is `False`), include only the leading-order
        term with the (0, 0, 0) multipole as an approximation.
    spline : int, optional
        Spline order for interpolation (default is :py:const:`SPLINE`).
        See :mod:`scipy.interpolate` for more details.

    Returns
    -------
    float
        Integral constraint.

    """
    # Unify multipole key type.
    window_multipoles = {
        Multipole(key): val
        for key, val in window_multipoles.items()
    }
    zeta = {
        Multipole(key): val
        for key, val in zeta.items()
    }
    monopole = Multipole((0, 0, 0))

    # Unify sample point type.
    if isinstance(window_sampts, dict):
        window_sampts = {
            Multipole(key): val
            for key, val in window_sampts.items()
        }
    else:
        window_sampts = {
            Multipole(key): window_sampts
            for key in window_multipoles.keys()
        }
    if isinstance(r, dict):
        r = {Multipole(key): val for key, val in r.items()}
    else:
        r = {Multipole(key): r for key in zeta.keys()}

    # Find common multipoles.
    if approx:
        multipoles_common = [monopole]
    else:
        multipoles_common = set(window_multipoles.keys()) & set(zeta.keys())
        if monopole not in multipoles_common:
            raise ValueError(
                "The (0, 0, 0) monopole is required for calculating "
                "the integral constraint."
            )
    if not set(multipoles_common).issubset(window_sampts.keys()):
        raise ValueError(
            "The window function separation sample points must contain "
            "all the multipoles required for the integral constraint."
        )
    if not set(multipoles_common).issubset(r.keys()):
        raise ValueError(
            "The input 3PCF separation sample points must contain all "
            "the multipoles required for the integral constraint."
        )

    # Set common separation sample points.
    if r_common is None:
        r_common = r
    elif not isinstance(r_common, dict):
        r_common = {multipole: r_common for multipole in multipoles_common}

    if not _check_conv_range(r_common, r):
        warnings.warn(
            "Convolution range is not fully covered by the sample points "
            "of the input 3PCF multipoles. "
            "Inaccurate extrapolation may occur."
        )
    if not _check_conv_range(r_common, window_sampts):
        warnings.warn(
            "Convolution range is not fully covered by the sample points "
            "of the window function multipoles. "
            "Inaccurate extrapolation may occur."
        )

    # Resample at common separation sample points to compute integrals.
    Q_in = {
        multipole: RectBivariateSpline(
            window_sampts[multipole], window_sampts[multipole],
            window_multipoles[multipole],
            kx=spline, ky=spline
        )(r_common[multipole], r_common[multipole])
        for multipole in multipoles_common
    }
    Z_in = {
        multipole: RectBivariateSpline(
            r[multipole], r[multipole],
            zeta[multipole].real,
            kx=spline, ky=spline
        )(r_common[multipole], r_common[multipole])
        for multipole in multipoles_common
    }

    ic_nom = 0.
    for multipole in multipoles_common:
        # `float` conversion necessary to avoid slow-down in computation.
        _H2 = float(_H_factor(*multipole.index)) ** 2
        _N = _N_prefactor(multipole=multipole)

        ic_nom += _N * _H2 * _integrate_2d_samples(
            r_common[multipole], r_common[multipole],
            Q_in[multipole] * Z_in[multipole]
        )

    ic_denom = _integrate_2d_samples(
        r_common[monopole], r_common[monopole], Q_in[monopole]
    )

    return ic_nom / ic_denom


class WinConvBase:
    """Generic window convolution for correlation functions (CF).

    Parameters
    ----------
    formulae : str, |formulae_dict_type| or |formulae_type|
        Window convolution formulae (see
        :class:`~triumvirate.winconv.WinConvFormulae`).  If a string,
        it is assumed to be a named formula (see
        :py:const:`NAMED_FORMULAE`).
    window_sampts : |sampts_type|
        Window function multipole separation sample points, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    window_multipoles : dict of {|multipole_type|: array of float}
        Window function multipole samples (each key is a multipole).

    Attributes
    ----------
    r_in : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array} \
           or None
        Separation sample points for the input CF multipoles.
    r_out : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array} \
            or None
        Separation sample points for the output CF multipoles.

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

        for _multipole in self._formulae._multipoles_Z_true:
            self._multipole_rep = _multipole
            break
        else:
            raise ValueError(
                "The window function multipoles must contain at least "
                "one of the unwindowed CF multipole indices."
            )

        # Store window multipoles.
        if isinstance(window_sampts, dict):
            self._rQ_in = {
                Multipole(key): val
                for key, val in window_sampts.items()
            }
        else:
            self._rQ_in = {
                Multipole(key): window_sampts
                for key in window_multipoles.keys()
            }
        self._Q_in = {
            Multipole(key): val
            for key, val in window_multipoles.items()
        }

        if not set(self._formulae.multipoles_Q).issubset(self._rQ_in.keys()):
            raise ValueError(
                "The window function separation sample points must contain "
                "all the multipole keys required by the convolution formulae."
            )
        if not set(self._formulae.multipoles_Q).issubset(self._Q_in.keys()):
            raise ValueError(
                "The window function multipole samples must contain all "
                "the multipole keys required by the convolution formulae."
            )

        # Set null sample points.
        self.r_in, self.r_out = None, None

    def initialise_sampts(self, r_in, r_out=None, enforce_coverage=False):
        """Initialise sample points :attr:`r_in` and :attr:`r_out`
        for the input and output CF multipoles.

        Parameters
        ----------
        r_in : |sampts_type| or |sampts_3pt_type|
            Separation sample points for the input CF multipoles, either
            the same for all multipoles or a dictionary of sample points
            (value) for each multipole (key).
        r_out : |sampts_type| or |sampts_3pt_type|
            Separation sample points for the output CF multipoles, either
            the same for all multipoles or a dictionary of sample points
            (value) for each multipole (key).  If `None` (default),
            it is set to `r_in`.
        enforce_coverage : bool, optional
            If `True` (default is `False`) and `r_out` is not provided,
            enforce that the output CF separation sample points are fully
            covered by the separation sample points of the window function
            and input CF multipoles, and that they are logarithmically
            spaced.

        """
        # Set input and output separation sample points.
        if isinstance(r_in, dict):
            # Coerce multipole type and check multipole matching.
            self.r_in = {
                Multipole(multipole): r_in[multipole]
                for multipole in r_in.keys()
            }
            if not set(self.r_in.keys()).issubset(self._formulae.multipoles_Z):
                raise ValueError(
                    "The input 2PCF separation sample points must contain all "
                    "the multipole keys required by the convolution formulae."
                )
        else:
            self.r_in = {
                multipole: r_in
                for multipole in self._formulae.multipoles_Z
            }

        if r_out is None:
            if enforce_coverage:
                self.r_out = {}
                for multipole in self._formulae.multipoles:
                    coords = []
                    for term in self._formulae[multipole]:
                        coords.append(self._rQ_in[term.multipole_Q])
                        if term.multipole_Z not in \
                                WinConvTerm._SPEC_INDICES:
                            coords.append(self.r_in[term.multipole_Z])
                    self.r_out[multipole] = _find_conv_range(*coords)
            else:
                # Check multipole matching for sample points.
                try:
                    self.r_out = {
                        multipole: self.r_in[multipole]
                        for multipole in self._formulae.multipoles
                    }
                except KeyError:
                    raise ValueError(
                        "Cannot set the output CF separation sample points "
                        "to the input CF separation sample points, since "
                        "there is no matching input CF multipole for "
                        "the output CF multipole. "
                        "Please provide the output CF separation "
                        "sample points explicitly by setting `r_out` "
                        "or set `enforce_coverage` to True."
                    )
        elif isinstance(r_out, dict):
            self.r_out = {
                Multipole(key): val
                for key, val in r_out.items()
            }
        else:
            self.r_out = {
                multipole: r_out
                for multipole in self._formulae.multipoles
            }

        # Check coverage of sample points.
        for multipole in self._formulae.multipoles:
            for term in self._formulae[multipole]:
                if not _check_conv_range(
                    self.r_out[multipole], self._rQ_in[term.multipole_Q]
                ):
                    msg = (
                        "pole: {}; conv: {}, {}; wind: {}, {}".format(
                            multipole,
                            self.r_out[multipole].min(),
                            self.r_out[multipole].max(),
                            self._rQ_in[term.multipole_Q].min(),
                            self._rQ_in[term.multipole_Q].max()
                        )
                    )
                    warnings.warn(
                        "The convolution range `r_out` is not fully covered "
                        "by the window function separation sample points. "
                        "Inaccurate extrapolation may occur." + msg,
                        category=ConvolutionRangeWarning
                    )
                if (
                    term.multipole_Z not in WinConvTerm._SPEC_INDICES
                    and not _check_conv_range(
                        self.r_out[multipole], self.r_in[term.multipole_Z]
                    )
                ):
                    warnings.warn(
                        "The convolution range `r_out` is not fully covered "
                        "by the input CF separation sample points. "
                        "Inaccurate extrapolation may occur.",
                        category=ConvolutionRangeWarning
                    )


class TwoPCFWinConv(WinConvBase):
    """Window convolution of two-point correlation function (2PCF).

    Parameters
    ----------
    formulae : str, |formulae_dict_type| or |formulae_type|
        Window convolution formulae (see
        :class:`~triumvirate.winconv.WinConvFormulae`).  If a string,
        it is assumed to be a named formula (see
        :py:const:`NAMED_FORMULAE`).
    window_sampts : |sampts_type|
        Window function multipole separation sample points, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    window_multipoles : dict of {|multipole_type|: 1-d array of float}
        Window function multipole samples (each key is a multipole).
    r_in : |sampts_type|
        Separation sample points for the input 2PCF multipoles, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    r_out : 1-d array of float, optional
        Separation sample points for the output 2PCF multipoles.
        If `None` (default), it is set to `r_in`.
    spline : int, optional
        Spline order for interpolation (default is :py:const:`SPLINE`).
        See :mod:`scipy.interpolate` for more details.

    Attributes
    ----------
    r_in : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
        Separation sample points for the input 2PCF multipoles.
    r_out : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
        Separation sample points for the output 2PCF multipoles.
    spline : int
        Spline order for interpolation (default is :py:const:`SPLINE`).


    .. attention::

        Integral constraint is ignored.

    .. seealso:: :class:`~triumvirate.winconv.WinConvBase`

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 r_in, r_out=None, spline=SPLINE):

        super().__init__(formulae, window_sampts, window_multipoles)

        self.initialise_sampts(r_in, r_out=r_out)
        self.spline = spline

    def convolve(self, xi_in):
        """Convolve 2PCF multipoles.

        Parameters
        ----------
        xi_in : dict of {|multipole_type|: 1-d array of float}
            Input 2PCF multipole samples (each key is a multipole)
            at sample points :attr:`r_in`.

        Returns
        -------
        xi_conv_out : dict of {:class:`~triumvirate.winconv.Multipole`: \
                      1-d :class:`numpy.ndarray`}
            Output windowed 2PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_out`.

        """
        xi_conv_out = {}
        for multipole in self._formulae.multipoles:
            summands = []
            for term in self._formulae[multipole]:
                Qpole_in = InterpolatedUnivariateSpline(
                    self._rQ_in[term.multipole_Q],
                    self._Q_in[term.multipole_Q],
                    k=self.spline, ext=_EXT
                )(self.r_out[multipole])
                Zpole_in = InterpolatedUnivariateSpline(
                    self.r_in[term.multipole_Z],
                    xi_in[term.multipole_Z],
                    k=self.spline
                )(self.r_out[multipole])
                summands.append(float(term.coeff) * Qpole_in * Zpole_in)

            xi_conv_out[multipole] = np.sum(summands, axis=0)

        return xi_conv_out


class PowspecWinConv(WinConvBase):
    """Window convolution of power spectrum.

    Parameters
    ----------
    formulae : str, |formulae_dict_type| or |formulae_type|
        Window convolution formulae (see
        :class:`~triumvirate.winconv.WinConvFormulae`).  If a string,
        it is assumed to be a named formula (see
        :py:const:`NAMED_FORMULAE`).
    window_sampts : |sampts_type|
        Window function multipole separation sample points, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).  If not logarithmically
        spaced, these are respaced with the same range and number
        of points; `window_multipoles` are resampled accordingly.
    window_multipoles : dict of {|multipole_type| : 1-d array of float}
        Window function multipole samples (each key is a multipole).
        If `window_sampts` needs to be logarithmically respaced, these are
        resampled accordingly.
    k_in : |sampts_type|
        Wavenumber sample points for the input power spectrum multipoles,
        either the same for all multipoles or a dictionary of
        sample points (value) for each multipole (key).  Must be
        logarithmically spaced.
    k_out : 1-d array of float, optional
        Wavenumber sample points for the output power spectrum multipoles.
        If `None` (default), it is automatically derived.
    r_common : 1-d array of float, optional
        Separation sample points for multipoles of both the
        window function and the intermediate 2PCF (transformed from
        `pk_in`).  If `None` (default), it is automatically derived;
        otherwise, `r_common` must be logarithmically spaced.
    transform_kwargs : dict, optional
        FFTLog transform keyword arguments to be passed to
        :class:`~triumvirate.transforms.SphericalBesselTransform`
        (default is `None`).
    spline : int or tuple of int, optional
        Spline order for interpolation (default is :py:const:`SPLINE`).
        If a tuple, the first element is for configuration-space
        interpolation and the second element is for Fourier-space
        interpolation.  See :mod:`scipy.interpolate` for more details.

    Attributes
    ----------
    k_in : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
        Logarithmically spaced wavenumber sample points for the
        input power spectrum multipoles.
    k_out : 1-d array of float
        Wavenumber sample points for the output power spectrum.
    r_common : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
        Logarithmically spaced separation sample points for multipoles
        of both the window function and the intermediate 2PCF
        (transformed from `P_in`).
    spline : tuple of int
        Spline order (default is :py:const:`SPLINE`) for interpolation
        in Fourier space and in configuration space.

    Raises
    ------
    :exc:`~triumvirate._arrayops.SpacingError`
        When `k_in` or `r_common` is not logarithmically spaced.


    .. attention::

        Integral constraint is ignored.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 k_in, k_out=None, r_common=None,
                 transform_kwargs=None, spline=SPLINE):

        super().__init__(formulae, window_sampts, window_multipoles)

        # Set default transform parameters.
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

        self.spline = (spline, spline) if isinstance(spline, int) else spline

        # Set input wavenumber sample points.
        if isinstance(k_in, dict):
            # Coerce multipole type and check multipole matching.
            self.k_in = {
                Multipole(key): val
                for key, val in k_in.items()
            }
            if not set(self.k_in.keys()).issubset(self._formulae.multipoles_Z):
                raise ValueError(
                    "The input power spectrum wavenumber sample points must "
                    "contain all the multipole keys required by the "
                    "convolution formulae."
                )
        else:
            self.k_in = {
                multipole_Z: k_in
                for multipole_Z in self._formulae.multipoles_Z
            }

        try:
            for k in self.k_in.values():
                _check_1d_array(k, check_loglin=True)
        except SpacingError:
            raise SpacingError(
                "Input power spectrum wavenumber sample points are not "
                "logarithmically spaced."
            )

        # Instantiate forward transforms for each unwindowed
        # power spectrum multipole.
        self._bsjt = {
            multipole_Z: SphericalBesselTransform(
                degree=multipole_Z[0],
                bias=0,
                sample_pts=self.k_in[multipole_Z],
                **self._transform_kwargs
            )
            for multipole_Z in self._formulae._multipoles_Z_true
        }

        # Set intermediate 2PCF separation sample points.
        if r_common is not None:
            try:
                _check_1d_array(r_common, check_loglin=True)
            except SpacingError:
                raise SpacingError(
                    "Common separation sample points are not "
                    "logarithmically spaced."
                )
            self.initialise_sampts(
                r_in={
                    multipole_Z: self._bsjt[multipole_Z]._post_sampts
                    for multipole_Z in self._formulae._multipoles_Z_true
                },
                r_out=r_common
            )
        else:
            self.initialise_sampts(
                r_in={
                    multipole_Z: self._bsjt[multipole_Z]._post_sampts
                    for multipole_Z in self._formulae._multipoles_Z_true
                },
                enforce_coverage=False
            )
        self.r_common = self.r_out

        # Instantiate backward transforms for each windowed 2PCF multipole.
        self._fsjt = {
            multipole: SphericalBesselTransform(
                degree=multipole[0],
                bias=0,
                sample_pts=self.r_common[multipole],
                **self._transform_kwargs
            )
            for multipole in self._formulae.multipoles
        }

        # Set output power spectrum wavenumber sample points.
        if k_out is None:
            self.k_out = _find_conv_range(
                *[
                    self._fsjt[multipole]._post_sampts
                    for multipole in self._formulae.multipoles
                ]
            )
        else:
            self.k_out = k_out

        # Set up 2PCF convolution.
        self._winconv = TwoPCFWinConv(
            self._formulae, self._rQ_in, self._Q_in,
            self.r_in, r_out=self.r_common, spline=self.spline[0]
        )

    def convolve(self, P_in, ret_xi=False):
        """Convolve power spectrum multipoles.

        Parameters
        ----------
        P_in : dict of {|multipole_type|: 1-d array of float}
            Input power spectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_in`.
        ret_xi : bool, optional
            Whether to return the intermediate 2PCF multipoles
            (default is `False`).

        Returns
        -------
        P_conv_out : dict of {:class:`~triumvirate.winconv.Multipole`: \
                      1-d :class:`numpy.ndarray`}
            Output windowed power spectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_out`.
        xi_conv : dict of {:class:`~triumvirate.winconv.Multipole`: \
                  1-d :class:`numpy.ndarray`}, optional
            Intermediate windowed 2PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_common`
            (only returned if `ret_xi` is `True`).

        """
        # Backward transform power spectrum multipoles to 2PCF multipoles.
        xi_in = {}
        for multipole_Z in self._formulae._multipoles_Z_true:
            _, xi_in[multipole_Z] = self._bsjt[multipole_Z].\
                transform_cosmo_multipoles(-1, P_in[multipole_Z])

        # Perform 2PCF convolution.
        xi_conv = self._winconv.convolve(xi_in)

        # Forward transform 2PCF multipoles to power spectrum multipoles.
        P_conv_out = {}
        for multipole in self._formulae.multipoles:
            k_conv, pk_conv_pole = self._fsjt[multipole].\
                transform_cosmo_multipoles(1, xi_conv[multipole])
            P_conv_out[multipole] = InterpolatedUnivariateSpline(
                k_conv, pk_conv_pole.real, k=self.spline[-1]
            )(self.k_out) + 1.j * InterpolatedUnivariateSpline(
                k_conv, pk_conv_pole.imag, k=self.spline[-1]
            )(self.k_out)

        if ret_xi:
            return P_conv_out, xi_conv
        return P_conv_out


class ThreePCFWinConv(WinConvBase):
    """Window convolution of three-point correlation function (3PCF).

    Parameters
    ----------
    formulae : str, |formulae_dict_type| or |formulae_type|
        Window convolution formulae (see
        :class:`~triumvirate.winconv.WinConvFormulae`).  If a string,
        it is assumed to be a named formula (see
        :py:const:`NAMED_FORMULAE`).
    window_sampts : |sampts_3pt_type|
        Window function multipole separation sample points, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    window_multipoles : dict of {|multipole_3pt_type|: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    r_in : |sampts_3pt_type|
        Separation sample points for the input 3PCF multipoles, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).
    r_out : 1-d array of float, optional
        Separation sample points for all output 3PCF multipoles.
        If `None` (default), it is set to `r_in`.
    spline : int, optional
        Spline order for interpolation (default is :py:const:`SPLINE`).
        See :mod:`scipy.interpolate` for more details.

    Attributes
    ----------
    r_in : dict of {str or tuple: 1-d array of float}
        Separation sample points for the input 3PCF.
    r_out : dict of {str or tuple: 1-d array of float}
        Separation sample points for the output 3PCF.
    spline : int
        Spline order for interpolation (default is :py:const:`SPLINE`).


    .. attention::

        All convolution assumes that the window function and
        3PCF multipole samples are square matrices, i.e. they are sampled
        at the same sample points in both dimensions.

    .. seealso:: :class:`~triumvirate.winconv.WinConvBase`

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 r_in, r_out=None, spline=SPLINE):

        super().__init__(formulae, window_sampts, window_multipoles)
        self.spline = spline

        self.initialise_sampts(r_in, r_out=r_out)

    def convolve(self, zeta_in, ic=None):
        r"""Convolve 3PCF multipoles.

        Parameters
        ----------
        zeta_in : dict of {|multipole_3pt_type|: 2-d array of float}
            Input 3PCF multipole samples (each key is a multipole)
            at sample points :attr:`r_in`.
        ic : float, optional
            Integral constraint :math:`\bar{\zeta}`.  If `None` (default)
            and required by the convolution formula, it is calculated
            from the input 3PCF and the window function.

        Returns
        -------
        zeta_conv_out : dict of {:class:`~triumvirate.winconv.Multipole`: \\
                        2-d :class:`numpy.ndarray`}
            Output windowed 3PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_out`.

        Raises
        ------
        ValueError
            When the dimensions of the input 3PCF multipole samples do not
            match the separation sample points of the input 3PCF.

        """
        # Enforce integral constraint.
        if 'ic' in self._formulae.multipoles_Z and ic is None:
            ic = calc_threept_ic(
                self._rQ_in, self._Q_in, self.r_in, zeta_in,
                r_common=None, approx=False
            )

        # Perform convolution term by term.
        zeta_conv_out = {}
        for multipole in self._formulae.multipoles:
            summands = []
            for term in self._formulae[multipole]:
                Qpole_in = RectBivariateSpline(
                    self._rQ_in[term.multipole_Q],
                    self._rQ_in[term.multipole_Q],
                    self._Q_in[term.multipole_Q],
                    kx=self.spline, ky=self.spline
                )(self.r_out[multipole], self.r_out[multipole])

                if term.multipole_Z == 'ic':
                    Zpole_in = ic
                else:
                    try:
                        Zpole_in = RectBivariateSpline(
                            self.r_in[term.multipole_Z],
                            self.r_in[term.multipole_Z],
                            zeta_in[term.multipole_Z].real,
                            kx=self.spline, ky=self.spline
                        )(self.r_out[multipole], self.r_out[multipole])
                    # Raised by :mod:`scipy.interpolate._fitpack2`.
                    except IndexError:
                        raise ValueError(
                            "The input 3PCF multipole samples must be "
                            "square matrices matching the dimension of "
                            "the input 3PCF separation sample points."
                        )

                summands.append(float(term.coeff) * Qpole_in * Zpole_in)

            zeta_conv_out[multipole] = np.sum(summands, axis=0)

        return zeta_conv_out

    def convolve_diag(self, zeta_diag_in, ic=None):
        r"""Convolve diagonal 3PCF multipoles.

        Parameters
        ----------
        zeta_diag_in : dict of {|multipole_3pt_type|: 1-d array of float}
            Input diagonal 3PCF multipole samples (each key
            is a multipole) at sample points :attr:`r_in`.
        ic : float, optional
            Integral constraint :math:`\bar{\zeta}`.  Cannot be `None`
            (default) if required by the convolution formula.

        Returns
        -------
        zeta_diag_conv_out : dict of {\\
                             :class:`~triumvirate.winconv.Multipole`: \\
                             1-d :class:`numpy.ndarray`}
            Output windowed diagonal 3PCF multipole samples (each key
            is a multipole) at sample points :attr:`r_out`.

        Raises
        ------
        ValueError
            When the integral constraint is required by the convolution
            formulae but `ic` is `None`.

        """
        # Enforce integral constraint.
        if 'ic' in self._formulae.multipoles_Z and ic is None:
            raise ValueError(
                "The integral constraint is required by the convolution "
                "formulae but not provided."
            )

        # Perform convolution term by term.
        zeta_diag_conv_out = {}
        for multipole in self._formulae.multipoles:
            summands = []
            for term in self._formulae[multipole]:
                Qpole_in = InterpolatedUnivariateSpline(
                    self._rQ_in[term.multipole_Q],
                    np.diag(self._Q_in[term.multipole_Q]),
                    k=self.spline, ext=_EXT
                )(self.r_out[multipole])

                if term.multipole_Z == 'ic':
                    Zpole_in = ic
                else:
                    try:
                        Zpole_in = InterpolatedUnivariateSpline(
                            self.r_in[term.multipole_Z],
                            zeta_diag_in[term.multipole_Z],
                            k=self.spline
                        )(self.r_out[multipole])
                    # Raised by :mod:`scipy.interpolate._fitpack2`.
                    except IndexError:
                        raise ValueError(
                            "The input diagonal 3PCF multipole samples "
                            "must be vectors matching the dimension of "
                            "the input 3PCF separation sample points."
                        )

                summands.append(float(term.coeff) * Qpole_in * Zpole_in)

            zeta_diag_conv_out[multipole] = np.sum(summands, axis=0)

        return zeta_diag_conv_out


class BispecWinConv(WinConvBase):
    """Window convolution of bispectrum.

    Parameters
    ----------
    formulae : str, |formulae_dict_type| or |formulae_type|
        Window convolution formulae (see
        :class:`~triumvirate.winconv.WinConvFormulae`).  If a string,
        it is assumed to be a named formula (see
        :py:const:`NAMED_FORMULAE`).
    window_sampts : |sampts_3pt_type|
        Window function separation sample points, either
        the same for all multipoles or a dictionary of sample points
        (value) for each multipole (key).  If not logarithmically
        spaced, these are respaced with the same range and number
        of points; `window_multipoles` are resampled accordingly.
    window_multipoles : dict of {|multipole_3pt_type|: 2-d array of float}
        Window function multipole samples (each key is a multipole).
        If `window_sampts` needs to be logarithmically respaced, these are
        resampled accordingly.
    k_in : |sampts_3pt_type|
        Wavenumber sample points for the input bispectrum multipoles,
        either the same for all multipoles or a dictionary of
        sample points (value) for each multipole (key).  Must be
        logarithmically spaced.
    k_out : 1-d array of float, optional
        Wavenumber sample points for the output bispectrum multipoles.
        If `None` (default), it is automatically derived.
    r_common : 1-d array of float, optional
        Separation sample points for multipoles of both the
        window function and the intermediate 3PCF (transformed from
        `B_in`).  If `None` (default), it is automatically derived;
        otherwise, `r_common` must be logarithmically spaced.
    transform_kwargs : dict, optional
        FFTLog transform keyword arguments to be passed to
        :class:`~triumvirate.transforms.DoubleSphericalBesselTransform`
        (default is `None`).
    spline : int or tuple of int, optional
        Spline order for interpolation (default is :py:const:`SPLINE`).
        If a tuple, the first element is for configuration-space
        interpolation and the second element is for Fourier-space
        interpolation. See :mod:`scipy.interpolate` for more details.

    Attributes
    ----------
    k_in : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
        Logarithmically spaced wavenumber sample points for the
        input bispectrum.
    k_out : 1-d array of float
        Wavenumber sample points for the output bispectrum.
    r_common : dict of {:class:`~triumvirate.winconv.Multipole`: 1-d array}
        Logarithmically spaced separation sample points for multipoles
        of both the window function and the intermediate 3PCF
        (transformed from `B_in`).
    spline : tuple of int
        Spline order (default is :py:const:`SPLINE`) for interpolation
        in Fourier space and in configuration space.

    Raises
    ------
    :exc:`~triumvirate._arrayops.SpacingError`
        When `k_in` or `r_common` is not logarithmically spaced.


    .. attention::

        All convolution assumes that the 2-d window function and
        bispectrum multipole samples are square matrices, i.e. they are
        sampled at the same sample points in both dimensions.  This also
        applies to any integral constraint components.

    """

    def __init__(self, formulae, window_sampts, window_multipoles,
                 k_in, k_out=None, r_common=None,
                 transform_kwargs=None, spline=SPLINE):

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

        self.spline = (spline, spline) if isinstance(spline, int) else spline

        # Set input bispectrum wavenumber sample points.
        if isinstance(k_in, dict):
            # Coerce multipole type and check multipole matching.
            self.k_in = {
                Multipole(key): val
                for key, val in k_in.items()
            }
            if not set(self.k_in.keys()).issubset(self._formulae.multipoles_Z):
                raise ValueError(
                    "The input bispectrum wavenumber sample points must "
                    "contain all the multipole keys required by the "
                    "convolution formulae."
                )
        else:
            self.k_in = {
                multipole_Z: k_in
                for multipole_Z in self._formulae.multipoles_Z
            }

        try:
            for k in self.k_in.values():
                _check_1d_array(k, check_loglin=True)
        except SpacingError:
            raise SpacingError(
                "Input bispectrum wavenumber sample points are not "
                "logarithmically spaced."
            )

        # Instantiate forward transforms for each unwindowed
        # bispectrum multipole.
        self._bdsjt = {
            multipole_Z: DoubleSphericalBesselTransform(
                degrees=(multipole_Z[0], multipole_Z[1]),
                biases=(0, 0),
                sample_pts=self.k_in[multipole_Z],
                **self._transform_kwargs
            )
            for multipole_Z in self._formulae._multipoles_Z_true
        }

        # Set intermediate 3PCF separation sample points.
        if r_common is not None:
            try:
                _check_1d_array(r_common, check_loglin=True)
            except SpacingError:
                raise SpacingError(
                    "Common separation sample points are not "
                    "logarithmically spaced."
                )
            self.initialise_sampts(
                r_in={
                    multipole_Z: self._bdsjt[multipole_Z]._post_sampts
                    for multipole_Z in self._formulae._multipoles_Z_true
                },
                r_out=r_common
            )
        else:
            self.initialise_sampts(
                r_in={
                    multipole_Z: self._bdsjt[multipole_Z]._post_sampts
                    for multipole_Z in self._formulae._multipoles_Z_true
                },
                enforce_coverage=False
            )
        self.r_common = self.r_out

        # Instantiate backward transforms for each windowed 3PCF multipole.
        self._fdsjt = {
            multipole: DoubleSphericalBesselTransform(
                degrees=(multipole[0], multipole[1]),
                biases=(0, 0),
                sample_pts=self.r_common[multipole],
                **self._transform_kwargs
            )
            for multipole in self._formulae.multipoles
        }

        # Set output bispectrum wavenumber sample points.
        if k_out is None:
            self.k_out = _find_conv_range(
                *[
                    self._fdsjt[multipole]._post_sampts
                    for multipole in self._formulae.multipoles
                ]
            )
        else:
            self.k_out = k_out

        # Set up 3PCF convolution.
        self._winconv = ThreePCFWinConv(
            self._formulae, self._rQ_in, self._Q_in,
            self.r_in, r_out=self.r_common, spline=self.spline[0]
        )

    def convolve(self, B_in, ic=None, ret_zeta=False):
        r"""Convolve bispectrum multipoles.

        Parameters
        ----------
        B_in : dict of {|multipole_3pt_type|: 2-d array of float}
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
        B_conv_out : dict of {:class:`~triumvirate.winconv.Multipole`: \\
                     2-d :class:`numpy.ndarray`}
            Output windowed bispectrum multipole samples (each key is
            a multipole) at sample points :attr:`k_out`.
        zeta_conv : dict of {:class:`~triumvirate.winconv.Multipole`: \\
                    2-d :class:`numpy.ndarray`}
            Intermediate windowed 3PCF multipole samples (each key is
            a multipole) at sample points :attr:`r_common`
            (only returned if `ret_zeta` is `True`).

        """
        # Backward transform bispectrum multipoles to 3PCF multipoles.
        zeta_in = {}
        for multipole_Z in self._formulae._multipoles_Z_true:
            _, zeta_in[multipole_Z] = self._bdsjt[multipole_Z].\
                transform_cosmo_multipoles(-1, B_in[multipole_Z])

        # Perform 3PCF convolution.
        zeta_conv = self._winconv.convolve(zeta_in, ic=ic)

        # Forward transform 3PCF multipoles to bispectrum multipoles.
        B_conv_out = {}
        for multipole in self._formulae.multipoles:
            transformer = self._fdsjt[multipole]

            k_conv = transformer._post_sampts
            _, bk_conv_pole = transformer.transform_cosmo_multipoles(
                1, zeta_conv[multipole]
            )

            B_conv_out[multipole] = RectBivariateSpline(
                k_conv, k_conv, bk_conv_pole.real,
                kx=self.spline[-1], ky=self.spline[-1]
            )(self.k_out, self.k_out) + 1j * RectBivariateSpline(
                k_conv, k_conv, bk_conv_pole.imag,
                kx=self.spline[-1], ky=self.spline[-1]
            )(self.k_out, self.k_out)

        if ret_zeta:
            return B_conv_out, zeta_conv
        return B_conv_out

    # def convolve_diag(self, B_diag_in, ic=None, ret_zeta_diag=False):
    #     r"""Convolve bispectrum multipoles assuming a diagonal form.

    #     Parameters
    #     ----------
    #     B_diag_in : dict of {|multipole_3pt_type|: 1-d array of float}
    #         Input bispectrum multipole samples (each key is a multipole)
    #         at sample points :attr:`k_in`.
    #     ic : float, optional
    #         Integral constraint :math:`\bar{\zeta}` (default is `None`).
    #         If required by the convolution formula, it must be provided.
    #     ret_zeta_diag : bool, optional
    #         Whether to return the intermediate diagonal 3PCF multipoles
    #         (default is `False`).

    #     Returns
    #     -------
    #     B_diag_conv_out : dict of \\
    #                       {:class:`~triumvirate.winconv.Multipole`: \\
    #                       1-d :class:`numpy.ndarray`}
    #         Output windowed diagonal bispectrum multipole samples
    #         (each key is a multipole) at sample points :attr:`k_out`.
    #     zeta_diag_conv : dict of \\
    #                      {:class:`~triumvirate.winconv.Multipole`: \\
    #                      1-d :class:`numpy.ndarray`}
    #         Intermediate windowed diagonal 3PCF multipole samples
    #         (each key is a multipole) at sample points :attr:`r_common`
    #         (only returned if `ret_zeta` is `True`).

    #     """
    #     # Backward transform bispectrum multipoles to 3PCF multipoles.
    #     zeta_in = {}
    #     for multipole_Z in self._formulae._multipoles_Z_true:
    #         _, zeta_in[multipole_Z] = self._bdsjt[multipole_Z].\
    #             transform_cosmo_multipoles(-1, np.diag(B_diag_in[multipole_Z]))

    #     # Perform 3PCF convolution.
    #     zeta_conv = self._winconv.convolve(zeta_in, ic=ic)
    #     zeta_diag_conv = {
    #         multipole: np.diag(Zpole_conv)
    #         for multipole, Zpole_conv in zeta_conv.items()
    #     }

    #     # Forward transform 3PCF multipoles to bispectrum multipoles.
    #     B_conv_out = {}
    #     for multipole in self._formulae.multipoles:
    #         transformer = self._fdsjt[multipole]

    #         k_conv = transformer._post_sampts
    #         _, bk_conv_pole = transformer.transform_cosmo_multipoles(
    #             1, zeta_conv[multipole]
    #         )

    #         B_conv_out[multipole] = RectBivariateSpline(
    #             k_conv, k_conv, bk_conv_pole.real,
    #             kx=self.spline[-1], ky=self.spline[-1]
    #         )(self.k_out, self.k_out) + 1j * RectBivariateSpline(
    #             k_conv, k_conv, bk_conv_pole.imag,
    #             kx=self.spline[-1], ky=self.spline[-1]
    #         )(self.k_out, self.k_out)
    #     B_diag_conv_out = {
    #         multipole: np.diag(Bpole_conv)
    #         for multipole, Bpole_conv in B_conv_out.items()
    #     }

    #     if ret_zeta_diag:
    #         return B_diag_conv_out, zeta_diag_conv
    #     return B_diag_conv_out
