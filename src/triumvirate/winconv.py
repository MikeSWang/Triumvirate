"""
Window Convolution (:mod:`~triumvirate.winconv`)
==========================================================================

.. versionadded:: 0.4.0


Perform window convolution of two- and three-point statistics.

.. autosummary::
    WinConvTerm
    WinConvFormulae
    TwoPointWinConv
    ThreePointWinConv

"""
import warnings
from collections import namedtuple
from fractions import Fraction

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline, RectBivariateSpline

from triumvirate.transforms import (
    DoubleSphericalBesselTransform,
    SphericalBesselTransform,
    resample_lglin,
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
    Q_\\mathrm{000} \\zeta_\\mathrm{000} \
    - Q_\\mathrm{000} \\zeta_\\mathrm{ic} \
    + \\frac{1}{3} Q_\\mathrm{110} \\zeta_\\mathrm{110}


    .. |formulae_dicttype| replace:: \
        dict of \
        {|pole_type|: sequence of tuples of \
        (|pole_type|, |pole_type|, |coeff_type|)}

    .. |formulae_type| replace:: \
        dict of \
        {|pole_type|: list of :class:`~triumvirate.winconv.WinConvTerm`}

    .. |pole_type| replace:: str or int or tuple of int

    .. |coeff_type| replace:: int, float or :class:`fractions.Fraction`

    """

    def __init__(self, formulae):

        self.formulae = {
            multipole: [WinConvTerm(*term) for term in formula]
            for multipole, formula in formulae.items()
        }

        self.multipoles = formulae.keys()
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


class TwoPointWinConv:
    pass


class ThreePointWinConv:
    """Window convolution of three-point statistics.

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

        All convolution assumes that the 2-d window function,
        three-point correlation function (3PCF) and the bispectrum
        multipole samples are square matrices, i.e. they are sampled
        at the same sample points in both dimensions.

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

    def conv_bispec(self):
        pass

    def conv_3pcf(self):
        pass

    def conv_3pcf_diag(self):
        pass


def wconv_3pcf_diag(formulae, Q, zeta, r_common, r_Q=None, r_zeta=None):
    r"""Convolve multipoles of the three-point correlation function (3PCF)
    and the window function in diagonal form using specified formulae.

    Parameters
    ----------
    formulae : dict of {str: sequence of (str, str, float)} or 'sugiyama+19'
        Window convolution formulae encoded as a dictionary where for
        each convolved-multipole degree (the key), the value is
        a sequence of tuples each corresponding to a window convolution
        term including the numerical coefficient.

        For example,

        .. code-block::

            formulae = {
                '000': [
                    ('000', '000', 1),
                    ('000', 'ic', -1),
                    ('110', '110', 1./3),
                ],
            }

        computes only the '000' monopole

        .. math::

            \tilde{\zeta}_{000} =
                1 \cdot Q_{000} \zeta_{000}
                + (-1) \cdot Q_{000} \zeta_{\mathrm{ic}}
                + \frac{1}{3} \cdot Q_{110} \zeta_{110} \,.

        If a named formula (e.g. 'sugiyama19') is passed, the window
        convolution formula follows that prescription
        (e.g. :cite:t:`Sugiyama:2019`).

    Q : dict of {str: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    zeta : dict of {str: 2-d array of float}
        Unwindowed 3PCF multipole samples (each key is a multipole).
    r_common : 1-d array of float
        Common separation sample points for both `Q` and `zeta`.
        If `r_Q` or `r_zeta` is provided, all multipole samples are
        re-interpolated over `r_common`; otherwise all multipole
        sample points are assumed to be `r_common`.
    r_Q : 1-d array of float, optional
        Separation sample points for `Q`.  If `None` (default),
        `r_common` is assumed instead.
    r_zeta : 1-d array of float, optional
        Separation sample points for `zeta`.  If `None` (default),
        `r_common` is assumed instead.

    Returns
    -------
    zeta_wconv : dict of {str: :class:`numpy.ndarray`}
        Window-convolved 3PCF samples.

    """
    # Fetch named formulae as a dictionary.
    if isinstance(formulae, str) and formulae.lower() == 'sugiyama+19':
        formulae = _FORMULAE_SUGIYAMA19

    # Convert dictionary to formular object.
    if not isinstance(formulae, WinConvFormulae):
        formulae = WinConvFormulae(formulae)

    # Coerce onto the same coordinate samples.
    if r_Q is None:
        Q_samples = Q
    else:
        r_Q = r_Q if r_Q is not None else r_common
        Q_samples = {
            multipole: interp1d(r_Q, Q[multipole])(r_common)
            for multipole in Q.keys()
        }

    if r_zeta is None:
        zeta_samples = zeta
    else:
        r_zeta = r_zeta if r_zeta is not None else r_common
        zeta_samples = {
            multipole: interp1d(r_zeta, zeta[multipole])(r_common)
            for multipole in zeta.keys()
        }

    # Perform convolution as multiplication.
    zeta_conv = {}
    for multipole in formulae.multipoles:
        zeta_tilde = 0.  # technically not a scalar
        for term in formulae[multipole]:
            zeta_tilde += term.coeff \
                * Q_samples[term.ind_Q] * zeta_samples[term.ind_Z]
        zeta_conv[multipole] = zeta_tilde

    return zeta_conv


def wconv_3pcf(formulae, Q, zeta, r1_common, r2_common,
               r1_Q=None, r2_Q=None, r1_zeta=None, r2_zeta=None):
    r"""Convolve multipoles of the three-point correlation function (3PCF)
    and the window function using specified formulae.

    Parameters
    ----------
    formulae : dict of {str: sequence of (str, str, float)} or 'sugiyama+19'
        Window convolution formulae encoded as a dictionary where for
        each convolved-multipole degree (the key), the value is
        a sequence of tuples each corresponding to a window convolution
        term including the numerical coefficient.

        For example,

        .. code-block::

            formulae = {
                '000': [
                    ('000', '000', 1),
                    ('000', 'ic', -1),
                    ('110', '110', 1./3),
                ],
            }

        computes only the '000' monopole

        .. math::

            \tilde{\zeta}_{000} =
                1 \cdot Q_{000} \zeta_{000}
                + (-1) \cdot Q_{000} \zeta_{\mathrm{ic}}
                + \frac{1}{3} \cdot Q_{110} \zeta_{110} \,.

        If a named formula (e.g. 'sugiyama19') is passed, the window
        convolution formula follows that prescription
        (e.g. :cite:t:`Sugiyama:2019`).

    Q : dict of {str: 2-d array of float}
        Window function multipole samples (each key is a multipole).
    zeta : dict of {str: 2-d array of float}
        Unwindowed 3PCF multipole samples (each key is a multipole).
    r1_common, r2_common : 1-d array of float
        Common separation sample points for both `Q` and `zeta`.
        If `r1_Q`/`r1_zeta` or `r2_Q`/`r2_zeta` is provided,
        all multipole samples are re-interpolated over `r1_common`
        and `r2_common`; otherwise all multipole sample points are assumed
        to be `r1_common` and `r2_common`.
    r1_Q, r2_Q : 1-d array of float, optional
        Separation sample points for `Q`.  If `None` (default),
        `r1_common` and `r2_common` are assumed instead.
    r1_zeta, r2_zeta : 1-d array of float, optional
        Separation sample points for `zeta`.  If `None` (default),
        `r1_common` and `r2_common` are assumed instead.

    Returns
    -------
    zeta_wconv : dict of {str: :class:`numpy.ndarray`}
        Window-convolved 3PCF samples

    """
    # Fetch named formulae as a dictionary.
    if isinstance(formulae, str) and formulae.lower() == 'sugiyama+19':
        formulae = _FORMULAE_SUGIYAMA19

    # Convert dictionary to formular object.
    if not isinstance(formulae, WinConvFormulae):
        formulae = WinConvFormulae(formulae)

    # Coerce onto the same coordinate samples.
    if r1_Q is None and r2_Q is None:
        Q_samples = Q
    else:
        r1_Q = r1_Q if r1_Q is not None else r1_common
        r2_Q = r2_Q if r2_Q is not None else r2_common
        Q_samples = {
            multipole: RectBivariateSpline(
                r1_Q, r2_Q, Q[multipole]
            )(r1_common, r2_common)
            for multipole in Q.keys()
        }

    if r1_zeta is None and r2_zeta is None:
        zeta_samples = zeta
    else:
        r1_zeta = r1_zeta if r1_zeta is not None else r1_common
        r2_zeta = r2_zeta if r2_zeta is not None else r2_common
        zeta_samples = {
            multipole: RectBivariateSpline(
                r1_zeta, r2_zeta, zeta[multipole]
            )(r1_common, r2_common)
            for multipole in zeta.keys()
        }

    # Perform convolution as multiplication.
    zeta_conv = {}
    for multipole in formulae.multipoles:
        zeta_tilde = 0.  # technically not a scalar
        for term in formulae[multipole]:
            zeta_tilde += term.coeff \
                * Q_samples[term.ind_Q] * zeta_samples[term.ind_Z]
        zeta_conv[multipole] = zeta_tilde

    return zeta_conv


def wconv_bispec(formulae, Q, bk, k, k_out, r_common, r_Q=None,
                 transform_kwargs=None):
    r"""Convolve multipoles of the bispectrum and the window function
    via configuration space using specified formulae.

    Parameters
    ----------
    formulae : dict of {str: sequence of (str, str, float)} or 'sugiyama+19'
        See :func:`~triumvirate.winconv.wconv_3pcf`.
    Q : dict of {str: 2-d array of float}
        Configuration-space window function multipole samples
        (each key is a multipole).
    bk : dict of {str: 2-d array of float}
        Unwindowed bispectrum multipole samples (each key is a multipole).
    k : 1-d array of float
        Unwindowed bispectrum wavenumber samples.
    k_out : 1-d array of float
        Output bispectrum wavenumber samples.
    r_common : 1-d array of float
        Common separation sample points for both `Q` and 3PCF :math:`zeta`
        (transformed from `bk`).  If `r_Q` is provided, all multipole
        samples are re-interpolated over `r_common`; otherwise all
        multipole sample points are assumed to be `r_common`.
    r_Q : 1-d array of float, optional
        Separation sample points for `Q`.  If `None` (default),
        `r_common` is assumed instead.
    transform_kwargs : dict, optional
        Keyword arguments to be passed to
        :class:`~triumvirate.bihankel.transform_bispec_to_3pcf` and
        :class:`~triumvirate.bihankel.transform_3pcf_to_bispec` for
        transforming to/from configuration space.  If `None` (default),
        the FFTLog sample number is set to 1024 and no extrapolation
        is performed.

    Returns
    -------
    bk_wconv : dict of {str: :class:`numpy.ndarray`}
        Window-convolved bispectrum samples.


    .. attention::

        This function assumes the wavenumber samples :math:`k_{1,2}` in
        the bispectrum :math:`B_{\ell_1 \ell_2 L}(k_1, k_2)` are the same.

    """
    if transform_kwargs is None:
        transform_kwargs = {}
    n_fftlog = transform_kwargs.get('n_fftlog', 1024)

    zeta = {}
    for multipole, bk_in in bk.items():
        ell1, ell2 = int(multipole[0]), int(multipole[1])
        zeta[multipole] = transform_bispec_to_3pcf(
            ell1, ell2, k, bk_in, r_common, n_fftlog,
            **transform_kwargs
        )['zeta']

    zeta_wconv = wconv_3pcf(
        formulae, Q, zeta, r_common, r_common, r1_Q=r_Q, r2_Q=r_Q
    )

    bk_wconv = {}
    for multipole, zeta_in in zeta_wconv.items():
        ell1, ell2 = int(multipole[0]), int(multipole[1])
        bk_wconv[multipole] = transform_3pcf_to_bispec(
            ell1, ell2, r_common, zeta_in, k_out, n_fftlog
        )

    return bk_wconv
