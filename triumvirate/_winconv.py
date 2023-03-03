"""
Window Convolution (:mod:`~triumvirate.winconv`)
==========================================================================

Perform window convolution of three-point statistics.

"""
from typing import NamedTuple

import numpy as np
from scipy.interpolate import interp1d, RectBivariateSpline

from triumvirate._bihankel import (
    transform_bispec_to_3pcf,
    transform_3pcf_to_bispec,
)


# Default convolution formula from Sugiyama et al. (2019) [arXiv:1803.02132].
_FORMULAE_SUGIYAMA19 = {
    '000': [
        ('000', '000', 1),
        ('000', 'ic', -1),
        ('110', '110', 1./3),
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
        ('110', '112', 1./3),
        ('112', '110', 1./3),
    ],
    '112': [
        ('000', '112', 1),
        ('112', '000', 1),
        ('112', 'ic', -1),
        ('022', '110', 2./5),
        ('202', '110', 2./5),
        ('110', '022', 2./5),
        ('110', '202', 2./5),
    ],
}


class _WConvTerm(NamedTuple):
    r"""Window convolution term as a tuple of its three factors.

    For example, ``term = _WConvTerm('000', 'ic', -1)``, with
    parameters/attributes/factors ``term.deg_Q = '000'``,
    ``term.deg_zeta = 'ic'`` and ``term.coeff = -1``, corresponds to the
    term :math:`- 1 \cdot Q_{000} \zeta_{\mathrm{ic}}`.

    Attributes
    ----------
    deg_Q : str
        Window function multipole degree.
    deg_zeta : str
        Unwindowed 3-point correlation funciton multipole degree.
    coeff : float
        Numerical factor for the convolution term.

    """
    deg_Q: str
    deg_zeta: str
    coeff: float


class _WConvFormulae:
    """Window convolution formulae.

    The full formulae are encoded as a dictionary, where each key
    corresponds to a windowed 3-point correlation function (3PCF)
    multipole, and each value is a sequence of (named) tuples each
    corresponding to a window convolution term.

    Parameters
    ----------
    formulae : dict of {str: sequence of tuples of (str, str, float)}
        Window convolution formulae.

    Attributes
    ----------
    multipoles_Q : list of str
        List of degrees of the required window function multipoles.
    multipoles_zeta : list of str
        List of degrees of the required 3PCF multipoles.

    """

    def __init__(self, formulae):

        self._formulae = {
            multipole: [_WConvTerm(*term) for term in formula]
            for multipole, formula in formulae.items()
        }

        self.multipoles = formulae.keys()

        self.multipoles_Q = np.unique([
            term.deg_Q
            for formula in self._formulae.values()
            for term in formula
        ])
        self.multipoles_zeta = np.unique([
            term.deg_zeta
            for formula in self._formulae.values()
            for term in formula
        ])

    def __getitem__(self, multipole):
        """Get the formula for the window convolution of a specific
        3PCF multipole.

        Parameters
        ----------
        multipole : str
            Degrees of the specified windowed 3PCF multipole.

        Returns
        -------
        list of :class:`~triumvirate.winconv._WConvTerm`
            List of window convolution terms for `multipole`.

        """
        return self._formulae[multipole]


def wconv_3pcf_diag(formulae, Q, zeta, r_common, r_Q=None, r_zeta=None):
    r"""Convolve multipoles of the 3-point correlation function (3PCF) and
    the window function in diagonal form using specified formulae.

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
    if not isinstance(formulae, _WConvFormulae):
        formulae = _WConvFormulae(formulae)

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
                * Q_samples[term.deg_Q] * zeta_samples[term.deg_zeta]
        zeta_conv[multipole] = zeta_tilde

    return zeta_conv


def wconv_3pcf(formulae, Q, zeta, r1_common, r2_common,
               r1_Q=None, r2_Q=None, r1_zeta=None, r2_zeta=None):
    r"""Convolve multipoles of the 3-point correlation function (3PCF) and
    the window function using specified formulae.

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
    if not isinstance(formulae, _WConvFormulae):
        formulae = _WConvFormulae(formulae)

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
                * Q_samples[term.deg_Q] * zeta_samples[term.deg_zeta]
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
