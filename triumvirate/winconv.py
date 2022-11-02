"""
Window Convolution (:mod:`~triumvirate.winconv`)
==========================================================================

Perform window convolution of three-point statistics.

"""
from typing import NamedTuple

import numpy as np
from scipy.interpolate import RectBivariateSpline

from triumvirate.bihankel import (
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
    """Window convolution term as a tuple of its three factors.

    For example, `term = _WConvTerm('000', 'ic', -1)`, with
    parameters/attributes/factors ``term.deg_Q = '000'``,
    ``term.deg_zeta = 'ic'` and ``term.coeff = -1`` corresponds to the
    term :math:`- Q_{000} \zeta_{\mathrm{ic}}`.

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

    def __getitem__(self, deg):
        """Get the formula for the window convolution of a specific
        3PCF multipole.

        Parameters
        ----------
        deg : str
            Degree of the specified windowed 3PCF multipole.

        Returns
        -------
        _type_
            _description_
        """
        return self._formulae[deg]


def wconv_3pcf(formulae, zeta, Q, r1_com, r2_com,
               r1_zeta=None, r2_zeta=None, r1_Q=None, r2_Q=None):
    """Convolve multipoles of the 3-point correlation function (3PCF) and
    the window function using specified formulae.

    Parameters
    ----------
    formulae : dict of {str: sequence of tuples of (str, str, float)} or 'sugiyama+19'
        Window convolution formulae encoded as a dictionary where for
        each convolved-multipole degree key, the value is a sequence of
        tuples each corresponding to a window convolution term's
        multipole factors and the numerical coefficient.

        Example:

        .. code-block::

            formulae = {
                '000': [
                    ('000', '000', 1),
                    ('000', 'ic', -1),
                    ('110', '110', 1./3),
                ],
            }

        computes only the multipole

        .. math::

            \tilde{\zeta}_{000} =
                Q_{000} (\zeta_{000} - \zeta_{\mathrm{ic}})
                + \frac{1}{3} Q_{110} \zeta_{110} \,.

        If a string (namely 'sugiyama19') is passed, the window
        convolution formula from Sugiyama et al. (2019)
        [`arXiv:1803.02132 <Sugiyama+19>`_].

    zeta : dict of {str: 2-d array of float}
        Unwindowed 3PCF multipole samples.
    Q : dict of {str: 2-d array of float}
        Window function multipole samples.
    r1_com, r2_com : 1-d array of float
        Common separation sample points for both `zeta` and `Q`.
        If `r1_zeta` and `r1_Q` (or `r2_zeta` and `r2_Q`) are provided,
        all multipole samples are interpolated over `r1_com` and `r2_com`;
        otherwise all multipole samples are assumed to correspond to
        `r1_com` and `r2_com`.
    r1_zeta, r2_zeta : 1-d array of float, optional
        Separation sample points for `zeta`.  If `None` (default),
        `r1_com` and `r2_com` are used instead.
    r1_Q, r2_Q : 1-d array of float, optional
        Separation sample points for `Q`.  If `None` (default),
        `r1_com` and `r2_com` are used instead.

    Returns
    -------
    zeta_wconv : dict of {str: :class:`numpy.ndarray`}
        Window-convolved 3PCF samples.

    .. _Sugiyama+19: https://arxiv.org/abs/1803.02132

    """
    if isinstance(formulae, str) and formulae.lower() == 'sugiyama+19':
        formulae = _FORMULAE_SUGIYAMA19
    if not isinstance(formulae, _WConvFormulae):
        formulae = _WConvFormulae(formulae)

    if r1_zeta is None and r2_zeta is None:
        zeta_samples = zeta
    else:
        r1_zeta = r1_zeta if r1_zeta is not None else r1_com
        r2_zeta = r2_zeta if r2_zeta is not None else r2_com
        zeta_samples = {
            deg: RectBivariateSpline(
                r1_zeta, r2_zeta, zeta[deg]
            )(r1_com, r2_com)
            for deg in zeta.keys()
        }

    if r1_Q is None and r2_Q is None:
        Q_samples = Q
    else:
        r1_Q = r1_Q if r1_Q is not None else r1_com
        r2_Q = r2_Q if r2_Q is not None else r2_com
        Q_samples = {
            deg: RectBivariateSpline(r1_Q, r2_Q, Q[deg])(r1_com, r2_com)
            for deg in Q.keys()
        }

    zeta_conv = {}
    for multipole in formulae.multipoles:
        zeta_tilde = 0.
        for term in formulae[multipole]:
            zeta_tilde += term.coeff \
                * Q_samples[term.deg_Q] * zeta_samples[term.deg_zeta]
        zeta_conv[multipole] = zeta_tilde

    return zeta_conv


def wconv_bispec(formulae, bk, k_in, k_out, Q, r_com, r1_Q=None, r2_Q=None,
                 transform_kwargs=None):
    """Convolve multipoles of the bispectrum and the window function
    via configuration space using specified formulae.

    Parameters
    ----------
    formulae : dict of {str: sequence of tuples of (str, str, float)} or 'sugiyama+19'
        See :func:`~winconv.wconv.wconv_3pcf`.
    bk : dict of {str: 2-d array of float}
        Unwindowed input bispectrum multipoles.
    k_in : 1-d array of float
        Input wavenumbers for which `bk` is evaluated.
    k_out : 1-d array of float
        Output wavenumbers for which the window-convolved bispectrum
        is evaluated.
    Q : dict of {str: 2-d array of float}
        Window function multipoles.
    r_com : 1-d array of float
        Common separation sample points for both :math:`zeta` (transformed
        from `bk`) and `Q`.  If `r1_Q` (or `r2_Q`) is provided, all
        multipole samples are interpolated over `r_com`; otherwise all
        multipole samples are assumed to correspond to `r_com`.
    r1_Q, r2_Q : 1-d array of float, optional
        See :func:`~winconv.wconv.wconv_3pcf`.
    transform_kwargs : dict, optional
        Keyword arguments to be passed to
        :class:`triumvirate.bihankel.transform_bispec_to_3pcf` and
        :class:`triumvirate.bihankel.transform_3pcf_to_bispec` for
        transforming to/from configuration space.  If `None` (default),
        the sample number for FFTLog transform is set to 1024 and no
        extrapolation is performed.

    Returns
    -------
    bk_wconv : dict of {str: :class:`numpy.ndarray`}
        Window-convolved bispectrum samples.

    """
    if transform_kwargs is None:
        transform_kwargs = {}
    n_fftlog = transform_kwargs.get('n_fftlog', 1024)

    zeta = {}
    for deg, bk_in in bk.items():
        zeta[deg] = transform_bispec_to_3pcf(
            int(deg[0]), int(deg[1]), bk_in, k_in, r_com, n_fftlog, **transform_kwargs
        )['zeta']

    zeta_wconv = wconv_3pcf(
        formulae, zeta, Q, r_com, r_com, r1_Q=r1_Q, r2_Q=r2_Q
    )

    bk_wconv = {}
    for deg, zeta_in in zeta_wconv.items():
        bk_wconv[deg] = transform_3pcf_to_bispec(
            int(deg[0]), int(deg[1]), zeta_in, r_com, k_out, n_fftlog
        )

    return bk_wconv
