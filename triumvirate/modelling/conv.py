"""
Window Convolution (:mod:`~triumvirate.modelling.conv`)
==========================================================================

Perform window convolution of (wide-angle corrected) bispectrum models.

"""
import numpy as np
from scipy.interpolate import RectBivariateSpline as Spline2d


def _interp_2d_grid(r1_sample, r2_sample, f_samples, r1_grid, r2_grid):
    r"""Interpolate samples of a bivariate function :math:`f(r_1, r_2)`
    over a 2-dimensional regular grid.

    Parameters
    ----------
    r1_sample : array of float
        Sample points corresponding to :math:`r_1`.
    r2_sample : array of float
        Sample points corresponding to :math:`r_2`.
    f_samples : array of float
        2-d function samples with each row a value of :math:`r_1` and
        each column a value of :math:`r_2`.
    r1_grid : array of float
        Grid points corresponding to :math:`r_1` over which the function
        is interpolated.
    r2_grid : array of float
        Grid points corresponding to :math:`r_2` over which the function
        is interpolated.

    Returns
    -------
    f_grid : array of float
        2-d function values interpolated over the grid with each row
        a value of :math:`r_1` and each column a value of :math:`r_2`.

    """
    f = Spline2d(r1_sample, r2_sample, f_samples)

    f_grid = f(r1_grid, r2_grid)

    return f_grid


def winconv_3pcf(order, xi_samples, r1_xi, r2_xi, Q_samples, r1_Q, r2_Q,
                 xi_ic, r1_grid=None, r2_grid=None):
    """Convolve multipoles of the 3-point correlation function with
    those of the window function.

    Currently this is only implemented for :math:`(\ell_1, \ell_2, L) \in
    {(0, 0, 0), (1, 1, 0), (1, 1, 2), (2, 0, 2)}`.

    Parameters
    ----------
    order : tuple of (int, int, int)
        Multipole order.
    xi_samples : dict of {tuple: :class:`numpy.ndarray`}
        Samples of the unwindowed 3-point correlation function.
    r1_xi, r2_xi : :class:`numpy.ndarray`
        Sample points corresponding to :math:`r_1` and :math:`r_2`
        for `xi_samples`.
    Q_samples : dict of {tuple: :class:`numpy.ndarray`}
        Samples of the window function multipoles.
    r1_Q, r2_Q : :class:`numpy.ndarray`
        Sample points corresponding to :math:`r_1` and :math:`r_2`
        for `Q_samples`.
    xi_ic : float
        Integral constraint value.
    r1_grid, r2_grid : :class:`numpy.ndarray`, optional
        Grid points corresponding to :math:`r_1` and :math:`r_2`
        over which the functions may be interpolated before convolution.
        If either is `None` (default), `xi_samples` and `Q_samples` are
        assumed to have the same sample points.

    Returns
    -------
    xi_tilde_samples : :class:`numpy.ndarray`
        Samples of the windowed 3-point correlation function.

    """
    if r1_grid is None:
        assert np.allclose(r1_xi, r1_Q), (
            "`r1_xi` and `r1_Q` are not close. "
            "Please provide `r1_grid` and `r2_grid`."
        )
    if r2_grid is None:
        assert np.allclose(r2_xi, r2_Q), (
            "`r2_xi` and `r2_Q` are not close. "
            "Please provide `r1_grid` and `r2_grid`."
        )

    if r1_grid is not None and r2_grid is not None:
        xi_samples_ = {}
        for ell_order in xi_samples.keys():
            xi_samples_[ell_order] = _interp_2d_grid(
                r1_xi, r2_xi, xi_samples[ell_order], r1_grid, r2_grid
            )

        Q_samples_ = {}
        for ell_order in Q_samples.keys():
            Q_samples_[ell_order] = _interp_2d_grid(
                r1_Q, r2_Q, Q_samples[ell_order], r1_grid, r2_grid
            )
    else:
        xi_samples_ = xi_samples
        Q_samples_ = Q_samples

    if order == (0, 0, 0):
        xi_tilde_samples = (
            Q_samples_[order] * (xi_samples_[order] - xi_ic)
            + 1./3 * Q_samples_[(1, 1, 0)] * xi_samples_[(1, 1, 0)]
        )
    elif order == (1, 1, 0):
        xi_tilde_samples = (
            Q_samples_[(0, 0, 0)] * xi_samples_[order]
            + Q_samples_[order] * (xi_samples_[(0, 0, 0)] - xi_ic)
        )
    elif order == (1, 1, 2):
        xi_tilde_samples = (
            Q_samples_[(0, 0, 0)] * xi_samples_[order]
            + Q_samples_[order] * (xi_samples_[(0, 0, 0)] - xi_ic)
            + 2./5 * Q_samples_[(1, 1, 0)]
                * (xi_samples_[(0, 2, 2)] + xi_samples_[(2, 0, 2)])
            + 2./5 * (Q_samples_[(0, 2, 2)] + Q_samples_[(2, 0, 2)])
                * xi_samples_[(1, 1, 0)]
        )
    elif order == (2, 0, 2):
        xi_tilde_samples = (
            Q_samples_[(0, 0, 0)] * xi_samples_[order]
            + Q_samples_[order] * (xi_samples_[(0, 0, 0)] - xi_ic)
            + 1./3 * Q_samples_[(1, 1, 0)] * xi_samples_[(1, 1, 2)]
            + 1./3 * Q_samples_[(1, 1, 2)] * xi_samples_[(1, 1, 0)]
        )
    else:
        raise RuntimeError(
            f"Input `order` {order} not supported by this function yet."
        )

    return xi_tilde_samples
