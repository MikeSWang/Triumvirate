"""
Array Operations (:mod:`~triumvirate._arrayops`)
==========================================================================

Perform array operations.

.. autosummary::
    ShapeError

1-d arrays:

.. autosummary::
    extrap_lin
    extrap_loglin
    extrap_pad

2-d arrays:

.. autosummary:
    extrap2d_lin
    extrap2d_loglin
    extrap2d_pad

"""
from collections.abc import Sequence

import numpy as np


class ShapeError(ValueError):
    """Exception raised when the array shape is incorrect.

    """
    pass


# ========================================================================
# Extrapolation
# ========================================================================

def _check_1d_array(a, check_sorted=False, check_sign=False,
                    check_lin=False, check_loglin=False):
    """Check the input array is a 1-d array, possibly satisfying
    additional conditions.

    Parameters
    ----------
    a : array of float
        Input array.
    check_sorted : bool, optional
        If `True` (default is `False`), check the input array is sorted.
    check_sign : bool, optional
        If `True` (default is `False`), check for mixed-sign or zero
        entries in the input array.
    check_lin : bool, optional
        If `True` (default is `False`), check the input array is
        linearly spaced.
    check_loglin : bool, optional
        If `True` (default is `False`), check the input array is
        log-linearly spaced (and all entries have the same sign).

    Returns
    -------
    a : 1-d array of float
        Squeezed 1-d array.

    Raises
    ------
    ShapeError
        When the input array is not equivalently 1-d.
    ValueError
        When the input array fails any of the checks specified.

    """
    a = np.squeeze(a)

    if a.ndim != 1:
        raise ShapeError("Input array is not 1-d.")

    def _check_sign(a):
        same_sgn = np.all(a > 0. if a[0] > 0. else a < 0.)
        if not same_sgn:
            raise ValueError("Input array contains mixed-sign entries.")

    if check_sorted:
        if (not np.all(a[:-1] <= a[1:])) or (not np.all(a[:-1] >= a[1:])):
            raise ValueError("Input array is not a sorted sequence.")
    if check_sign:
        _check_sign(a)
    if check_lin:
        spacing = np.mean(a[1:] - a[:-1])
        if not np.allclose(a[1:] - a[:-1], spacing):
            raise ValueError("Input array is not linearly spaced.")
    if check_loglin:
        _check_sign(a)
        log_spacing = np.mean(np.log(a[1:] / a[:-1]))
        if not np.allclose(np.log(a[1:] / a[:-1]), log_spacing):
            raise ValueError("Input array is not log-linearly spaced.")

    return a


def _check_2d_array(a, check_sorted=False, check_sign=False,
                    check_lin=False, check_loglin=False):
    """Check the input array is a sorted 2-d array, possibly satisfying
    additional conditions.

    Parameters
    ----------
    a : array of float
        Input array.
    check_sorted : bool, optional
        If `True` (default is `False`), check the input array is sorted.
    check_sign : bool, optional
        If `True` (default is `False`), check for mixed-sign or zero
        entries in the input array.
    check_lin : bool, optional
        If `True` (default is `False`), check the input array is
        linearly spaced.
    check_loglin : bool, optional
        If `True` (default is `False`), check the input array is
        log-linearly spaced (and all entries have the same sign).

    Returns
    -------
    a : 2-d array of float
        Squeezed 2-d array.

    Raises
    ------
    ShapeError
        When the input array is not equivalently 2-d.
    ValueError
        When the input array fails any of the checks specified.

    """
    a = np.squeeze(a)

    if a.ndim != 2:
        raise ShapeError("Input array is not 2-d.")

    def _check_sign(a):
        same_sgn = np.all(a > 0. if a[0, 0] > 0. else a < 0.)
        if not same_sgn:
            raise ValueError("Input array contains mixed-sign entries.")

    if check_sorted:
        for a_, direction in zip([a, a.transpose()], ['row', 'column']):
            if (
                (
                    not np.all([
                        np.all(a__row[:-1] <= a__row[1:]) for a__row in a_
                    ])
                ) or (
                    not np.all([
                        np.all(a__row[:-1] >= a__row[1:]) for a__row in a_
                    ])
                )
            ):
                raise ValueError(
                    "Input array is not a sorted sequence "
                    f"along the {direction}s."
                )
    if check_sign:
        _check_sign(a)
    if check_lin:
        for a_, direction in zip([a, a.transpose()], ['row', 'column']):
            spacing = np.mean(a_[:, 1:] - a_[:, :-1], axis=1)
            if not np.allclose(a_[:, 1:] - a_[:, :-1], spacing[:, None]):
                raise ValueError(
                    "Input array is not linearly spaced "
                    f"along the {direction}s."
                )
    if check_loglin:
        _check_sign(a)
        for a_, direction in zip([a, a.transpose()], ['row', 'column']):
            log_spacing = np.mean(np.log(a_[:, 1:] / a_[:, :-1]), axis=1)
            if not np.allclose(
                np.log(a_[:, 1:] / a_[:, :-1]), log_spacing[:, None]
            ):
                raise ValueError(
                    "Input array is not log-linearly spaced "
                    f"along the {direction}s."
                )

    return a


def extrap_lin(a, n_ext):
    """Extrapolate a 1-d array linearly.

    The linear spacing used is the distance between the two points
    at either end separately.

    Parameters
    ----------
    a : array of float
        Input 1-d array.
    n_ext : int
        Number of extra elements on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a = _check_1d_array(a)

    a_l = a[0] + (a[1] - a[0]) * np.arange(-n_ext, 0)
    a_r = a[-1] + (a[-1] - a[-2]) * np.arange(1, n_ext + 1)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def extrap_loglin(a, n_ext):
    """Extrapolate a 1-d array log-linearly.

    The log-linear spacing used is the log-distance between the two points
    at either end separately.

    Parameters
    ----------
    a : array of float
        Input 1-d array.
    n_ext : int
        Number of extra elements on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a = _check_1d_array(a, check_sign=True)

    a_l = a[0] * (a[1] / a[0]) ** np.arange(-n_ext, 0)
    a_r = a[-1] * (a[-1] / a[-2]) ** np.arange(1, n_ext + 1)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def extrap_pad(a, n_ext, c_lower, c_upper):
    """Extrapolate a 1-d array by constant padding.

    Parameters
    ----------
    a : array of float
        Input 1-d array.
    n_ext : int
        Number of extra elements on either side.
    c_lower, c_upper : float
        Constant value used for padding on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a = _check_1d_array(a)

    a_l = np.full(n_ext, c_lower)
    a_r = np.full(n_ext, c_upper)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def _extrap2d_row_lin(a, ncol_ext):
    """Extrapolate a 2-d array linearly along each row
    (i.e. horizontally).

    Parameters
    ----------
    a : 2-d :class:`numpy.ndarray`
        Input 2-d array.
    ncol_ext : int
        Number of extra columns on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    """
    a_l = a[:, [0]] + np.arange(-ncol_ext, 0) * (a[:, [1]] - a[:, [0]])
    a_r = a[:, [-1]] + np.arange(1, ncol_ext + 1) * (a[:, [-1]] - a[:, [-2]])

    a_out = np.c_[a_l, a, a_r]

    return a_out


def _extrap2d_row_loglin(a, ncol_ext):
    """Extrapolate a 2-d array log-linearly along each row
    (i.e. horizontally).

    Parameters
    ----------
    a : 2-d :class:`numpy.ndarray`
        Input 2-d array.
    ncol_ext : int
        Number of extra columns on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    """
    a_l = a[:, [0]] * (a[:, [1]] / a[:, [0]]) ** np.arange(-ncol_ext, 0)
    a_r = a[:, [-1]] * (a[:, [-1]] / a[:, [-2]]) ** np.arange(1, ncol_ext + 1)

    a_out = np.c_[a_l, a, a_r]

    return a_out


def extrap2d_lin(a, n_ext):
    """Extrapolate a 2-d array bilinearly.

    Parameters
    ----------
    a : array of float
        Input 2-d array.
    n_ext : int or tuple of int
        Number of extra elements on either side of each row and/or column.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    """
    a = _check_2d_array(a)

    nrow_ext, ncol_ext = \
        n_ext if isinstance(n_ext, Sequence) else (n_ext, n_ext)

    # Extrapolate horizontally then vertically.
    a_ = _extrap2d_row_lin(a, ncol_ext)
    a_out = _extrap2d_row_lin(a_.transpose(), nrow_ext).transpose()

    return a_out


def extrap2d_loglin(a, n_ext):
    """Extrapolate a 2-d array log-bilinearly.

    Parameters
    ----------
    a : array of float
        Input 2-d array.
    n_ext : int or tuple of int
        Number of extra elements on either side of each row and/or column.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    """
    a = _check_2d_array(a, check_sign=True)

    nrow_ext, ncol_ext = \
        n_ext if isinstance(n_ext, Sequence) else (n_ext, n_ext)

    # Extrapolate horizontally then vertically.
    a_ = _extrap2d_row_loglin(a, ncol_ext)
    a_out = _extrap2d_row_loglin(a_.transpose(), nrow_ext).transpose()

    return a_out


def extrap2d_pad(a, n_ext, c_lower, c_upper):
    """Extrapolate a 2-d array by constant padding bi-directionally.

    Parameters
    ----------
    a : array of float
        Input 2-d array.
    n_ext : int or tuple of int
        Number of extra elements on either side of each row and/or column.
    c_lower, c_upper : float or tuple of float
        Constant value used for padding on either side of
        each row and/or column.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Notes
    -----
    The padded array will have the following structure schematically:

    .. code-block:: none

        c_lower[1]  c_lower[1]  ...  c_lower[1]   c_lower[1]
        c_lower[0]  a[0, 0]     ...  a[0, -1]     c_upper[0]
        ...         ...         ...  ...          ...
        c_lower[0]  a[-1, 0]    ...  a[-1, -1]    c_upper[0]
        c_upper[1]  c_upper[1]  ...  c_upper[-1]  c_upper[1]

    """
    a = _check_2d_array(a)

    nrow_ext, ncol_ext = \
        n_ext if isinstance(n_ext, Sequence) else (n_ext, n_ext)
    crow_lower, ccol_lower = \
        c_lower if isinstance(c_lower, Sequence) else (c_lower, c_lower)
    crow_upper, ccol_upper = \
        c_upper if isinstance(c_upper, Sequence) else (c_upper, c_upper)

    a_l = np.full((a.shape[0], ncol_ext), crow_lower)
    a_r = np.full((a.shape[0], ncol_ext), crow_upper)
    a_ = np.c_[a_l, a, a_r]

    a__u = np.full((nrow_ext, a_.shape[-1]), ccol_lower)
    a__d = np.full((nrow_ext, a_.shape[-1]), ccol_upper)
    a_out = np.r_[a__u, a_, a__d]

    return a_out
