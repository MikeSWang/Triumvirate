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


class OrderError(ValueError):
    """Exception raised when the array order is incorrect.

    """
    pass


class MixedSignError(ValueError):
    """Exception raised when the array signs are mixed.

    """
    pass


class SpacingError(ValueError):
    """Exception raised when the array spacing is incorrect.

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
            raise MixedSignError("Input array contains mixed-sign entries.")

    if check_sorted:
        if (not np.all(a[:-1] <= a[1:])) or (not np.all(a[:-1] >= a[1:])):
            raise OrderError("Input array is not a sorted sequence.")
    if check_sign:
        _check_sign(a)
    if check_lin:
        spacing = np.mean(a[1:] - a[:-1])
        if not np.allclose(a[1:] - a[:-1], spacing):
            raise SpacingError("Input array is not linearly spaced.")
    if check_loglin:
        _check_sign(a)
        log_spacing = np.mean(np.log(a[1:] / a[:-1]))
        if not np.allclose(np.log(a[1:] / a[:-1]), log_spacing):
            raise SpacingError("Input array is not log-linearly spaced.")

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
            raise MixedSignError("Input array contains mixed-sign entries.")

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
                raise OrderError(
                    "Input array is not a sorted sequence "
                    f"along the {direction}s."
                )
    if check_sign:
        _check_sign(a)
    if check_lin:
        for a_, direction in zip([a, a.transpose()], ['row', 'column']):
            spacing = np.mean(a_[:, 1:] - a_[:, :-1], axis=1)
            if not np.allclose(a_[:, 1:] - a_[:, :-1], spacing[:, None]):
                raise SpacingError(
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
                raise SpacingError(
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


def extrap2d_pad(a, n_ext, crow_lower, crow_upper, ccol_lower, ccol_upper):
    """Extrapolate a 2-d array by constant padding bi-directionally.

    Parameters
    ----------
    a : array of float
        Input 2-d array.
    n_ext : int or tuple of int
        Number of extra elements on either side of each row and/or column.
    crow_lower, crow_upper, ccol_lower, ccol_upper : (sequence of) float
        Constant values used for padding on either side of
        each row and/or column.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Notes
    -----
    The padded array will have the following structure schematically:

    .. code-block:: none

        ... ...             ...            ...  ...             ...             ...
        ... ...             ccol_lower[0]  ...  ccol_lower[-1]  ...             ...
        ... crow_lower[0]   a[0, 0]        ...  a[0, -1]        crow_upper[0]   ...
        ... ...             ...            ...  ...             ...             ...
        ... crow_lower[-1]  a[-1, 0]       ...  a[-1, -1]       crow_upper[-1]  ...
        ... ...             ccol_upper[0]  ...  ccol_upper[-1]  ...             ...
        ... ...             ...            ...  ...             ...             ...

    If any of `crow_lower`, `crow_upper`, `ccol_lower` and `ccol_upper`
    is a scalar, it is repeated to form a 1-d array of the
    appropriate length.

    For consistency, it is required that
    ``crow_lower[0] == ccol_lower[0]``,
    ``crow_upper[0] == ccol_lower[-1]``,
    ``crow_lower[-1] == ccol_upper[0]``, and
    ``crow_upper[-1] == ccol_upper[-1]``.

    """  # noqa: E501
    a = _check_2d_array(a)

    if np.isscalar(crow_lower):
        crow_lower = np.full(a.shape[0], crow_lower)
    if np.isscalar(crow_upper):
        crow_upper = np.full(a.shape[0], crow_upper)
    if np.isscalar(ccol_lower):
        ccol_lower = np.full(a.shape[-1], ccol_lower)
    if np.isscalar(ccol_upper):
        ccol_upper = np.full(a.shape[-1], ccol_upper)

    if crow_lower[0] != ccol_lower[0]:
        raise ValueError(
            "Inconsistent lower padding values for rows and columns."
        )
    if crow_upper[0] != ccol_lower[-1]:
        raise ValueError(
            "Inconsistent lower and upper padding values for rows and columns."
        )
    if crow_lower[-1] != ccol_upper[0]:
        raise ValueError(
            "Inconsistent upper and lower padding values for rows and columns."
        )
    if crow_upper[-1] != ccol_upper[-1]:
        raise ValueError(
            "Inconsistent upper padding values for rows and columns."
        )

    nrow_ext, ncol_ext = \
        n_ext if isinstance(n_ext, Sequence) else (n_ext, n_ext)

    a_l = np.tile(crow_lower, (nrow_ext, 1)).T
    a_r = np.tile(crow_upper, (nrow_ext, 1)).T
    a_ = np.c_[a_l, a, a_r]

    ccol_lower = np.concatenate([a_l[0], ccol_lower, a_r[0]])
    ccol_upper = np.concatenate([a_l[-1], ccol_upper, a_r[-1]])

    a__u = np.tile(ccol_lower, (ncol_ext, 1))
    a__d = np.tile(ccol_upper, (ncol_ext, 1))
    a_out = np.r_[a__u, a_, a__d]

    return a_out
