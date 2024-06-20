"""
Array Operations (:mod:`~triumvirate._arrayops`)
==========================================================================

Perform array operations.

.. autosummary::
    ShapeError
    OrderError
    MixedSignError
    SpacingError
    DivergenceWarning

1-d arrays:

.. autosummary::
    extrap_pad
    extrap_lin
    extrap_loglin
    extrap_loglin_oscil

2-d arrays:

.. autosummary:
    extrap2d_lin
    extrap2d_loglin
    extrap2d_pad
    extrap2d_loglin_oscil
    reshape_threept_datatab

"""
import warnings
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


class DivergenceWarning(RuntimeWarning):
    """Warning issued when divergent behaviour is detected in an array.

    """
    pass


# ========================================================================
# Extrapolation
# ========================================================================

def _check_1d_array(a, check_sorted=False, check_epsign=False,
                    check_lin=False, check_loglin=False):
    """Check the input array is a 1-d array, possibly satisfying
    additional conditions.

    Parameters
    ----------
    a : array of float
        Input array.
    check_sorted : bool, optional
        If `True` (default is `False`), check the input array is sorted.
    check_epsign : bool, optional
        If `True` (default is `False`), check for mixed sign entries at
        the endpoints in the input array, and for nonzero endpoints with
        zero next-to-endpoint entries.
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

    if a.ndim != 1 or a.size < 2:
        raise ShapeError("Input array is not 1-d.")

    def _check_epsign(a, strict=False):
        if strict:
            same_sign_l = a[0] * a[1] > 0.
            same_sign_r = a[-1] * a[-2] > 0.
        else:
            same_sign_l = (a[0] * a[1] > 0.) or (a[0] == 0.)
            same_sign_r = (a[-1] * a[-2] > 0.) or (a[-1] == 0.)
        if not (same_sign_l and same_sign_r):
            raise MixedSignError(
                "Input array contains mixed-sign entries at the endpoints."
            )

    if check_sorted:
        if (not np.all(a[:-1] <= a[1:])) or (not np.all(a[:-1] >= a[1:])):
            raise OrderError("Input array is not a sorted sequence.")
    if check_epsign:
        _check_epsign(a)
    if check_lin:
        spacing = np.mean(a[1:] - a[:-1])
        if not np.allclose(a[1:] - a[:-1], spacing):
            raise SpacingError("Input array is not linearly spaced.")
    if check_loglin:
        _check_epsign(a, strict=True)
        log_spacing = np.mean(np.log(a[1:] / a[:-1]))
        if not np.allclose(np.log(a[1:] / a[:-1]), log_spacing):
            raise SpacingError("Input array is not log-linearly spaced.")

    return a


def _check_2d_array(a, check_sorted=False, check_epsign=False,
                    check_lin=False, check_loglin=False):
    """Check the input array is a sorted 2-d array, possibly satisfying
    additional conditions.

    Parameters
    ----------
    a : array of float
        Input array.
    check_sorted : bool, optional
        If `True` (default is `False`), check the input array is sorted.
    check_epsign : bool, optional
        If `True` (default is `False`), check for mixed sign entries at
        the endpoints in the input array, and for nonzero endpoints with
        zero next-to-endpoint entries.
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

    if a.ndim != 2 or a.size < 4:
        raise ShapeError("Input array is not 2-d.")

    def _check_epsign(a, strict=False):
        if strict:
            same_sign_l = np.all(a[:, 0] * a[:, 1] > 0.)
            same_sign_r = np.all(a[:, -1] * a[:, -2] > 0.)
            same_sign_u = np.all(a[0, :] * a[1, :] > 0.)
            same_sign_d = np.all(a[-1, :] * a[-2, :] > 0.)
        else:
            same_sign_l = np.all(np.logical_or(
                a[:, 0] * a[:, 1] > 0., a[:, 0] == 0.
            ))
            same_sign_r = np.all(np.logical_or(
                a[:, -1] * a[:, -2] > 0., a[:, -1] == 0.
            ))
            same_sign_u = np.all(np.logical_or(
                a[0, :] * a[1, :] > 0., a[0, :] == 0.
            ))
            same_sign_d = np.all(np.logical_or(
                a[-1, :] * a[-2, :] > 0., a[-1, :] == 0.
            ))

        if not (same_sign_l and same_sign_r and same_sign_u and same_sign_d):
            raise MixedSignError(
                "Input array contains mixed-sign entries at the endpoints."
            )

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
    if check_epsign:
        _check_epsign(a)
    if check_lin:
        for a_, direction in zip([a, a.transpose()], ['row', 'column']):
            spacing = np.mean(a_[:, 1:] - a_[:, :-1], axis=1)
            if not np.allclose(a_[:, 1:] - a_[:, :-1], spacing[:, None]):
                raise SpacingError(
                    "Input array is not linearly spaced "
                    f"along the {direction}s."
                )
    if check_loglin:
        _check_epsign(a, strict=True)
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
    a = _check_1d_array(a, check_epsign=True)

    a_l = a[0] * (a[0] / a[1]) ** np.arange(n_ext, 0, -1)
    a_r = a[-1] * (a[-1] / a[-2]) ** np.arange(1, n_ext + 1, 1)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def extrap_loglin_oscil(a, n_ext):
    """Extrapolate a 1-d array with oscillatory behaviour log-linearly.

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

    if len(a) < 8:
        warnings.warn(
            "Array length is less than 8. Extrapolation may not be accurate.",
            UserWarning
        )

    # Check for sign changes and restrict extrapolation to a range
    # with at least 2 points on either 'bank' of the 'zero river'.
    signdiff = np.diff(np.sign(a))
    signchange = np.nonzero(signdiff)[0]
    if len(signchange) > 0:
        idx_lbank = max(signchange[0] + 1, 3)
        idx_rbank = min(signchange[-1] + 1, a.size - 4)
    else:
        idx_lbank = a.size - 1
        idx_rbank = 0

    # Check for oscillations and further restrict extrapolation to a range
    # excluding the innermost peak/trough with at least 2 points on either
    # 'bank' of the 'zero river'.
    lngrad_lbank = np.nan_to_num(np.gradient(np.log(np.abs(a[:idx_lbank]))))
    signchange_lngrad_lbank = np.nonzero(np.diff(np.sign(lngrad_lbank)))[0]
    if len(signchange_lngrad_lbank) > 0:
        idx_ldbank = max(signchange_lngrad_lbank[-1] + 1, 3)
    else:
        idx_ldbank = a.size - 1

    lngrad_rbank = np.nan_to_num(np.gradient(np.log(np.abs(a[idx_rbank:]))))
    signchange_lngrad_rbank = np.nonzero(np.diff(np.sign(lngrad_rbank)))[0]
    if len(signchange_lngrad_rbank) > 0:
        idx_rdbank = min(idx_rbank + signchange_lngrad_rbank[0], a.size - 4)
    else:
        idx_rdbank = 0

    # Divide the array into two 'halves' either side of the peak with
    # each 'half' halved into 'quarters'.
    idx_peak = np.argmax(np.abs(a))
    if idx_peak == 0 or idx_peak == a.size - 1:
        idx_peak = a.size // 2
    idx_lquater = max((idx_peak + 1) // 2 - 1, 3)
    idx_rquater = min(idx_peak + (a.size - 1 - idx_peak) // 2, a.size - 4)

    # Define the tails used as the extrapolation range.
    idx_ltail = min(idx_lbank, idx_ldbank, idx_lquater)
    idx_rtail = max(idx_rbank, idx_rdbank, idx_rquater)

    # Extrapolate the left and right tails.
    lngrad_ltail = np.nan_to_num(np.gradient(np.log(np.abs(a[:idx_ltail]))))
    # Oscillatory behaviour
    if np.diff(np.sign(lngrad_ltail)).any():
        # Use mean gradient.
        if (mean_lngrad_ltail := np.mean(lngrad_ltail)) >= 0.:
            powlaw_ltail = - mean_lngrad_ltail
            # print("Left tail extrap: trend")
        # Use constant padding when mean gradient is divergent.
        else:
            powlaw_ltail = 0.
            warnings.warn(
                "Divergent behaviour detected at the left tail. "
                "Constant padding is used instead.",
                DivergenceWarning
            )
            # print("Left tail extrap: trend -> const")
    # No oscillations
    else:
        # Use endpoint gradient.
        if (
            lngrad_lend := np.nan_to_num(
                np.log(np.abs(a[1])) - np.log(np.abs(a[0]))
            )
        ) >= 0.:
            powlaw_ltail = - lngrad_lend
            # print("Left tail extrap: endpoint")
        # Endpoint gradient is divergent possibly with sign changes.
        else:
            mean_lngrad_ltail = np.nan_to_num(np.mean(
                np.gradient(np.log(np.abs(a[:idx_lquater])))
            ))
            # Use mean gradient.
            if mean_lngrad_ltail >= 0.:
                powlaw_ltail = - mean_lngrad_ltail
                # print("Left tail extrap: endpoint -> trend")
            # Use constant padding when mean gradient is divergent.
            else:
                powlaw_ltail = 0.
                warnings.warn(
                    "Divergent behaviour detected at the left endpoints. "
                    "Constant padding is used instead.",
                    DivergenceWarning
                )
                # print("Left tail extrap: endpoint -> const")

    lngrad_rtail = np.nan_to_num(np.gradient(np.log(np.abs(a[idx_rtail:]))))
    if np.diff(np.sign(lngrad_rtail)).any():
        if (mean_lngrad_rtail := np.mean(lngrad_rtail)) <= 0.:
            powlaw_rtail = mean_lngrad_rtail
            # print("Right tail extrap: trend")
        else:
            powlaw_rtail = 0.
            warnings.warn(
                "Divergent behaviour detected at the right tail. "
                "Constant padding is used instead.",
                DivergenceWarning
            )
            # print("Right tail extrap: trend -> const")
    else:
        if (
            lngrad_rend := np.nan_to_num(
                np.log(np.abs(a[-1])) - np.log(np.abs(a[-2]))
            )
        ) <= 0.:
            powlaw_rtail = lngrad_rend
            # print("Right tail extrap: endpoint")
        else:
            mean_lngrad_rtail = np.nan_to_num(np.mean(
                np.gradient(np.log(np.abs(a[idx_rquater:])))
            ))
            if mean_lngrad_rtail <= 0.:
                powlaw_rtail = mean_lngrad_rtail
                # print("Right tail extrap: endpoint -> trend")
            else:
                powlaw_rtail = 0.
                warnings.warn(
                    "Divergent behaviour detected at the right endpoints. "
                    "Constant padding is used instead.",
                    DivergenceWarning
                )
                # print("Right tail extrap: endpoint -> const")

    a_l = a[0] * np.exp(powlaw_ltail * np.arange(n_ext, 0, -1))
    a_r = a[-1] * np.exp(powlaw_rtail * np.arange(1, n_ext + 1, 1))
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
    a_l = a[:, [0]] * (a[:, [0]]/a[:, [1]]) ** np.arange(ncol_ext, 0, -1)
    a_r = a[:, [-1]] * (a[:, [-1]]/a[:, [-2]]) ** np.arange(1, ncol_ext + 1, 1)

    a_out = np.c_[a_l, a, a_r]

    return a_out


def _extrap2d_row_loglin_oscil(a, ncol_ext):
    """Extrapolate a 2-d array with oscillatory behaviour log-linearly
    along each row (i.e. horizontally).

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
    a_out = []
    for a_row in a:
        a_out.append(extrap_loglin_oscil(a_row, ncol_ext))
    a_out = np.array(a_out)

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
    a = _check_2d_array(a, check_epsign=True)

    nrow_ext, ncol_ext = \
        n_ext if isinstance(n_ext, Sequence) else (n_ext, n_ext)

    # Extrapolate horizontally then vertically.
    a_ = _extrap2d_row_loglin(a, ncol_ext)
    a_out = _extrap2d_row_loglin(a_.transpose(), nrow_ext).transpose()

    return a_out


def extrap2d_loglin_oscil(a, n_ext):
    """Extrapolate a 2-d array with oscillatory behaviour log-bilinearly.

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
    a_ = _extrap2d_row_loglin_oscil(a, ncol_ext)
    a_out = _extrap2d_row_loglin_oscil(a_.transpose(), nrow_ext).transpose()

    return a_out


# ========================================================================
# Reshaping
# ========================================================================

def reshape_threept_datatab(paired_coords, flat_multipole, shape='triu'):
    """Reshape flattened three-point clustering measurements into 2-d
    arrays.

    Parameters
    ----------
    paired_coords : 2-d array of float
        Paired coordinate points, e.g. rows of :math:`(k_1, k_2)` values.
    flat_multipole : 1-d array of float
        Flattened multipole values corresponding to the
        paired coordinates, assuming the two coordinate variables
        are sampled at the same points.
    shape : {'triu', 'full'}, optional
        Shape of the input array without flattening.  If 'triu' (default),
        it corresponds to an upper triangular matrix and the output
        2-d array is symmetric; if 'full', it corresponds to a
        full matrix.

    Returns
    -------
    common_coords : 1-d array of float
        Common coordinate points as a vector.
    multipole : 2-d array of float
        Multipole values as a symmetric matrix.

    Raises
    ------
    ValueError
        If `paired_coords` is not a 2-d array.
    ValueError
        If `paired_coords` and `flat_multipole` have different lengths.
    ValueError
        If `paired_coords` and `flat_multipole` do not correspond to an
        upper triangular matrix.

    """
    # Check dimensions and reshape.
    if np.array(paired_coords).ndim != 2:
        raise ValueError("Paired coordinate points must be a 2-d array.")
    if (flat_len := len(paired_coords)) != len(flat_multipole):
        raise ValueError(
            "Paired coordinate points and flattened multipole "
            "array must have the same length."
        )

    if shape == 'triu':
        # Determine data dimensions.
        nbins = (np.sqrt(8*flat_len + 1) - 1) / 2.
        if not np.isclose(nbins, int(nbins)):
            raise ValueError(
                "Input arrays do not correspond to an upper triangular matrix."
            )
        nbins = int(nbins)

        # Reshape the flattened multipole array into a symmetric matrix.
        common_coords = paired_coords[:nbins, 1]

        triu_indices = np.triu_indices(nbins)

        multipole = np.zeros((nbins, nbins))
        multipole[triu_indices] = flat_multipole
        multipole.T[triu_indices] = flat_multipole
    elif shape == 'full':
        # Determine data dimensions.
        nbins = np.sqrt(flat_len)
        if not np.isclose(nbins, int(nbins)):
            raise ValueError(
                "Input arrays do not correspond to a square matrix."
            )
        nbins = int(nbins)

        # Reshape the flattened multipole array into a full matrix.
        common_coords = paired_coords[:nbins, 1]
        multipole = flat_multipole.reshape((nbins, nbins))
    else:
        raise ValueError(
            f"Invalid shape specified: '{shape}' not 'triu' or 'full'."
        )

    return common_coords, multipole
