"""
Double Hankel Transforms (:mod:`~triumvirate.bihankel`)
==========================================================================

Perform Hankel-related double transforms using the FFTLog algorithm.

Array operations
----------------

.. autosummary::
    ShapeError

1-d arrays:

.. autosummary::
    extrap_lin
    extrap_loglin
    extrap_pad

2-d arrays:

.. autosummary:
    extrap2d_bilin
    extrap2d_biloglin
    extrap2d_bipad

Transforms
----------

.. autosummary::
    transform_bispec_to_3pcf
    transform_3pcf_to_bispec

"""
import numpy as np
from scipy import interpolate

import triumvirate._fftlog as _fftlog


class ShapeError(ValueError):
    """Exception raised when the array shape is incorrect.

    """
    pass


# ========================================================================
# Extrapolation
# ========================================================================

def _check_1d_array(a):
    """Check the input array is a sorted 1-d array.

    Parameters
    ----------
    a : array of float
        Input 1-d array.

    Returns
    -------
    a : array of float
        Input 1-d array.

    Raises
    ------
    ShapeError
        When the input array is not 1-d.

    """
    a = np.squeeze(a)

    if a.ndim != 1:
        raise ShapeError("Input array is not 1-d.")
    # if not np.all(a[:-1] <= a[1:]):
    #     raise ValueError("Input array is not an increasing sequence.")

    return a


def _check_2d_array(a):
    """Check the input array is a sorted 2-d array.

    Parameters
    ----------
    a : array of float
        Input 2-d array.

    Returns
    -------
    a : array of float
        Input 2-d array.

    Raises
    ------
    ShapeError
        When the input array is not 2-d.

    """
    a = np.squeeze(a)
    if a.ndim != 2:
        raise ShapeError("Input array is not 1-d.")
    # if (
    #     not np.all([np.all(a_row[:-1] <= a_row[1:]) for a_row in a])
    # ) or (
    #     not np.all([
    #         np.all(a_col[:-1] <= a_col[1:]) for a_col in a.transpose()
    #     ])
    # ):
    #     raise ValueError("Input array is not an increasing sequence.")

    return a


def extrap_lin(a, n_ext):
    """Extrapolate a 1-d array linearly.

    The linear spacing used is the distance between the two points
    at either edge separately.

    Parameters
    ----------
    a : 1-d array of float
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

    Parameters
    ----------
    a : 1-d array of float
        Inout 1-d array.
    n_ext : int
        Number of extra elements on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    Raises
    ------
    ValueError
        When the input array contains non-positive entries.

    """
    a = _check_1d_array(a)
    if np.any(a <= 0.):
        raise ValueError("Input array contains non-positive entries.")

    a_l = a[0] * (a[1] / a[0]) ** np.arange(-n_ext, 0)
    a_r = a[-1] * (a[-1] / a[-2]) ** np.arange(1, n_ext + 1)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def extrap_pad(a, n_ext, a_const=0.):
    """Extrapolate a 1-d array by constant padding.

    Parameters
    ----------
    a : 1-d array of float
        Input 1-d array.
    n_ext : int
        Number of extra elements on either side.
    a_const : float, optional
        Constant value (default is 0.) used for padding.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a = _check_1d_array(a)

    a_l = a_r = np.full(n_ext, a_const)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def _extrap2d_lin(a, ncol_ext):
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


def _extrap2d_loglin(a, ncol_ext):
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


def extrap2d_bilin(a, ncol_ext, nrow_ext=None):
    """Extrapolate a 2-d array bilinearly.

    Parameters
    ----------
    a : 2-d array of float
        Input 2-d array.
    ncol_ext : int
        Number of extra columns on either side of each row.
    nrow_ext : int, optional
        Number of extra rows on either side of each column.  If `None`
        (default), this is set to `ncol_ext`.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Raises
    ------
    ValueError
        If the input array is not equivalently 2-d.

    """
    a = _check_2d_array(a)

    # Extrapolate horizontally.
    a_ = _extrap2d_lin(a, ncol_ext)

    # Extrapolate vertically.  Use the same algorithm for the horizontal
    # extrapolation by transposing.
    nrow_ext = ncol_ext if nrow_ext is None else nrow_ext

    a_out = _extrap2d_lin(a_.transpose(), nrow_ext).transpose()

    return a_out


def extrap2d_logbilin(a, ncol_ext, nrow_ext=None):
    """Extrapolate a 2-d array log-bilinearly.

    Parameters
    ----------
    a : 2-d array of float
        Input 2-d array.
    ncol_ext : int
        Number of extra columns on either side of each row.
    nrow_ext : int, optional
        Number of extra rows on either side of each column.  If `None`
        (default), this is set to `ncol_ext`.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Raises
    ------
    ValueError
        When the input array contains non-positive entries.

    """
    a = _check_2d_array(a)
    if np.any(a <= 0.):
        raise ValueError("Input array contains non-positive entries.")

    # Extrapolate horizontally.
    a_ = _extrap2d_loglin(a, ncol_ext)

    # Extrapolate vertically.  Use the same algorithm for the horizontal
    # extrapolation by transposing.
    nrow_ext = ncol_ext if nrow_ext is None else nrow_ext

    a_out = _extrap2d_loglin(a_.transpose(), nrow_ext).transpose()

    return a_out


def extrap2d_bipad(a, ncol_ext, nrow_ext=None,
                   a_const_col=0., a_const_row=None):
    """Extrapolate a 2-d array by constant padding.

    Parameters
    ----------
    a : 2-d array of float
        Input 2-d array.
    ncol_ext : int
        Number of extra columns on either side of each row.
    nrow_ext : int, optional
        Number of extra rows on either side of each column.  If `None`
        (default), this is set to `ncol_ext`.
    a_const_col : float, optional
        Constant value (default is 0.) used for padding in each row.
    a_const_row : float, optional
        Constant value used for padding in each column.  If `None`
        (default), this is set to `a_const_col`.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Raises
    ------
    ValueError
        If the input array is not equivalently 2-d.

    """
    a = _check_2d_array(a)

    nrow_ext = ncol_ext if nrow_ext is None else nrow_ext
    a_const_row = a_const_col if a_const_row is None else a_const_row

    a_l = a_r = np.full((a.shape[0], ncol_ext), a_const_col)
    a_ = np.c_[a_l, a, a_r]

    a__u = a__d = np.full((nrow_ext, a_.shape[-1]), a_const_row)
    a_out = np.r_[a__u, a_, a__d]

    return a_out


# ========================================================================
# Transforms
# ========================================================================

def transform_bispec_to_3pcf(ell1, ell2, k_in, bk_in, r_out,
                             n_fftlog, extrap=None, n_extrap=None):
    r"""Transform bispectrum samples to three-point correlation function
    (3PCF) samples using double 1-d FFTLog operations.

    Parameters
    ----------
    ell1, ell2 : int
        Multipole degree of the bispectrum.
    k_in : 1-d array of float
        Input wavenumbers.
    bk_in : 2-d array of float
        Input bispectrum.
    r_out : 1-d array of float
        Output separations.
    n_fftlog : int
        FFTLog sample number.
    extrap : {'lin', 'loglin', 'zero'}, optional
        If not `None` (default), set one of the following options---

        - 'lin': input bispectrum is extrapolated linearly;
        - 'loglin': input bispectrum is extrapolated log-linearly;
        - 'zero': input bispectrum is zero padded.

        In each case, `n_extrap` must be set, and the input wavenumbers
        `k_in` are extrapolated log-linearly.
    n_extrap : int, optional
        Extrapolation sample number (one-sided) (default is `None`).

    Returns
    -------
    zeta_dict : dict
        Output 3PCF including its attributes---

        - 'r', 'zeta': post-transform separation and 3PCF samples;
        - 'r_fftlog', 'zeta_fftlog': FFTLog separation and 3PCF samples;
        - 'r1_mesh', 'r2_mesh': post-transform separation mesh grid;
        - 'n_fftlog': FFTLog sample number;
        - 'ell1', 'ell2': multipole degrees.

    Raises
    ------
    ValueError
        When `n_extrap` is not set but `extrap` is.


    .. attention::

        This function assumes the wavenumber samples :math:`k_{1,2}` in
        the bispectrum :math:`B_{\ell_1 \ell_2 L}(k_1, k_2)` are the same.

    """
    # Prepare samples.
    if extrap is not None:
        if n_extrap is None:
            raise ValueError("`n_extrap` must be provided if `extrap` is set.")
        k_sample = extrap_loglin(k_in, n_extrap)

        if extrap == 'lin':
            bk_sample = extrap2d_bilin(bk_in, n_extrap)
        elif extrap == 'loglin':
            bk_sample = extrap2d_logbilin(bk_in, n_extrap)
        elif extrap == 'zero':
            bk_sample = extrap2d_bipad(bk_in, n_extrap)
        else:
            raise ValueError(f"Unknown `extrap` option: {extrap}.")
    else:
        k_sample, bk_sample = k_in, bk_in

    n_sample = len(k_sample)

    # Perform FFTLog transform.
    # CAVEAT: The choice of interpolator is discretionary but should
    # have no effect.
    k_fftlog = np.logspace(*np.log(k_sample[[0, -1]]), n_fftlog, base=np.e)

    zeta_mesh_partial = np.zeros((n_sample, n_fftlog))
    for i in range(n_sample):
        interpolator_bk = interpolate.interp1d(
            k_sample, bk_sample[i, :],
            fill_value='extrapolate', kind='cubic'
        )
        bk_fftlog = interpolator_bk(k_fftlog)

        r_fftlog, zeta_fftlog = _fftlog.sj_transform(
            ell2, 2, n_fftlog, k_fftlog, bk_fftlog
        )

        zeta_mesh_partial[i, :] = zeta_fftlog

    zeta_mesh_full = np.zeros((n_fftlog, n_fftlog))
    for j in range(n_fftlog):
        interpolator_bk = interpolate.interp1d(
            k_sample, zeta_mesh_partial[:, j],
            fill_value='extrapolate', kind='cubic'
        )
        bk_fftlog = interpolator_bk(k_fftlog)

        r_fftlog, zeta_fftlog = _fftlog.sj_transform(
            ell1, 2, n_fftlog, k_fftlog, bk_fftlog
        )

        zeta_mesh_full[:, j] = zeta_fftlog

    # An overall complex `parity` is needed here as it is not included
    # in the backend transforms to avoid redundant computation.
    parity = np.real(1j ** (ell1 + ell2))
    zeta_fftlog = parity * zeta_mesh_full

    # Prepare output.
    interpolator2d_zeta = interpolate.RectBivariateSpline(
        r_fftlog, r_fftlog, zeta_fftlog
    )

    zeta_out = interpolator2d_zeta(r_out, r_out)

    (r2_mesh, r1_mesh) = np.meshgrid(r_out, r_out)

    zeta_dict = {
        'r': r_out,
        'zeta': zeta_out,
        'r_fftlog': r_fftlog,
        'zeta_fftlog': zeta_fftlog,
        'r1_mesh': r1_mesh,
        'r2_mesh': r2_mesh,
        'n_fftlog': n_fftlog,
        'ell1': ell1,
        'ell2': ell2,
    }

    return zeta_dict


def transform_3pcf_to_bispec(ell1, ell2, r_in, zeta_in, k_out,
                             n_fftlog, extrap=None, n_extrap=None):
    r"""Transform three-point correlation function (3PCF) samples to
    bispectrum samples using double 1-d FFTLog operations.

    Parameters
    ----------
    ell1, ell2 : int
        Multipole degree of the 3PCF.
    r_in : 1-d array of float
        Input separations.
    zeta_in : 2-d array of float
        Input 3PCF.
    k_out : 1-d array of float
        Output wavenumbers.
    n_fftlog : int
        FFTLog sample number.
    extrap : {'lin', 'loglin', 'zero'}, optional
        If not `None` (default), set one of the following options:

        - 'lin': input 3PCF is extrapolated linearly;
        - 'loglin': input 3PCF is extrapolated log-linearly;
        - 'zero': input 3PCF is zero padded.

        In each case, `n_extrap` must be set, and the input separations
        `r_in` are extrapolated log-linearly.
    n_extrap : int, optional
        Extrapolation sample number (one-sided) (default is `None`).

    Returns
    -------
    bk_dict : dict
        Output bispectrum including its attributes:

        - 'k', 'bk': post-transform wavenumber and bispectrum samples;
        - 'k_fftlog', 'bk_fftlog': FFTLog wavenumber and
          bispectrum samples;
        - 'k1_mesh', 'k2_mesh': post-transform wavenumber mesh grid;
        - 'n_fftlog': FFTLog sample number;
        - 'ell1', 'ell2': multipole degrees.

    Raises
    ------
    ValueError
        When `n_extrap` is not set but `extrap` is.


    .. attention::

        This function assumes the separation samples :math:`r_{1,2}` in
        the 3PCF :math:`\zeta_{\ell_1 \ell_2 L}(r_1, r_2)` are the same.

    """
    # Prepare samples.
    if extrap is not None:
        if n_extrap is None:
            raise ValueError("`n_extrap` must be provided if `extrap` is set.")
        r_sample = extrap_loglin(r_in, n_extrap)

        if extrap == 'lin':
            zeta_sample = extrap2d_bilin(zeta_in, n_extrap)
        elif extrap == 'loglin':
            zeta_sample = extrap2d_logbilin(zeta_in, n_extrap)
        elif extrap == 'zero':
            zeta_sample = extrap2d_bipad(zeta_in, n_extrap)
        else:
            raise ValueError(f"Unknown `extrap` option: {extrap}.")
    else:
        r_sample, zeta_sample = r_in, zeta_in

    n_sample = len(r_sample)

    # Perform FFTLog transform.
    # CAVEAT: The choice of interpolator is discretionary but should
    # have no effect.
    r_fftlog = np.logspace(*np.log(r_sample[[0, -1]]), n_fftlog, base=np.e)

    bk_mesh_partial = np.zeros((n_sample, n_fftlog))
    for i in range(n_sample):
        interpolator_zeta = interpolate.interp1d(
            r_sample, zeta_sample[i, :],
            fill_value='extrapolate', kind='cubic'
        )
        zeta_fftlog = interpolator_zeta(r_fftlog)

        k_fftlog, bk_fftlog = _fftlog.sj_transform(
            ell2, 2, n_fftlog, r_fftlog, zeta_fftlog
        )

        bk_mesh_partial[i, :] = bk_fftlog

    bk_mesh_full = np.zeros((n_fftlog, n_fftlog))
    for j in range(n_fftlog):
        interpolator_zeta = interpolate.interp1d(
            r_sample, bk_mesh_partial[:, j],
            fill_value='extrapolate', kind='cubic'
        )
        zeta_fftlog = interpolator_zeta(r_fftlog)

        k_fftlog, bk_fftlog = _fftlog.sj_transform(
            ell1, 2, n_fftlog, r_fftlog, zeta_fftlog
        )

        bk_mesh_full[:, j] = bk_fftlog

    # An overall complex `parity` is needed here as it is not included
    # in the backend transforms to avoid redundant computation.
    # Factors of pi are needed here because the forward transform is used
    # in the backend as the backward transform.
    parity = np.real((-1j) ** (ell1 + ell2))
    bk_fftlog = (2*np.pi) ** 6 * parity * bk_mesh_full

    # Prepare output.
    interpolator2d_bk = interpolate.RectBivariateSpline(
        k_fftlog, k_fftlog, bk_fftlog
    )

    bk_out = interpolator2d_bk(k_out, k_out)

    (k2_mesh, k1_mesh) = np.meshgrid(k_out, k_out)

    bk_dict = {
        'k': k_out,
        'bk': bk_out,
        'k_fftlog': k_fftlog,
        'bk_fftlog': bk_fftlog,
        'k1_mesh': k1_mesh,
        'k2_mesh': k2_mesh,
        'n_fftlog': n_fftlog,
        'ell1': ell1,
        'ell2': ell2,
    }

    return bk_dict
