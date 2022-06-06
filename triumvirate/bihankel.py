"""
Double Hankel Transforms (:mod:`~triumvirate.bihankel`)
==========================================================================

Perform Hankel-related double transforms using the FFTLog algorithm.

"""
import numpy as np
from scipy import interpolate

import _fftlog


# ########################################################################
# Extrapolation
# ########################################################################

def extrap_lin(a, n_ext):
    """Extrapolate a 1-d array linearly.

    Parameters
    ----------
    a : array_like
        Inout 1-d array.
    n_ext : int
        Extra element number on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a_l = a[0] + (a[1] - a[0]) * np.arange(-n_ext, 0)
    a_r = a[-1] + (a[-1] - a[-2]) * np.arange(1, n_ext + 1)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def extrap_loglin(a, n_ext):
    """Extrapolate a 1-d array log-linearly.

    Parameters
    ----------
    a : array_like
        Inout 1-d array.
    n_ext : int
        Extra element number on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a_l = a[0] * (a[1] / a[0]) ** np.arange(-n_ext, 0)
    a_r = a[-1] * (a[-1] / a[-2]) ** np.arange(1, n_ext + 1)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def extrap_padding(a, n_ext, a_const=0.):
    """Extrapolate a 1-d array by constant padding.

    Parameters
    ----------
    a : array_like
        Inout 1-d array.
    n_ext : int
        Extra element number on either side.
    a_const : float, optional
        Constant value (default is 0.) used for padding.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 1-d array.

    """
    a_l = a_r = np.full(n_ext, a_const)
    a_out = np.concatenate([a_l, a, a_r])

    return a_out


def _extrap2d_lin(a, n_ext):
    """Extrapolate a 2-d array linearly along each row
    (i.e. horizontally).

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        Input 2-d array.
    n_ext : int
        Extra column number on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    """
    da_l = np.diff(a[:, [0, 1]], axis=1)
    da_r = np.diff(a[:, [-2, -1]], axis=1)

    a_l = a[:, [0]] + np.arange(-n_ext, 0) * da_l
    a_r = a[:, [-1]] + np.arange(1, n_ext + 1) * da_r

    a_out = np.c_[a_l, a, a_r]

    return a_out


def _extrap2d_loglin(a, n_ext):
    """Extrapolate a 2-d array log-linearly along each row
    (i.e. horizontally).

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        Input 2-d array.
    n_ext : int
        Extra column number on either side.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    """
    a_l = a[:, [0]] * (a[:, [1]] / a[:, [0]]) ** np.arange(-n_ext, 0)
    a_r = a[:, [-1]] * (a[:, [-1]] / a[:, [-2]]) ** np.arange(1, n_ext + 1)

    a_out = np.c_[a_l, a, a_r]

    return a_out


def extrap2d_bilin(a, n_ext, n_ext_col=None):
    """Extrapolate a 2-d array bilinearly.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        Input 2-d array.
    n_ext : int
        Extra-column number on either side of each row.
    n_ext_col : int, optional
        Extra-row number on either side of each column.  If `None`
        (default), this is set to `n_ext`.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Raises
    ------
    ValueError
        If the input array is not equivalently 2-d.

    """
    a = np.squeeze(a)
    if a.ndim != 2:
        raise ValueError("Input array is not 2-d.")

    if n_ext_col is None:
        n_ext_col = n_ext

    # Extrapolate horizontally.
    a_ = _extrap2d_lin(a, n_ext)

    # Extrapolate vertically.  Use the same algorithm for the horizontal
    # extrapolation by transposing.
    a_ = np.transpose(a_)
    a_out = _extrap2d_lin(a_, n_ext_col)
    a_out = np.transpose(a_out)

    return a_out


def extrap2d_logbilin(a, n_ext, n_ext_col=None):
    """Extrapolate a 2-d array log-bilinearly.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        Input 2-d array.
    n_ext : int
        Extra-column number on either side of each row.
    n_ext_col : int, optional
        Extra-row number on either side of each column.  If `None`
        (default), this is set to `n_ext`.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Raises
    ------
    ValueError
        If the input array is not equivalently 2-d, or contain
        non-positive entries.

    """
    a = np.squeeze(a)
    if a.ndim != 2:
        raise ValueError("Input array is not 2-d.")
    if np.any(a <= 0.):
        raise ValueError("Input array contains non-positive entries.")

    if n_ext_col is None:
        n_ext_col = n_ext

    # Extrapolate horizontally.
    a_ = _extrap2d_loglin(a, n_ext)

    # Extrapolate vertically.  Use the same algorithm for the horizontal
    # extrapolation by transposing.
    a_ = np.transpose(a_)
    a_out = _extrap2d_loglin(a_, n_ext_col)
    a_out = np.transpose(a_out)

    return a_out


def extrap2d_bipad(a, n_ext, n_ext_col=None, a_const=0., a_const_col=None):
    """Extrapolate a 2-d array by constant padding.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        Input 2-d array.
    n_ext : int
        Extra-column number on either side of each row.
    n_ext_col : int, optional
        Extra-row number on either side of each column.  If `None`
        (default), this is set to `n_ext`.
    a_const : float, optional
        Constant value (default is 0.) used for padding in each row.
    a_const_col : float, optional
        Constant value used for padding in each column.  If `None`
        (default), this is set to `a_const`.

    Returns
    -------
    a_out : :class:`numpy.ndarray`
        Extrapolated 2-d array.

    Raises
    ------
    ValueError
        If the input array is not equivalently 2-d.

    """
    a = np.squeeze(a)
    if a.ndim != 2:
        raise ValueError("Input array is not 2-d.")

    if n_ext_col is None:
        n_ext_col = n_ext
    if a_const_col is None:
        a_const_col = a_const

    a_l = a_r = np.full((a.shape[0], n_ext), a_const)
    a_ = np.c_[a_l, a, a_r]
    a__u = a__d = np.full((n_ext_col, a_.shape[-1]), a_const_col)
    a_out = np.r_[a__u, a_, a__d]

    return a_out


# ########################################################################
# Transforms
# ########################################################################

def transform_bispec_to_3pcf(ell1, ell2, bk_in, k_in, r_out,
                             n_fftlog, n_extrap=None, extrap=None):
    """Transform bispectrum samples to 3-point correlation function (3PCF)
    samples using 1-d FFTLog.

    Parameters
    ----------
    ell1, ell2 : int
        Multipole degree of the bispectrum.
    bk_in, k_in : :class:`numpy.ndarray`
        Input bispectrum and wavenumbers.
    r_out : array_like
        Output separations.
    n_fftlog : int
        FFTLog sample number.
    n_extrap : int, optional
        Extrapolation sample number (one-sided) (default is `None`).
    extrap : {'lin', 'loglin', '0-pad'}, optional
        If not `None` (default), set one of the following options:
        * 'lin' -- input bispectrum is extrapolated linearly;
        * 'loglin' -- input bispectrum is extrapolated log-linearly;
        * '0-pad' -- input bispectrum is zero padded.
        In each case, `n_extrap` must be set, and the input wavenumbers
        `k_in` are extrapolated log-linearly.

    Returns
    -------
    zeta_dict : dict
        Output 3PCF including its attributes:
        * 'zeta', 'r' -- transformed 3PCF and separations;
        * 'r1_mesh', 'r2_mesh' -- transformed separation coordinate array;
        * 'zeta_fftlog', 'r_fftlog' -- FFTLog-sampled 3PCF
          and separations;
        * 'ell1', 'ell2' -- multipole degrees;
        * 'n_fftlog' -- FFTLog sample number used.

    Raises
    ------
    ValueError
        If `n_extrap` is not set when `extrap` is.

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
        elif extrap == '0-pad':
            bk_sample = extrap2d_bipad(bk_in, n_extrap)
        else:
            raise ValueError(f"Unknown option for `extrap`: {extrap}.")
    else:
        k_sample, bk_sample = k_in, bk_in

    n_sample = len(k_sample)

    # Perform FFTLog transform.
    k_fftlog = np.logspace(
        np.log(k_sample[0]), np.log(k_sample[-1]), n_fftlog, base=np.e
    )

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

    # An overall complex 'sign' is needed here as it is not included
    # in the backend transforms to avoid redundant computation.
    sign = np.real(1j ** (ell1 + ell2))
    zeta_fftlog = sign * zeta_mesh_full

    # Prepare output.
    interpolator2d_zeta = interpolate.RectBivariateSpline(
        r_fftlog, r_fftlog, zeta_fftlog
    )

    zeta_out = interpolator2d_zeta(r_out, r_out)

    (r2_mesh, r1_mesh) = np.meshgrid(r_out, r_out)

    zeta_dict = {
        'zeta': zeta_out,
        'r': r_out,
        'r1_mesh': r1_mesh,
        'r2_mesh': r2_mesh,
        'zeta_fftlog': zeta_fftlog,
        'r_fftlog': r_fftlog,
        'ell1': ell1,
        'ell2': ell2,
        'n_fftlog': n_fftlog,
    }

    return zeta_dict


def transform_3pcf_to_bispec(ell1, ell2, zeta_in, r_in, k_out,
                             n_fftlog, n_extrap=None, extrap=None):
    """Transform 3-point correlation function (3PCF) samples to bispectrum
    samples using 1-d FFTLog.

    Parameters
    ----------
    ell1, ell2 : int
        Multipole degree of the bispectrum.
    zeta_in, r_in : :class:`numpy.ndarray`
        Input 3PCF and separations.
    k_out : array_like
        Output wavenumbers.
    n_fftlog : int
        FFTLog sample number.
    n_extrap : int, optional
        Extrapolation sample number (one-sided) (default is `None`).
    extrap : {'lin', 'loglin', '0-pad'}, optional
        If not `None` (default), set one of the following options:
        * 'lin' -- input 3PCF is extrapolated linearly;
        * 'loglin' -- input 3PCF is extrapolated log-linearly;
        * '0-pad' -- input 3PCF is zero padded.
        In each case, `n_extrap` must be set, and the input separations
        `r_in` are extrapolated log-linearly.

    Returns
    -------
    bk_dict : dict
        Output bispectrum including its attributes:
        * 'bk', 'k' -- transformed bispectrum and wavenumbers;
        * 'k1_mesh', 'k2_mesh' -- transformed wavenumber coordinate array;
        * 'bk_fftlog', 'k_fftlog' -- FFTLog-sampled bispectrum
          and wavenumbers;
        * 'ell1', 'ell2' -- multipole degrees;
        * 'n_fftlog' -- FFTLog sample number used.

    Raises
    ------
    ValueError
        If `n_extrap` is not set when `extrap` is.

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
        elif extrap == '0-pad':
            zeta_sample = extrap2d_bipad(zeta_in, n_extrap)
        else:
            raise ValueError(f"Unknown option for `extrap`: {extrap}.")
    else:
        r_sample, zeta_sample = r_in, zeta_in

    n_sample = len(r_sample)

    # Perform FFTLog transform.
    r_fftlog = np.logspace(
        np.log(r_sample[0]), np.log(r_sample[-1]), n_fftlog, base=np.e
    )

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

    # An overall complex 'sign' is needed here as it is not included
    # in the backend transforms to avoid redundant computation.
    # Factors of pi are needed here because the forward transform is used
    # in the backend as the backward transform.
    sign = (2*np.pi) ** 6 * np.real((-1j) ** (ell1 + ell2))
    bk_fftlog = sign * bk_mesh_full

    # Prepare output.
    interpolator2d_bk = interpolate.RectBivariateSpline(
        k_fftlog, k_fftlog, bk_fftlog
    )

    bk_out = interpolator2d_bk(k_out, k_out)

    (k2_mesh, k1_mesh) = np.meshgrid(k_out, k_out)

    bk_dict = {
        'bk': bk_out,
        'k': k_out,
        'k1_mesh': k1_mesh,
        'k2_mesh': k2_mesh,
        'bk_fftlog': bk_fftlog,
        'k_fftlog': k_fftlog,
        'ell1': ell1,
        'ell2': ell2,
        'n_fftlog': n_fftlog,
    }

    return bk_dict
