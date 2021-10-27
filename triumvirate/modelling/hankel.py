"""
Hankel Transform (:mod:`~triumvirate.modelling.hankel`)
==========================================================================

Perform cosmological Hankel transform (multi-fold spherical
Bessel transform).

"""
import numpy as np
from scipy import interpolate

import hankel_fftlog


def transform_bispec_to_3pcf(bk_dict, rbin, N_fftlog=1024):
    """Double Hankel transform bispectrum samples to 3-point
    correlation function (3PCF) samples using FFTLog.

    Parameters
    ----------
    bk_dict : dict
        Input bispectrum including its attributes.  Mandatory keys are:
        * multipole degrees 'ell1', 'ell2', 'ELL';
        * log-spaced sample points 'kbin1_fftlog', 'kbin2_fftlog' (values
          of which are 2-d arrays) for FFTLog;
        * bispectrum samples 'bk_fftlog' (value of which is a 2-d array).
    rbin : array of float
        Output separation bins.
    N_fftlog : int, optional
        FFTLog sample number (default is 1024).

    Returns
    -------
    zeta_dict : dict
        Output 3PCF including its attributes.  Included keys are:
        * multipole degrees 'ell1', 'ell2', 'ELL';
        * input sample points 'rbin1', 'rbin2' (values of which are
          2-d arrays);
        * output 3PCF samples 'zeta' (value of which is a 2-d array)
          correponding to 'rbin1' and 'rbin2';
        * log-spaced sample points 'rbin1_fftlog', 'rbin2_fftlog' (values
          of which are 2-d arrays) generated by FFTLog;
        * 3PCF samples 'zeta_fftlog' (value of which is a 2-d array)
          correponding to 'rbin1_fftlog' and 'rbin2_fftlog';
        * FFTLog sampling number 'N_fftlog'.

    """
    # Set multipole degrees.
    ell1, ell2, ELL = bk_dict['ell1'], bk_dict['ell2'], bk_dict['ELL']

    # Set input arrays.
    k_in, bk_in = bk_dict['kbin1_fftlog'][:, 0], bk_dict['bk_fftlog']

    # Compute 3PCF.
    k_sample = np.logspace(
        np.log(k_in[0]), np.log(k_in[-1]), N_fftlog, base=np.e
    )

    zeta_mesh_partial = np.zeros((len(k_in), N_fftlog))
    for i in range(len(k_in)):
        interpolator_bk = interpolate.interp1d(
            k_in, bk_in[i, :],
            fill_value='extrapolate', kind='cubic'
        )
        bk_sample = interpolator_bk(k_sample)

        r_fftlog, zeta_fftlog = np.zeros(N_fftlog), np.zeros(N_fftlog)
        hankel_fftlog.hankel_transform(
            ell2, 2, N_fftlog, k_sample, bk_sample, r_fftlog, zeta_fftlog
        )

        interpolator_zeta = interpolate.interp1d(
            r_fftlog, zeta_fftlog,
            fill_value='extrapolate', kind='cubic'
        )
        zeta_mesh_partial[i, :] = interpolator_zeta(r_fftlog[:])

    zeta_mesh_full = np.zeros((N_fftlog, N_fftlog))
    for j in range(N_fftlog):
        interpolator_bk = interpolate.interp1d(
            k_in, zeta_mesh_partial[:, j],
            fill_value='extrapolate', kind='cubic'
        )
        bk_sample = interpolator_bk(k_sample)

        r_fftlog, zeta_fftlog = np.zeros(N_fftlog), np.zeros(N_fftlog)
        hankel_fftlog.hankel_transform(
            ell1, 2, N_fftlog, k_sample, bk_sample, r_fftlog, zeta_fftlog
        )

        interpolator_zeta = interpolate.interp1d(
            r_fftlog, zeta_fftlog,
            fill_value='extrapolate', kind='cubic'
        )
        zeta_mesh_full[:, j] = interpolator_zeta(r_fftlog[:])

    # An overall complex 'sign' is needed here as it is not included
    # in the backend transforms to avoid redundant computation.
    sign = np.real(1j ** (ell1 + ell2))
    zeta_fftlog = sign * zeta_mesh_full

    # Prepare output.
    interpolator2d_zeta = interpolate.interp2d(
        r_fftlog, r_fftlog, zeta_fftlog, kind='cubic'
    )

    zeta_out = interpolator2d_zeta(rbin[:], rbin[:])

    (rbin2_out, rbin1_out) = np.meshgrid(rbin, rbin)
    (rbin2_fftlog, rbin1_fftlog) = np.meshgrid(r_fftlog, r_fftlog)

    zeta_dict = {
        'zeta': zeta_out,
        'rbin1': rbin1_out,
        'rbin2': rbin2_out,
        'zeta_fftlog': zeta_fftlog,
        'rbin1_fftlog': rbin1_fftlog,
        'rbin2_fftlog': rbin2_fftlog,
        'ell1': ell1,
        'ell2': ell2,
        'ELL': ELL,
        'N_fftlog': N_fftlog,
    }

    return zeta_dict


def transform_3pcf_to_bispec(zeta_dict, kbin, N_fftlog=1024):
    """Double Hankel transform 3-point correlation function (3PCF)
    samples to bispectrum samples using FFTLog.

    Parameters
    ----------
    zeta_dict : dict
        Input 3PCF including its attributes.  Mandatory keys are:
        * multipole degrees 'ell1', 'ell2', 'ELL';
        * log-spaced sample points 'rbin1_fftlog', 'rbin2_fftlog' (values
          of which are 2-d arrays) for FFTLog;
        * bispectrum samples 'zeta_fftlog' (value of which is a 2-d array).
    kbin : array of float
        Output wavenumber bins.
    N_fftlog : int, optional
        FFTLog sample number (default is 1024).

    Returns
    -------
    bk_dict : dict
        Output bispectrum including its attributes.  Included keys are:
        * multipole degrees 'ell1', 'ell2', 'ELL';
        * input sample points 'kbin1', 'kbin2' (values of which are
          2-d arrays);
        * output 3PCF samples 'bk' (value of which is a 2-d array)
          correponding to 'kbin1' and 'kbin2';
        * log-spaced sample points 'kbin1_fftlog', 'kbin2_fftlog' (values
          of which are 2-d arrays) generated by FFTLog;
        * 3PCF samples 'bk_fftlog' (value of which is a 2-d array)
          correponding to 'kbin1_fftlog' and 'kbin2_fftlog';
        * FFTLog sampling number 'N_fftlog'.

    """
    # Set multipole degrees.
    ell1, ell2, ELL = zeta_dict['ell1'], zeta_dict['ell2'], zeta_dict['ELL']

    # Set input arrays.
    r_in, zeta_in = zeta_dict['rbin1_fftlog'][:, 0], zeta_dict['zeta_fftlog']

    # Compute bispectrum.
    r_sample = np.logspace(
        np.log(r_in[0]), np.log(r_in[-1]), N_fftlog, base=np.e
    )

    bk_mesh_partial = np.zeros((len(r_in), N_fftlog))
    for i in range(len(r_in)):
        interpolator_zeta = interpolate.interp1d(
            r_in, zeta_in[i, :],
            fill_value='extrapolate', kind='cubic'
        )
        zeta_sample = interpolator_zeta(r_sample)

        k_fftlog, bk_fftlog = np.zeros(N_fftlog), np.zeros(N_fftlog)
        hankel_fftlog.hankel_transform(
            ell2, 2, N_fftlog, r_sample, zeta_sample, k_fftlog, bk_fftlog
        )

        interpolato_bk = interpolate.interp1d(
            k_fftlog, bk_fftlog,
            fill_value='extrapolate', kind='cubic'
        )
        bk_mesh_partial[i, :] = interpolato_bk(k_fftlog[:])

    bk_mesh_full = np.zeros((N_fftlog, N_fftlog))
    for j in range(N_fftlog):
        interpolator_zeta = interpolate.interp1d(
            r_in, bk_mesh_partial[:, j],
            fill_value='extrapolate', kind='cubic'
        )
        zeta_sample = interpolator_zeta(r_sample)

        k_fftlog, bk_fftlog = np.zeros(N_fftlog), np.zeros(N_fftlog)
        hankel_fftlog.hankel_transform(
            ell1, 2, N_fftlog, r_sample, zeta_sample, k_fftlog, bk_fftlog
        )

        interpolato_bk = interpolate.interp1d(
            k_fftlog, bk_fftlog,
            fill_value='extrapolate', kind='cubic'
        )
        bk_mesh_full[:, j] = interpolato_bk(k_fftlog[:])

    # An overall complex 'sign' is needed here as it is not included
    # in the backend transforms to avoid redundant computation.
    # Factors of pi are needed here because the forward transform is used
    # in the backend as the backward transform.
    sign = (2*np.pi) ** 6 * np.real((-1j) ** (ell1 + ell2))
    bk_fftlog = sign * bk_mesh_full

    # Prepare output.
    interpolator2d_bk = interpolate.interp2d(
        k_fftlog, k_fftlog, bk_fftlog, kind='cubic'
    )

    bk_out = interpolator2d_bk(kbin[:], kbin[:])

    (kbin2_out, kbin1_out) = np.meshgrid(kbin, kbin)
    (kbin2_fftlog, kbin1_fftlog) = np.meshgrid(k_fftlog, k_fftlog)

    bk_dict = {
        'bk': bk_out,
        'kbin1': kbin1_out,
        'kbin2': kbin2_out,
        'bk_fftlog': bk_fftlog,
        'kbin1_fftlog': kbin1_fftlog,
        'kbin2_fftlog': kbin2_fftlog,
        'ell1': ell1,
        'ell2': ell2,
        'ELL': ELL,
        'N_fftlog': N_fftlog,
    }

    return bk_dict
