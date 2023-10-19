"""
Double Hankel Transforms (:mod:`~triumvirate.bihankel`)
==========================================================================

Perform Hankel-related double transforms using the FFTLog algorithm.

.. autosummary::
    transform_bispec_to_3pcf
    transform_3pcf_to_bispec

"""
import numpy as np
from scipy import interpolate

import triumvirate._fftlog as _fftlog


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
            bk_sample = extrap2d_lin(bk_in, n_extrap)
        elif extrap == 'loglin':
            bk_sample = extrap2d_loglin(bk_in, n_extrap)
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
            zeta_sample = extrap2d_lin(zeta_in, n_extrap)
        elif extrap == 'loglin':
            zeta_sample = extrap2d_loglin(zeta_in, n_extrap)
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
