import numpy as np
import pytest

try:
    import triumvirate.bihankel as hankel
except (ImportError, ModuleNotFoundError):
    import os, sys

    # Add to Python search path.
    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../"
    ))
    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../../"
    ))

    import triumvirate.bihankel as hankel


@pytest.mark.parametrize(
    "ell1,ell2,ELL,n_fftlog,n_extrap,extrap,case",
    [
        (0, 0, 0, 1024, 32, 'loglin', ''),
    ]
)
def test_transforms(ell1, ell2, ELL, n_fftlog, n_extrap, extrap, case):

    # Load test data.
    bk_data = np.load(
        f"triumvirate/tests/test_input/test_bk{ell1}{ell2}{ELL}{case}.npy",
        allow_pickle=True
    ).item()

    # Set up binning.
    k_out = bk_data['k']

    r_out = np.logspace(
        np.log(1. * 1e-1), np.log(1./3 * 1e4),
        num=n_fftlog, base=np.e
    )
    # r_out = np.linspace(0.5, 200, 400)

    # Perform transformation.
    zeta_data = hankel.transform_bispec_to_3pcf(
        ell1, ell2,
        bk_data['bk'], bk_data['k'], r_out,
        n_fftlog, n_extrap=n_extrap, extrap=extrap
    )

    bk_data_reverse = hankel.transform_3pcf_to_bispec(
        ell1, ell2,
        zeta_data['zeta_fftlog'], zeta_data['r_fftlog'], k_out,
        n_fftlog
    )

    # Verify results.
    num_kbins = len(bk_data['k'])
    bk_in = np.array([
        bk_data['k1_mesh'].reshape(num_kbins ** 2),
        bk_data['k2_mesh'].reshape(num_kbins ** 2),
        bk_data['bk'].reshape(num_kbins ** 2)
    ]).T
    np.savetxt(
        f"triumvirate/tests/test_output/bk{ell1}{ell2}{ELL}{case}_ref.csv",
        bk_in,
        fmt='%.9e'
    )

    num_rbins = len(zeta_data['r'])
    zeta = np.array([
        zeta_data['r1_mesh'].reshape(num_rbins ** 2),
        zeta_data['r2_mesh'].reshape(num_rbins ** 2),
        zeta_data['zeta'].reshape(num_rbins ** 2)
    ]).T
    np.savetxt(
        f"triumvirate/tests/test_output/zeta{ell1}{ell2}{ELL}{case}.csv",
        zeta,
        fmt='%.9e'
    )

    num_kbins = len(bk_data_reverse['k'])
    bk_out = np.array([
        bk_data_reverse['k1_mesh'].reshape(num_kbins ** 2),
        bk_data_reverse['k2_mesh'].reshape(num_kbins ** 2),
        bk_data_reverse['bk'].reshape(num_kbins ** 2)
    ]).T
    np.savetxt(
        f"triumvirate/tests/test_output/bk{ell1}{ell2}{ELL}{case}_reverse.csv",
        bk_out,
        fmt='%.9e'
    )

    assert np.allclose(bk_out, bk_in, rtol=0.001), \
        "Inverse Hankel transform does not recover input bispectrum!"
