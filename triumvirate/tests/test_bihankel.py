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
    "ell1,ell2,ELL,N_fftlog,case",
    [
        (0, 0, 0, 1024, ''),
    ]
)
def test_transforms(ell1, ell2, ELL, N_fftlog, case):

    # Load test data.
    bk_dict = np.load(
        f"triumvirate/tests/test_input/test_bk{ell1}{ell2}{ELL}{case}.npy",
        allow_pickle=True
    ).item()

    kbin = bk_dict['kbin2'][0]

    # Set up binning.
    rbin = np.logspace(
        np.log(1. * 1e-1), np.log(1./3 * 1e4),
        num=N_fftlog, base=np.e
    )
    # rbin = np.linspace(0.5, 200, 400)

    # Perform transformation.
    zeta_dict = hankel.transform_bispec_to_3pcf(
        bk_dict, rbin, N_fftlog=N_fftlog
    )

    bk_dict_reverse = hankel.transform_3pcf_to_bispec(
        zeta_dict, kbin, N_fftlog=N_fftlog
    )

    # Verify results.
    num_kbins = len(bk_dict['kbin1'])
    bk_in = np.array([
        bk_dict['kbin1'].reshape(num_kbins ** 2),
        bk_dict['kbin2'].reshape(num_kbins ** 2),
        bk_dict['bk'].reshape(num_kbins ** 2)
    ]).T
    np.savetxt(
        f"triumvirate/tests/test_output/bk{ell1}{ell2}{ELL}{case}_ref.csv",
        bk_in,
        fmt='%.9e'
    )

    num_rbins = len(zeta_dict['rbin1'])
    zeta = np.array([
        zeta_dict['rbin1'].reshape(num_rbins ** 2),
        zeta_dict['rbin2'].reshape(num_rbins ** 2),
        zeta_dict['zeta'].reshape(num_rbins ** 2)
    ]).T
    np.savetxt(
        f"triumvirate/tests/test_output/zeta{ell1}{ell2}{ELL}{case}.csv",
        zeta,
        fmt='%.9e'
    )

    num_kbins = len(bk_dict_reverse['kbin1'])
    bk_out = np.array([
        bk_dict_reverse['kbin1'].reshape(num_kbins ** 2),
        bk_dict_reverse['kbin2'].reshape(num_kbins ** 2),
        bk_dict_reverse['bk'].reshape(num_kbins ** 2)
    ]).T
    np.savetxt(
        f"triumvirate/tests/test_output/bk{ell1}{ell2}{ELL}{case}_reverse.csv",
        bk_out,
        fmt='%.9e'
    )

    assert np.allclose(bk_out, bk_in, rtol=0.001), \
        "Inverse Hankel transform does not recover input bispectrum!"
