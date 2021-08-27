import numpy as np
import pytest

try:
    import triumvirate.modelling.hankel as hankel
except ImportError:
    import os, sys

    # Add to Python search path.
    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../../"
    ))
    sys.path.insert(0, os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "../../../"
    ))

    import triumvirate.modelling.hankel as hankel

@pytest.mark.parametrize(
    "ell1,ell2,ELL,N_fftlog,case",
    [
        (0, 0, 0, 256, 'tree'),
    ]
)
def test_transforms(ell1, ell2, ELL, N_fftlog, case):

    # Set up binning.
    kbin = np.linspace(0.01, 0.2, 20)
    rbin = np.linspace(0.5, 200, 400)

    # Perform transformations.
    bk_dict = np.load(
        f"test_input/test_bk_{case}.npy", allow_pickle=True
    ).item()
    zeta_dict = hankel.transform_bispec_to_3pcf(
        bk_dict, rbin, N_fftlog=N_fftlog
    )
    bk_dict_reverse = hankel.transform_3pcf_to_bispec(zeta_dict, kbin)

    # Verify results.
    num_kbins = len(bk_dict['kbin1'])
    bk_in = np.array([
        bk_dict['kbin1'].reshape(num_kbins ** 2),
        bk_dict['kbin2'].reshape(num_kbins ** 2),
        bk_dict['bk'].reshape(num_kbins ** 2)
    ]).T
    np.savetxt(
        "test_output/modelling/bk{:d}{:d}{:d}_{}.dat".format(
            ell1, ell2, ELL, case
        ),
        bk_in, fmt='%.7e'
    )

    num_rbins = len(rbin)
    zeta_mid = np.array([
        zeta_dict['rbin1'].reshape(num_rbins ** 2),
        zeta_dict['rbin2'].reshape(num_rbins ** 2),
        zeta_dict['zeta'].reshape(num_rbins ** 2)
    ]).T
    np.savetxt(
        "test_output/modelling/zeta{:d}{:d}{:d}_{}.dat".format(
            ell1, ell2, ELL, case
        ),
        zeta_mid, fmt='%.7e'
    )

    bk_out = np.array([
        bk_dict_reverse['kbin1'].reshape(num_kbins ** 2),
        bk_dict_reverse['kbin2'].reshape(num_kbins ** 2),
        bk_dict_reverse['bk'].reshape(num_kbins ** 2)
    ]).T
    np.savetxt(
        "test_output/modelling/bk{:d}{:d}{:d}_{}_reverse.dat".format(
            ell1, ell2, ELL, case
        ),
        bk_out, fmt='%.7e'
    )

    assert np.allclose(bk_out, bk_in, rtol=0.001), \
        "Inverse Hankel transform does not recover input bispectrum!"
