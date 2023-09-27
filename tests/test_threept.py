"""Test :mod:`~triumvirate.threept`.

"""
import numpy as np
import pytest

from triumvirate.threept import (
    compute_3pcf,
    compute_3pcf_in_gpp_box,
    compute_3pcf_window,
    compute_bispec,
    compute_bispec_in_gpp_box
)


@pytest.mark.slow
@pytest.mark.parametrize(
    "degrees, form, idx_bin",
    [
        ((0, 0, 0), 'diag', None),
        ((2, 0, 2), 'diag', None),
        ((0, 0, 0), 'row', 0),
    ]
)
def test_compute_bispec(degrees, form, idx_bin,
                        test_data_catalogue, test_rand_catalogue,
                        test_binning_fourier,
                        test_paramset,
                        test_logger,
                        test_stats_dir):

    degree_tag = ''.join(map(str, degrees))
    if form == 'off-diag':
        form_tag = 'offdiag{:d}'.format(idx_bin)
    elif form == 'row':
        form_tag = 'row{:d}'.format(idx_bin)
    else:
        form_tag = form

    measurements = compute_bispec(
        test_data_catalogue, test_rand_catalogue,
        degrees=degrees,
        binning=test_binning_fourier,
        form=form,
        idx_bin=idx_bin,
        paramset=test_paramset,
        logger=test_logger
    )

    measurements_ext = np.loadtxt(
        test_stats_dir/f"bk{degree_tag}_{form_tag}_lpp.txt", unpack=True
    )

    assert np.allclose(measurements['k1_bin'], measurements_ext[0]), \
        "Measurement bins do not match."
    assert np.allclose(measurements['k1_eff'], measurements_ext[1]), \
        "Measured coordinates do not match."
    assert np.allclose(measurements['nmodes_1'], measurements_ext[2]), \
        "Measured mode counts do not match."
    assert np.allclose(measurements['k2_bin'], measurements_ext[3]), \
        "Measurement bins do not match."
    assert np.allclose(measurements['k2_eff'], measurements_ext[4]), \
        "Measured coordinates do not match."
    assert np.allclose(measurements['nmodes_2'], measurements_ext[5]), \
        "Measured mode counts do not match."
    assert np.allclose(
        measurements['bk_raw'],
        measurements_ext[-4] + 1j * measurements_ext[-3]
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['bk_shot'],
        measurements_ext[-2] + 1j * measurements_ext[-1],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."


@pytest.mark.slow
@pytest.mark.parametrize(
    "degrees, form, idx_bin",
    [
        ((0, 0, 0), 'diag', None),
        ((2, 0, 2), 'diag', None),
        ((0, 0, 0), 'row', 0),
    ]
)
def test_compute_3pcf(degrees, form, idx_bin,
                      test_data_catalogue, test_rand_catalogue,
                      test_binning_config,
                      test_paramset,
                      test_logger,
                      test_stats_dir):

    degree_tag = ''.join(map(str, degrees))
    if form == 'off-diag':
        form_tag = 'offdiag{:d}'.format(idx_bin)
    elif form == 'row':
        form_tag = 'row{:d}'.format(idx_bin)
    else:
        form_tag = form

    measurements = compute_3pcf(
        test_data_catalogue, test_rand_catalogue,
        degrees=degrees,
        binning=test_binning_config,
        form=form,
        idx_bin=idx_bin,
        paramset=test_paramset,
        logger=test_logger
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"zeta{degree_tag}_{form_tag}_lpp.txt", unpack=True
    )

    assert np.allclose(
        measurements['r1_bin'], measurements_ext[0]
    ), "Measurement bins do not match."
    assert np.allclose(
        measurements['r1_eff'], measurements_ext[1]
    ), "Measured coordinates do not match."
    assert np.allclose(
        measurements['npairs_1'], measurements_ext[2]
    ), "Measured pair counts do not match."
    assert np.allclose(
        measurements['r2_bin'], measurements_ext[3]
    ), "Measurement bins do not match."
    assert np.allclose(
        measurements['r2_eff'], measurements_ext[4]
    ), "Measured coordinates do not match."
    assert np.allclose(
        measurements['npairs_2'], measurements_ext[5]
    ), "Measured pair counts do not match."
    assert np.allclose(
        measurements['zeta_raw'],
        measurements_ext[-4] + 1j * measurements_ext[-3]
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['zeta_shot'],
        measurements_ext[-2] + 1j * measurements_ext[-1],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."


@pytest.mark.slow
@pytest.mark.parametrize(
    "degrees, form, idx_bin",
    [
        ((0, 0, 0), 'diag', None),
        ((2, 0, 2), 'diag', None),
        ((0, 0, 0), 'row', 0),
    ]
)
def test_compute_bispec_in_gpp_box(degrees, form, idx_bin,
                                   test_data_catalogue,
                                   test_binning_fourier,
                                   test_paramset,
                                   test_logger,
                                   test_stats_dir):

    degree_tag = ''.join(map(str, degrees))
    if form == 'off-diag':
        form_tag = 'offdiag{:d}'.format(idx_bin)
    elif form == 'row':
        form_tag = 'row{:d}'.format(idx_bin)
    else:
        form_tag = form

    measurements = compute_bispec_in_gpp_box(
        test_data_catalogue,
        degrees=degrees,
        binning=test_binning_fourier,
        form=form,
        idx_bin=idx_bin,
        paramset=test_paramset,
        logger=test_logger
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"bk{degree_tag}_{form_tag}_gpp.txt", unpack=True
    )

    assert np.allclose(measurements['k1_bin'], measurements_ext[0]), \
        "Measurement bins do not match."
    assert np.allclose(measurements['k1_eff'], measurements_ext[1]), \
        "Measured coordinates do not match."
    assert np.allclose(measurements['nmodes_1'], measurements_ext[2]), \
        "Measured mode counts do not match."
    assert np.allclose(measurements['k2_bin'], measurements_ext[3]), \
        "Measurement bins do not match."
    assert np.allclose(measurements['k2_eff'], measurements_ext[4]), \
        "Measured coordinates do not match."
    assert np.allclose(measurements['nmodes_2'], measurements_ext[5]), \
        "Measured mode counts do not match."
    assert np.allclose(
        measurements['bk_raw'],
        measurements_ext[-4] + 1j * measurements_ext[-3]
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['bk_shot'],
        measurements_ext[-2] + 1j * measurements_ext[-1],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."


@pytest.mark.slow
@pytest.mark.parametrize(
    "degrees, form, idx_bin",
    [
        ((0, 0, 0), 'diag', None),
        ((2, 0, 2), 'diag', None),
        ((0, 0, 0), 'row', 0),
    ]
)
def test_compute_3pcf_in_gpp_box(degrees, form, idx_bin,
                                 test_data_catalogue,
                                 test_binning_config,
                                 test_paramset,
                                 test_logger,
                                 test_stats_dir):

    degree_tag = ''.join(map(str, degrees))
    if form == 'off-diag':
        form_tag = 'offdiag{:d}'.format(idx_bin)
    elif form == 'row':
        form_tag = 'row{:d}'.format(idx_bin)
    else:
        form_tag = form

    measurements = compute_3pcf_in_gpp_box(
        test_data_catalogue,
        degrees=degrees,
        binning=test_binning_config,
        form=form,
        idx_bin=idx_bin,
        paramset=test_paramset,
        logger=test_logger
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"zeta{degree_tag}_{form_tag}_gpp.txt", unpack=True
    )

    assert np.allclose(
        measurements['r1_bin'], measurements_ext[0]
    ), "Measurement bins do not match."
    assert np.allclose(
        measurements['r1_eff'], measurements_ext[1]
    ), "Measured coordinates do not match."
    assert np.allclose(
        measurements['npairs_1'], measurements_ext[2]
    ), "Measured pair counts do not match."
    assert np.allclose(
        measurements['r2_bin'], measurements_ext[3]
    ), "Measurement bins do not match."
    assert np.allclose(
        measurements['r2_eff'], measurements_ext[4]
    ), "Measured coordinates do not match."
    assert np.allclose(
        measurements['npairs_2'], measurements_ext[5]
    ), "Measured pair counts do not match."
    assert np.allclose(
        measurements['zeta_raw'],
        measurements_ext[-4] + 1j * measurements_ext[-3]
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['zeta_shot'],
        measurements_ext[-2] + 1j * measurements_ext[-1],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."


@pytest.mark.slow
@pytest.mark.parametrize(
    "degrees, form, idx_bin",
    [
        ((0, 0, 0), 'diag', None),
        ((2, 0, 2), 'diag', None),
        ((0, 0, 0), 'row', 0),
    ]
)
def test_compute_3pcf_window(degrees, form, idx_bin,
                             test_rand_catalogue,
                             test_binning_config,
                             test_paramset,
                             test_logger,
                             test_stats_dir):

    degree_tag = ''.join(map(str, degrees))
    if form == 'off-diag':
        form_tag = 'offdiag{:d}'.format(idx_bin)
    elif form == 'row':
        form_tag = 'row{:d}'.format(idx_bin)
    else:
        form_tag = form

    measurements = compute_3pcf_window(
        test_rand_catalogue,
        degrees=degrees,
        binning=test_binning_config,
        form=form,
        idx_bin=idx_bin,
        paramset=test_paramset,
        logger=test_logger
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"zetaw{degree_tag}_{form_tag}.txt", unpack=True
    )

    assert np.allclose(
        measurements['r1_bin'], measurements_ext[0]
    ), "Measurement bins do not match."
    assert np.allclose(
        measurements['r1_eff'], measurements_ext[1]
    ), "Measured coordinates do not match."
    assert np.allclose(
        measurements['npairs_1'], measurements_ext[2]
    ), "Measured pair counts do not match."
    assert np.allclose(
        measurements['r2_bin'], measurements_ext[3]
    ), "Measurement bins do not match."
    assert np.allclose(
        measurements['r2_eff'], measurements_ext[4]
    ), "Measured coordinates do not match."
    assert np.allclose(
        measurements['npairs_2'], measurements_ext[5]
    ), "Measured pair counts do not match."
    assert np.allclose(
        measurements['zeta_raw'],
        measurements_ext[-4] + 1j * measurements_ext[-3]
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['zeta_shot'],
        measurements_ext[-2] + 1j * measurements_ext[-1],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."
