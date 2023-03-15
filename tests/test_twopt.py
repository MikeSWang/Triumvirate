"""Test :mod:`~triumvirate.twopt`.

"""
import numpy as np
import pytest

from triumvirate.twopt import (
    compute_corrfunc,
    compute_corrfunc_in_gpp_box,
    compute_corrfunc_window,
    compute_powspec,
    compute_powspec_in_gpp_box
)


@pytest.mark.slow
@pytest.mark.parametrize(
    "degree",
    [0, 2,]  # noqa: E231
)
def test_compute_powspec_in_gpp_box(degree, test_data_catalogue,
                                    test_paramset, test_logger,
                                    test_stats_dir):
    measurements = compute_powspec_in_gpp_box(
        test_data_catalogue, degree=degree, paramset=test_paramset,
        logger=test_logger
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"pk{degree}_gpp_nbk.txt", unpack=True
    )
    assert np.allclose(measurements['kbin'], measurements_ext[0]), \
        "Measurement bins do not match."
    assert np.allclose(measurements['keff'], measurements_ext[1]), \
        "Measured coordinates do not match."
    assert np.allclose(measurements['nmodes'], measurements_ext[2]), \
        "Measured mode counts do not match."
    assert np.allclose(
        measurements['pk_raw'],
        measurements_ext[3] + 1j * measurements_ext[4],
        rtol=1.e-3
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['pk_shot'],
        measurements_ext[5] + 1j * measurements_ext[6],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."


@pytest.mark.slow
@pytest.mark.parametrize(
    "degree",
    [0, 2,]  # noqa: E231
)
def test_compute_powspec(degree, test_data_catalogue, test_rand_catalogue,
                         test_paramset, test_logger, test_stats_dir):
    measurements = compute_powspec(
        test_data_catalogue, test_rand_catalogue,
        degree=degree, paramset=test_paramset,
        logger=test_logger
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"pk{degree}_lpp_nbk.txt", unpack=True
    )
    assert np.allclose(measurements['kbin'], measurements_ext[0]), \
        "Measurement bins do not match."
    assert np.allclose(measurements['keff'], measurements_ext[1]), \
        "Measured coordinates do not match."
    assert np.allclose(measurements['nmodes'], measurements_ext[2]), \
        "Measured mode counts do not match."
    assert np.allclose(
        measurements['pk_raw'],
        measurements_ext[3] + 1j * measurements_ext[4]
    ), "Measured raw statistics do not match."
    assert np.allclose(
        measurements['pk_shot'],
        measurements_ext[5] + 1j * measurements_ext[6],
        atol=1.e-6
    ), "Measured shot noise contributions do not match."
