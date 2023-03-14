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


@pytest.mark.parametrize(
    "degree",
    [0, 2,]  # noqa: E231
)
def test_compute_powspec_in_gpp_box(degree, test_data_catalogue,
                                    test_paramset, test_stats_dir):
    measurements = compute_powspec_in_gpp_box(
        test_data_catalogue, degree=degree, paramset=test_paramset
    )
    measurements_ext = np.loadtxt(
        test_stats_dir/f"pk{degree}_sim_nbk.txt", unpack=True
    )
    assert np.allclose(measurements['kbin'], measurements_ext[0])
    assert np.allclose(measurements['keff'], measurements_ext[1])
    assert np.allclose(measurements['nmodes'], measurements_ext[2])
    assert np.allclose(
        measurements['pk_raw'],
        measurements_ext[3] + 1j * measurements_ext[4]
    )
    assert np.allclose(
        measurements['pk_shot'],
        measurements_ext[5] + 1j * measurements_ext[6],
        atol=1.e-6
    )
