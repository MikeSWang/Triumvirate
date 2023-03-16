"""Test :mod:`~triumvirate.dataobjs`.

"""
import numpy as np
import pytest

from triumvirate.dataobjs import Binning


@pytest.fixture
def default_binning(valid_paramset):
    binning = Binning.from_parameter_set(valid_paramset)
    return binning


@pytest.mark.parametrize(
    "space, scheme, bin_min, bin_max, num_bins",
    [
        ('config', 'logpad', 1.e-1, 1.e4, 20),
        ('fourier', 'lin', None, None, None),
    ]
)
def test_Binning___cinit__(space, scheme, bin_min, bin_max, num_bins):
    Binning(space, scheme, bin_min=bin_min, bin_max=bin_max, num_bins=num_bins)


def test_Binning_from_parameter_set(valid_paramset):
    Binning.from_parameter_set(valid_paramset)


@pytest.mark.parametrize(
    "bin_min, bin_max, num_bins",
    [
        (1.e-1, 1.e4, 20),
        (5.e-3, 1.e-1+5.e-3, 10),
    ]
)
def test_Binning_set_bins(bin_min, bin_max, num_bins, default_binning):
    default_binning.set_bins(bin_min, bin_max, num_bins)
    assert np.allclose(
        default_binning.bin_edges,
        np.linspace(bin_min, bin_max, num=num_bins+1)
    ), "Bin resetting is not correct."


@pytest.mark.parametrize(
    "boxsize, ngrid",
    [
        (1000., [16,]*3),  # noqa: E231
        ([1000.,]*3, 16),  # noqa: E231
    ]
)
def test_Binning_set_grid_based_bins(boxsize, ngrid, default_binning):

    default_binning.set_grid_based_bins(boxsize, ngrid)

    bin_width = 2*np.pi / np.max(boxsize)

    bin_min = 0.
    bin_max = np.pi * np.min(ngrid) / np.max(boxsize) + bin_width/2.
    num_bins = np.min(ngrid) // 2

    assert np.allclose(
        default_binning.bin_edges,
        np.linspace(bin_min, bin_max, num=num_bins+1)
    ), "Grid-based bin resetting is not correct."


@pytest.mark.parametrize(
    "bin_edges",
    [
        np.sort(np.random.rand(10)),
    ]
)
def test_Binning_set_custom_bins(bin_edges, default_binning):

    default_binning.set_custom_bins(bin_edges)

    bin_centres = (bin_edges[1:] + bin_edges[:-1]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]

    assert np.allclose(
        default_binning.bin_edges,
        bin_edges
    ), "Customised bin resetting is not correct in bin edges."
    assert np.allclose(
        default_binning.bin_centres,
        bin_centres
    ), "Customised bin resetting is not correct in bin centres."
    assert np.allclose(
        default_binning.bin_widths,
        bin_widths
    ), "Customised bin resetting is not correct in bin widths."
