"""Test :mod:`~triumvirate._arrayops`.

"""
import warnings

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal, assert_raises

from triumvirate._arrayops import (
    _check_1d_array,
    _check_2d_array,
    ShapeError,
    OrderError,
    MixedSignError,
    SpacingError,
    extrap_pad,
    extrap_lin,
    extrap_loglin,
    extrap_loglin_oscil,
    extrap2d_lin,
    extrap2d_loglin,
    extrap2d_pad,
    extrap2d_loglin_oscil,
    reshape_threept_datatab,
)


@pytest.mark.parametrize(
    "a, check_sorted, check_epsign, check_lin, check_loglin, exception",
    [
        # Test case 1: check sorted array
        (
            np.array([1., 3., 2.]), True, False, False, False, OrderError,
        ),
        # Test case 2: check mixed sign entries at endpoints
        (
            np.array([1., -2., 3.]), False, True, False, False, MixedSignError,
        ),
        # Test case 3: check linearly spaced array
        (
            np.array([1., 2., 4.]), False, False, True, False, SpacingError,
        ),
        # Test case 4: check log-linearly spaced array
        (
            np.array([1., 2., 3.]), False, False, False, True, SpacingError,
        ),
        # Test case 5: check array dimension
        (
            np.array([[1., 2., 3.]]*2), False, False, False, False, ShapeError,
        ),
    ],
)
def test__check_1d_array(a,
                         check_sorted, check_epsign, check_lin, check_loglin,
                         exception):
    if exception is not None:
        with pytest.raises(exception):
            a_out = _check_1d_array(
                a, check_sorted, check_epsign, check_lin, check_loglin
            )
        return

    a_out = _check_1d_array(
        a, check_sorted, check_epsign, check_lin, check_loglin
    )
    assert_allclose(a_out, a)


@pytest.mark.parametrize(
    "a, check_sorted, check_epsign, check_lin, check_loglin, exception",
    [
        # Test case 1: check sorted array
        (
            np.array([[1., 2., 3.], [2., 3., 2.]]),
            True, False, False, False, OrderError,
        ),
        # Test case 2: check mixed sign entries at endpoints
        (
            np.array([[1., 2., 3.], [2., -3., 2.]]),
            False, True, False, False, MixedSignError,
        ),
        # Test case 3: check linearly spaced array
        (
            np.array([[1., 2., 3.], [1., 2., 4.]]),
            False, False, True, False, SpacingError,
        ),
        # Test case 4: check log-linearly spaced array
        (
            np.array([[1., 2., 3.], [1., 2., 4.]]),
            False, False, False, True, SpacingError,
        ),
        # Test case 5: check array dimension
        (
            np.array([1, 2, 3]), False, False, False, False, ShapeError,
        ),
    ],
)
def test__check_2d_array(a,
                         check_sorted, check_epsign, check_lin, check_loglin,
                         exception):
    if exception is not None:
        with pytest.raises(exception):
            a_out = _check_2d_array(
                a, check_sorted, check_epsign, check_lin, check_loglin
            )
        return

    a_out = _check_2d_array(
        a, check_sorted, check_epsign, check_lin, check_loglin
    )
    assert_allclose(a_out, a)


@pytest.mark.parametrize(
    "a, n_ext, c_lower, c_upper, expected_output",
    [
        # Test case 1: extrapolate with positive padding
        (
            np.array([1, 2, 3,]), 2, -1, 1,
            np.array([-1, -1, 1, 2, 3, 1, 1,]),
        ),
        # Test case 2: extrapolate with no padding
        (
            np.array([1, 2, 3,]), 0, 0, 0,
            np.array([1, 2, 3,]),
        ),
    ],
)
def test_extrap_pad(a, n_ext, c_lower, c_upper, expected_output):
    a_out = extrap_pad(a, n_ext, c_lower, c_upper)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext, expected_output",
    [
        # Test case 1: extrapolate
        (
            np.array([1, 2, 2,]), 2,
            np.array([-1, 0, 1, 2, 2, 2, 2,]),
        ),
        # Test case 2: extrapolate not
        (
            np.array([1, 2, 3,]), 0,
            np.array([1, 2, 3,]),
        ),
    ],
)
def test_extrap_lin(a, n_ext, expected_output):
    a_out = extrap_lin(a, n_ext)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext, expected_output",
    [
        # Test case 1: extrapolate log-linearly
        (
            np.array([1., 2., 3.,]), 1,
            np.array([.5, 1., 2., 3., 4.5,]),
        ),
        # Test case 2: extrapolate log-linearly with negative values
        (
            np.array([-1., -2., -3.,]), 1,
            np.array([-.5, -1., -2., -3., -4.5,]),
        ),
        # Test case 3: extrapolate log-linearly with zero values
        (
            np.array([0, 1, 2,]), 2,
            np.array([0, 0, 0, 1, 2, 4, 8,]),
        ),
    ],
)
def test_extrap_loglin(a, n_ext, expected_output):
    a_out = extrap_loglin(a, n_ext)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext, expected_output",
    [
        # Test case 1: extrapolate log-linearly
        (
            np.array([1., 2., 4., 8., 16., 8., 4., 2., 1.,]), 1,
            np.array([.5, 1., 2., 4., 8., 16., 8., 4., 2., 1., .5,]),
        ),
        # Test case 2: extrapolate log-linearly with negative values
        (
            - np.array([1., 2., 4., 8., 16., 8., 4., 2., 1.,]), 1,
            - np.array([.5, 1., 2., 4., 8., 16., 8., 4., 2., 1., .5,]),
        ),
    ],
)
def test_extrap_loglin_oscil(a, n_ext, expected_output):
    a_out = extrap_loglin_oscil(a, n_ext)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext,"
    "crow_lower, crow_upper, ccol_lower, ccol_upper, expected_output",
    [
        # Test case 1: extrapolate with positive padding
        (
            np.array([[1, 2, 3], [4, 5, 6]]), 1,
            [-1, 1], [-1, -1], [-1, -2, -1], [1,  2, -1],
            np.array([
                [-1, -1, -2, -1, -1,],
                [-1,  1,  2,  3, -1,],
                [ 1,  4,  5,  6, -1,],  # noqa: E201
                [ 1,  1,  2, -1, -1,],  # noqa: E201
            ]),
        ),
    ],
)
def test_extrap2d_pad(a, n_ext, crow_lower, crow_upper, ccol_lower, ccol_upper,
                      expected_output):
    a_out = extrap2d_pad(
        a, n_ext, crow_lower, crow_upper, ccol_lower, ccol_upper
    )
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext, expected_output",
    [
        # Test case 1: extrapolate bilinearly
        (
            np.array([[1, 2, 3], [4, 5, 6]]), 1,
            np.array([
                [-3, -2, -1,  0,  1,],
                [ 0,  1,  2,  3,  4,],  # noqa: E201
                [ 3,  4,  5,  6,  7,],  # noqa: E201
                [ 6,  7,  8,  9, 10,],  # noqa: E201
            ]),
        ),
    ],
)
def test_extrap2d_lin(a, n_ext, expected_output):
    a_out = extrap2d_lin(a, n_ext)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext, expected_output",
    [
        # Test case 1: Extrapolate log-linearly
        (
            np.array([[1., 2., 4.,], [2., 4., 8.,]]), 1,
            np.array([
                [ .25,  .5, 1.,  2.,  4.,],  # noqa: E201
                [ .5 , 1. , 2.,  4.,  8.,],  # noqa: E201, E203
                [1.  , 2. , 4.,  8., 16.,],  # noqa: E201, E203
                [2.  , 4. , 8., 16., 32.,],  # noqa: E201, E203
            ]),
        ),
    ],
)
def test_extrap2d_loglin(a, n_ext, expected_output):
    a_out = extrap2d_loglin(a, n_ext)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "a, n_ext, expected_output",
    [
        # Test case 1: Extrapolate log-linearly with oscillatory behavior
        (
            np.array([
                [1. ,  2.,  4.,  8.,  4.,  2., 1.],  # noqa: E203
                [2. ,  4.,  8., 16.,  8.,  4., 2.],  # noqa: E203
                [4. ,  8., 16., 32., 16.,  8., 4.],  # noqa: E203
                [8. , 16., 32., 64., 32., 16., 8.],  # noqa: E203
                [4. ,  8., 16., 32., 16.,  8., 4.],  # noqa: E203
                [2. ,  4.,  8., 16.,  8.,  4., 2.],  # noqa: E203
                [1. ,  2.,  4.,  8.,  4.,  2., 1.],  # noqa: E203
            ]),
            1,
            np.array([
                [ .25,  .5,  1.,  2.,  4.,  2.,  1.,  .5,  .25,],  # noqa: E201
                [ .5 , 1. ,  2.,  4.,  8.,  4.,  2., 1. ,  .5 ,],  # noqa: E201, E203, E501
                [1.  , 2. ,  4.,  8., 16.,  8.,  4., 2. , 1.  ,],  # noqa: E201, E203, E501
                [2.  , 4. ,  8., 16., 32., 16.,  8., 4. , 2.  ,],  # noqa: E201, E203, E501
                [4.  , 8. , 16., 32., 64., 32., 16., 8. , 4.  ,],  # noqa: E201, E203, E501
                [2.  , 4. ,  8., 16., 32., 16.,  8., 4. , 2.  ,],  # noqa: E201, E203, E501
                [1.  , 2. ,  4.,  8., 16.,  8.,  4., 2. , 1.  ,],  # noqa: E201, E203, E501
                [ .5 , 1. ,  2.,  4.,  8.,  4.,  2., 1. ,  .5 ,],  # noqa: E201, E203, E501
                [ .25,  .5,  1.,  2.,  4.,  2.,  1.,  .5,  .25,],  # noqa: E201
            ])
        )
    ]
)
def test_extrap2d_loglin_oscil(a, n_ext, expected_output):
    warnings.filterwarnings(
        "ignore", category=UserWarning,
        message="Array length is less than 8. "
                "Extrapolation may not be accurate."
    )
    a_out = extrap2d_loglin_oscil(a, n_ext)
    assert_allclose(a_out, expected_output)


@pytest.mark.parametrize(
    "paired_coords, flat_multipole, shape,"
    "expected_common_coords, expected_multipole",
    [
        # Test case 1: Valid input
        (
            np.array([[0, 0], [0, 1], [1, 0]]),
            np.array([1, 2, 3]),
            'triu',
            np.array([0, 1]),
            np.array([[1, 2,], [2, 3],]),
        ),
        # Test case 2: Invalid input - paired_coords is not a 2-d array
        (
            np.array([1, 2, 3]),
            np.array([1, 2, 3]),
            None,
            None,
            None,
        ),
        # Test case 3: Invalid input - unmatching dimensions
        (
            np.array([[1, 2], [2, 3], [3, 4]]),
            np.array([1, 2, 3, 4]),
            None,
            None,
            None,
        ),
        # Test case 4: Invalid input - not an upper triangular matrix
        (
            np.array([[1, 2], [2, 3], [3, 4]]),
            np.array([1, 2, 3, 4, 5, 6]),
            None,
            None,
            None,
        ),
        # Test case 5: Valid input
        (
            np.array([[0, 0], [0, 1], [1, 0], [1, 1]]),
            np.array([1, 2, 3, 4]),
            'full',
            np.array([0, 1]),
            np.array([[1, 2,], [3, 4],])
        ),
    ]
)
def test_reshape_threept_datatab(paired_coords, flat_multipole, shape,
                                 expected_common_coords, expected_multipole):
    if expected_common_coords is not None and expected_multipole is not None:
        common_coords, multipole = reshape_threept_datatab(
            paired_coords, flat_multipole, shape=shape
        )
        assert_array_equal(common_coords, expected_common_coords)
        assert_array_equal(multipole, expected_multipole)
    else:
        assert_raises(
            ValueError, reshape_threept_datatab, paired_coords, flat_multipole
        )
