"""Test :mod:`~triumvirate.transforms`.

"""
import numpy as np
import pytest
from scipy.special import gamma

from triumvirate._fftlog import HankelTransform
from triumvirate.transforms import (
    DoubleSphericalBesselTransform,
    SphericalBesselTransform,
    resample_lglin,
    resample_lin,
)


@pytest.fixture
def hankel_func_pair(request):
    """Analytical function pair for testing the Hankel transform.

    """
    if request.param == 'sym':
        return (
            lambda r, mu=0., **kwargs: r ** (mu + 1.) * np.exp(-r*r/2.),
            lambda k, mu=0., **kwargs: k ** (mu + 1.) * np.exp(-k*k/2.),
        )
    if request.param == 'asym':
        return (
            lambda r, mu=0., **kwargs: r * (1. + r**2)**(-3./2),
            lambda k, mu=0., **kwargs: k * np.exp(-k),
        )
    raise ValueError("Unknown analytical function pair.")


@pytest.fixture
def sj_func_pair(request):
    """Analytical function pair for testing the spherical
    Bessel transform.

    """
    if request.param == 'sym':
        return (
            lambda r, ell=0, **kwargs:
                r ** ell * np.exp(-r*r/2.),
            lambda k, ell=0, **kwargs:
                k ** ell * np.exp(-k*k/2.) * (2*np.pi)**(3./2),
        )
    if request.param == 'asym':
        return (
            lambda r, s=1, **kwargs: r**(s - 3),
            lambda k, s=1, ell=0, **kwargs:
                np.pi**(3./2)
                * gamma(s/2. + ell/2.) / gamma(3./2 - s/2. + ell/2.)
                * (2./k)**s,
        )


@pytest.mark.parametrize(
    "order, bias, nsamp, lgrange_samp, extrap, hankel_func_pair, lgrange_test",
    [
        (2, 0, 2**10, [-5., 5.], 0, 'sym', [-3., .5]),
        (0, 0, 2**11, [-5., 5.], 4, 'asym', [-3., .5]),
    ],
    indirect=['hankel_func_pair',]
)
def test_hankeltransform(order, bias, nsamp, lgrange_samp, extrap,
                         hankel_func_pair, lgrange_test):

    f, g = hankel_func_pair

    pre_sampts = np.logspace(*lgrange_samp, nsamp, base=10, endpoint=False)
    transformer = HankelTransform(
        order, bias, pre_sampts, kr_c=1., lowring=True, extrap=extrap
    )
    post_sampts = np.asarray(transformer._post_sampts)

    pre_samples = f(pre_sampts, mu=order)
    _, post_samples = transformer.transform(pre_samples)
    post_results = g(post_sampts, mu=order)

    test_range = (lgrange_test[0] < np.log10(post_sampts)) \
        & (np.log10(post_sampts) < lgrange_test[-1])
    assert np.allclose(
        post_samples[test_range].real, post_results[test_range],
        rtol=1.e-3
    ), "Transform not accurate in the expected range."


@pytest.mark.parametrize(
    "degree, bias, nsamp, lgrange_samp, extrap, sj_func_pair, lgrange_test",
    [
        (0, 0, 2**10, [-5., 5.], 0, 'sym', [-3.5, .5]),
        (2, 0, 2**11, [-5., 5.], 4, 'asym', [-1., 3.]),
    ],
    indirect=['sj_func_pair',]
)
def test_sphericalbesseltransform(degree, bias, nsamp, lgrange_samp, extrap,
                                  sj_func_pair, lgrange_test):

    f, g = sj_func_pair

    pre_sampts = np.logspace(*lgrange_samp, nsamp, base=10, endpoint=False)
    pre_samples = f(pre_sampts, ell=degree)

    transformer = SphericalBesselTransform(
        degree, bias, pre_sampts, lowring=True, extrap=extrap
    )

    # Test only
    # :meth:`~triumvirate.transforms.SphericalBesselTransform.transform`
    # as a proxy for
    # :class:`~triumvirate.transforms.SphericalBesselTransform.transform_cosmo_multipoles`.
    post_sampts, post_samples = transformer.transform(pre_samples)
    post_results = g(post_sampts, ell=degree)

    test_range = (lgrange_test[0] < np.log10(post_sampts)) \
        & (np.log10(post_sampts) < lgrange_test[-1])
    assert np.allclose(
        post_samples[test_range].real, post_results[test_range],
        rtol=1.e-3
    ), "Transform not accurate in the expected range."


@pytest.mark.parametrize(
    "degrees, biases, nsamp, lgrange_samp, extrap, sj_func_pair, lgrange_test",
    [
        ((0, 0), (0, 0), 2**10, [-5., 5.], 0, 'sym', [-3., 0.]),
        ((2, 2), (0, 0), 2**11, [-5., 5.], 4, 'asym', [-.5, 3.]),
    ],
    indirect=['sj_func_pair',]
)
def test_doublesphericalbesseltransform(degrees, biases, nsamp, lgrange_samp,
                                        extrap, sj_func_pair, lgrange_test):

    f, g = sj_func_pair

    pre_sampts = np.logspace(*lgrange_samp, nsamp, base=10, endpoint=False)
    pre_samples_1, pre_samples_2 = (
        f(pre_sampts, ell=ell) for ell in degrees
    )
    pre_samples = pre_samples_1[:, None] * pre_samples_2[None, :]

    transformer = DoubleSphericalBesselTransform(
        degrees, biases, pre_sampts, lowring=True, extrap=extrap
    )

    # Test only
    # :meth:`~triumvirate.transforms.DoubleSphericalBesselTransform.transform`
    # as a proxy for
    # :class:`~triumvirate.transforms.DoubleSphericalBesselTransform.transform_cosmo_multipoles`.
    post_sampts, post_samples = transformer.transform(pre_samples)
    post_results = np.multiply(*[
        g(post_sampts_i, ell=ell)
        for (post_sampts_i, ell) in zip(post_sampts, degrees)
    ])

    test_range = np.logical_and.reduce([
        np.log10(post_sampts[0]) > lgrange_test[0],
        np.log10(post_sampts[0]) < lgrange_test[-1],
        np.log10(post_sampts[-1]) > lgrange_test[0],
        np.log10(post_sampts[-1]) < lgrange_test[-1],
    ])
    assert np.allclose(
        post_samples[test_range].real, post_results[test_range],
        rtol=1.e-3
    ), "Transform not accurate in the expected range."


@pytest.mark.parametrize(
    "sampts, size, func, ndim",
    [
        (np.arange(1., 11.), None, lambda x: x**2, 1),
        (
            np.arange(1., 51.) * 1.e-3, 10,
            lambda x, y: np.exp(-x*np.sqrt(y)),
            2
        ),
    ]
)
def test_resample_lglin(sampts, size, func, ndim):

    if ndim == 1:
        samples = func(sampts)
    elif ndim == 2:
        samples = func(*np.meshgrid(sampts, sampts, indexing='ij'))
    else:
        raise ValueError("Resampling only supports 1- or 2-d arrays.")

    resampts, resamples = resample_lglin(sampts, samples, size=size)

    if ndim == 1:
        resamples_exact = func(resampts)
    if ndim == 2:
        resamples_exact = func(
            *np.meshgrid(resampts[0], resampts[-1], indexing='ij')
        )

    assert np.allclose(resamples, resamples_exact), "Resampling not accurate."


@pytest.mark.parametrize(
    "sampts, size, func, ndim",
    [
        (np.arange(1., 11.), None, lambda x: x**2, 1),
        (np.arange(1., 51.) * 1.e-3, 10, lambda x, y: x + y, 2),
    ]
)
def test_resample_lin(sampts, size, func, ndim):

    if ndim == 1:
        samples = func(sampts)
    elif ndim == 2:
        samples = func(*np.meshgrid(sampts, sampts, indexing='ij'))
    else:
        raise ValueError("Resampling only supports 1- or 2-d arrays.")

    resampts, resamples = resample_lin(sampts, samples, size=size)

    if ndim == 1:
        resamples_exact = func(resampts)
    if ndim == 2:
        resamples_exact = func(
            *np.meshgrid(resampts[0], resampts[-1], indexing='ij')
        )

    assert np.allclose(resamples, resamples_exact), "Resampling not accurate."
