"""Test :mod:`~triumvirate.transforms`.

"""
import numpy as np
import pytest

from triumvirate._fftlog import HankelTransform
from triumvirate.transforms import SphericalBesselTransform


def fa_analytic(r, mu=0.):
    return r ** (mu + 1.) * np.exp(-r*r/2.)


def ga_analytic(k, mu=0.):
    return k ** (mu + 1.) * np.exp(-k*k/2.)


def fb_analytic(r):
    return r * (1 + r**2)**(-3/2)


def gb_analytic(k):
    return k * np.exp(-k)


def fsph_analytic(r, ell=0):
    return r ** ell * np.exp(-r*r/2.)


def gsph_analytic(k, ell=0):
    return (2*np.pi)**(3./2) * k ** ell * np.exp(-k*k/2.)


@pytest.mark.parametrize(
    "order, bias, nsamp, lg_range, test_case",
    [
        (2, 0, 2**10, [-5., 5.], 'a'),
        (0, 0, 2**11, [-5., 5.], 'b'),
    ]
)
def test_hankeltransform(order, bias, nsamp, lg_range, test_case):

    if test_case not in ['a', 'b']:
        raise ValueError("Invalid test case.")

    pre_sampts = np.logspace(*lg_range, nsamp, base=10, endpoint=False)
    transformer = HankelTransform(
        order, bias, pre_sampts, kr_c=1., lowring=True
    )
    post_sampts = np.asarray(transformer._post_sampts)

    if test_case == 'a':
        pre_samples = fa_analytic(pre_sampts, mu=order)
    if test_case == 'b':
        pre_samples = fb_analytic(pre_sampts)

    _, post_samples = transformer.transform(pre_samples)

    if test_case == 'a':
        post_results = ga_analytic(post_sampts, mu=order)
    if test_case == 'b':
        post_results = gb_analytic(post_sampts)

    accu_range = (-3. < np.log10(post_sampts)) & (np.log10(post_sampts) < 0.5)
    assert np.allclose(
        post_samples[accu_range], post_results[accu_range],
        rtol=1.e-3
    ), "Transform not accurate in the expected range."


@pytest.mark.parametrize(
    "degree, bias, nsamp, lg_range",
    [
        (0, 0, 2**10, [-5., 5.]),
        (2, 0, 2**11, [-5., 5.]),
    ]
)
def test_sphericalbesseltransform(degree, bias, nsamp, lg_range):

    pre_sampts = np.logspace(*lg_range, nsamp, base=10, endpoint=False)
    pre_samples = fsph_analytic(pre_sampts, ell=degree)

    transformer = SphericalBesselTransform(
        degree, bias, pre_sampts, lowring=True
    )
    assert transformer.size == nsamp, "Sample size not as expected."

    # Test only
    # :meth:`~triumvirate.transforms.SphericalBesselTransform.transform`
    # as a proxy for
    # :class:`~triumvirate.transforms.SphericalBesselTransform.transform_cosmo_multipoles`.
    post_sampts, post_samples = transformer.transform(pre_samples)
    post_results = gsph_analytic(post_sampts, ell=degree)

    accu_range = (-3.5 < np.log10(post_sampts)) & (np.log10(post_sampts) < 0.5)
    assert np.allclose(
        post_samples[accu_range], post_results[accu_range],
        rtol=1.e-3
    ), "Transform not accurate in the expected range."
