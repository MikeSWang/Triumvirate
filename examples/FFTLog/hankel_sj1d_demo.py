"""Demonstrate the use of the FFTLog algorithm for Hankel-like transforms.

.. attention::

    Requires Python >=3.11 for the following language features:

    - Match statement (PEP 634);
    - Self type (PEP 673).

"""
import argparse
import os.path as osp
import sys
import warnings
from collections.abc import Callable
from functools import wraps
from time import time
from typing import Literal, Self

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.special import gamma

from triumvirate._fftlog import HankelTransform
from triumvirate.transforms import SphericalBesselTransform

try:
    import mcfit
except ImportError:
    mcfit = None
    warnings.warn(
        "Could not import `mcfit` package. "
        "Some of the test cases will not be available."
    )

try:
    import hankl
except ImportError:
    hankl = None
    warnings.warn(
        "Could not import `hankl` package. "
        "Some of the test cases will not be available."
    )


def parse_parameters() -> argparse.Namespace:
    """Initialise the command-line argument parser.

    Returns
    -------
    pars
        Parsed command-line parameters.

    """
    parser = argparse.ArgumentParser(
        description="FFTLog Hankel transform demo."
    )

    parser.add_argument(
        '-t', "--test-case", type=str,
        choices=['hankel-sym', 'hankel-asym', 'sj-sym', 'sj-asym', 'sj-cosmo'],
        required=True,
        help="test case"
    )
    parser.add_argument(
        '--order', type=int, default=0,
        help="order of the Hankel transform (default is %(default)d)"
    )
    parser.add_argument(
        '--degree', type=int, default=0,
        help="degree of the spherical Bessel transform "
             "(default is %(default)d)"
    )
    parser.add_argument(
        '--extrap', type=int, nargs='?', const=4, default=0,
        help="extrapolation option for smoothness "
             "(default is %(default)d if disabled or %(const)d if enabled)"
    )
    parser.add_argument(
        '--lgcentre', type=float, default=0.,
        help="logarithmic centre for the sample points "
             "(default is %(default).1f)"
    )
    parser.add_argument(
        '--shorten', type=float, default=1.,
        help="shorten the range of sample points by this factor "
             "(default is %(default).1f)"
    )
    parser.add_argument(
        '--show-diff', action='store_true',
        help="show the difference between the transform results"
    )
    parser.add_argument(
        '--time', action='store_true',
        help="time the transform without showing a plotted figure"
    )
    parser.add_argument(
        '--threaded', action='store_true',
        help="use multi-threaded FFTLog transform"
    )

    pars = parser.parse_args()

    # HACK: Override the degree to be 0 to match the external samples.
    if pars.test_case == 'sj-cosmo' and pars.degree != 0:
        pars.degree = 0
        warnings.warn(
            "Spherical Bessel transform degree overridden to 0 "
            "when external cosmological multipole samples are referenced "
            "as they are computed for the monopole only."
        )

    return pars


def get_sampts(centre: float, shorten: float) -> np.ndarray:
    """Get the log-range of sample points.

    Parameters
    ----------
    centre
        Logarithmic centre for the sample points.
    shorten
        Shorten the symmetric range of sample points by this factor.

    Returns
    -------
    Tuple[float, float]
        Log-range of sample points.

    """
    return np.logspace(
        *((val - centre) / shorten for val in LGRANGE),
        round(NSAMP / shorten),
        base=10
    )


def get_analy_func_pair(
    testcase: Literal['hankel-sym', 'hankel-asym', 'sj-sym', 'sj-asym']
) -> tuple[tuple[Callable, str], tuple[Callable, str]]:
    """Get pre- and post-transform analytical functions for the given
    test case.

    Parameters
    ----------
    testcase
        Test cases:

        - 'hankel-sym' for symmetric-form functions for the
          Hankel transform;
        - 'hankel-asym' for asymmetric-form functions for the
          Hankel transform;
        - 'sj-sym' for symmetric-form functions for the spherical
          Bessel transform;
        - 'sj-asym' for asymmetric-form functions for the spherical
          Bessel transform.

    Returns
    -------
    Pre- and post-transform functions with their TeX strings.

    """  # numpydoc ignore=RT03
    match testcase:
        case 'hankel-sym':
            return ((
                lambda r, mu=0.: r**(mu + 1.) * np.exp(-r*r/2.),
                r"$f(r) = r^{\mu + 1} \mathrm{e}^{-r^2 / 2}$"
            ), (
                lambda k, mu=0.: k**(mu + 1.) * np.exp(-k*k/2.),
                r"$g(k) = k^{\mu + 1} \mathrm{e}^{-k^2 / 2}$"
            ))
        case 'hankel-asym':
            return ((
                lambda r, **kwargs: r * (1 + r**2)**(-3/2),
                r"$f(r) = r (1 + r^2)^{-3/2}$",
            ), (
                lambda k, **kwargs: k * np.exp(-k),
                r"$g(k) = k \mathrm{e}^{-k}$"
            ))
        case 'sj-sym':
            return ((
                lambda r, ell=0: r**ell * np.exp(-r*r/2.),
                r"$f(r) = r^\ell \mathrm{e}^{-r^2 / 2}$"
            ), (
                lambda k, ell=0: (2*np.pi)**(3./2) * k**ell * np.exp(-k*k/2.),
                r"$g(k) = (2\pi)^{3/2} k^\ell \mathrm{e}^{-k^2 / 2}$"
            ))
        case 'sj-asym':
            return ((
                lambda r, s=1, **kwargs: r**(s - 3),
                r"$f(r) = r^{s - 3}$"
            ), (
                lambda k, s=1, ell=0, **kwargs:
                    np.pi**(3./2)
                    * gamma(s/2. + ell/2.) / gamma(3./2 - s/2. + ell/2.)
                    * (2./k)**s,
                r"$g(k) = \pi^{3/2} "
                r"\frac{\Gamma(s/2 + \ell/2)}{\Gamma(3/2 - s/2 + \ell/2)} "
                r"\left(\frac{2}{k}\right)^s$"
            ))
        case _:
            raise ValueError(
                f"No analytical functions for test case: {testcase}."
            )


def get_external_samples(
    presamp_file: str,
    postsamp_file: str,
    samp_dir: str = "examples/FFTLog/storage/input/samps"
) -> tuple[tuple[np.ndarray, np.ndarray], tuple[np.ndarray, np.ndarray]]:
    """Get pre- and post-transform samples from external files.

    Parameters
    ----------
    presamp_file
        File name for pre-transform samples.
    postsamp_file
        File name for post-transform samples.
    samp_dir
        Directory containing the sample files.

    Returns
    -------
    presamples, postsamples
        Pre- and post-transform samples.

    """
    presamples = np.loadtxt(osp.join(samp_dir, presamp_file), unpack=True)
    postsamples = np.loadtxt(osp.join(samp_dir, postsamp_file), unpack=True)

    return presamples, postsamples


def translate_extrap_mcfit(extrap_opt: int) -> bool | Literal['const']:
    """Translate the extrapolation option for `mcfit`.

    Parameters
    ----------
    extrap_opt : int
        Extrapolation option.

    Returns
    -------
    bool or str
        Extrapolation option for `mcfit`.

    Raises
    ------
    ValueError
        When `extrap_opt` is not supported by `mcfit`.

    """
    if extrap_opt == 0:
        return False
    if extrap_opt == 1:
        raise ValueError("Zero padding not available for `mcfit`.")
    if extrap_opt == 1:
        return 'const'
    if extrap_opt == 3:
        raise ValueError("Linear extrapolation not available for `mcfit`.")
    if extrap_opt == 4:
        return True
    if extrap_opt == 5:
        raise ValueError(
            "Log-linear extrapolation with oscillations "
            "not available for `mcfit`."
        )
    raise ValueError(f"Unrecognised extrapolation option: {extrap_opt}.")


class comparison_plot:
    """Comparison plot context manager.

    Parameters
    ----------
    comp_type
        Type of comparison plot:

        - 'analytic' for analytical functions;
        - 'samples' for external samples.

    plot_diff
        If `True` (default), plot the difference between the
        transform results.
    xlim, ylim
        Limits for the x- and y-axes.
    xlabel, ylabel
        Labels for the x- and y-axes.
    title
        Title for the plot.
    timed
        If `True`, time the transform and skip showing the figure.

    """

    def __init__(
        self,
        comp_type: Literal['analytic', 'samples'],
        plot_diff: bool = True,
        xlim: tuple[float, float] | None = None,
        ylim: tuple[float, float] | None = None,
        xlabel: str | None = None,
        ylabel: str | None = None,
        title: str | None = None,
        timed: bool = False,
    ) -> None:
        # Store figure attributes.
        self._comp_type = comp_type
        self._flag_diff = plot_diff
        self._xlim = xlim
        self._ylim = ylim
        self._xlabel = xlabel
        self._ylabel = ylabel
        self._title = title
        self._timed = timed

    def __enter__(self) -> Self:
        # Set up the figure.
        self._fig = plt.figure()

        if self._flag_diff:
            self._ax_comp = plt.subplot2grid((4, 1), (0, 0), rowspan=3)
            self._ax_diff = plt.subplot2grid(
                (4, 1), (3, 0), sharex=self._ax_comp
            )
        else:
            self._ax_comp = plt.subplot2grid((1, 1), (0, 0))
            self._ax_diff = None

        return self

    def __exit__(self, exc_type, exc_value, traceback):
        # Apply figure attributes.
        if self._xlim:
            self._ax_comp.set_xlim(*self._xlim)
        if self._ylim:
            self._ax_comp.set_ylim(*self._ylim)
        if self._flag_diff:  # in %
            self._ax_diff.axhline(-1.e-1, c='k', alpha=.5, ls='--')
            self._ax_diff.axhline(1.e-1, c='k', alpha=.5, ls='--')
            self._ax_diff.set_ylim(-1.75e-1, 1.75e-1)

        self._ax_comp.set_xscale('log')
        if self._comp_type == 'analytic':
            self._ax_comp.set_yscale('log')

        if not self._xlabel:
            if self._comp_type == 'analytic':
                self._xlabel = r"$k$"
            if self._comp_type == 'samples':
                self._xlabel = r"$r$"
        if self._flag_diff:
            self._ax_comp.tick_params(labelbottom=False)
            self._ax_diff.set_xlabel(self._xlabel)
        else:
            self._ax_comp.set_xlabel(self._xlabel)

        self._ax_comp.set_ylabel(self._ylabel)
        if self._flag_diff:
            self._ax_diff.set_ylabel("diff. (%)")

        self._ax_comp.legend()
        self._ax_comp.set_title(self._title)

        plt.subplots_adjust(wspace=0, hspace=0)
        plt.show(block=not self._timed)

    def __call__(self, runner: Callable) -> Callable:
        @wraps(runner)
        def context_wrapper(*args, **kwargs):
            with self:
                return runner(*args, canvas=self, **kwargs)
        return context_wrapper


def get_testcase(pars: argparse.Namespace) -> None:
    """Get the test case.

    Parameters
    ----------
    pars
        Parsed command-line parameters.

    """  # numpydoc ignore=RT01
    testcase = pars.test_case
    showdiff = pars.show_diff
    timed = pars.time

    @comparison_plot(
        comp_type='analytic',
        plot_diff=showdiff,
        xlim=(1.e-5, 1.e1),
        ylim=(1.e-6, 1.e0),
        ylabel=r"$g(k)$",
        timed=timed,
    )
    def run_hankel(canvas: comparison_plot) -> None:
        # Get analytical functions.
        (f, fstr), (g, gstr) = get_analy_func_pair(testcase)

        # Set up samples.
        r = get_sampts(pars.lgcentre, pars.shorten)
        fr = f(r, mu=pars.order)

        # Compute the Hankel transforms.
        t_start = time()
        k_fftlog, gk_fftlog = HankelTransform(
            pars.order, BIAS, r, PIVOT, LOWRING,
            extrap=pars.extrap, extrap_exp=EXPAND,
            threaded=pars.threaded
        ).transform(fr + 0.j)
        t_end = time()
        if timed:
            print(f"FFTLog transform time: {t_end - t_start:.3f} s")

        k_fftlog_analy, gk_fftlog_analy = k_fftlog, g(k_fftlog, mu=pars.order)

        if (_mcfit := mcfit) is not None:
            try:
                extrap_mcfit = translate_extrap_mcfit(pars.extrap)
            except ValueError as err:
                _mcfit = None
                warnings.warn(str(err))
            else:
                k_mcfit, gk_mcfit = mcfit.Hankel(
                    r, lowring=LOWRING,
                    N=EXPAND*1j if extrap_mcfit else NSAMP
                )(fr / r, extrap=extrap_mcfit)
                gk_mcfit *= k_mcfit
                gk_mcfit_analy = g(k_mcfit, mu=pars.order)
        if hankl is not None:
            k_hankl, gk_hankl = hankl.FFTLog(
                r, fr, mu=pars.order, q=BIAS, xy=PIVOT,
                lowring=LOWRING, ext=pars.extrap
            )
            gk_hankl_analy = g(k_hankl, mu=pars.order)

        # Calculate the differences.
        if pars.show_diff:
            with np.errstate(divide='ignore', invalid='ignore'):
                dgk_fftlog = gk_fftlog.real / gk_fftlog_analy - 1.
                if _mcfit is not None:
                    dgk_mcfit = gk_mcfit / gk_mcfit_analy - 1.
                if hankl is not None:
                    dgk_hankl = gk_hankl / gk_hankl_analy - 1.

        # Plot the results.
        plot_fftlog = canvas._ax_comp.plot(
            k_fftlog, gk_fftlog.real, ls='-', label='FFTLog'
        )
        if _mcfit is not None:
            plot_mcfit = canvas._ax_comp.plot(
                k_mcfit, gk_mcfit, ls='--', label='mcfit'
            )
        if hankl is not None:
            plot_hankl = canvas._ax_comp.plot(
                k_hankl, gk_hankl, ls='--', label='hankl'
            )
        _ = canvas._ax_comp.plot(
            k_fftlog_analy, gk_fftlog_analy, ls=':', label='analytical'
        )

        if pars.show_diff:
            canvas._ax_diff.plot(
                k_fftlog, 100 * dgk_fftlog,
                c=plot_fftlog[0].get_color(), ls='-', label='FFTLog'
            )
            if _mcfit is not None:
                canvas._ax_diff.plot(
                    k_mcfit, 100 * dgk_mcfit,
                    c=plot_mcfit[0].get_color(), ls='--', label='mcfit'
                )
            if hankl is not None:
                canvas._ax_diff.plot(
                    k_hankl, 100 * dgk_hankl,
                    c=plot_hankl[0].get_color(), ls='--', label='hankl'
                )

        canvas._title = (
            fstr + r"$\,$, " + gstr + fr" ($\mu = {pars.order}$)"
        )

    @comparison_plot(
        comp_type='analytic',
        plot_diff=showdiff,
        xlim={
            'sj-sym': (1.e-5, 1.e1),
            'sj-asym': (1.e-1, 1.e3),
        }.get(testcase),
        ylim={
            'sj-sym': (1.e-10, 1.e2),
            'sj-asym': (1.e-0, 1.e5),
        }.get(testcase),
        ylabel=r"$k^2 g(k)$",
        timed=timed,
    )
    def run_sj(canvas: comparison_plot) -> None:
        # Get analytical functions.
        (f, fstr), (g, gstr) = get_analy_func_pair(testcase)

        # Set up samples.
        r = get_sampts(pars.lgcentre, pars.shorten)
        fr = f(r, ell=pars.degree)

        # Compute the Hankel transforms.
        t_start = time()
        k_fftlog, gk_fftlog = SphericalBesselTransform(
            pars.degree, BIAS, r, PIVOT, LOWRING,
            extrap=pars.extrap, extrap_exp=EXPAND,
            threaded=pars.threaded
        ).transform(fr + 0.j)
        t_end = time()
        if timed:
            print(f"FFTLog transform time: {t_end - t_start:.3f} s")

        k_fftlog_analy, gk_fftlog_analy = \
            k_fftlog, g(k_fftlog, ell=pars.degree)

        if (_mcfit := mcfit) is not None:
            try:
                extrap_mcfit = translate_extrap_mcfit(pars.extrap)
            except ValueError as err:
                _mcfit = None
                warnings.warn(str(err))
            else:
                k_mcfit, gk_mcfit = mcfit.SphericalBessel(
                    r, lowring=LOWRING,
                    N=EXPAND*1j if extrap_mcfit else NSAMP
                )(fr * (2*np.pi)**(3./2), extrap=extrap_mcfit)
                gk_mcfit_analy = g(k_mcfit, ell=pars.degree)

        # Calculate the differences.
        if pars.show_diff:
            with np.errstate(divide='ignore', invalid='ignore'):
                dgk_fftlog = gk_fftlog.real / gk_fftlog_analy - 1.
                if _mcfit is not None:
                    dgk_mcfit = gk_mcfit / gk_mcfit_analy - 1.

        # Plot the results.
        plot_fftlog = canvas._ax_comp.plot(
            k_fftlog, k_fftlog**2 * gk_fftlog.real, ls='-', label='FFTLog'
        )
        if _mcfit is not None:
            plot_mcfit = canvas._ax_comp.plot(
                k_mcfit, k_mcfit**2 * gk_mcfit, ls='--', label='mcfit'
            )
        _ = canvas._ax_comp.plot(
            k_fftlog_analy, k_fftlog_analy**2 * gk_fftlog_analy,
            ls=':', label='analytical'
        )

        if pars.show_diff:
            canvas._ax_diff.plot(
                k_fftlog, 100 * dgk_fftlog,
                c=plot_fftlog[0].get_color(), ls='-', label='FFTLog'
            )
            if _mcfit is not None:
                canvas._ax_diff.plot(
                    k_mcfit, 100 * dgk_mcfit,
                    c=plot_mcfit[0].get_color(), ls='--', label='mcfit'
                )

        canvas._title = (
            fstr + r"$\,$, " + gstr + fr" ($\ell = {pars.degree}$)"
        )

    @comparison_plot(
        comp_type='samples',
        plot_diff=showdiff,
        xlim=(2.5e0, 2.5e2),
        ylim=(-.25e-0, 2.25e0),
        ylabel=r"$r \xi(r)$",
        timed=timed,
    )
    def run_sj_cosmo(canvas: comparison_plot) -> None:
        # Set up samples.
        (k, pk), (r, xi) = get_external_samples(
            presamp_file=PK_PRE_SAMP_FILE,
            postsamp_file=XIR_POST_SAMP_FILE.format(
                '_extrap' if pars.extrap else ''
            )
        )

        # Compute the FFTLog transform.
        t_start = time()
        r_fftlog, xi_fftlog = SphericalBesselTransform(
            pars.degree, BIAS, k, lowring=LOWRING, extrap=pars.extrap,
            threaded=pars.threaded
        ).transform_cosmo_multipoles(-1, pk)
        t_end = time()
        if timed:
            print(f"FFTLog transform time: {t_end - t_start:.3f} s")

        if hankl is not None:
            r_hankl, xi_hankl = hankl.P2xi(
                k, pk, pars.degree, n=BIAS,
                lowring=LOWRING, ext=pars.extrap
            )

        xi_ext_fftlog = InterpolatedUnivariateSpline(r, xi)(r_fftlog)
        if hankl is not None:
            xi_ext_hankl = InterpolatedUnivariateSpline(r, xi)(r_hankl)

        # Calculate the differences.
        if pars.show_diff:
            with np.errstate(divide='ignore'):
                dxi_fftlog = xi_fftlog.real / xi_ext_fftlog - 1.
                if hankl is not None:
                    dxi_hankl = xi_hankl.real / xi_ext_hankl - 1.

        # Plot the results.
        plot_fftlog = canvas._ax_comp.plot(
            r_fftlog, r_fftlog * xi_fftlog.real, ls='-', label='FFTLog'
        )
        if hankl is not None:
            plot_hankl = canvas._ax_comp.plot(
                r_hankl, r_hankl * xi_hankl.real, ls='--', label='hankl'
            )
        _ = canvas._ax_comp.plot(
            r, r * xi, ls=':', label='extern (mcfit)'
        )

        if pars.show_diff:
            canvas._ax_diff.plot(
                r_fftlog, 100 * dxi_fftlog,
                c=plot_fftlog[0].get_color(), ls='-', label='FFTLog'
            )
            if hankl is not None:
                canvas._ax_diff.plot(
                    r_hankl, 100 * dxi_hankl,
                    c=plot_hankl[0].get_color(), ls='--', label='hankl'
                )

    match testcase:
        case 'hankel-sym' | 'hankel-asym':
            return run_hankel
        case 'sj-sym' | 'sj-asym':
            return run_sj
        case 'sj-cosmo':
            return run_sj_cosmo
        case _:
            raise ValueError(f"Undefined test case: {testcase}.")


def main() -> Literal[0]:
    """Main function.

    Returns
    -------
    0
        Exit code (0 on success).

    """  # numpydoc ignore=RT01
    pars = parse_parameters()
    runner = get_testcase(pars)
    runner()

    return 0


# Global parameters
BIAS: float = 0.
PIVOT: float = 1.
LOWRING: bool = True
NSAMP: int = 768
LGRANGE: tuple[float, float] = (-5., 5.)
EXPAND: float = 1.25
PK_PRE_SAMP_FILE: str = "pk_lg_presamps.dat"
XIR_POST_SAMP_FILE: str = "xir_lg_postsamps{}.dat"

if __name__ == '__main__':
    sys.exit(main())
