"""Demonstrate the use of the FFTLog algorithm for double spherical
Bessel transforms.

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
from scipy.interpolate import RectBivariateSpline
from scipy.special import gamma

from triumvirate.transforms import DoubleSphericalBesselTransform


def parse_parameters() -> argparse.Namespace:
    """Initialise the command-line argument parser.

    Returns
    -------
    pars
        Parsed command-line parameters.

    """
    parser = argparse.ArgumentParser(
        description="Double spherical Bessel transform demo."
    )

    parser.add_argument(
        '-t', "--test-case", type=str,
        choices=['sym', 'asym', 'cosmo'],
        required=True,
        help="test case"
    )
    parser.add_argument(
        '--degrees', type=int, nargs=2, default=(0, 0),
        help="degrees of the double spherical Bessel transform "
             "(default is %(default)s)"
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
    if pars.test_case == 'cosmo' and pars.degrees != (0, 0):
        pars.degrees = (0, 0)
        warnings.warn(
            "Double spherical Bessel transform degree overridden to (0, 0) "
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
    testcase: Literal['sym', 'asym']
) -> tuple[tuple[Callable, str], tuple[Callable, str]]:
    """Get pre- and post-transform analytical functions for the given
    test case.

    Parameters
    ----------
    testcase
        Test cases:

        - 'sym' for symmetric-form functions for the double spherical
          Bessel transform;
        - 'asym' for asymmetric-form functions for the double spherical
          Bessel transform.

    Returns
    -------
    Pre- and post-transform functions with their TeX strings.

    """  # numpydoc ignore=RT03
    match testcase:
        case 'sym':
            return ((
                lambda r, ell=0: r**ell * np.exp(-r*r/2.),
                r"$f(r) = r^\ell \mathrm{e}^{-r^2 / 2}$"
            ), (
                lambda k, ell=0: (2*np.pi)**(3./2) * k**ell * np.exp(-k*k/2.),
                r"$g(k) = (2\pi)^{3/2} k^\ell \mathrm{e}^{-k^2 / 2}$"
            ))
        case 'asym':
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
                self._ax_diff.axhline(-1.e-3, c='k', alpha=.5, ls='--')
                self._ax_diff.axhline(1.e-3, c='k', alpha=.5, ls='--')
                self._ax_diff.set_ylim(-1.75e-3, 1.75e-3)

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
        xlim={
            'sym': (1.e-5, 1.e1),
            'asym': (1.e-1, 1.e3),
        }.get(testcase),
        ylim={
            'sym': (1.e-18, 1.e3),
            'asym': (1.e-0, 1.e9),
        }.get(testcase),
        xlabel=r"$k_1 = k_2$",
        ylabel=r"$k_1^2 k_2^2 g(k_1) g(k_2)$",
        timed=timed,
    )
    def run_analy(canvas: comparison_plot) -> None:
        # Get analytical functions.
        (f, fstr), (g, gstr) = get_analy_func_pair(testcase)

        # Set up samples.
        r = get_sampts(pars.lgcentre, pars.shorten)
        fr_1, fr_2 = (f(r, ell=ell) for ell in pars.degrees)

        # Compute the Hankel transforms.
        frr = fr_1[:, None] * fr_2[None, :]

        t_start = time()
        kk_fftlog, gkk_fftlog = DoubleSphericalBesselTransform(
            pars.degrees, BIASES, r, PIVOT, LOWRING,
            extrap=pars.extrap, extrap_exp=EXPAND,
            threaded=pars.threaded
        ).transform(frr + 0.j)
        t_end = time()
        if timed:
            print(f"FFTLog transform time: {t_end - t_start:.3f} s")

        gkk_analy = np.multiply(*[
            g(ki_fftlog, ell=ell)
            for (ki_fftlog, ell) in zip(kk_fftlog, pars.degrees)
        ])

        k_diag = kk_fftlog[0][np.equal(*kk_fftlog)]

        # Calculate the differences.
        if pars.show_diff:
            with np.errstate(divide='ignore', invalid='ignore'):
                dgkk_fftlog = gkk_fftlog.real / gkk_analy - 1.

        # Plot the results.
        plot_fftlog = canvas._ax_comp.plot(
            k_diag, k_diag**4 * np.diag(gkk_fftlog).real,
            ls='-', label='FFTLog'
        )
        _ = canvas._ax_comp.plot(
            k_diag, k_diag**4 * np.diag(gkk_analy),
            ls=':', label='analytical'
        )

        if pars.show_diff:
            canvas._ax_diff.plot(
                k_diag, 100 * np.diag(dgkk_fftlog),
                c=plot_fftlog[0].get_color(), ls='-', label='FFTLog'
            )

        canvas._title = (
            fstr + r"$\,$, " + gstr + fr" ($\ell = {pars.degrees}$)"
        )

    @comparison_plot(
        comp_type='samples',
        plot_diff=showdiff,
        xlim=(2.5e0, 2.5e2),
        ylim=(-.25e-0, 2.25e0),
        xlabel=r"$r_1 = r_2$",
        ylabel=r"$[r_1 r_2 \xi(r_1) \xi(r_2)]^{1/2}$",
        timed=timed,
    )
    def run_samples(canvas: comparison_plot) -> None:
        # Set up samples.
        (k, pk), (r, xi) = get_external_samples(
            presamp_file=PK_PRE_SAMP_FILE,
            postsamp_file=XIR_POST_SAMP_FILE.format(
                '_extrap' if pars.extrap else ''
            )
        )

        # Compute the FFTLog transform.
        bkk = pk[:, None] * pk[None, :]

        t_start = time()
        rr_fftlog, xirr_fftlog = DoubleSphericalBesselTransform(
            pars.degrees, BIASES, k, lowring=LOWRING, extrap=pars.extrap,
            threaded=pars.threaded
        ).transform_cosmo_multipoles(-1, bkk)
        t_end = time()
        if timed:
            print(f"FFTLog transform time: {t_end - t_start:.3f} s")

        xirr = xi[:, None] * xi[None, :]
        xirr_ext_fftlog = RectBivariateSpline(r, r, xirr)(
            *rr_fftlog, grid=False
        )

        r_diag = rr_fftlog[0][np.equal(*rr_fftlog)]
        # r_row = rr_fftlog[0][:, 0]

        # Calculate the differences.
        if pars.show_diff:
            with np.errstate(divide='ignore'):
                dxirr_fftlog = xirr_fftlog.real / xirr_ext_fftlog - 1.

        # Plot the results.
        plot_fftlog = canvas._ax_comp.plot(
            r_diag, r_diag * np.sqrt(np.diag(xirr_ext_fftlog).real),
            ls='-', label='FFTLog'
        )
        _ = canvas._ax_comp.plot(
            r, r * np.sqrt(np.diag(xirr)), ls=':', label='extern (mcfit)'
        )

        if pars.show_diff:
            canvas._ax_diff.plot(
                r_diag, 100 * np.diag(dxirr_fftlog),
                c=plot_fftlog[0].get_color(), ls='-', label='FFTLog'
            )

    match testcase:
        case 'sym' | 'asym':
            return run_analy
        case 'cosmo':
            return run_samples
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
BIASES: tuple[int, int] = (0, 0)
PIVOT: float = 1.
LOWRING: bool = True
NSAMP: int = 768
LGRANGE: tuple[float, float] = (-5., 5.)
EXPAND: float = 1.25
PK_PRE_SAMP_FILE: str = "pk_lg_presamps.dat"
XIR_POST_SAMP_FILE: str = "xir_lg_postsamps{}.dat"

if __name__ == '__main__':
    sys.exit(main())
