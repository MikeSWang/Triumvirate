"""Benchmark performance of `triumvirate`.

"""
import inspect
import os
import warnings
from argparse import ArgumentParser
from functools import partial
from pathlib import Path
from timeit import Timer

import numpy as np


root_dir = Path(os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename
))).parent.parent


def configure():
    """Configure benchmarking."""

    config = ArgumentParser(
        description="Benchmark performance of `Triumvirate`."
    )

    config.add_argument(
        '-bl', dest='run_bispec_lpp', action='store_true',
        help="Benchmark local plane-parallel bispectrum performance."
    )
    config.add_argument(
        '-bg', dest='run_bispec_gpp', action='store_true',
        help="Benchmark global plane-parallel bispectrum performance."
    )
    config.add_argument(
        '-pl', dest='run_powspec_lpp', action='store_true',
        help="Benchmark local plane-parallel power spectrum performance."
    )
    config.add_argument(
        '-pg', dest='run_powspec_gpp', action='store_true',
        help="Benchmark global plane-parallel power spectrum performance."
    )
    config.add_argument(
        '--aniso', action='store_true',
        help="Benchmark parameter: if used, compute the anisotropic quadrupole."
    )
    config.add_argument(
        '-g', '--ngrid', dest='ngrid', type=int, nargs='+',
        help="Benchmark parameter: grid number per dimension."
    )
    config.add_argument(
        '-n', '--niter', dest='niter', type=int,
        help="Benchmark parameter: number of iterations."
    )

    return config.parse_args()


def setup_benchmarker(algo, case, multipole):
    """Set up benchmarker.

    Parameters
    ----------
    algo : str, {'bispec', 'powspec'}
        Benchmarked algorithm, either bispectrum ('bispec') or power
        spectrum ('powspec').
    case : str, {'lpp', 'gpp'}
        Benchmarked case, either local plane-parallel ('lpp') or
        global plane-parallel ('gpp').
    multipole : (sequence of) int
        Multipole degree(s).

    Returns
    -------
    callable
        The algorithmic function to be benchmarked.

    """
    if algo not in {'bispec', 'powspec'}:
        raise ValueError(f"Unrecognised algorithm to be benchmarked: '{algo}'.")
    if case not in {'lpp', 'gpp'}:
        raise ValueError(f"Unrecognised case to be benchmarked: '{case}'.")

    # Set catalogue paths.
    test_catalogue_dir = root_dir/"triumvirate"/"tests/test_input"/"catalogues"

    test_catalogue_paths = []
    test_catalogue_fields = None
    if case == 'lpp':
        test_catalogue_paths.append(test_catalogue_dir/"test_catalogue_data.dat")
        test_catalogue_paths.append(test_catalogue_dir/"test_catalogue_rand.dat")
        test_catalogue_fields = ['x', 'y', 'z', 'nz', 'ws']
    if case == 'gpp':
        test_catalogue_paths.append(test_catalogue_dir/"test_catalogue_sim.dat")
        test_catalogue_fields = ['x', 'y', 'z', 'ws']

    # Set parameter defaults.
    BOXSIZE = 1000.
    ASSIGNMENT = 'tsc'
    KMIN = 0.005
    KMAX = 0.205
    NBIN = 20

    # Create parameter set.
    from triumvirate.parameters import (
        fetch_paramset_template,
        _modify_sampling_parameters,
    )

    paramset = _modify_sampling_parameters(
        fetch_paramset_template('dict'),
        params_sampling={
            'boxsize': [BOXSIZE,] * 3,
            'assignment': ASSIGNMENT,
        }
    )

    paramset['verbose'] = 50

    # Create binning scheme.
    from triumvirate.dataobjs import Binning

    binning = Binning(
        'fourier', 'lin', bin_min=KMIN, bin_max=KMAX, num_bins=NBIN
    )

    # Load test catalogue(s).
    from triumvirate.catalogue import ParticleCatalogue

    test_catalogues = []
    with warnings.catch_warnings():
        warnings.filterwarnings(
            action='ignore', message=".*field is not provided.*"
        )
        for tcltg_path in test_catalogue_paths:
            test_catalogues.append(
                ParticleCatalogue.read_from_file(
                    tcltg_path, names=test_catalogue_fields
                )
            )

    # Declare test function.
    if algo == 'bispec':
        from triumvirate.threept import compute_bispec, compute_bispec_in_gpp_box
        if case == 'lpp':
            def run_triumvirate(ngrid):
                compute_bispec(
                    *test_catalogues, degrees=multipole, binning=binning,
                    sampling_params={'ngrid': [ngrid,]*3},
                    paramset=paramset
                )
        if case == 'gpp':
            def run_triumvirate(ngrid):
                compute_bispec_in_gpp_box(
                    *test_catalogues, degrees=multipole, binning=binning,
                    sampling_params={'ngrid': [ngrid,]*3},
                    paramset=paramset
                )
    if algo == 'powspec':
        from triumvirate.twopt import compute_powspec, compute_powspec_in_gpp_box
        if case == 'lpp':
            def run_triumvirate(ngrid):
                compute_powspec(
                    *test_catalogues, degree=multipole, binning=binning,
                    sampling_params={'ngrid': [ngrid,]*3},
                    paramset=paramset
                )
        if case == 'gpp':
            def run_triumvirate(ngrid):
                compute_powspec_in_gpp_box(
                    *test_catalogues, degree=multipole, binning=binning,
                    sampling_params={'ngrid': [ngrid,]*3},
                    paramset=paramset
                )

    return run_triumvirate


if __name__ == '__main__':

    config = configure()

    # Determine benchmarked cases.
    benchmarkers = {}
    if config.run_bispec_lpp:
        multipole = (0, 0, 0) if not config.aniso else (2, 0, 2)
        benchmarkers['bispec_lpp'] = setup_benchmarker('bispec', 'lpp', multipole)
    if config.run_bispec_gpp:
        multipole = (0, 0, 0) if not config.aniso else (2, 0, 2)
        benchmarkers['bispec_gpp'] = setup_benchmarker('bispec', 'gpp', multipole)
    if config.run_powspec_lpp:
        multipole = 0 if not config.aniso else 2
        benchmarkers['powspec_lpp'] = setup_benchmarker('powspec', 'lpp', multipole)
    if config.run_powspec_gpp:
        multipole = 0 if not config.aniso else 2
        benchmarkers['powspec_gpp'] = setup_benchmarker('powspec', 'gpp', multipole)

    # Perform benchmarking runs.
    results = {}
    for algo, func in benchmarkers.items():
        runtimes = {}
        for ngrid in config.ngrid:
            timer = Timer(partial(func, ngrid))
            runtimes[ngrid] = timer.repeat(repeat=config.niter, number=1)
        results[algo] = runtimes

    # Export benchmark results.
    if results:
        np.save("benchmark_results.npy", results, allow_pickle=True)
    else:
        warnings.warn("No benchmarking runs.")
