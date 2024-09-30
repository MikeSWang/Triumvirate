"""Benchmark performance of ``triumvirate``.

"""
import inspect
import os
import warnings
from argparse import ArgumentParser
from functools import partial
from pathlib import Path
from timeit import Timer

import numpy as np

from triumvirate.threept import (
    compute_3pcf,
    compute_3pcf_in_gpp_box,
    compute_bispec,
    compute_bispec_in_gpp_box,
)
from triumvirate.twopt import (
    compute_corrfunc,
    compute_corrfunc_in_gpp_box,
    compute_powspec,
    compute_powspec_in_gpp_box,
)


ALGOS = {
    'bispec': {
        'lpp': {
            'fourier': compute_bispec,
            'config': compute_3pcf,
        },
        'gpp': {
            'fourier': compute_bispec_in_gpp_box,
            'config': compute_3pcf_in_gpp_box,
        }
    },
    'powspec': {
        'lpp': {
            'fourier': compute_powspec,
            'config': compute_corrfunc,
        },
        'gpp': {
            'fourier': compute_powspec_in_gpp_box,
            'config': compute_corrfunc_in_gpp_box,
        }
    },
}

root_dir = Path(os.path.dirname(os.path.abspath(
    inspect.getframeinfo(inspect.currentframe()).filename
))).parent.parent


def configure():
    """Configure benchmarking.

    Returns
    -------
    argparse.Namespace
        Parsed command-line arguments.

    """
    cfg = ArgumentParser(
        description="Benchmark performance of `Triumvirate`."
    )

    cfg.add_argument(
        '-blf', dest='run_bispec_lpp', action='store_true',
        help="benchmark local plane-parallel bispectrum performance"
    )
    cfg.add_argument(
        '-bgf', dest='run_bispec_gpp', action='store_true',
        help="benchmark global plane-parallel bispectrum performance"
    )
    cfg.add_argument(
        '-blc', dest='run_3pcf_lpp', action='store_true',
        help="benchmark local plane-parallel 3PCF performance"
    )
    cfg.add_argument(
        '-bgc', dest='run_3pcf_gpp', action='store_true',
        help="benchmark global plane-parallel 3PCF performance"
    )

    cfg.add_argument(
        '-plf', dest='run_powspec_lpp', action='store_true',
        help="benchmark local plane-parallel power spectrum performance"
    )
    cfg.add_argument(
        '-pgf', dest='run_powspec_gpp', action='store_true',
        help="benchmark global plane-parallel power spectrum performance"
    )
    cfg.add_argument(
        '-plc', dest='run_2pcf_lpp', action='store_true',
        help="benchmark local plane-parallel 2PCF performance"
    )
    cfg.add_argument(
        '-pgc', dest='run_2pcf_gpp', action='store_true',
        help="benchmark global plane-parallel 2PCF performance"
    )

    cfg.add_argument(
        '-a', '--aniso', action='store_true',
        help="benchmark parameter: "
             "if used, compute the anisotropic quadrupole"
    )
    cfg.add_argument(
        '-g', '--ngrid', dest='ngrid', type=int, nargs='+',
        help="benchmark parameter: mesh grid number per dimension"
    )
    cfg.add_argument(
        '-n', '--niter', dest='niter', type=int,
        help="benchmark parameter: number of iterations"
    )
    cfg.add_argument(
        '--output-tag', default='',
        help="output file tag for saved benchmarking results"
    )
    cfg.add_argument(
        '--output-dir', default="",
        help="output directory to save benchmarking results to"
    )

    return cfg.parse_args()


def setup_benchmarker(algo, case, space, multipole):
    """Set up benchmarker.

    Parameters
    ----------
    algo : str, {'bispec', 'powspec'}
        Benchmarked algorithm, either bispectrum ('bispec') or power
        spectrum ('powspec').
    case : str, {'lpp', 'gpp'}
        Benchmarked case, either local plane-parallel ('lpp') or
        global plane-parallel ('gpp').
    space : str, {'fourier', 'config'}
        Coordinate space, either Fourier ('fourier') or
        configuration ('config').
    multipole : (sequence of) int
        Multipole degree(s).

    Returns
    -------
    callable
        The algorithmic function to be benchmarked.

    """
    if algo not in {'bispec', 'powspec'}:
        raise ValueError(
            f"Unrecognised algorithm to be benchmarked: '{algo}'."
        )
    if case not in {'lpp', 'gpp'}:
        raise ValueError(f"Unrecognised case to be benchmarked: '{case}'.")
    if space not in {'fourier', 'config'}:
        raise ValueError(f"Unrecognised space to be benchmarked: '{space}'.")

    # Set catalogue paths.
    test_catalogue_dir = root_dir/"triumvirate"/"tests/test_input"/"catalogues"

    test_catalogue_paths = []
    test_catalogue_fields = None
    if case == 'lpp':
        test_catalogue_paths.append(
            test_catalogue_dir/"test_catalogue_data.dat"
        )
        test_catalogue_paths.append(
            test_catalogue_dir/"test_catalogue_rand.dat"
        )
        test_catalogue_fields = ['x', 'y', 'z', 'nz', 'ws']
    if case == 'gpp':
        test_catalogue_paths.append(
            test_catalogue_dir/"test_catalogue_sim.dat"
        )
        test_catalogue_fields = ['x', 'y', 'z', 'ws']

    # Set parameter defaults.
    BOXSIZE = 1000.
    ASSIGNMENT = 'tsc'

    NBIN = 20
    if space == 'fourier':
        BIN_MIN, BIN_MAX = 0.005, 0.205
    if space == 'config':
        BIN_MIN, BIN_MAX = 5., 205.

    # Create parameter set.
    from triumvirate.parameters import (
        fetch_paramset_template,
        _modify_sampling_parameters,
    )

    paramset = _modify_sampling_parameters(
        fetch_paramset_template('dict'),
        params_sampling={
            'boxsize': [BOXSIZE,] * 3,  # noqa: E231
            'assignment': ASSIGNMENT,
        }
    )

    paramset['verbose'] = 50

    # Create binning scheme.
    from triumvirate.dataobjs import Binning

    binning = Binning(
        space, 'lin', bin_min=BIN_MIN, bin_max=BIN_MAX, num_bins=NBIN
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
        def run_triumvirate(ngrid):
            ALGOS[algo][case][space](
                *test_catalogues, degrees=multipole, binning=binning,
                sampling_params={'ngrid': [ngrid,]*3},  # noqa: E231
                paramset=paramset
            )
    elif algo == 'powspec':
        def run_triumvirate(ngrid):
            ALGOS[algo][case][space](
                *test_catalogues, degree=multipole, binning=binning,
                sampling_params={'ngrid': [ngrid,]*3},  # noqa: E231
                paramset=paramset
            )

    return run_triumvirate


if __name__ == '__main__':

    cfg = configure()

    # Determine benchmarked cases.
    benchmarkers = {}
    if cfg.run_bispec_lpp:
        multipole = (0, 0, 0) if not cfg.aniso else (2, 0, 2)
        benchmarkers['bispec_lpp'] = setup_benchmarker(
            'bispec', 'lpp', 'fourier', multipole
        )
    if cfg.run_bispec_gpp:
        multipole = (0, 0, 0) if not cfg.aniso else (2, 0, 2)
        benchmarkers['bispec_gpp'] = setup_benchmarker(
            'bispec', 'gpp', 'fourier', multipole
        )
    if cfg.run_3pcf_lpp:
        multipole = (0, 0, 0) if not cfg.aniso else (2, 0, 2)
        benchmarkers['3pcf_lpp'] = setup_benchmarker(
            'bispec', 'lpp', 'config', multipole
        )
    if cfg.run_3pcf_gpp:
        multipole = (0, 0, 0) if not cfg.aniso else (2, 0, 2)
        benchmarkers['3pcf_gpp'] = setup_benchmarker(
            'bispec', 'gpp', 'config', multipole
        )
    if cfg.run_powspec_lpp:
        multipole = 0 if not cfg.aniso else 2
        benchmarkers['powspec_lpp'] = setup_benchmarker(
            'powspec', 'lpp', 'fourier', multipole
        )
    if cfg.run_powspec_gpp:
        multipole = 0 if not cfg.aniso else 2
        benchmarkers['powspec_gpp'] = setup_benchmarker(
            'powspec', 'gpp', 'fourier', multipole
        )
    if cfg.run_2pcf_lpp:
        multipole = 0 if not cfg.aniso else 2
        benchmarkers['2pcf_lpp'] = setup_benchmarker(
            'powspec', 'lpp', 'config', multipole
        )
    if cfg.run_2pcf_gpp:
        multipole = 0 if not cfg.aniso else 2
        benchmarkers['2pcf_gpp'] = setup_benchmarker(
            'powspec', 'gpp', 'config', multipole
        )

    # Set outputs.
    output_path = Path(
        cfg.output_dir, f"benchmark_results{cfg.output_tag}.npy"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # Perform benchmarking runs.
    results = {}
    for algo, func in benchmarkers.items():
        runtimes = {}
        for ngrid in cfg.ngrid:
            timer = Timer(partial(func, ngrid))
            print(
                f"Benchmarking: {algo=}, aniso={cfg.aniso}, "
                f"{ngrid=}, niter={cfg.niter}"
            )
            runtimes[ngrid] = timer.repeat(repeat=cfg.niter, number=1)

            # Periodically export the results.
            results[algo] = runtimes
            if results:
                np.save(output_path, results, allow_pickle=True)

    if not results:
        warnings.warn("No benchmarking runs.")
