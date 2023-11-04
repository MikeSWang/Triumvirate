"""Measure power spectrum using `pypower`.

"""
import os.path as osp
from argparse import ArgumentParser
from itertools import repeat
from pathlib import Path
from pprint import pformat
from tempfile import TemporaryDirectory

import numpy as np
from pypower import CatalogFFTPower, normalization_from_nbar


def configure():
    """Parse command-line arguments to runtime parameters.

    Returns
    -------
    :class:`argparse.Namespace`
        Parameter namespace.

    """
    parser = ArgumentParser(
        description="Measure power spectrum using `pypower`."
    )

    parser.add_argument(
        '--input-dir', type=str, default="./",
        help="input (catalogue) directory"
    )
    parser.add_argument(
        '--output-dir', type=str, default="./",
        help="output (measurement) directory"
    )
    parser.add_argument(
        '--data-ctlgfile', type=str, required=True,
        help="data-source catalogue file"
    )
    parser.add_argument(
        '--rand-ctlgfile', type=str, default=None,
        help="random-source catalogue file"
    )
    parser.add_argument(
        '--ctlg-fields', type=str, nargs='+', required=True,
        help="catalogue fields"
    )
    parser.add_argument(
        '--output-file', type=str, required=True,
        help="output measurement file name"
    )

    parser.add_argument(
        '--boxsize', type=float, nargs='+', required=True,
        help="box size for all dimensions"
    )
    parser.add_argument(
        '--ngrid', type=int, nargs='+', required=True,
        help="mesh grid number"
    )
    parser.add_argument(
        '--scheme', type=str, default='pcs',
        help="mesh assignment scheme"
    )
    parser.add_argument(
        '--interlace', action='store_true',
        help="interlacing switch"
    )

    parser.add_argument(
        '--multipole', type=int, required=True,
        help="power spectrum multipole degree"
    )
    parser.add_argument(
        '--range', type=float, nargs=2, required=True,
        help="wavenumber range (in h/Mpc)"
    )
    parser.add_argument(
        '--num-bins', type=int, required=True,
        help="number of bins"
    )

    parser.add_argument(
        '--norm-only', action='store_true',
        help="normalisation-only switch"
    )
    parser.add_argument(
        '--no-norm', action='store_true',
        help="no-normalisation switch"
    )

    config = parser.parse_args()

    if len(config.boxsize) == 1:
        config.boxsize = [config.boxsize] * 3
    if len(config.ngrid) == 1:
        config.ngrid = [config.ngrid] * 3

    return config


def run_pypower(config):
    """Run `pypower` to measure the power spectrum.

    Parameters
    ----------
    config : :class:`argparse.Namespace`
        Parameters.

    Returns
    -------
    :class:`pypower.CatalogFFTPower`
        Power spectrum measurements.

    """
    def get_dtype():
        return list(zip(config.ctlg_fields, repeat(float)))

    # Load catalogues.
    data_ctlgpath = osp.join(config.input_dir, config.data_ctlgfile)
    if config.rand_ctlgfile is None:
        rand_ctlgpath = None
    else:
        rand_ctlgpath = osp.join(config.input_dir, config.rand_ctlgfile)

    print(
        "Loading data and random catalogues: {}; {}.".format(
            data_ctlgpath, rand_ctlgpath
        )
    )

    gcat = np.loadtxt(data_ctlgpath, dtype=get_dtype())
    if rand_ctlgpath is None:
        rcat = None
    else:
        rcat = np.loadtxt(rand_ctlgpath, dtype=get_dtype())

    # Process catalogues.
    data_positions1 = np.column_stack((gcat['x'], gcat['y'], gcat['z']))
    if rcat is None:
        randoms_positions1 = None
    else:
        randoms_positions1 = np.column_stack((rcat['x'], rcat['y'], rcat['z']))

    if rcat is None:
        data_weights1 = gcat['ws']
        randoms_weights1 = None
    else:
        data_weights1 = gcat['ws'] * gcat['wc']
        randoms_weights1 = rcat['ws'] * rcat['wc']

    # Calculate normalisation.
    if config.norm_only:
        data_norm = normalization_from_nbar(gcat['nz'], weights=data_weights1)
        rand_norm = normalization_from_nbar(
            rcat['nz'], weights=randoms_weights1, data_weights=data_weights1
        )
        print(
            "Normalisation from particles: {} (data), {} (rand)"
            .format(1./data_norm, 1./rand_norm)
        )
        exit(0)

    # Make measurements.
    print(
        "Measuring power spectrum: "
        "multipole={}, range={}, num_bins={}, "
        "ngrid={}, boxsize={}, scheme={}, interlace={}"
        .format(
            config.multipole, config.range, config.num_bins,
            config.ngrid, config.boxsize, config.scheme, config.interlace
        )
    )

    bin_edges = np.linspace(*config.range, num=config.num_bins+1)

    results = CatalogFFTPower(
        data_positions1,
        randoms_positions1=randoms_positions1,
        data_weights1=data_weights1,
        randoms_weights1=randoms_weights1,
        position_type='pos', weight_type='auto', los='endpoint',
        ells=[config.multipole], edges=bin_edges,
        nmesh=config.ngrid, boxsize=config.boxsize,
        resampler=config.scheme, interlacing=config.interlace
    )

    print("Measured power spectrum.")

    return results


def save_results(results, config):
    """Save `pypower` measurements to file.

    Parameters
    ----------
    results : :class:`pypower.CatalogFFTPower`
        `pypower` measurements.
    config : :class:`argparse.Namespace`
        Parameters.

    """
    # Get data table.
    norm = 1. / results.poles.wnorm if not config.no_norm else 1.

    k_cen = 1./2 * np.add(
        results.poles.edges[0][:-1], results.poles.edges[0][1:]
    )
    k_eff = results.poles.modes[0]
    nmodes = results.poles.nmodes
    pk0 = np.multiply(norm, results.poles.power_nonorm[0])
    pk0_shot = np.full_like(
        pk0, np.multiply(norm, results.poles.shotnoise_nonorm)
    )

    # Get other attributes as header.
    with TemporaryDirectory() as tmpdir:
        tmpfile = osp.join(tmpdir, config.output_file + '.pypower.npy')
        results.save(tmpfile)
        results = np.load(tmpfile, allow_pickle=True).item()

    # norm = 1. / results['poles']['wnorm']

    # k_cen = 1./2 * np.add(
    #     results['poles']['edges'][0][:-1], results['poles']['edges'][0][1:]
    # )
    # k_eff = results['poles']['modes'][0]
    # nmodes = results['poles']['nmodes']
    # pk0 = np.multiply(norm, results['poles']['power_nonorm'][0])
    # pk0_shot = np.full_like(
    #     pk0, np.multiply(norm, results['poles']['shotnoise_nonorm'])
    # )

    header = '\n'.join([
        pformat(results['attrs']),
        "'powerzero_nonorm': {}".format(results['poles']['power_zero_nonorm']),
        "'shotnoise_nonorm': {}".format(results['poles']['shotnoise_nonorm']),
        "'wnorm': {}".format(results['poles']['wnorm']),
        ", ".join([
            "[{:d}] {}".format(colidx, colname)
            for colidx, colname in enumerate([
                "k_cen", "k_eff", "nmodes",
                "Re{{pk{:d}_raw}}".format(config.multipole),
                "Im{{pk{:d}_raw}}".format(config.multipole),
                "Re{{pk{:d}_shot}}".format(config.multipole),
                "Im{{pk{:d}_shot}}".format(config.multipole)
            ])
        ])
    ])
    fmt = '\t'.join(['%.9e'] * 2 + ['%10d'] + ['% .9e'] * 4)

    # Save to file.
    outtable = np.column_stack([
        k_cen, k_eff, nmodes, pk0.real, pk0.imag, pk0_shot.real, pk0_shot.imag
    ])
    outpath = osp.join(config.output_dir, config.output_file)

    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(outpath, outtable, header=header, fmt=fmt)

    print("Saved power spectrum to {}.".format(outpath))


if __name__ == '__main__':
    config = configure()
    results = run_pypower(config)
    save_results(results, config)
