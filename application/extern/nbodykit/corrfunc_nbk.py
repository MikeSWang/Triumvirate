"""Measure correlation function using `nbodykit`.

"""
from argparse import ArgumentParser
import os.path as osp
from pathlib import Path

import numpy as np
from nbodykit import setup_logging
from nbodykit.lab import CSVCatalog
from nbodykit.lab import FFTCorr


def configure():
    """Parse command-line arguments to runtime parameters.

    Returns
    -------
    :class:`argparse.Namespace`
        Parameter namespace.

    """
    parser = ArgumentParser(
        description="Measure flat-sky correlation function using `nbodykit`."
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
        '--ctlg-file', type=str, required=True,
        help="catalogue file"
    )
    parser.add_argument(
        '--ctlg-fields', type=str, nargs='+',
        help="catalogue fields (must contain sample-weight column 'Weight')"
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
        '--scheme', type=str,
        help="mesh assignment scheme"
    )
    parser.add_argument(
        '--interlace', action='store_true',
        help="enable interlacing"
    )

    parser.add_argument(
        '--multipole', type=int,
        help="correlation function multipole degree"
    )
    parser.add_argument(
        '--range', type=float, nargs=2,
        help="separation range (in Mpc/h)"
    )
    parser.add_argument(
        '--num-bins', type=int,
        help="number of bins"
    )

    config = parser.parse_args()

    if len(config.boxsize) == 1:
        config.boxsize = [config.boxsize] * 3
    if len(config.ngrid) == 1:
        config.ngrid = [config.ngrid] * 3

    return config


def add_position_field_to_catalogue(catalogue):
    """Add the 'Position' column from 'x', 'y', 'z' columns to a catalogue.

    Parameters
    ----------
    catalogue : :class:`nbodykit.base.CatalogueSource`
        Catalogue.  Must contain Cartesian coordinate columns 'x', 'y'
        and 'z'.

    """
    catalogue['Position'] = \
        catalogue['x'][:, None] * [1., 0., 0.] \
        + catalogue['y'][:, None] * [0., 1., 0.] \
        + catalogue['z'][:, None] * [0., 0., 1.]


def run_nbodykit(config):
    """Run `nbodykit` to measure correlation function.

    Parameters
    ----------
    config : :class:`argparse.Namespace`
        Parameters.

    Returns
    -------
    :class:`nbodykit.algorithms.fftcorr.FFTCorr`
        Correlation function measurements.

    """
    setup_logging()

    # Load catalogue.
    data_ctlgpath = osp.join(config.input_dir, config.data_ctlgfile)

    print("Loading catalogue: {}.".format(data_ctlgpath))

    ctlg_data = CSVCatalog(data_ctlgpath, config.ctlg_fields)

    # Process catalogue.
    add_position_field_to_catalogue(ctlg_data)

    # Make measurements.
    print(
        "Measuring correlation function: "
        "multipole={}, range={}, num_bins={}, "
        "ngrid={}, boxsize={}, scheme={}, interlace={}"
        .format(
            config.multipole, config.range, config.num_bins,
            config.ngrid, config.boxsize, config.scheme, config.interlace
        )
    )

    mesh = ctlg_data.to_mesh(
        BoxSize=config.boxsize,
        Nmesh=config.ngrid,
        resampler=config.scheme,
        compensated=True,
        interlaced=config.interlace
    )

    dr = np.diff(config.range) / config.num_bins

    results = FFTCorr(
        mesh, mode='1d', poles=[config.multipole,],  # noqa: E231
        rmin=config.range[0], rmax=config.range[1]+dr/100., dr=dr
    )

    return results


def save_results(results, config):
    """Save `nbodykit` measurements to file.

    Parameters
    ----------
    results : :class:`nbodykit.algorithms.fftcorr.FFTCorr`
        Correlation function measurements.
    config : :class:`argparse.Namespace`
        Parameters.

    """
    header = ", ".join([
        "[{:d}] {}".format(colidx, colname)
        for colidx, colname in enumerate([
            "r_cen", "r_eff", "npairs",
            "Re{{xi{:d}_raw}}".format(config.multipole),
            "Im{{xi{:d}_raw}}".format(config.multipole),
        ])
    ])
    fmt = '\t'.join(['%.9e'] * 2 + ['%10d'] + ['% .9e'] * 2)

    outtable = np.transpose([
        (results.poles.edges['r'][:-1] + results.poles.edges['r'][1:]) / 2.,
        results.poles['r'],
        results.poles['modes'],
        results.poles[f'corr_{config.multipole}'].real,
        results.poles[f'corr_{config.multipole}'].imag,
    ])
    outpath = osp.join(config.output_dir, config.output_file)

    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(outpath, outtable, fmt=fmt, header=header)

    print("Saved correlation function to {}.".format(outpath))


if __name__ == '__main__':
    config = configure()
    results = run_nbodykit(config)
    save_results(results, config)
