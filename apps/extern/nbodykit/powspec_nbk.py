"""Measure power spectrum using `nbodykit`.

"""
import os.path as osp
from argparse import ArgumentParser
from pathlib import Path

import numpy as np
from nbodykit import setup_logging
from nbodykit.lab import CSVCatalog, FKPCatalog
from nbodykit.lab import ConvolvedFFTPower, FFTPower


def configure():
    """Parse command-line arguments to runtime parameters.

    Returns
    -------
    :class:`argparse.Namespace`
        Parameter namespace.

    """
    parser = ArgumentParser(
        description="Measure power spectrum using `nbodykit`."
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
        '--data-ctlgfile', type=str,
        help="data-source catalogue file"
    )
    parser.add_argument(
        '--rand-ctlgfile', type=str,
        help="random-source catalogue file"
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
        help="power spectrum multipole degree"
    )
    parser.add_argument(
        '--range', type=float, nargs=2,
        help="wavenumber range (in h/Mpc)"
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
    """Run `nbodykit` to measure power spectrum.

    Parameters
    ----------
    config : :class:`argparse.Namespace`
        Parameters.

    Returns
    -------
    :class:`nbodykit.algorithms.convpower.fkp.ConvolvedFFTPower`
        Power spectrum measurements.

    """
    setup_logging()

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

    ctlg_data = CSVCatalog(data_ctlgpath, config.ctlg_fields)
    if rand_ctlgpath is None:
        ctlg_rand = None
    else:
        ctlg_rand = CSVCatalog(rand_ctlgpath, config.ctlg_fields)

    # Process catalogues.
    add_position_field_to_catalogue(ctlg_data)
    if ctlg_rand:
        add_position_field_to_catalogue(ctlg_rand)

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

    if ctlg_rand:
        mesh = FKPCatalog(ctlg_data, ctlg_rand).to_mesh(
            BoxSize=config.boxsize,
            Nmesh=config.ngrid,
            resampler=config.scheme,
            compensated=True,
            interlaced=config.interlace
        )
    else:
        mesh = ctlg_data.to_mesh(
            BoxSize=config.boxsize,
            Nmesh=config.ngrid,
            resampler=config.scheme,
            compensated=True,
            interlaced=config.interlace
        )

    dk = np.diff(config.range) / config.num_bins

    if config.rand_ctlgfile:
        results = ConvolvedFFTPower(
            mesh, poles=[config.multipole,],  # noqa: E231
            kmin=config.range[0], kmax=config.range[1]+dk/100., dk=dk
        )
    else:
        results = FFTPower(
            mesh, mode='1d', poles=[config.multipole,],  # noqa: E231
            kmin=config.range[0], kmax=config.range[1]+dk/100., dk=dk
        )

    return results


def save_results(results, config):
    """Save `nbodykit` measurements to file.

    Parameters
    ----------
    results : :class:`nbodykit.algorithms.convpower.fkp.ConvolvedFFTPower`
        Power spectrum measurements.
    config : :class:`argparse.Namespace`
        Parameters.

    """
    header = ", ".join([
        "[{:d}] {}".format(colidx, colname)
        for colidx, colname in enumerate([
            "k_cen", "k_eff", "nmodes",
            "Re{{pk{:d}_raw}}".format(config.multipole),
            "Im{{pk{:d}_raw}}".format(config.multipole),
            "Re{{pk{:d}_shot}}".format(config.multipole),
            "Im{{pk{:d}_shot}}".format(config.multipole)
        ])
    ])
    fmt = '\t'.join(['%.9e'] * 2 + ['%10d'] + ['% .9e'] * 4)

    outtable = np.transpose([
        (results.poles.edges['k'][:-1] + results.poles.edges['k'][1:]) / 2.,
        results.poles['k'],
        results.poles['modes'],
        results.poles[f'power_{config.multipole}'].real,
        results.poles[f'power_{config.multipole}'].imag,
        results.poles.attrs['shotnoise'] * np.full_like(
            results.poles['k'], config.multipole == 0
        ),
        results.poles.attrs['shotnoise'] * np.zeros_like(results.poles['k'])
    ])
    outpath = osp.join(config.output_dir, config.output_file)

    Path(outpath).parent.mkdir(parents=True, exist_ok=True)
    np.savetxt(outpath, outtable, fmt=fmt, header=header)

    print(f"Saved power spectrum to: {outpath}")


if __name__ == '__main__':
    config = configure()
    results = run_nbodykit(config)
    save_results(results, config)
