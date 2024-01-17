"""Generate integration test dataset.

"""
import os.path
from argparse import ArgumentParser

from triumvirate.catalogue import ParticleCatalogue
from triumvirate.logger import setup_logger
from triumvirate.parameters import ParameterSet, fetch_paramset_template
from triumvirate.threept import (
    compute_3pcf,
    compute_3pcf_in_gpp_box,
    compute_3pcf_window,
    compute_bispec,
    compute_bispec_in_gpp_box,
)
from triumvirate.twopt import (
    compute_corrfunc,
    compute_corrfunc_in_gpp_box,
    compute_corrfunc_window,
    compute_powspec,
    compute_powspec_in_gpp_box,
)


def setup_config():
    """Set up configuration for integration testing.

    Returns
    -------
    :class:`argparse.Namespace`
        Integration testing configuration.

    """
    parser = ArgumentParser(
        description="Configure integration test dataset generation."
    )

    parser.add_argument(
        '--input-dir', type=str, default="tests/test_input/ctlgs",
        help="input (catalogue) directory"
    )
    parser.add_argument(
        '--output-dir', type=str,
        default="storage/deploy-test/output/stats",
        help="output (measurement) directory"
    )
    parser.add_argument(
        '--data-ctlgfile', type=str, default="test_data_catalogue.txt",
        help="data-source catalogue file"
    )
    parser.add_argument(
        '--rand-ctlgfile', type=str, default="test_rand_catalogue.txt",
        help="random-source catalogue file"
    )
    parser.add_argument(
        '--ctlg-fields', type=str, nargs='+', default=['x', 'y', 'z', 'nz'],
        help="catalogue fields (must contain sample-weight column 'Weight')"
    )
    parser.add_argument(
        '--output-tag', type=str,
        help="output file tag"
    )

    parser.add_argument(
        '--boxsize', type=float, default=1000.,
        help="box size for all dimensions"
    )
    parser.add_argument(
        '--ngrid', type=int, default=64,
        help="mesh grid number"
    )
    parser.add_argument(
        '--scheme', type=str, default='tsc',
        help="mesh assignment scheme"
    )
    parser.add_argument(
        '--interlace', action='store_true',
        help="enable interlacing"
    )

    parser.add_argument(
        '--case', type=str, choices=['survey', 'sim', 'random'],
        help="catalogue type/case"
    )
    parser.add_argument(
        '--statistic', type=str, choices=[
            'powspec', '2pcf',
            'bispec', '3pcf',
        ],
        help="statistic type"
    )
    parser.add_argument(
        '--multipole', type=str,
        help="multipole degree(s)"
    )
    parser.add_argument(
        '--range', type=float, nargs=2,
        help="wavenumber or separation range (in h/Mpc or Mpc/h)"
    )
    parser.add_argument(
        '--num-bins', type=int, default=4,
        help="number of bins"
    )
    parser.add_argument(
        '--form', type=str, choices=['diag', 'off-diag', 'row', 'full'],
        default='diag',
        help="3-point statistic form"
    )

    config = parser.parse_args()

    return config


def setup_params(config):
    """Set up parameter set.

    Parameters
    ----------
    config : :class:`argparse.namespace`
        Integration testing configuration.

    Returns
    -------
    :class:`triumvirate.parameters.ParameterSet`
        Integration test parameter set.

    """
    param_dict = fetch_paramset_template('dict')

    param_dict['directories']['measurements'] = config.output_dir
    param_dict['tags']['output'] = config.output_tag

    for ax_name in ['x', 'y', 'z']:
        param_dict['boxsize'][ax_name] = config.boxsize
        param_dict['ngrid'][ax_name] = config.ngrid

    param_dict.update({
        'scheme': config.scheme,
        'interlace': config.interlace,
        'catalogue_type': config.case,
        'statistic_type': config.statistic,
        'form': config.form,
        'range': config.range,
        'num_bins': config.num_bins,
        'idx_bin': 0,
    })

    if len(config.multipole) == 1:
        param_dict['degrees'] = {
            'ell1': None, 'ell2': None, 'ELL': int(config.multipole),
        }
    elif len(config.multipole) == 3:
        param_dict['degrees'] = {
            'ell1': int(config.multipole[0]),
            'ell2': int(config.multipole[1]),
            'ELL': int(config.multipole[2]),
        }
    else:
        raise ValueError("Invalid multipole: {}.".format(config.multipole))

    return ParameterSet(param_dict=param_dict, logger=logger)


if __name__ == '__main__':
    # Set up test.
    config = setup_config()
    logger = setup_logger()
    params = setup_params(config)

    # Determine and run test case.
    if params['catalogue_type'] == 'survey':
        ctlg_data = ParticleCatalogue.read_from_file(
            os.path.join(config.input_dir, config.data_ctlgfile),
            names=config.ctlg_fields,
            logger=logger
        )
        ctlg_rand = ParticleCatalogue.read_from_file(
            os.path.join(config.input_dir, config.rand_ctlgfile),
            names=config.ctlg_fields,
            logger=logger
        )
    elif params['catalogue_type'] == 'sim':
        ctlg_data = ParticleCatalogue.read_from_file(
            os.path.join(config.input_dir, config.data_ctlgfile),
            names=config.ctlg_fields,
            logger=logger
        )
    elif params['catalogue_type'] == 'random':
        ctlg_rand = ParticleCatalogue.read_from_file(
            os.path.join(config.input_dir, config.rand_ctlgfile),
            names=config.ctlg_fields,
            logger=logger
        )
    else:
        raise ValueError("Invalid case: {}.".format(config.case))

    if params['statistic_type'] == 'powspec' \
            and params['catalogue_type'] == 'survey':
        compute_powspec(
            ctlg_data, ctlg_rand, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == 'powspec' \
            and params['catalogue_type'] == 'sim':
        compute_powspec_in_gpp_box(
            ctlg_data, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == '2pcf' \
            and params['catalogue_type'] == 'survey':
        compute_corrfunc(
            ctlg_data, ctlg_rand, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == '2pcf' \
            and params['catalogue_type'] == 'sim':
        compute_corrfunc_in_gpp_box(
            ctlg_data, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == '2pcf' \
            and params['catalogue_type'] == 'random':
        compute_corrfunc_window(
            ctlg_rand, paramset=params, save='.txt', logger=logger
        )

    if params['statistic_type'] == 'bispec' \
            and params['catalogue_type'] == 'survey':
        compute_bispec(
            ctlg_data, ctlg_rand, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == 'bispec' \
            and params['catalogue_type'] == 'sim':
        compute_bispec_in_gpp_box(
            ctlg_data, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == '3pcf' \
            and params['catalogue_type'] == 'survey':
        compute_3pcf(
            ctlg_data, ctlg_rand, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == '3pcf' \
            and params['catalogue_type'] == 'sim':
        compute_3pcf_in_gpp_box(
            ctlg_data, paramset=params, save='.txt', logger=logger
        )
    if params['statistic_type'] == '3pcf' \
            and params['catalogue_type'] == 'random':
        compute_3pcf_window(
            ctlg_rand, paramset=params, save='.txt', logger=logger
        )
