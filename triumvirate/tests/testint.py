"""Integration test.

Test cases include any supported clustering statistics with results saved.

"""
from triumvirate.logger import setup_logger
from triumvirate.parameters import ParameterSet
from triumvirate.catalogue import ParticleCatalogue
from triumvirate.twopt import (
    compute_corrfunc, compute_corrfunc_in_gpp_box, compute_corrfunc_window,
    compute_powspec, compute_powspec_in_gpp_box,
)
from triumvirate.threept import (
    compute_3pcf, compute_3pcf_in_gpp_box, compute_3pcf_window,
    compute_bispec, compute_bispec_in_gpp_box,
)


if __name__ == '__main__':
    # Set up logger.
    logger = setup_logger()

    # Set up parameters.
    pars = ParameterSet(
        "triumvirate/tests/test_input/params/test_params.yml", logger=logger
    )

    # Load test catalogues.
    if pars['catalogue_type'] == 'sim':
        cat_data = ParticleCatalogue.read_from_file(
            "{}/{}".format(
                pars['directories']['catalogues'],
                pars['files']['data_catalogue']
            ),
            names=['x', 'y', 'z', 'ws'],
            logger=logger
        )
    if pars['catalogue_type'] == 'survey':
        cat_data = ParticleCatalogue.read_from_file(
            "{}/{}".format(
                pars['directories']['catalogues'],
                pars['files']['data_catalogue']
            ),
            names=['x', 'y', 'z', 'nz', 'ws'],
            logger=logger
        )
    if pars['catalogue_type'] == 'survey' or pars['catalogue_type'] == 'random':
        cat_rand = ParticleCatalogue.read_from_file(
            "{}/{}".format(
                pars['directories']['catalogues'],
                pars['files']['rand_catalogue']
            ),
            names=['x', 'y', 'z', 'nz', 'ws'],
            logger=logger
        )

    # Run test case.
    if (pars['catalogue_type'], pars['statistic_type']) == ('sim', 'powspec'):
        measurements = compute_powspec_in_gpp_box(
            cat_data, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
    if (pars['catalogue_type'], pars['statistic_type']) == ('sim', '2pcf'):
        measurements = compute_corrfunc_in_gpp_box(
            cat_data, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
    if (pars['catalogue_type'], pars['statistic_type']) == ('sim', 'bispec'):
        measurements = compute_bispec_in_gpp_box(
            cat_data, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
    if (pars['catalogue_type'], pars['statistic_type']) == ('sim', '3pcf'):
        measurements = compute_3pcf_in_gpp_box(
            cat_data, paramset=pars, save='.txt', logger=logger
        )
        exit(0)

    if (pars['catalogue_type'], pars['statistic_type']) == ('survey', 'powspec'):
        measurements = compute_powspec(
            cat_data, cat_rand, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
    if (pars['catalogue_type'], pars['statistic_type']) == ('survey', '2pcf'):
        measurements = compute_corrfunc(
            cat_data, cat_rand, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
    if (pars['catalogue_type'], pars['statistic_type']) == ('survey', 'bispec'):
        measurements = compute_bispec(
            cat_data, cat_rand, paramset=pars, save='.txt', logger=logger
        )
    if (pars['catalogue_type'], pars['statistic_type']) == ('survey', '3pcf'):
        measurements = compute_3pcf(
            cat_data, cat_rand, paramset=pars, save='.txt', logger=logger
        )
        exit(0)

    if (pars['catalogue_type'], pars['statistic_type']) == ('random', '2pcf-win'):
        measurements = compute_corrfunc_window(
            cat_rand, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
    if (pars['catalogue_type'], pars['statistic_type']) == ('random', '3pcf-win'):
        measurements = compute_3pcf_window(
            cat_rand, paramset=pars, save='.txt', logger=logger
        )
        exit(0)
