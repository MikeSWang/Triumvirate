"""Configure `pytest`.

"""
import inspect
import warnings
from pathlib import Path

import numpy as np
import pytest

from triumvirate.catalogue import ParticleCatalogue
from triumvirate.dataobjs import Binning
from triumvirate.logger import setup_logger
from triumvirate.parameters import ParameterSet, fetch_paramset_template


# ========================================================================
# Paths
# ========================================================================

conftest_file = inspect.getframeinfo(inspect.currentframe()).filename


# ========================================================================
# Configuration
# ========================================================================

def pytest_addoption(parser):
    # Add '--runslow' option for slow test functions.
    parser.addoption(
        '--runslow', action='store_true', default=False, help="run slow tests"
    )


def pytest_configure(config):
    # Add 'slow' marker for tests to run only when running slow tests
    # is enabled.
    config.addinivalue_line('markers', "slow: mark test as slow to run")

    # Add 'no_capture' marker for tests to run only when output capturing
    # is disabled.
    config.addinivalue_line(
        'markers',
        "no_capture: mark test to run only when output capturing is disabled"
    )


def pytest_collection_modifyitems(config, items):
    # Skip tests marked with 'no_capture' when capturing is not disabled.
    if not config.getoption('--runslow'):
        skip_slow = pytest.mark.skip(
            reason="need '--runslow' option to run slow tests"
        )
        for item in items:
            if 'slow' in item.keywords:
                item.add_marker(skip_slow)
    if config.getoption('--capture') != 'no':
        skip_captured = pytest.mark.skip(
            reason="need '--capture=no' or '-s' option to run"
        )
        for item in items:
            if 'no_capture' in item.keywords:
                item.add_marker(skip_captured)


# ========================================================================
# Fixtures
# ========================================================================

@pytest.fixture(scope='session')
def test_dir():
    return Path(conftest_file).parent.resolve()


@pytest.fixture(scope='session')
def test_input_dir(test_dir):
    return test_dir/"test_input"


@pytest.fixture(scope='session')
def test_output_dir(test_dir):
    return test_dir/"test_output"


@pytest.fixture(scope='session')
def test_param_dir(test_input_dir):
    return test_input_dir/"params"


@pytest.fixture(scope='session')
def test_ctlg_dir(test_input_dir):
    return test_input_dir/"ctlgs"


@pytest.fixture(scope='session')
def test_stats_dir(test_input_dir):
    return test_input_dir/"stats"


@pytest.fixture
def test_logger():
    return setup_logger()


@pytest.fixture
def valid_paramset():
    param_dict = fetch_paramset_template('dict')

    for ax_name in ['x', 'y', 'z']:
        param_dict['boxsize'][ax_name] = 1000.
        param_dict['ngrid'][ax_name] = 64.

    param_dict.update({
        'statistic_type': 'bispec',
        'degrees': {'ell1': 0, 'ell2': 0, 'ELL': 0},
        'range': [0.005, 0.105],
        'num_bins': 10,
        'verbose': 60,
    })

    return ParameterSet(param_dict=param_dict)


@pytest.fixture
def test_paramset(test_param_dir):
    return ParameterSet(param_filepath=test_param_dir/"test_params.yml")


@pytest.fixture
def test_binning_fourier():
    return Binning('fourier', 'lin', bin_min=0.005, bin_max=0.105, num_bins=4)


@pytest.fixture
def test_binning_config():
    return Binning('config', 'lin', bin_min=50., bin_max=150., num_bins=4)


@pytest.fixture
def test_catalogue_properties():
    return {
        'separation': 100.,
        'boxsize': 1.e3,
        'nparticle': 3,
        'contrast': 1.e4,
    }


@pytest.fixture
def test_data_catalogue(test_ctlg_dir, test_logger, test_catalogue_properties):

    warnings.filterwarnings('ignore', message=".*field is not provided.*")

    if (test_ctlg_dir/"test_data_catalogue.txt").exists():
        return ParticleCatalogue.read_from_file(
            test_ctlg_dir/"test_data_catalogue.txt",
            names=['x', 'y', 'z', 'nz'],
            logger=test_logger
        )

    R = test_catalogue_properties['separation']
    L = test_catalogue_properties['boxsize']
    N = test_catalogue_properties['nparticle']

    x = [-R/2., 0., R/2.]
    y = [-R/2./np.sqrt(3), R/np.sqrt(3), -R/2./np.sqrt(3)]
    z = [0., 0., 0.]
    nz = N / L**3

    return ParticleCatalogue(x, y, z, nz=nz, logger=test_logger)


@pytest.fixture
def test_rand_catalogue(test_ctlg_dir, test_logger, test_catalogue_properties):

    warnings.filterwarnings('ignore', message=".*field is not provided.*")

    if (test_ctlg_dir/"test_rand_catalogue.txt").exists():
        return ParticleCatalogue.read_from_file(
            test_ctlg_dir/"test_rand_catalogue.txt",
            names=['x', 'y', 'z', 'nz'],
            logger=test_logger
        )

    alpha = test_catalogue_properties['contrast']
    L = test_catalogue_properties['boxsize']
    N = test_catalogue_properties['nparticle']

    generator = np.random.default_rng(seed=42)

    x, y, z = generator.uniform(-L/2., L/2., size=(3, int(alpha*N)))
    nz = N / L**3

    return ParticleCatalogue(x, y, z, nz=nz, logger=test_logger)


@pytest.fixture
def test_uniform_catalogue(test_ctlg_dir, test_logger,
                           test_catalogue_properties):

    warnings.filterwarnings('ignore', message=".*field is not provided.*")

    if (test_ctlg_dir/"test_uniform_catalogue.txt").exists():
        return ParticleCatalogue.read_from_file(
            test_ctlg_dir/"test_uniform_catalogue.txt",
            names=['x', 'y', 'z', 'nz'],
            logger=test_logger
        )

    alpha = test_catalogue_properties['contrast'] / 10.
    L = test_catalogue_properties['boxsize']
    N = test_catalogue_properties['nparticle']

    generator = np.random.default_rng(seed=69)

    x, y, z = generator.uniform(-L/2., L/2., size=(3, int(alpha*N)))
    nz = N / L**3

    return ParticleCatalogue(x, y, z, nz=nz, logger=test_logger)
