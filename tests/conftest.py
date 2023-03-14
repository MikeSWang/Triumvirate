"""Configure `pytest`.

"""
import inspect
from pathlib import Path

import numpy as np
import pytest

from triumvirate.catalogue import ParticleCatalogue
from triumvirate.logger import setup_logger
from triumvirate.parameters import ParameterSet


# ========================================================================
# Paths
# ========================================================================

conftest_file = inspect.getframeinfo(inspect.currentframe()).filename

test_dir = Path(conftest_file).parent.resolve()


# ========================================================================
# Configuration
# ========================================================================

def pytest_configure(config):
    # Add 'no_capture' marker for tests to run only when output capturing
    # is disabled.
    config.addinivalue_line(
        'markers',
        "no_capture: mark test to run only when output capturing is disabled"
    )


def pytest_collection_modifyitems(config, items):
    # Skip tests marked with 'no_capture' when capturing is not disabled.
    if config.getoption('--capture') != 'no':
        skip_captured = pytest.mark.skip(
            reason="need '--capture=no' or '-s' option to run"
        )
        for item in items:
            if "no_capture" in item.keywords:
                item.add_marker(skip_captured)


# ========================================================================
# Fixtures
# ========================================================================

@pytest.fixture(scope='session')
def test_input_dir():
    return test_dir/"test_input"


@pytest.fixture(scope='session')
def test_output_dir():
    return test_dir/"test_output"


@pytest.fixture(scope='session')
def test_param_dir(test_input_dir):
    return test_input_dir/"params"


@pytest.fixture(scope='session')
def test_ctlg_dir(test_input_dir):
    return test_input_dir/"catalogues"


@pytest.fixture(scope='session')
def test_stats_dir(test_input_dir):
    return test_input_dir/"statistics"


@pytest.fixture(scope='session')
def logger():
    return setup_logger()


@pytest.fixture(scope='session')
def valid_paramset(test_param_dir):
    return ParameterSet(param_filepath=test_param_dir/"valid_params.yml")


@pytest.fixture(scope='session')
def test_catalogue_properties():
    return {
        'separation': 100.,
        'boxsize': 1.e3,
        'nparticle': 3,
        'contrast': 1.e4,
    }


@pytest.fixture(scope='session')
def test_data_catalogue(test_catalogue_properties):

    R = test_catalogue_properties['separation']
    L = test_catalogue_properties['boxsize']
    N = test_catalogue_properties['nparticle']

    x = [-R/2., 0., R/2.]
    y = [-R/2./np.sqrt(3), R/np.sqrt(3), -R/2./np.sqrt(3)]
    z = [0., 0., 0.]
    nz = N / L**3

    return ParticleCatalogue(x, y, z, nz=nz)


def test_rand_catalogue(test_catalogue_properties):

    alpha = test_catalogue_properties['contrast']
    L = test_catalogue_properties['boxsize']
    N = test_catalogue_properties['nparticle']

    generator = np.random.default_rng(seed=42)

    x, y, z = generator.uniform(-L/2., L/2., size=(3, int(alpha*N)))
    nz = N / L**3

    return ParticleCatalogue(x, y, z, nz=nz)
