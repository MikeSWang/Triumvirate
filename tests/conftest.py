"""Configure `pytest`.

"""
import inspect
from pathlib import Path

import pytest

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
def logger():
    return setup_logger()


@pytest.fixture(scope='session')
def valid_paramset(test_param_dir):
    return ParameterSet(param_filepath=test_param_dir/"valid_params.yml")
