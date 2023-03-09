"""Configure `pytest`.

"""
import pytest


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
        return
