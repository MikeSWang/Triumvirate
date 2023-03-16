"""Test :mod:`~triumvirate.logger`.

"""
import re

import pytest

from triumvirate.logger import setup_logger

# Define the RegEx expression for logging format.
_LOG_FORMAT = (
    r"\[\d{4}-\d{2}-\d{2}\s\d{2}:\d{2}:\d{2}\s\(\+\d+:\d{2}:\d{2}\)\s\w{4}\]"
    r"\s.*\(in C\+\+\)"
)


@pytest.fixture
def default_logger():
    return setup_logger()


@pytest.mark.no_capture
def test_setup_logger_formatter(default_logger, capfd):
    # Test logging formatter with C++ code state indication.
    default_logger.info("Test log entry.", cpp_state=True)
    assert re.match(_LOG_FORMAT, capfd.readouterr().out) is not None, \
        "Misformatted logging."


def test_setup_logger_adapter(default_logger, caplog):
    # Test logger C++ code state indication only.
    default_logger.info("Test log entry with C++ indication.", cpp_state=True)
    assert re.match(r'.*\(in C\+\+\)', caplog.text) is not None, \
        "No C++ indication shown by logger."
