import re

import pytest

from triumvirate.logger import setup_logger

# _LOG_FORMAT = (
#     r"\[\d{4}-\d{2}-\d{2}\s\d{2}:\d{2}:\d{2}\s\(\+\d+:\d{2}:\d{2}\)\s\w{4}\]"
#     r"\s.*\(in C\+\+\)"
# )


@pytest.fixture(scope='module')
def logger():
    return setup_logger()


def test_setup_logger(logger, caplog):

    # Test C++ code state indication.
    logger.info("Test log entry with C++ indication.", cpp_state=True)

    assert re.match(r'.*\(in C\+\+\)', caplog.text) is not None, \
        "No C++ indication shown by logger."
