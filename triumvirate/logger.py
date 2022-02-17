"""
Program Logger (:mod:`~triumvirate.logger`)
===========================================================================

Provide the program logger.

"""
import logging
import sys
import time

# TODO: Add C++ context to logger.


class _LogFormatter(logging.Formatter):
    """Customised logging formatter.

    """

    _start_time = time.time()

    def format(self, record):
        """Modify the default logging record by adding elapsed time in
        hours, minutes and seconds.

        Parameters
        ----------
        record : :class:`Logging.LogRecord`
            Default logging record object.

        Returns
        -------
        str
            Modified record message with elapsed time.

        """
        elapsed_time = record.created - self._start_time
        h, remainder_time = divmod(elapsed_time, 3600)
        m, s = divmod(remainder_time, 60)

        record.elapsed = "(+{}:{:02d}:{:02d})".format(int(h), int(m), int(s))

        return logging.Formatter.format(self, record)


def setup_logger():
    """Return a system-stream logger.

    Returns
    -------
    logger : :class:`logging.Logger`
        Logger.

    """
    # Clear logger handlers.
    logger = logging.getLogger()
    if logger.hasHandlers():
        logger.handlers.clear()

    # Set logger handler and its formatter.
    logging_formatter = _LogFormatter(
        fmt='[%(asctime)s %(elapsed)s %(levelname)s]%(incpp_state)s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging_handler = logging.StreamHandler(sys.stdout)
    logging_handler.setFormatter(logging_formatter)

    logger.addHandler(logging_handler)

    return logger
