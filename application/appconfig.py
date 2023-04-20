"""Configure application programs.

"""
import logging
import sys
import time


class _LoggerFormatter(logging.Formatter):
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

        record.elapsed = \
            "(+{:02d}:{:02d}:{:02d})".format(int(h), int(m), int(s))

        return logging.Formatter.format(self, record)


def setup_logger(level=logging.INFO):
    """Return the root logger formatted with elapsed time and piped
    to ``stdout``.

    Parameters
    ----------
    level : int, optional
        Logging level (default is :const:`logging.INFO`).

    Returns
    -------
    logger : :class:`logging.Logger`
        Formatted root logger.

    """
    logger = logging.getLogger()

    logging_handler = logging.StreamHandler(sys.stdout)
    logging_formatter = _LoggerFormatter(
        fmt='[%(asctime)s %(elapsed)s %(levelname).4s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    logging_handler.setFormatter(logging_formatter)
    logger.addHandler(logging_handler)

    logger.setLevel(level)

    logging.captureWarnings(True)

    return logger
