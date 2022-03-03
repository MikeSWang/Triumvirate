"""
Program Logger (:mod:`~triumvirate.logger`)
===========================================================================

Configure the program logger.

"""
import logging
import sys
import time

from numpy import isin


class _ElapsedLogFormatter(logging.Formatter):
    """Elapsed-time logging formatter.

    """

    _start_time = time.time()

    def format(self, record):
        """Add elapsed time in hours, minutes and seconds to the log
        record.

        Parameters
        ----------
        record : :class:`Logging.LogRecord`
            See :class:`logging.LogRecord`.

        Returns
        -------
        str
            Modified log record with elapsed time.

        """
        elapsed_time = record.created - self._start_time
        h, remainder_time = divmod(elapsed_time, 3600)
        m, s = divmod(remainder_time, 60)

        record.elapsed = "(+{}:{:02d}:{:02d})".format(int(h), int(m), int(s))

        return logging.Formatter.format(self, record)


class _CppLogAdapter(logging.LoggerAdapter):
    """C++ logging adapter.

    """

    def process(self, msg, kwargs):
        """Adapt logger message.

        Parameters
        ----------
        msg : str
            See :class:`logging.LoggerAdapter`.
        kwargs : dict
            See :class:`logging.LoggerAdapter`.

        Returns
        -------
        str
            Adapted log message with C++ state indication.

        """
        # Extract passed state variable or resort to default from `extra`.
        cpp_state = kwargs.pop('cpp_state', self.extra['cpp_state'])

        if isinstance(cpp_state, str):
            if cpp_state.lower() == 'start':
                return "(C++ start) %s" % msg, kwargs
            if cpp_state.lower() == 'end':
                return "(C++ end) %s" % msg, kwargs

        if cpp_state:
            return "(C++) %s" % msg, kwargs
        return "%s" % msg, kwargs


def setup_logger():
    """Set up and return a customised logger with elapsed time and
    C++ state indication.

    Returns
    -------
    logger : :class:`logging.LoggerAdapter`
        Customised logger.

    """
    # Set formatter.
    formatter = _ElapsedLogFormatter(
        fmt='[%(asctime)s %(elapsed)s %(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

    # Set handler.
    handler = logging.StreamHandler(sys.stdout)
    handler.setFormatter(formatter)

    # Instantiate logger.
    logger_ = logging.getLogger()

    if logger_.hasHandlers():
        logger_.handlers.clear()
    logger_.addHandler(handler)

    logger_.setLevel(logging.INFO)

    # Adapt logger for C++ code indication.
    logger = _CppLogAdapter(logger_, {'cpp_state': False})

    return logger
