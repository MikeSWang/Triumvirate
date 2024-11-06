"""
Program Logger (:mod:`~triumvirate.logger`)
===========================================

Configure the program logger.

.. autosummary::
    setup_logger

"""  # numpydoc ignore=SS01
import logging
import os
import sys
import time
import warnings
from copy import copy


class _ElapsedColourLogFormatter(logging.Formatter):
    """Elapsed-time colourised logging formatter.

    """
    # Log level colours
    _IN_BOLD = 1  # default to: 0 for non-bold, 1 for bold
    _COLOURS = {
        'NOTSET': "\033[0m",                    # reset
        'DEBUG': f"\033[{_IN_BOLD};36m",        # cyan
        'INFO': f"\033[{_IN_BOLD};32m",         # green (blue)
        'WARNING': f"\033[{_IN_BOLD};33m",      # yellow
        'ERROR': f"\033[{_IN_BOLD};31m",        # red
        'CRITICAL': f"\033[{_IN_BOLD};37;41m",  # white on red
    }

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
        fmt_record = copy(record)

        elapsed_time = record.created - self._start_time
        h, remainder_time = divmod(elapsed_time, 3600)
        m, s = divmod(remainder_time, 60)

        fmt_record.elapsed = f"(+{int(h):02d}:{int(m):02d}:{int(s):02d})"

        reset = self._COLOURS.get('NOTSET', "\033[0m")
        colour = self._COLOURS.get(record.levelname, reset)

        fmt_record.levelname = f"{record.levelname:.4s}"\
            .replace('DEBU', 'DBUG')\
            .replace('CRIT', 'FATAL')
        if 'color' in os.getenv('TERM', '') \
                and os.getenv('TRV_INTERACTIVE') is not None:
            fmt_record.levelname = colour + fmt_record.levelname + reset

        return logging.Formatter.format(self, fmt_record)


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
            Adapted log message with C++ runtime indication.

        """
        # Extract passed state variable or resort to default from `extra`.
        cpp_state = kwargs.pop('cpp_state', self.extra['cpp_state'])

        if isinstance(cpp_state, str):
            if cpp_state.lower() == 'start':
                return "%s (entering C++)" % msg, kwargs
            if cpp_state.lower() == 'end':
                return "%s (exited C++)" % msg, kwargs

        if cpp_state:  # merely indicating in CPP state, no extra text
            return "%s (in C++)" % msg, kwargs
        return "%s" % msg, kwargs


# Modify `warnings.formatwarning`.
def _format_warning(message, category, filename, lineno, line=None):
    """Modify formatter for warning messages in logged record.

    See Also
    --------
    :func:`warnings.formatwarning`
        For implementation details.

    """  # numpydoc ignore=PR01,RT01
    msg = warnings.WarningMessage(
        message, category, filename, lineno, None, line
    )

    msg_txt = "{} ({}:{}:{}{{}})".format(
        str(msg.message).replace('{', '{{').replace('}', '}}'),
        msg.filename, msg.lineno, msg.category.__name__
    )

    if msg.line is None:
        try:
            import linecache
            linesource = linecache.getline(msg.filename, msg.lineno)
        except Exception:
            linesource = None
    else:
        linesource = msg.line

    if linesource:
        linesource = linesource.strip()
        msg_txt = msg_txt.format(" --> \"" + linesource + "\"")
    else:
        msg_txt = msg_txt.format("")

    return msg_txt


def setup_logger(log_level=logging.INFO):
    """Set up and return a customised logger with elapsed time,
    C++ runtime indication and formatted warning messages.

    Parameters
    ----------
    log_level : int, optional
        Logging devel (default is `logging.INFO`).

    Returns
    -------
    logger : :class:`logging.LoggerAdapter`
        Customised logger.

    See Also
    --------
    :mod:`logging`
        For more details of the Python logging facility.

    """
    # Set formatter.
    formatter = _ElapsedColourLogFormatter(
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

    logger_.setLevel(log_level)

    # Adapt logger for C++ code indication.
    logger = _CppLogAdapter(logger_, {'cpp_state': False})

    # Adapt logger to capture warnings.
    warnings.formatwarning = _format_warning
    logging.captureWarnings(True)

    return logger
