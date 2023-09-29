"""
Program Logger (:mod:`~triumvirate.logger`)
==========================================================================

Configure the program logger.

.. autosummary::
    setup_logger

"""
import logging
import sys
import time
import warnings


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

        record.elapsed = \
            "(+{:02d}:{:02d}:{:02d})".format(int(h), int(m), int(s))

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

    See :func:`warnings.formatwarning` for implementation details.

    """
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
    formatter = _ElapsedLogFormatter(
        fmt='[%(asctime)s %(elapsed)s %(levelname).4s] %(message)s',
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
