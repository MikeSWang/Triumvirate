# cython: c_string_type=unicode, c_string_encoding=utf8
"""
Parameter Set Configuration (:mod:`~triumvirate.config.parameters`)
===========================================================================

Configure program parameter set(s).

"""
from os.path import abspath
from pprint import pformat

import yaml


cdef class ParameterSet:
    """Program parameters.

    Parameters
    ----------
    filepath : str or :class:`pathlib.Path`
        Parameter file path.

    """

    def __cinit__(self, filepath, *args, **kwargs):

        self._source = abspath(filepath)
        self._status = 'original'.encode('utf-8')

        with open(filepath, 'r') as file_stream:
            self._params = yaml.load(file_stream, Loader=yaml.Loader)

    def __str__(self):
        return "\n".join([
            f"Source: {self._source}",
            f"Status: {self._status}",
            f"Values: {pformat(self._params, sort_dicts=False)}",
        ])

    def __getitem__(self, key):
        return self._params[key]

    def __setitem__(self, key, val):
        self._params.update({key: val})
        self._status = 'modified'.encode('utf-8')
