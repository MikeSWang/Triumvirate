# cython: c_string_type=unicode, c_string_encoding=utf8
"""
Parameter Set (:mod:`~triumvirate.parameter_set`)
===========================================================================

Configure program parameter set.

"""
from os.path import abspath
from pprint import pformat

import yaml


cdef class ParameterSet:
    """Program parameter set.

    Parameters
    ----------
    filepath : str or :class:`pathlib.Path`
        Parameter file path.

    """

    # Add `args` and `kwargs` for Cython subclass 'init'-compatibility.
    def __cinit__(self, filepath, *args, **kwargs):

        self._source = abspath(filepath)
        self._status = 'original'.encode('utf-8')

        with open(filepath, 'r') as filestream:
            self._params = yaml.load(filestream, Loader=yaml.Loader)

    def __str__(self):
        return "\n".join([
            f"Source: {self._source}",
            f"Status: {self._status}",
            f"Params: {pformat(self._params, sort_dicts=False)}",
        ])

    def __getitem__(self, key):
        return self._params[key]

    def __setitem__(self, key, val):
        self._params.update({key: val})
        self._status = 'modified'.encode('utf-8')

    def print_to_file(self, filepath):
        """Print validated parameters to a file.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`
            Output file path.

        """
        with open(filepath, 'w') as outfile:
            yaml.dump(
                self._params, outfile,
                sort_keys=False, default_flow_style=False
            )
