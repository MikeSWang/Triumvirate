# cython: c_string_type=unicode, c_string_encoding=utf8
"""
Parameter Set (:mod:`~triumvirate.parameters`)
===========================================================================

Configure program parameter set.

"""
from os.path import abspath
from pprint import pformat

import numpy as np
import yaml


cdef class ParameterSet:
    """Program parameter set.

    The value of a parameter is accessed through its key in the same
    way as for :type:`dict`.

    Parameters
    ----------
    filepath : str or :class:`pathlib.Path`
        Parameter file path.
    logger : :class:`logging.Logger`, optional
        Program logger (default is `None`).

    """

    # Add `args` and `kwargs` for Cython subclass 'init'-compatibility.
    def __cinit__(self, filepath, logger=None, *args, **kwargs):

        self.thisptr = new CppParameterSet()

        self._logger = logger

        self._source = abspath(filepath)
        self._status = 'original'.encode('utf-8')

        with open(filepath, 'r') as filestream:
            self._params = yaml.load(filestream, Loader=yaml.Loader)

        self._parse_attrs()

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        return "\n".join([
            f"Source: {self._source}",
            f"Status: {self._status}",
            f"Params: {pformat(self._params, sort_dicts=False)}",
        ])

    def __len__(self):
        return len(self._params)

    def __getitem__(self, key):
        return self._params[key]

    def __setitem__(self, key, val):
        self._params.update({key: val})
        self._status = 'modified'.encode('utf-8')

    def printout(self, filepath=None):
        """Print out validated parameters to a file.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`, optional
            Printout file path.  If `None` (default), parameters are
            printed out to a default file in measurement directory.

        """
        if filepath:
            with open(filepath, 'w') as out_file:
                yaml.dump(
                    self._params, out_file,
                    sort_keys=False, default_flow_style=False
                )

            if self._logger:
                self._logger.info("Printed out parameters to %s.", filepath)
        else:
            if self._logger:
                try:
                    self._logger.info("", cpp_state='start')
                except:
                    self._logger.info("Entering C++ run...")

            if self.thisptr.printout() != 0:
                raise RuntimeError(
                    "Failed to print out extracted parameters to file."
                )

            if self._logger:
                try:
                    self._logger.info("", cpp_state='end')
                except:
                    self._logger.info("... exited C++ run.")

                self._logger.info(
                    "Printed out parameters to measurement directory."
                )

    def _parse_attrs(self):
        """Parse the parameter set into wrapped C++ parameter class.

        """
        # RFE: Implement serialisation of catalogue files for
        # e.g. suite of simulations.
        self.thisptr.catalogue_dir = \
            self._params['directories']['catalogues']
        self.thisptr.measurement_dir = \
            self._params['directories']['measurements']
        self.thisptr.data_catalogue_file = \
            self.thisptr.catalogue_dir \
            + self._params['filenames']['data_catalogue']
        self.thisptr.rand_catalogue_file = \
            self.thisptr.catalogue_dir \
            + self._params['filenames']['rand_catalogue']
        self.thisptr.output_tag = \
            self._params['tags']['output']

        self.thisptr.boxsize = [
            self._params['boxsize']['x'],
            self._params['boxsize']['y'],
            self._params['boxsize']['z'],
        ]
        self.thisptr.nmesh = [
            self._params['nmesh']['x'],
            self._params['nmesh']['y'],
            self._params['nmesh']['z'],
        ]

        self.thisptr.assignment = self._params['assignment']
        self.thisptr.norm_convention = self._params['norm_convention']

        self.thisptr.catalogue_type = self._params['catalogue_type']
        self.thisptr.measurement_type = self._params['measurement_type']

        self.thisptr.ell1 = self._params['degrees']['ell1']
        self.thisptr.ell2 = self._params['degrees']['ell2']
        self.thisptr.ELL = self._params['degrees']['ELL']

        self.thisptr.i_wa = self._params['wa_orders']['i']
        self.thisptr.j_wa = self._params['wa_orders']['j']

        self.thisptr.binning = self._params['binning']
        self.thisptr.form = self._params['form']

        if 'spec' in str(self.thisptr.measurement_type):
            self.thisptr.kmin = self._params['range'][0]
            self.thisptr.kmax = self._params['range'][1]
            self.thisptr.num_kbin = self._params['dim']
            self.thisptr.num_rbin = self._params['dim']  # placeholder
            self.thisptr.ith_kbin = self._params['index']
        if 'pcf' in str(self.thisptr.measurement_type):
            self.thisptr.rmin = self._params['range'][0]
            self.thisptr.rmax = self._params['range'][1]
            self.thisptr.num_rbin = self._params['dim']
            self.thisptr.num_kbin = self._params['dim']  # placeholder
            self.thisptr.ith_rbin = self._params['index']

        self.thisptr.volume = np.prod(list(self._params['boxsize'].values()))
        self.thisptr.nmesh_tot = np.prod(list(self._params['nmesh'].values()))

        self._validate()

    def _validate(self):
        """Validate extracted parameters.

        """
        if self._logger:
            self._logger.info("Validating parameters...")
            try:
                self._logger.info("", cpp_state='start')
            except:
                self._logger.info("Entering C++ run...")

        if self.thisptr.validate() != 0:
            raise ValueError("Invalid measurement parameters.")

        if self._logger:
            try:
                self._logger.info("", cpp_state='end')
            except:
                self._logger.info("... exited C++ run.")
            self._logger.info("... validated parameters.")
