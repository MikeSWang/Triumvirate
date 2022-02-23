# cython: c_string_type=unicode, c_string_encoding=utf8
"""
Measurement Parameters (:mod:`~triumvirate.parameters`)
===========================================================================

Load measurement parameters.

"""
from cython.operator import dereference as deref

import numpy as np
import yaml


cdef class Parameters(ParameterSet):
    """Measurement parameters.

    Parameters
    ----------
    filepath : str or :class:`pathlib.Path`
        Parameter file path.

    """

    def __cinit__(self, filepath):

        super().__init__()

        self.thisptr = new CppParameters()

        self._parse_attrs()

    def printout(self):
        """Print out validated parameters to a file in measurement
        directory.

        """
        exit_status = self.thisptr.printout()
        if exit_status != 0:
            raise RuntimeError(
                "Failed to print out extracted parameters to file."
            )

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

    def _parse_attrs(self):

        self.thisptr.catalogue_dir = \
            self._params['directories']['catalogues']
        self.thisptr.measurement_dir = \
            self._params['directories']['measurements']
        self.thisptr.data_catalogue_file = \
            self._params['filenames']['data_catalogue']
        self.thisptr.rand_catalogue_file = \
            self._params['filenames']['rand_catalogue']
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
            self.thisptr.num_rbin = self._params['dim']
            self.thisptr.ith_kbin = self._params['index']
        if 'pcf' in str(self.thisptr.measurement_type):
            self.thisptr.rmin = self._params['range'][0]
            self.thisptr.rmax = self._params['range'][1]
            self.thisptr.num_rbin = self._params['dim']
            self.thisptr.num_kbin = self._params['dim']
            self.thisptr.ith_rbin = self._params['index']

        self.thisptr.volume = np.prod(list(self._params['boxsize'].values()))
        self.thisptr.nmesh_tot = np.prod(list(self._params['nmesh'].values()))

        self._validate()

    def _validate(self):
        """Validate extracted parameters.

        """
        exit_status = self.thisptr.validate()
        if exit_status != 0:
            raise ValueError("Invalid measurement parameters.")
