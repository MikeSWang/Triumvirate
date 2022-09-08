# cython: c_string_type=unicode, c_string_encoding=utf8
"""
Parameter Set (:mod:`~triumvirate.parameters`)
===========================================================================

Configure program parameter set.

"""
from os.path import abspath
from pathlib import Path
from pprint import pformat

import numpy as np
import yaml


class InvalidParameter(ValueError):
    """:exc:`ValueError` raised when a program parameter is invalid.

    """
    pass


cdef class ParameterSet:
    """Program parameter set.

    Parameters
    ----------
    param_file : str or :class:`pathlib.Path`
        Parameter file path.
    logger : :class:`logging.Logger`, optional
        Program logger (default is `None`).

    Notes
    -----
    `param_file` should be point to a YAML file.  See
    :func:`~triumvirate.parameters.show_param_template()` to
    generate an example file.

    """

    def __cinit__(self, param_file, logger=None):

        self.thisptr = new CppParameterSet()

        self._logger = logger

        self._source = abspath(param_file)
        self._status = 'original'.encode('utf-8')
        self._validity = 'unvalidated'.encode('utf-8')

        with open(param_file, 'r') as pfile:
            self._params = yaml.load(pfile, Loader=yaml.Loader)

        self._parse_attrs()

    @classmethod
    def from_dict(cls, param_dict, logger=None):
        """Parameter set from a dictionary.

        Parameters
        ----------
        param_dict : dict
            Parameter dictionary.
        logger : :class:`logging.Logger`, optional
            Program logger (default is `None`).

        Notes
        -----
        `param_dict` is a nested dictionary.  See
        :func:`~triumvirate.parameters.show_param_template()` to
        generate an example `param_dict`.

        """
        self = object.__new__(cls)

        self._logger = logger

        self._source = 'dict'
        self._status = 'original'.encode('utf-8')
        self._validity = 'unvalidated'.encode('utf-8')

        self._params = param_dict

        self._parse_attrs()

    def __dealloc__(self):
        del self.thisptr

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
        self._validity = 'unvalidated'.encode('utf-8')

    def print_to_file(self, filepath=None):
        """Print out validated parameters to a YAML file.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`, optional
            Printout file path.  If `None` (default), parameters are
            printed out to a default file in measurement directory.

        """
        if filepath is None:
            filepath = (
                Path(self._params['directories']['measurements'])/
                "parameters_used{}".format(self._params['tags']['output'])
            )

        with open(filepath, 'w') as out_file:
            yaml.dump(
                self._params, out_file,
                sort_keys=False, default_flow_style=False
            )

        if self._logger:
            self._logger.info("Printed out parameters to %s.", filepath)

    def _parse_attrs(self):
        """Parse the parameter set into wrapped C++ parameter class.

        """
        # ----------------------------------------------------------------
        # I/O
        # ----------------------------------------------------------------

        try:
            if self._params['directories']['catalogues'] is None:
                self.thisptr.catalogue_dir = "".encode('utf-8')
            else:
                self.thisptr.catalogue_dir = \
                    self._params['directories']['catalogues']
        except KeyError:
            self.thisptr.catalogue_dir = "".encode('utf-8')

        try:
            if self._params['directories']['measurements'] is None:
                self.thisptr.measurement_dir = "".encode('utf-8')
            else:
                self.thisptr.measurement_dir = \
                    self._params['directories']['measurements']
        except KeyError:
            self.thisptr.measurement_dir = "".encode('utf-8')

        try:
            if self._params['filenames']['data_catalogue'] is None:
                self.thisptr.data_catalogue_file = "".encode('utf-8')
            else:
                self.thisptr.data_catalogue_file = \
                    self._params['filenames']['data_catalogue']
        except KeyError:
            self.thisptr.data_catalogue_file = "".encode('utf-8')

        try:
            if self._params['filenames']['rand_catalogue'] is None:
                self.thisptr.rand_catalogue_file = "".encode('utf-8')
            else:
                self.thisptr.rand_catalogue_file = \
                    self._params['filenames']['rand_catalogue']
        except KeyError:
            self.thisptr.rand_catalogue_file = "".encode('utf-8')

        try:
            if self._params['tags']['output'] is None:
                self.thisptr.output_tag = ''.encode('utf-8')
            else:
                self.thisptr.output_tag = self._params['tags']['output']
        except KeyError:
            self.thisptr.output_tag = ''.encode('utf-8')

        # ----------------------------------------------------------------
        # Mesh sampling
        # ----------------------------------------------------------------

        self.thisptr.boxsize = [
            float(self._params['boxsize']['x']),
            float(self._params['boxsize']['y']),
            float(self._params['boxsize']['z']),
        ]
        self.thisptr.ngrid = [
            self._params['ngrid']['x'],
            self._params['ngrid']['y'],
            self._params['ngrid']['z'],
        ]

        # Derived parameters.
        self.thisptr.volume = np.prod(list(self._params['boxsize'].values()))
        self.thisptr.nmesh = np.prod(list(self._params['ngrid'].values()))

        if self._params['alignment'] is not None:
            self.thisptr.alignment = self._params['alignment'].lower()
        if self._params['padscale'] is not None:
            self.thisptr.padscale = self._params['padscale'].lower()
        if self._params['padfactor'] is not None:
            self.thisptr.padfactor = self._params['padfactor']

        if self._params['assignment'] is not None:
            self.thisptr.assignment = self._params['assignment'].lower()
        if self._params['interlace'] is not None:
            self.thisptr.interlace = str(self._params['interlace']).lower()

        # ----------------------------------------------------------------
        # Measurement
        # ----------------------------------------------------------------

        if self._params['catalogue_type'] is not None:
            self.thisptr.catalogue_type = \
                self._params['catalogue_type'].lower()
        else:
            raise InvalidParameter("`catalogue_type` parameter must be set.")
        if self._params['measurement_type'] is not None:
            self.thisptr.measurement_type = \
                self._params['measurement_type'].lower()
        else:
            raise InvalidParameter("`measurement_type` parameter must be set.")

        if self._params['norm_convention'] is not None:
            self.thisptr.norm_convention = \
                self._params['norm_convention'].lower()
        if self._params['shotnoise_convention'] is not None:
            self.thisptr.shotnoise_convention = \
                self._params['shotnoise_convention'].lower()

        if self._params['binning'] is not None:
            self.thisptr.binning = self._params['binning'].lower()
        if self._params['form'] is not None:
            self.thisptr.form = self._params['form'].lower()

        if self._params['degrees']['ell1'] is not None:
            self.thisptr.ell1 = self._params['degrees']['ell1']
        if self._params['degrees']['ell2'] is not None:
            self.thisptr.ell2 = self._params['degrees']['ell2']
        if self._params['degrees']['ELL'] is not None:
            self.thisptr.ELL = self._params['degrees']['ELL']
        else:
            raise InvalidParameter("`ELL` parameter must be set.")

        if self._params['wa_orders']['i'] is not None:
            self.thisptr.i_wa = self._params['wa_orders']['i']
        if self._params['wa_orders']['j'] is not None:
            self.thisptr.j_wa = self._params['wa_orders']['j']

        self.thisptr.bin_min = float(self._params['range'][0])
        self.thisptr.bin_max = float(self._params['range'][1])
        if self._params['num_bins'] is not None:
            self.thisptr.num_bins = self._params['num_bins']
        else:
            raise InvalidParameter("`num_bins` parameter must be set.")
        if self._params['idx_bin'] is not None:
            self.thisptr.idx_bin = self._params['idx_bin']

        # ----------------------------------------------------------------
        # Misc
        # ----------------------------------------------------------------

        if self._params['verbose'] is not None:
            self.thisptr.verbose = self._params['verbose']

        self._validate()

    def _validate(self):
        """Validate extracted parameters.

        """
        try:
            self._logger.info("Validating parameters...", cpp_state='start')
        except (AttributeError, TypeError):
            pass

        if self.thisptr.validate() != 0:
            raise InvalidParameter("Invalid measurement parameters.")

        # Fetch validated parameters that have been derived or transmuted.
        # Transmuted I/O paths are not fetched.
        self._params['npoint'] = self.thisptr.npoint
        self._params['space'] = self.thisptr.space
        self._params['interlace'] = self.thisptr.interlace

        self._validity = 'validated'.encode('utf-8')

        try:
            self._logger.info("... validated parameters.", cpp_state='end')
        except (AttributeError, TypeError):
            pass


def show_param_template(format):
    """Show a parameter set template, either as text file contents or
    as a dictionary.

    Parameters
    ----------
    format : {'yaml', 'dict'}
        Template format, either file contents ('yaml')
        or dictionary ('dict').

    Returns
    -------
    template : str or dict
        Parameter set template.

    Raises
    ------
    ValueError
        If `format` is neither 'yaml' nor 'dict'.

    """
    pkg_root_dir = Path(__file__).resolve()
    tml_filepath = pkg_root_dir/"resources"/"params_example.yml"

    if format.lower() == 'yaml':
        template = Path(tml_filepath).read_text()
    elif format.lower() == 'dict':
        template = ParameterSet(tml_filepath)._params

    return template
