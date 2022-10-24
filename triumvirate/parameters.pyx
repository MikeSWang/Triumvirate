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
    """Parameter set.

    This reads parameters from a file, stores and prints out the extracted
    parameters, and validates the parameters.

    Parameters
    ----------
    param_filepath : str or :class:`pathlib.Path`
        Parameter file path.  This should point to a YAML file.
    logger : :class:`logging.Logger`, optional
        Program logger (default is `None`).

    Notes
    -----
    Use :func:`~triumvirate.parameters.show_paramset_template()` to
    generate a template parameter file (with explanatory text).

    """

    def __cinit__(self, param_filepath, logger=None):

        self.thisptr = new CppParameterSet()

        self._logger = logger

        self._source = abspath(param_file)
        self._original = True
        self._validity = False

        with open(param_filepath, 'r') as param_file:
            self._params = yaml.load(param_file, Loader=yaml.Loader)

        self._parse_params()  # validate

    @classmethod
    def from_dict(cls, param_dict, logger=None):
        """Create p arameter set from a dictionary.

        Parameters
        ----------
        param_dict : dict
            Parameter dictionary (nested).
        logger : :class:`logging.Logger`, optional
            Program logger (default is `None`).

        Notes
        -----
        Use :func:`~triumvirate.parameters.show_paramset_template()` to
        generate a template parameter dictionary.

        """
        self = object.__new__(cls)

        self._logger = logger

        self._source = 'dict'
        self._original = True
        self._validity = False

        self._params = param_dict

        self._parse_params()  # validate

    def __dealloc__(self):
        del self.thisptr

    def __str__(self):
        return "\n".join([
            f"ParameterSet(source={self._source}, original={self._original})",
            f"{pformat(self._params, sort_dicts=False)}",
        ])

    def __getitem__(self, key):
        return self._params[key]

    def __setitem__(self, key, val):
        self._params.update({key: val})
        self._original = False
        self._validity = False

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

        with open(filepath, 'w') as printout_file:
            yaml.dump(
                self._params, printout_file,
                sort_keys=False,
                default_flow_style=False  # CAVEAT: discretionary
            )

        if self._logger:
            self._logger.info("Printed out parameters to %s.", filepath)

    def _parse_params(self):
        """Parse the parameter set into the corresponding
        C++ parameter class.

        """
        # ----------------------------------------------------------------
        # Parsing
        # ----------------------------------------------------------------

        # -- I/O ---------------------------------------------------------

        # Directories, files and tags default to empty if unset.
        # Full path concatenation is handled by C++.

        # Attribute string parameters.
        try:
            if self._params['directories']['catalogues'] is None:
                self.thisptr.catalogue_dir = ""#.encode('utf-8')
            else:
                self.thisptr.catalogue_dir = \
                    self._params['directories']['catalogues']
        except KeyError:
            self.thisptr.catalogue_dir = ""#.encode('utf-8')

        try:
            if self._params['directories']['measurements'] is None:
                self.thisptr.measurement_dir = ""#.encode('utf-8')
            else:
                self.thisptr.measurement_dir = \
                    self._params['directories']['measurements']
        except KeyError:
            self.thisptr.measurement_dir = ""#.encode('utf-8')

        try:
            if self._params['files']['data_catalogue'] is None:
                self.thisptr.data_catalogue_file = ""#.encode('utf-8')
            else:
                self.thisptr.data_catalogue_file = \
                    self._params['files']['data_catalogue']
        except KeyError:
            self.thisptr.data_catalogue_file = ""#.encode('utf-8')

        try:
            if self._params['files']['rand_catalogue'] is None:
                self.thisptr.rand_catalogue_file = ""#.encode('utf-8')
            else:
                self.thisptr.rand_catalogue_file = \
                    self._params['files']['rand_catalogue']
        except KeyError:
            self.thisptr.rand_catalogue_file = ""#.encode('utf-8')

        try:
            if self._params['tags']['output'] is None:
                self.thisptr.output_tag = ''#.encode('utf-8')
            else:
                self.thisptr.output_tag = self._params['tags']['output']
        except KeyError:
            self.thisptr.output_tag = ''#.encode('utf-8')

        # -- Mesh sampling -----------------------------------------------

        # Attribute numerical parameters.
        if None in self._params['boxsize'].values():
            raise InvalidParameter("`boxsize` parameters must be set.")
        self.thisptr.boxsize = [
            float(self._params['boxsize']['x']),
            float(self._params['boxsize']['y']),
            float(self._params['boxsize']['z']),
        ]
        if None in self._params['ngrid'].values():
            raise InvalidParameter("`ngrid` parameters must be set.")
        self.thisptr.ngrid = [
            self._params['ngrid']['x'],
            self._params['ngrid']['y'],
            self._params['ngrid']['z'],
        ]

        if self._params['padfactor'] is not None:
            self.thisptr.padfactor = self._params['padfactor']

        # Attribute string parameters.
        if self._params['alignment'] is not None:
            self.thisptr.alignment = self._params['alignment'].lower()
        if self._params['padscale'] is not None:
            self.thisptr.padscale = self._params['padscale'].lower()
        if self._params['assignment'] is not None:
            self.thisptr.assignment = self._params['assignment'].lower()
        if self._params['interlace'] is not None:  # possibly convert from bool
            self.thisptr.interlace = str(self._params['interlace']).lower()

        # Attribute derived parameters.
        self.thisptr.volume = np.prod(list(self._params['boxsize'].values()))
        self.thisptr.nmesh = np.prod(list(self._params['ngrid'].values()))

        # -- Measurement -------------------------------------------------

        # Attribute numerical parameters.
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

        if None in self._params['range'].values():
            raise InvalidParameter("`range` parameters must be set.")
        self.thisptr.bin_min = float(self._params['range'][0])
        self.thisptr.bin_max = float(self._params['range'][1])

        if self._params['num_bins'] is not None:
            self.thisptr.num_bins = self._params['num_bins']
        else:
            raise InvalidParameter("`num_bins` parameter must be set.")
        if self._params['idx_bin'] is not None:
            self.thisptr.idx_bin = self._params['idx_bin']

        # Attribute string parameters.
        if self._params['catalogue_type'] is not None:
            self.thisptr.catalogue_type = \
                self._params['catalogue_type'].lower()
        else:
            raise InvalidParameter("`catalogue_type` parameter must be set.")
        if self._params['statistic_type'] is not None:
            self.thisptr.statistic_type = \
                self._params['statistic_type'].lower()
        else:
            raise InvalidParameter("`statistic_type` parameter must be set.")

        if self._params['norm_convention'] is not None:
            self.thisptr.norm_convention = \
                self._params['norm_convention'].lower()

        if self._params['binning'] is not None:
            self.thisptr.binning = self._params['binning'].lower()
        if self._params['form'] is not None:
            self.thisptr.form = self._params['form'].lower()

        # -- Misc --------------------------------------------------------

        if self._params['verbose'] is not None:
            self.thisptr.verbose = self._params['verbose']

        # ----------------------------------------------------------------
        # Validation
        # ----------------------------------------------------------------

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

        self._validity = True

        try:
            self._logger.info("... validated parameters.", cpp_state='end')
        except (AttributeError, TypeError):
            pass


def show_paramset_template(format):
    """Show a parameter set template, either as a YAML file
    (with explanatory text) or as a dictionary.

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
