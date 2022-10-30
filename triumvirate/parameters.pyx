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

from parameters cimport CppParameterSet


# NOTE: Nested entries must not have a non-None default value.
_DEFAULT_PARAM_DICT = {
    'directories': {
        'catalogues': None,
        'measurements': None,
    },
    'files': {
        'data_catalogue': None,
        'rand_catalogue': None,
    },
    'tags': {
        'output': None,
    },
    'boxsize': {'x': None, 'y': None, 'z': None},
    'ngrid': {'x': None, 'y': None, 'z': None},
    'alignment': 'centre',
    'padscale': 'box',
    'padfactor': None,
    'assignment': 'tsc',
    'interlace': False,
    'catalogue_type': None,
    'statistic_type': None,
    'norm_convention': 'particle',
    'binning': 'lin',
    'form': 'diag',
    'degrees': {'ell1': None, 'ell2': None, 'ELL': None},
    'wa_orders': {'i': None, 'j': None},
    'range': [None, None],
    'num_bins': None,
    'idx_bin': None,
    'verbose': 20,
}


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
    param_filepath : str or :class:`pathlib.Path`, optional
        Parameter file path.  This should point to a YAML file.
        If :keyword:`None` (default), :param:`param_dict` should be
        provided; otherwise :param:`param_dict` should be :keyword:`None`.
    param_dict : dict, optional
        Parameter dictionary (nested).  If :keyword:`None` (default),
        :param:`param_filepath` should be provided; otherwise
        :param:`param_filepath` should be :keyword:`None`.
    logger : :class:`logging.Logger`, optional
        Program logger (default is :keyword:`None`).

    Raises
    ------
    ValueError
        If neither or both of :param:`param_filepath` and
        :param:`param_dict` is/are :keyword:`None`.

    Notes
    -----
    Use :func:`~triumvirate.parameters.fetch_paramset_template()` to
    generate a template parameter file (with explanatory text).

    """

    def __cinit__(self, param_filepath=None, param_dict=None, logger=None):

        self.thisptr = new CppParameterSet()

        self._logger = logger

        if (param_filepath is not None and param_dict is not None) \
                or (param_filepath is None and param_dict is None):
            raise ValueError(
                "One and only one of `param_filepath` and `param_dict` "
                "should be set."
            )

        self._params = {}
        if param_filepath is not None:
            self._source = abspath(param_filepath)
            with open(param_filepath, 'r') as param_file:
                self._params = yaml.load(param_file, Loader=yaml.Loader)
        if param_dict is not None:
            self._source = 'dict'
            self._params = param_dict

        self._original = True
        self._validity = False

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

    def items(self):
        """Return the set of entries like a :class:`dict`.

        Returns
        -------
        :class:`dict_items`
            Parameters as dictionary items.

        """
        return self._params.items()

    def get(self, key):
        """Return a possibly non-existent entry like a :class:`dict`.

        Parameters
        ----------
        key : str
            Parameter name, possibly non-existent.

        Returns
        -------
        object
            Parameter value, :keyword:`None` if non-existent.

        """
        return self._params.get(key)

    def print_to_file(self, filepath=None):
        """Print out validated parameters to a YAML file.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`, optional
            Printout file path.  If :keyword:`None` (default), parameters
            are printed out to a default file in measurement directory.

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
                self.thisptr.catalogue_dir = "".encode('utf-8')
            else:
                self.thisptr.catalogue_dir = \
                    self._params['directories']['catalogues'].encode('utf-8')
        except KeyError:
            self.thisptr.catalogue_dir = "".encode('utf-8')

        try:
            if self._params['directories']['measurements'] is None:
                self.thisptr.measurement_dir = "".encode('utf-8')
            else:
                self.thisptr.measurement_dir = \
                    self._params['directories']['measurements'].encode('utf-8')
        except KeyError:
            self.thisptr.measurement_dir = "".encode('utf-8')

        try:
            if self._params['files']['data_catalogue'] is None:
                self.thisptr.data_catalogue_file = "".encode('utf-8')
            else:
                self.thisptr.data_catalogue_file = \
                    self._params['files']['data_catalogue'].encode('utf-8')
        except KeyError:
            self.thisptr.data_catalogue_file = "".encode('utf-8')

        try:
            if self._params['files']['rand_catalogue'] is None:
                self.thisptr.rand_catalogue_file = "".encode('utf-8')
            else:
                self.thisptr.rand_catalogue_file = \
                    self._params['files']['rand_catalogue'].encode('utf-8')
        except KeyError:
            self.thisptr.rand_catalogue_file = "".encode('utf-8')

        try:
            if self._params['tags']['output'] is None:
                self.thisptr.output_tag = ''.encode('utf-8')
            else:
                self.thisptr.output_tag = \
                    self._params['tags']['output'].encode('utf-8')
        except KeyError:
            self.thisptr.output_tag = ''.encode('utf-8')

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
            self.thisptr.alignment = \
                self._params['alignment'].lower().encode('utf-8')
        if self._params['padscale'] is not None:
            self.thisptr.padscale = \
                self._params['padscale'].lower().encode('utf-8')
        if self._params['assignment'] is not None:
            self.thisptr.assignment = \
                self._params['assignment'].lower().encode('utf-8')
        if self._params['interlace'] is not None:  # possibly convert from bool
            self.thisptr.interlace = \
                str(self._params['interlace']).lower().encode('utf-8')

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

        if None in self._params['range']:
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
                self._params['catalogue_type'].lower().encode('utf-8')
        else:
            raise InvalidParameter("`catalogue_type` parameter must be set.")
        if self._params['statistic_type'] is not None:
            self.thisptr.statistic_type = \
                self._params['statistic_type'].lower().encode('utf-8')
        else:
            raise InvalidParameter("`statistic_type` parameter must be set.")

        if self._params['norm_convention'] is not None:
            self.thisptr.norm_convention = \
                self._params['norm_convention'].lower().encode('utf-8')

        if self._params['binning'] is not None:
            self.thisptr.binning = \
                self._params['binning'].lower().encode('utf-8')
        if self._params['form'] is not None:
            self.thisptr.form = \
                self._params['form'].lower().encode('utf-8')

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
        self._params['npoint'] = self.thisptr.npoint.decode('utf-8')
        self._params['space'] = self.thisptr.space.decode('utf-8')
        self._params['interlace'] = self.thisptr.interlace.decode('utf-8')

        self._validity = True

        try:
            self._logger.info("... validated parameters.", cpp_state='end')
        except (AttributeError, TypeError):
            pass


def fetch_paramset_template(format, ret_defaults=False, params_sampling=None):
    """Fetch a parameter set template, either as a YAML file
    (with explanatory text) or as a dictionary.

    The returned template parameters are not validated as some mandatory
    values are unset.

    Parameters
    ----------
    format : {'text', 'dict'}
        Template format, either file contents ('text')
        or dictionary ('dict').
    ret_defaults : bool, optional
        If :keyword:`True` (default is :keyword:`False`), return a list
        of non-:keyword:`None` default parameters in the template.
        Only used when :param:`format` is 'dict'.

    Returns
    -------
    template : str or dict
        Parameter set template as text or a dictionary.
    defaults : dict
        Parameter set default entries.  Only returned when :param:`format`
        is 'dict' and :param:`ret_defaults` is :keyword:`True`.

    Raises
    ------
    ValueError
        If :param:`format` is neither 'yaml' nor 'dict'.

    """
    pkg_root_dir = Path(__file__).parent
    tml_filepath = pkg_root_dir/"resources"/"params_example.yml"

    text_template = Path(tml_filepath).read_text()
    dict_template = yaml.load(text_template, loader=yaml.Loader)

    if not ret_defaults:
        if format.lower() == 'text':
            return text_template
        if format.lower() == 'dict':
            return dict_template

    defaults = {}
    for key, val in dict_template.items():
        if isinstance(val, (tuple, list, set, dict)):
            continue
        if val is not None:
            defaults.update({key: val})

    if format.lower() == 'text':
        return text_template, defaults
    if format.lower() == 'dict':
        return dict_template, defaults

    raise ValueError("`format` must be either 'text' or 'dict'.")


def _modify_sampling_parameters(paramset, params_sampling=None,
                                params_default=None, ret_defaults=False):
    """Modify sampling parameters in a parameter set.

    Parameters
    ----------
    paramset : dict-like
        Parameter set to be updated.
    params_sampling : dict, optional
        Dictionary containing a subset of the following entries
        for sampling parameters:
            * 'boxalign': {'centre', 'pad'};
            * 'boxsize': [float, float, float];
            * 'ngrid': [int, int, int];
            * 'assignment': {'ngp', 'cic', 'tsc', 'pcs'};
            * 'interlace': bool;
        and one and only one of the following when 'boxalign' is 'pad':
            * 'boxpad': float;
            * 'gridpad': float.
        This will override corresponding parameters in the template.
        If :keyword:`None` (default), no modification happens and
        the original :param:`paramset` is returned.
    params_default : dict, optional
        If not :keyword:`None` (default), this is a sub-dictionary of
        :param:`paramset` and any of its entries unmodified by
        :param:`params_sampling` will be registered as default values.
    ret_defaults : bool, optional
        If :keyword:`True` (default is :keyword:`False`), return a list
        of non-:keyword:`None` default parameters which have not been
        modified.  Only used if :param:`params_default` is
        not :keyword:`None`.

    Returns
    -------
    paramset : dict-like
        Modified parameter set.
    params_default : dict, optional
        Only returned when :param:`ret_defaults` is :keyword:`True` and
        there are unmodified default parameters from
        :param:`params_default`.

    Raises
    ------
    ValueError
        When :param:`params_default` is not a sub-dictionary of
        :param:`paramset`.

    """
    if params_sampling is None:
        return paramset

    if params_default is not None:
        if not (params_default.items() <= paramset.items()):
            raise ValueError(
                "`params_default` should be a subdictionary of `paramset`."
            )

    if 'boxalign' in params_sampling.keys():
        paramset['alignment'] = params_sampling['boxalign']
        if params_default: params_default.pop('alignment')

    if paramset.get('alignment') == 'pad':
        if 'boxpad' in params_sampling.keys() \
                and 'gridpad' in params_sampling.keys():
            raise ValueError(
                "Cannot set both 'boxpad` and 'gridpad' in `params_sampling`."
            )
    if 'boxpad' in params_sampling.keys():
        paramset['padscale'] = 'box'
        paramset['padfactor'] = params_sampling['boxpad']
        if params_default:
            params_default.pop('padscale')
            params_default.pop('padfactor')
    if 'gridpad' in params_sampling.keys():
        paramset['padscale'] = 'grid'
        paramset['padfactor'] = params_sampling['gridpad']
        if params_default:
            params_default.pop('padscale')
            params_default.pop('padfactor')

    if 'boxsize' in params_sampling.keys():
        paramset['boxsize'] = dict(zip(
            ['x', 'y', 'z'], params_sampling['boxsize']
        ))
        if params_default: params_default.pop('boxsize')
    if 'ngrid' in params_sampling.keys():
        paramset['ngrid'] = dict(zip(
            ['x', 'y', 'z'], params_sampling['ngrid']
        ))
        if params_default: params_default.pop('ngrid')
    if 'assignment' in params_sampling.keys():
        paramset['assignment'] = params_sampling['assignment']
        if params_default: params_default.pop('assignment')
    if 'interlace' in params_sampling.keys():
        paramset['interlace'] = params_sampling['interlace']
        if params_default: params_default.pop('interlace')

    if ret_defaults:
        return paramset, params_default
    return paramset


def _modify_measurement_parameters(paramset, params_measure=None,
                                   params_default=None, ret_defaults=False):
    """Modify measurement parameters in a parameter set.

    Parameters
    ----------
    paramset : dict-like
        Parameter dictionary to be updated.
    params_measure : dict, optional
        Dictionary containing a subset of the following entries
        for measurement parameters:
            * 'ell1', 'ell2', 'ELL': int;
            * 'i_wa', 'j_wa': int;
            * 'idx_bin': int;
            * 'form': {'diag', 'full'}.
        This will override corresponding parameters in the template.
        If :keyword:`None` (default), no modification happens and
        the original :param:`paramset` is returned.
    params_default : dict, optional
        If not :keyword:`None` (default), this is a sub-dictionary of
        :param:`paramset` and any of its entries unmodified by
        :param:`params_measure` will be registered as default values.
    ret_defaults : bool, optional
        If :keyword:`True` (default is :keyword:`False`), return a list
        of non-:keyword:`None` default parameters which have not been
        modified.  Only used if :param:`params_default` is
        not :keyword:`None`.

    Returns
    -------
    paramset : dict-like
        Modified parameter dictionary.
    params_default : dict, optional
        Only returned when :param:`ret_defaults` is :keyword:`True` and
        there are unmodified default parameters from
        :param:`params_default`.

    Raises
    ------
    ValueError
        When :param:`params_default` is not a subdictionary of
        :param:`paramset`.

    """
    if params_measure is None:
        return paramset

    if params_default is not None:
        if not (params_default.items() <= paramset.items()):
            raise ValueError(
                "`params_default` should be a subdictionary of `paramset`."
            )

    for deg_label in ['ell1', 'ell2', 'ELL']:
        if deg_label in params_measure.keys():
            paramset['degrees'][deg_label] = params_measure[deg_label]

    for waorder_label in ['i_wa', 'j_wa']:
        if waorder_label in params_measure.keys():
            paramset['wa_orders'][waorder_label.rstrip('_wa')] = \
                params_measure[waorder_label]

    if 'form' in params_measure.keys():
        paramset['form'] = params_measure['form']
        if params_default: params_default.pop('form')

    if 'idx_bin' in params_measure.keys():
        paramset['idx_bin'] = params_measure['idx_bin']
        if params_default: params_default.pop('idx_bin')

    if ret_defaults:
        return paramset, params_default
    return paramset
