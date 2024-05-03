"""Test :mod:`~triumvirate.parameters`.

"""
import os.path
import warnings
from collections.abc import Collection, Sequence
from copy import deepcopy
from pprint import pformat

import pytest
import yaml

from triumvirate.parameters import (
    InvalidParameterError,
    ParameterSet,
    fetch_paramset_template,
    _TMPL_PARAM_DICT,
)


def make_template_parameters_valid(tmpl_params):
    """Construct valid parameters from template parameters.

    Parameters
    ----------
    tmpl_params : str or dict
        Template parameters.

    Returns
    -------
    str or dict
        Valid parameters constructed from the template parameters.

    Raises
    ------
    TypeError
        If `template_parameters` type is not :class:`str`
        or :class:`dict`.

    """
    if isinstance(tmpl_params, str):
        def replace_parameter(_text, str_original, str_new):
            return _text.replace(str_original, str_new)

        param_text = tmpl_params
        parameter_replacements = [
            (
                "boxsize:\n  x:\n  y:\n  z:\n",
                "boxsize:\n  x: 1000.\n  y: 1000.\n  z: 1000.\n"
            ),
            (
                "ngrid:\n  x:\n  y:\n  z:\n",
                "ngrid:\n  x: 64\n  y: 64\n  z: 64\n"
            ),
            (
                "catalogue_type:",
                "catalogue_type: sim"
            ),
            (
                "statistic_type:",
                "statistic_type: bispec"
            ),
            (
                "ell1:",
                "ell1: 0"
            ),
            (
                "ell2:",
                "ell2: 0"
            ),
            (
                "ELL:",
                "ELL: 0"
            ),
            (
                "range: [~, ~]",
                "range: [0.005, 0.105]"
            ),
            (
                "num_bins:",
                "num_bins: 10"
            ),
            # Not necessary, but suppresses logging.
            (
                "verbose: 20",
                "verbose: 60"
            ),
        ]

        for _replacement in parameter_replacements:
            param_text = replace_parameter(param_text, *_replacement)

        return param_text

    if isinstance(tmpl_params, dict):
        param_dict = deepcopy(tmpl_params)

        for ax_name in ['x', 'y', 'z']:
            param_dict['boxsize'][ax_name] = 1000.
            param_dict['ngrid'][ax_name] = 64.

        param_dict.update({
            'catalogue_type': 'sim',
            'statistic_type': 'bispec',
            'degrees': {'ell1': 0, 'ell2': 0, 'ELL': 0},
            'range': [0.005, 0.105],
            'num_bins': 10,
            'verbose': 60,  # not necessary but suppresses logging
        })

        return param_dict

    raise TypeError("`template_parameters` is not str or dict.")


# Returns the subset of default parameters from the template.
@pytest.fixture
def default_parameters():
    return {
        'alignment': 'centre',
        'padscale': 'box',
        'assignment': 'tsc',
        'interlace': False,
        'form': 'diag',
        'norm_convention': 'particle',
        'binning': 'lin',
        'fftw_scheme': 'measure',
        'use_fftw_wisdom': False,
        'save_binned_vectors': False,
        'verbose': 20,
    }


@pytest.fixture
def template_parameters_text():
    return fetch_paramset_template('text')


@pytest.fixture
def template_parameters_dict():
    return fetch_paramset_template('dict')


@pytest.mark.parametrize(
    "source, ret_defaults",
    [
        ('text', False),
        ('dict', True),
    ]
)
def test_fetch_paramset_template(source, ret_defaults, default_parameters):

    fetched = fetch_paramset_template(source, ret_defaults=ret_defaults)
    if ret_defaults:
        template, defaults = fetched
    else:
        template = fetched

    # Check for 'text' source and exit early; else assume 'dict' source.
    if source == 'text':
        return
    if source != 'dict':
        raise ValueError("Invalid `pytest.mark.parametrize` parameters.")

    # Check `template` agrees with an independent, private template
    # parameter dictionary.
    assert template == _TMPL_PARAM_DICT, \
        "Fetched template parameters do not match internal records."

    # Check `defaults`` are pre-set default parameters.
    assert defaults == default_parameters, (
        "Fetched template default parameters do not match "
        "pre-set default parameters."
    )

    # Check `template` minus `defaults` are unset.
    non_defaults_keys = set(template.keys()) - set(defaults.keys())
    for param_name in non_defaults_keys:
        param_val = template[param_name]
        if isinstance(param_val, Sequence):
            assert all(
                param_val_element is None for param_val_element in param_val
            ), (
                "Non-NoneType element found in "
                "sequence-type non-default parameter value "
                "in parameter template."
            )
        elif isinstance(param_val, Collection):
            assert all(
                param_val_val is None for param_val_val in param_val.values()
            ), (
                "Non-NoneType value found in "
                "collection-type non-default parameter value "
                "in parameter template."
            )
        else:
            assert param_val is None, (
                "Non-NoneType non-default parameter value "
                "found in parameter template."
            )


@pytest.mark.parametrize(
    "param_filepath, param_dict",
    [
        ("tmpl_params.yml", None),
        (None, 'template_parameters_dict'),
        (
            "tmpl_params.yml",
            'template_parameters_dict',
        ),
    ]
)
def test_ParameterSet___cinit__(param_filepath, param_dict, test_param_dir,
                                request, tmp_path):

    # Patch paths for wheel testing.
    if param_filepath is not None:
        param_filepath = os.path.join(test_param_dir, param_filepath)

    # Convert parametrized parameters to fixtures.
    if param_dict is not None:
        param_dict = request.getfixturevalue(param_dict)

    # Check initialisation fails (and exit early if needed).
    if param_filepath is not None and param_dict is not None:
        with pytest.raises(ValueError):
            ParameterSet(
                param_filepath=param_filepath, param_dict=param_dict
            )
        return

    with pytest.raises(InvalidParameterError):
        ParameterSet(
            param_filepath=param_filepath, param_dict=param_dict
        )

    # Modify failed initialisation and check it succeeds.
    if param_filepath is None and param_dict is not None:
        param_dict = make_template_parameters_valid(param_dict)
        ParameterSet(param_dict=param_dict)
    elif param_filepath is not None and param_dict is None:
        with open(param_filepath, 'r') as tmpl_parameter_file:
            param_text = make_template_parameters_valid(
                tmpl_parameter_file.read()
            )

        param_filepath = tmp_path/"params.yml"
        with open(param_filepath, 'w') as param_file:
            param_file.write(param_text)

        ParameterSet(param_filepath=param_filepath)
    else:
        raise ValueError("Invalid `pytest.mark.parametrize` parameters.")


@pytest.mark.parametrize(
    "template_parameter_source",
    [
        'template_parameters_text',
        'template_parameters_dict',
    ]
)
def test_ParameterSet___str__(template_parameter_source, request, tmp_path):

    valid_parameter_source = make_template_parameters_valid(
        request.getfixturevalue(template_parameter_source)
    )

    if isinstance(valid_parameter_source, str):
        param_filename = "params.yml"
        with open(tmp_path/param_filename, 'w') as param_file:
            param_file.write(valid_parameter_source)

        paramset = ParameterSet(param_filepath=tmp_path/param_filename)
        assert param_filename in str(paramset)
        assert 'original=True' in str(paramset)
        return

    if isinstance(valid_parameter_source, dict):
        paramset = ParameterSet(param_dict=valid_parameter_source)
        assert 'dict' in str(paramset)
        assert 'original=True' in str(paramset)
        return

    raise ValueError(
        "Invalid `valid_parameter_source` fixture parametrisation."
    )


# Use `default_parameters` fixture to test the valid parameter set.
def test_ParameterSet___getitem__(valid_paramset, default_parameters):

    # Modify 'verbose' to suppress logging.
    default_parameters.update(verbose=60)

    for key, val in default_parameters.items():
        # This is three tests in one: TypeError, KeyError and value comparison.
        if valid_paramset[key] != val:
            warnings.warn(
                "`valid_paramset` fixture does not match default parameters "
                "for key '{}': {} != {}."
                .format(key, valid_paramset[key], val),
                category=RuntimeWarning
            )


@pytest.mark.parametrize(
    "param_name, param_value",
    [
        ('degrees', {'ell1': 2, 'ell2': 0, 'ELL': 2}),
        ('interlace', True),
    ]
)
def test_ParameterSet___setitem__(param_name, param_value, valid_paramset):

    # Parameter setting should not work when enabling interlacing
    # for three-point statistics.
    if valid_paramset.npoint == '3pt' and param_name == 'interlace':
        valid_paramset[param_name] = param_value
        with pytest.raises(AssertionError):
            assert valid_paramset[param_name] == param_value, \
                "Parameter setting overriden."
        valid_paramset['statistic_type'] = 'powspec'

    # Otherwise, parameter setting should work.
    valid_paramset[param_name] = param_value
    assert valid_paramset[param_name] == param_value, \
        "Parameter set value setting failed."


# Use `default_parameters` fixture to test the valid parameter set.
def test_ParameterSet__getattr__(valid_paramset, default_parameters):

    # Modify 'verbose' to suppress logging.
    default_parameters.update(verbose=60)

    for attr, val in default_parameters.items():
        # This is three tests in one: TypeError, KeyError and value comparison.
        if getattr(valid_paramset, attr) != val:
            warnings.warn(
                "`valid_paramset` fixture does not match default parameters "
                "for attribute '{}': {} != {}."
                .format(attr, getattr(valid_paramset, attr), val),
                category=RuntimeWarning
            )


@pytest.mark.parametrize(
    "param_name, param_value",
    [
        ('degrees', {'ell1': 2, 'ell2': 0, 'ELL': 2}),
        ('interlace', True),
    ]
)
def test_ParameterSet___setattr__(param_name, param_value, valid_paramset):

    # Parameter setting should not work when enabling interlacing
    # for three-point statistics.
    if valid_paramset.npoint == '3pt' and param_name == 'interlace':
        setattr(valid_paramset, param_name, param_value)
        with pytest.raises(AssertionError):
            assert getattr(valid_paramset, param_name) == param_value, \
                "Parameter setting overriden."
        valid_paramset.statistic_type = 'powspec'

    # Otherwise, parameter setting should work.
    setattr(valid_paramset, param_name, param_value)
    assert getattr(valid_paramset, param_name) == param_value, \
        "Parameter set value setting failed."


def test_ParameterSet_names(valid_paramset):
    # There are three derived parameters after validation:
    # 'assignment_order', 'npoint', 'space'.
    assert set(valid_paramset.names()) - set(_TMPL_PARAM_DICT.keys()) == {
        'assignment_order', 'npoint', 'space',
        'fftw_planner_flag', 'fftw_wisdom_file_f', 'fftw_wisdom_file_b'
    }, "Parameter set names do not match internal records."


def test_ParameterSet_items(valid_paramset):
    # There are two derived parameters after validation, 'npoint' and 'space'.
    assert valid_paramset.items(), "Parameter set items not returned."


@pytest.mark.parametrize(
    "param_name, param_value",
    [
        ('npoint', '3pt'),
        ('free_lunch', None),
    ]
)
def test_ParameterSet_get(param_name, param_value, valid_paramset):
    assert valid_paramset.get(param_name) == param_value, \
        "Parameter set parameter not returned correctly."


@pytest.mark.parametrize(
    "update_args, update_kwargs",
    [
        ([dict(verbose=60)], {'catalogue_type': 'random'}),
    ]
)
def test_ParameterSet_update(update_args, update_kwargs, valid_paramset):

    valid_paramset.update(*update_args, **update_kwargs)

    for update_arg in update_args:
        for _name, _value in update_arg.items():
            assert valid_paramset.get(_name) == _value, \
                "Parameter set update by positional argument failed."
    for _name, _value in update_kwargs.items():
        assert valid_paramset.get(_name) == _value, \
            "Parameter set update by keyword argument failed."


def test_ParameterSet_print(valid_paramset, capsys):
    valid_paramset.print()
    assert pformat(dict(valid_paramset.items()), sort_dicts=False) \
        in capsys.readouterr().out, "Parameter set misprinted."


def test_ParameterSet_save(valid_paramset, tmp_path):

    param_filepath = tmp_path/"params_saved.yml"

    valid_paramset.save(filepath=param_filepath)

    with open(param_filepath, 'r') as param_file:
        valid_param_dict = yaml.load(param_file, Loader=yaml.Loader)

    assert dict(valid_paramset.items()) == valid_param_dict, \
        "Loaded parameter set does not match the saved one."
