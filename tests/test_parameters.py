"""Test :mod:`~triumvirate.parameters`.

"""
from collections.abc import Collection, Sequence

import pytest
import yaml

from triumvirate.parameters import InvalidParameterError, NonParameterError
from triumvirate.parameters import ParameterSet
from triumvirate.parameters import fetch_paramset_template
from triumvirate.parameters import _TMPL_PARAM_DICT


@pytest.fixture(scope='module')
def template_parameters():
    return fetch_paramset_template('dict')


@pytest.fixture(scope='module')
def default_parameters():
    return {
        'alignment': 'centre',
        'padscale': 'box',
        'assignment': 'tsc',
        'interlace': False,
        'form': 'diag',
        'norm_convention': 'particle',
        'binning': 'lin',
        'verbose': 20,
    }


def test_fetch_paramset_template_text():

    return


@pytest.mark.parametrize(
    "case, ret_defaults",
    [
        ('text', False),
        ('dict', True),
    ]
)
def test_fetch_paramset_template_dict(case, ret_defaults, default_parameters):

    fetched = fetch_paramset_template(case, ret_defaults=ret_defaults)
    if ret_defaults:
        template, defaults = fetched
    else:
        template = fetched

    # Check for 'text' case and exit early.
    if case == 'text':
        return
    assert case == 'dict', f"Unrecognised argument value: {case}."

    # Check `template` agrees with an independent, private template
    # parameter dictionary.
    assert template == _TMPL_PARAM_DICT, \
        "Fetched template parameters do not match records."

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
    'param_filepath, param_dict',
    [
        ("triumvirate/resources/params_template.yml", None),
        (None, 'template_parameters'),
    ]
)
def test_ParameterSet___init__(param_filepath, param_dict, request):

    if param_dict is not None:
        param_dict = request.getfixturevalue(param_dict)

    # Check initialisation fails.
    with pytest.raises(InvalidParameterError):
        ParameterSet(
            param_filepath=param_filepath, param_dict=param_dict
        )

    # Modify failed initialisation and check it succeeds.
    if param_filepath is None and param_dict is not None:
        for ax_name in ['x', 'y', 'z']:
            param_dict['boxsize'][ax_name] = 1000.
            param_dict['ngrid'][ax_name] = 64.

        # NOTE: 'ell1' and 'ell2' are only mandatory for
        # three-point statistic algorithms.
        param_dict.update({
            'catalogue_type': 'sim',
            'statistic_type': 'bispec',
            'degrees': {'ell1': 0, 'ell2': 0, 'ELL': 0},
            'range': [0.005, 0.105],
            'num_bins': 10,
        })

        paramset = ParameterSet(
            param_filepath=param_filepath, param_dict=param_dict
        )
    elif param_filepath is not None and param_dict is None:

        with open(param_filepath, 'r') as tmpl_parameter_file:
            parameter_template = tmpl_parameter_file.read()

        def replace_parameter(parameter_template, str_original, str_new):
            return parameter_template.replace(str_original, str_new)

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
        ]

        for replacement in parameter_replacements:
            replace_parameter(parameter_template, *replacement)

        yaml.

        paramset = ParameterSet(
            param_filepath=param_filepath, param_dict=param_dict
        )
    else:
        raise ValueError("Invalid pytest.mark.parametrize parameters.")
