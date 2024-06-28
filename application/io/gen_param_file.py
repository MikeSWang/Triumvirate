r"""Generate parameter files from a template.

For C++, the INI parameter files use '#' as the comment delimiter, and
each parameter is specified by a single-line assignment in the form
``<param_name> = <param_value>``.

For Python, YAML parameter files are used.

Examples
--------
To generate a new parameter file ``<param-file>`` by modifying the
template parameter file ``<template-path>`` with the modified
parameter--value pairs ``<param-name1>`` and ``<param-value1>``, and
``<param-name2>`` and ``<param-value2>``, run:

.. code-block:: console

    $ python gen_param_file.py <template-path> -o <param-file> \
    >   -p <param-name1> <param-value1> -p <param-name2> <param-value2>

"""
import argparse
import warnings
from ast import literal_eval
from copy import deepcopy
from functools import reduce
from operator import getitem
from pathlib import Path

import yaml


class ParentKeyError(KeyError):
    pass


class ChildKeyError(KeyError):
    pass


def configure():
    """Generate modified-parameter configuration.

    Returns
    -------
    config : :class:`argparse.Namespace`
        Modified-parameter configuration.

    """
    parser = argparse.ArgumentParser(
        description="Generate parameter files by modifying a template."
    )

    parser.add_argument(
        'template_path', type=str,
        help="template parameter file path"
    )

    parser.add_argument(
        '--fmt', type=str, choices=['ini', 'yml', 'yaml'], default=None,
        help="parameter file format (default inferred from extension)"
    )

    parser.add_argument(
        '-o', '--output', type=str, default=None,
        help="output parameter file path"
    )

    parser.add_argument(
        '-p', '--modify', action='append', nargs=2,
        metavar=('PARAM_NAME', 'PARAM_VALUE'),
        help="pairwise modified parameter name(s) and value(s)"
    )

    config = parser.parse_args()

    if config.fmt is None:
        config.file_ext = Path(config.template_path).suffix
        config.fmt = config.file_ext.lstrip('.')
    else:
        config.file_ext = f".{config.fmt}"

    config.output = config.output or f"{config.file_stem}_new{config.file_ext}"

    return config


def setitem_nested_dict(dictionary, keychain, value):
    """Set the value of an item in a nested dictionary.

    Parameters
    ----------
    dictionary : dict
        Input dictionary.
    keychain : list
        List of keys in increasing depth.
    value
        Value to set to.

    Returns
    -------
    dictionary: dict
        Output dictionary.

    """
    key_parents, key_child = keychain[:-1], keychain[-1]

    try:
        base_dict = reduce(getitem, key_parents, dictionary)
    except KeyError as err:
        raise ParentKeyError(getattr(err, 'message', str(err))) from err

    try:
        if isinstance(base_dict[key_child], list):
            base_dict[key_child] = literal_eval(value)
        else:
            base_dict[key_child] = type(base_dict[key_child])(value)
    except KeyError as err:
        base_dict[key_child] = value
        raise ChildKeyError(getattr(err, 'message', str(err))) from err

    return dictionary


def gen_cpp_params(params_template, params_to_modify):
    """Generate C++ parameters as a string block.

    Parameters
    ----------
    params_template : list of str
        Parameter template as a list of strings each being a line
        in the template file.
    params_to_modify : list of list of str
        Name and value pairs of modified parameters.

    Returns
    -------
    params : str
        Modified parameters as a string block.

    """
    params = params_template.copy()
    for lineno, line in enumerate(params_template):
        pname_line = line.split('#')[0].split('=')[0].strip()
        if not pname_line:
            continue  # skip comment lines
        for pidx, (pname, pval) in enumerate(params_to_modify):
            if pname_line == pname:
                params[lineno] = f"{pname} = {pval}\n"
                del params_to_modify[pidx]  # remove once replaced

    if params_to_modify:  # check for unmodified parameters
        warnings.warn(
            "The following parameters are added: {}. "
            "Check for possible typos as they are missing in the template."
            .format(params_to_modify),
            category=RuntimeWarning
        )
        params.append('\n')
        params.extend([' = '.join(param) + '\n' for param in params_to_modify])

    return "".join(params)


def gen_python_params(params_template, params_to_modify):
    """Generate Python parameters as a nested dictionary.

    Parameters
    ----------
    params_template : dict
        Parameter template as a nested dictionary.
    params_to_modify : list of list of str
        Name and value pairs of modified parameters.

    Returns
    -------
    params : dict
        Modified parameters as a nested dictionary.

    """
    params = deepcopy(params_template)
    for pname, pval in params_to_modify:
        pkeychain = pname.split('.')
        try:
            setitem_nested_dict(params, pkeychain, pval)
        except ParentKeyError as err:
            raise RuntimeError(
                "Check the template for possible typos as "
                "the following parameter could not be modified: {} "
                "(missing key {})."
                .format(pname, getattr(err, 'message', str(err)).strip('"'))
            ) from err
        except ChildKeyError as err:
            warnings.warn(
                "The following parameter is added: {} (missing key {}). "
                "Check the template for possible typos."
                .format(pname, getattr(err, 'message', str(err)).strip('"')),
                category=RuntimeWarning
            )

    return params


def read_param_file(input_path, fmt):
    """Read parameters from a file.

    Parameters
    ----------
    input_path : str
        Input file path.
    fmt : {'ini', 'yml', 'yaml'}
        Input file format.

    Returns
    -------
    params : list[str] or dict
        Parameters read from the file.

    """
    if fmt == 'ini':
        with open(input_path, 'r') as input_file:
            params = input_file.readlines()
    elif fmt == 'yml' or fmt == 'yaml':
        with open(input_path, 'r') as input_file:
            params = yaml.load(input_file, Loader=yaml.SafeLoader)
    else:
        raise ValueError(f"Unsupported file format: {fmt}.")

    print("Read parameter template {}.".format(input_path))

    return params


def write_param_file(params, output_path, fmt):
    """Write parameters to a file.

    Parameters
    ----------
    params : str or dict
        Parameters to write.
    output_path : str
        Output file path.
    fmt : str
        Output file format.

    """
    # Ensure the output directory exists.
    Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    if fmt == 'ini':
        with open(output_path, 'w') as output_file:
            output_file.write(params)
    elif fmt == 'yml' or fmt == 'yaml':
        with open(output_path, 'w') as output_file:
            yaml.dump(
                params, output_file, sort_keys=False, default_flow_style=False
            )
    else:
        raise ValueError(f"Unsupported file format: {fmt}.")

    print("Generated parameter file {}.".format(output_path))


if __name__ == '__main__':

    # Configure.
    config = configure()

    # Load, generate and write anew.
    pars_template = read_param_file(config.template_path, config.fmt)

    if config.fmt == 'ini':
        pars_modified = gen_cpp_params(pars_template, config.modify)
    if config.fmt == 'yml' or config.fmt == 'yaml':
        pars_modified = gen_python_params(pars_template, config.modify)

    write_param_file(pars_modified, config.output, config.fmt)
