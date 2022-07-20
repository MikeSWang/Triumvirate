"""Generate parameter files from a template.

For C++, the parameter files are assumed to use '%' as the comment
delimiter, and parameter assignment is in the form
'<param_name> = <param_value>' each on a single line.

For Python, YAML files are used.

"""
import argparse
import warnings
from copy import deepcopy
from functools import reduce
from operator import getitem
from os.path import splitext
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
        prog='gen_param_file',
        description="Generate parameter files by modifying a template."
    )

    parser.add_argument(
        'template', type=str,
        help="(path to the) template parameter file"
    )

    parser.add_argument(
        '-i', '--cpp', action='store_true',
        help="if not in conflict, generate INI parameter file for C++"
    )

    parser.add_argument(
        '-y', '--python', action='store_true',
        help="if not in conflict, generate YAML parameter file for Python"
    )

    parser.add_argument(
        '-o', '--output', type=str,
        help="(path to the) output parameter file"
    )

    parser.add_argument(
        '-p', '--modify', action='append', nargs=2,
        metavar=('param_name', 'param_value'),
        help="pairwise modified parameter names and values"
    )

    config = parser.parse_args()

    if (config.cpp and config.python) or (not config.cpp and not config.python):
        raise ValueError(
            "One and only one of '-i' (or '--cpp') and '-y' (or '--python') "
            "flags must be set."
        )

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
        raise ParentKeyError(str(err))

    try:
        base_dict[key_child] = type(base_dict[key_child])(value)
    except KeyError as err:
        base_dict[key_child] = value
        raise ChildKeyError(str(err))

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
        pname_line = line.split('%')[0].split('=')[0].strip()
        if not pname_line: continue  # skip comment lines
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
                .format(pname, str(err).strip('"'))
            )
        except ChildKeyError as err:
            warnings.warn(
                "The following parameter is added: {} (missing key {}). "
                "Check the template for possible typos."
                .format(pname, str(err).strip('"')),
                category=RuntimeWarning
            )

    return params


if __name__ == '__main__':

    # Configure.
    config = configure()

    if config.output:
        # Use specified output path.
        output_path = config.output
    else:
        # Use the template path with '_new' tag.
        template_file_base, template_file_ext = splitext(config.template)
        output_path = "".join([template_file_base, '_new', template_file_ext])
    # Ensure the output directory exists.
    Path(output_path).parent.mkdir(exist_ok=True, parents=True)

    # Load and generate.
    if config.cpp:
        with open(config.template, 'r') as template_file:
            pars_template = template_file.readlines()
            pars = gen_cpp_params(pars_template, config.modify)
        with open(output_path, 'w') as output_file:
            output_file.write(pars)

    if config.python:
        with open(config.template, 'r') as template_file:
            pars_template = yaml.load(template_file, Loader=yaml.Loader)
            pars = gen_python_params(pars_template, config.modify)
        with open(output_path, 'w') as output_file:
            yaml.dump(
                pars, output_file, sort_keys=False, default_flow_style=False
            )
