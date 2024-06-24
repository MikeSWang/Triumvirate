r"""Exchange coordinate axes of 2-d three-point clustering statistics.

Examples
--------
To exchange the coordinate axes of a 2-d three-point clustering
measurement file with path ``<input_path>`` and save the transposed
measurement file to the output directory ``<output-dir>`` with the tag
``<tag>``, run:

.. code-block:: console

    $ python exchange_coord_axes.py <input_path> \
    >   --output-dir <output-dir> --output-tag <tag>

"""
import argparse
import re
from pathlib import Path

import numpy as np


def configure():
    """Configure parameters from command-line arguments.

    Returns
    -------
    config : argparse.Namespace
        Command-line arguments.

    """
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(
        'input_path', type=str,
        help="input measurement file path"
    )

    parser.add_argument(
        '--output-dir', type=str,
        help="output directory (defaults to the input directory)"
    )

    parser.add_argument(
        '--output-tag', type=str,
        help="output measurement file tag"
    )

    config = parser.parse_args()

    return config


def transpose_stat_str(pathlike, ret_multipole=False):
    """Transpose the multipole indices in a 2-d three-point clustering
    measurement path-like object.

    Parameters
    ----------
    pathlike : path-like
        Path-like containing the multipole indices in the first of
        multiple fields delimited by underscores.
    ret_multipole : bool, optional
        If `True` (default is `False`), return the original
        multipole indices.

    Returns
    -------
    filenamelike_ex : :class:`pathlib.Path`
        The transposed filename.
    stat_type : str
        The statistic type, e.g. 'bk'.
    multipole_ex : str
        The transposed multipole indices.
    multipole : str, optional
        The original multipole indices.

    """
    pathlike = Path(pathlike)

    fields = pathlike.stem.split('_')
    stat_field, nonstat_fields = fields[0], fields[1:]

    stat_type = re.match(r'\D+', stat_field).group()
    stat_indices = re.findall(r'\d+', stat_field)
    index_sep = ',' if ',' in stat_field else ''

    if len(stat_indices) == 1:
        stat_indices = list(stat_indices.pop())
    if len(stat_indices) != 3:
        raise ValueError(
            "The filename must start with the statistic type and multipole "
            "indices of length 3 before any underscores, "
            "e.g. 'bk000_' or 'bk0,0,0_'."
        )
    if stat_indices[0] == stat_indices[1]:
        raise ValueError("The first two multipole indices are the same.")
    stat_indices_ex = [stat_indices[1], stat_indices[0]] + stat_indices[2:]

    multipole = index_sep.join(stat_indices)
    multipole_ex = index_sep.join(stat_indices_ex)

    stat_field_ex = stat_type + multipole_ex

    pathlike_stem_ex = '_'.join([stat_field_ex] + nonstat_fields)
    filenamelike_ex = Path(pathlike_stem_ex).with_suffix(pathlike.suffix)

    if ret_multipole:
        return filenamelike_ex, stat_type, multipole_ex, multipole
    return filenamelike_ex, stat_type, multipole_ex


def transpose_datatab(datatab):
    """Transpose the data table of a 2-d three-point clustering
    measurement file.

    Parameters
    ----------
    datatab : :class:`numpy.ndarray`
        The input measurement data table.

    Returns
    -------
    datatab_ex : :class:`numpy.ndarray`
        The transposed data table.

    """
    datatab_ex = datatab.copy()

    # Exchange coordinate columns.
    coord_colspan = 3
    coord_colidx_a = 0
    coord_colidx_b = 3

    coord_cols_a = slice(coord_colidx_a, coord_colidx_a + coord_colspan)
    coord_cols_b = slice(coord_colidx_b, coord_colidx_b + coord_colspan)

    datatab_ex[:, coord_cols_a] = datatab[:, coord_cols_b]
    datatab_ex[:, coord_cols_b] = datatab[:, coord_cols_a]

    # Resort the data table.
    coords_ex_primary = datatab_ex[:, coord_colidx_a]
    coords_ex_secondary = datatab_ex[:, coord_colidx_b]

    datatab_ex = datatab_ex[np.lexsort(
        (coords_ex_secondary, coords_ex_primary)
    )]

    return datatab_ex


def main():
    """Main function.

    """
    cfg = configure()

    # Transpose the measurement filename.
    try:
        output_filename, stat_type, multipole_ex, multipole = \
            transpose_stat_str(cfg.input_path, ret_multipole=True)
    except ValueError as err:
        if "The first two multipole indices the same." in str(err):
            print(
                "No transposition is performed as the first two multipole "
                "indices are the same."
            )
            return
        # STYLE: Un-Pythonic 'else' clause after 'return' statement
        # but retained for clarity.
        else:
            raise err
    if not cfg.output_dir:
        output_dir = Path(cfg.input_path).parent
    else:
        output_dir = Path(cfg.output_dir)
    if cfg.output_tag:
        if output_filename.suffix not in ('.dat', '.txt', '.csv', '.tsv'):
            output_path = output_dir/(str(output_filename) + cfg.output_tag)
        else:
            output_path = output_dir/(
                output_filename.stem + cfg.output_tag + output_filename.suffix
            )
    else:
        output_path = output_dir/output_filename

    # Transpose the measurement file header.
    headerlines = []
    with open(cfg.input_path, 'r') as measurement_file:
        line = measurement_file.readline()
        while line.startswith('#'):
            headerlines.append(line)
            line = measurement_file.readline()
    headerlines[-1] = headerlines[-1].replace(
        stat_type + multipole, stat_type + multipole_ex
    )
    header = ''.join(headerlines).rstrip()

    # Transpose the data table.
    datatab = np.loadtxt(cfg.input_path)
    datatab_ex = transpose_datatab(datatab)

    # Save the transposed data table.
    fmt = '\t'.join(
        (['%.9e'] * 2 + ['%10d']) * 2 +
        ['% .9e'] * (datatab.shape[-1] - 6)
    )
    np.savetxt(output_path, datatab_ex, fmt=fmt, header=header, comments='')
    print(f"Transpose-exchanged measurement file saved to: {output_path}")


if __name__ == '__main__':
    main()
