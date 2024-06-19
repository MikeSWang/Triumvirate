"""Exchange coordinate axes of two-dimensional three-point clustering
statistics.

Examples
--------

.. code-block:: console

    $ python exchange_coord_axes.py <input_path> --output-dir <output-dir>

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


def transpose_filename(filepath):
    """Transpose the filename of a two-dimensional three-point clustering
    measurement file.

    Parameters
    ----------
    filepath : str
        The file path of the input file.

    Returns
    -------
    filename_ex : :class:`pathlib.Path`
        The transposed filename.

    """
    filepath = Path(filepath)
    filename = filepath.stem
    extension = filepath.suffix

    filename_fields = filename.split('_')
    stat_field = filename_fields[0]
    nonstat_fields = filename_fields[1:]

    stat_type = re.match(r'\D+', stat_field).group()
    stat_indices = re.findall(r'\d+', stat_field)

    if len(stat_indices) == 1:
        stat_indices = list(stat_indices.pop())
    if len(stat_indices) != 3:
        raise ValueError(
            "The filename must start with the statistic type and multipole "
            "indices of length 3 before any underscores, "
            "e.g. 'bk000_' or 'bk0,0,0_'."
        )
    stat_indices_ex = stat_indices[1:] + [stat_indices[0]]

    if ',' in stat_field:
        index_sep = ','
    else:
        index_sep = ''
    stat_field_ex = stat_type + index_sep.join(stat_indices_ex)

    filename_stem = '_'.join([stat_field_ex] + nonstat_fields)
    filename_ex = Path(filename_stem).with_suffix(extension)

    return filename_ex


def transpose_datatab(datatab):
    """Transpose the data table of a two-dimensional three-point
    clustering measurement file.

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

    coord_cols_a = slice(coord_colidx_a, coord_colidx_a + coord_colspan - 1)
    coord_cols_b = slice(coord_colidx_b, coord_colidx_b + coord_colspan - 1)

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
    output_filename = transpose_filename(cfg.input_path)
    if not cfg.output_dir:
        output_dir = Path(cfg.input_path).parent
    else:
        output_dir = Path(cfg.output_dir)
    if cfg.output_tag:
        if output_filename.suffix not in ('.dat', '.txt', '.csv', '.tsv'):
            output_path = output_dir/(output_filename + cfg.output_tag)
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
    header = ''.join(headerlines)

    # Transpose the data table.
    datatab = np.loadtxt(cfg.input_path)
    datatab_ex = transpose_datatab(datatab)

    # Save the transposed data table.
    np.savetxt(output_path, datatab_ex, fmt='%.9e', header=header)
    print(f"Transpose-exchanged measurement file saved to: {output_path}")


if __name__ == '__main__':
    main()
