"""Combine two-point and three-point clustering multipole data vectors.

.. caution::

    Data vectors are assumed to be bona fide 1-d, i.e.
    three-point clustering statistics are assumed to be in diagonal form.


Examples
--------
To combine the two-point clustering multipole data vectors ``pk0``,
``pk2`` and ``pk4`` over realisations 1 to 256 from the input directory
``<input-dir>`` with the regular expression tag
(e.g.) ``'_DR12SGC-COMPSAM_V6C_(.*)_'`` and save the combined data vector
to the output directory ``<output-dir>``, run:

.. code-block:: console

    $ python comb_data_vectors.py
    > --tag _DR12SGC-COMPSAM_V6C_(.*)_0.50z0.75
    > --range 1 256
    > --components pk0 pk2 pk4
    > --input-dir <input-dir>
    > --output-dir <output-dir>

"""
import argparse
import os
import re
import warnings

import numpy as np


def configure():
    """Configure data vector combination.

    Returns
    -------
    config : :class:`argparse.Namespace`
        Data vector combination configuration.

    """
    parser = argparse.ArgumentParser(
        description="Combine clustering multipole data vectors."
    )

    parser.add_argument(
        '--tag', type=str, default="(.*)",
        help="data vector file common tag including a single RegEx pattern "
             "matching a serial number, e.g. '_BOSS_PatchyMocks_(.*)' "
             "where '(.*)' should match a serial number"
    )
    parser.add_argument(
        '--range', type=int, nargs='*', default=None,
        help="range of the serial number (if any) in the common tag"
    )

    parser.add_argument(
        '-s', '--components', type=str, nargs='+', required=True,
        help="ordered list of statistic components "
             "with the following prefixes: {pk, xi, bk, zeta} "
             "(e.g. '--statistics pk0 pk2 pk4 bk000') "
    )

    parser.add_argument(
        '--input-dir', type=str, default="./",
        help="input directory of data vectors to be combined"
    )
    parser.add_argument(
        '--output-dir', type=str, default="./",
        help="output directory of the combined data vector"
    )

    config = parser.parse_args()

    return config


def split_component(component):
    """Split statistic component into the statistic tag and the
    multipole degree(s).

    Parameters
    ----------
    component : str
        Statistic component.  Must start with 'pk', 'xi', 'bk' or 'zeta'
        as a valid statistic tag and end with one or three single digits
        as (a) valid multipole degree(s).

    Returns
    -------
    stat : str
        Statistic tag.
    pole : str
        Multipole degree(s).

    Raises
    ------
    ValueError
        If the component is an unrecognised statistic or an invalid
        multipole.

    """
    stat, pole = re.findall(r'\d+|\D+', component)

    if stat.startswith(('pk', 'xi')):
        if len(pole) != 1:
            raise ValueError(
                "Multipole degree of two-point statistics "
                "must be a single digit."
            )
    elif stat.startswith(('bk', 'zeta')):
        if len(pole) != 3:
            raise ValueError(
                "Multipole degree of three-point statistics "
                "must be three single digits."
            )
    else:
        raise ValueError(f"Unrecognised statistic: {stat}.")

    return stat, pole


def read_measurements(stat, directory, fnpattern, sn_range=None):
    """Read measurement files matching the specification.

    Parameters
    ----------
    stat : {'pk', 'xi', bk', 'zeta'}
        Statictic.
    directory : path-like
        File directory.
    fnpattern : str
        Filename pattern.
    sn_range : sequence of int, optional
        Filename serial number range (default is `None`).

    Returns
    -------
    tuple of :class:`numpy.ndarray` of length 5
        Bin centres, bin effective coordinates, bin contribution counts,
        stacked raw measurements and stacked shot noise measurements.

    Raises
    ------
    ValueError
        If `sn_range` is not `None` or a sequence of length 2.

    """
    USECOLS = {
        'pk': list(range(7)),
        'xi': list(range(5)),
        'bk': list(range(2, 9)),
        'zeta': list(range(2, 9)),
    }

    # Check arguments.
    if sn_range is not None and len(sn_range) != 2:
        raise ValueError("`sn_range` must be a sequence of length 2.")

    # Match files.
    paths_found = []
    for _, _, filenames in os.walk(directory):
        for filename in sorted(filenames):
            if re.match(fnpattern, filename):
                paths_found.append(os.path.join(directory, filename))

    if len(paths_found) == 0:
        raise IOError(
            f"No filename of pattern {fnpattern} found in {directory}."
        )

    paths = []
    for path in paths_found:
        serial_num = int(re.search(fnpattern, path).group(1))
        if sn_range is None or (sn_range[0] <= serial_num <= sn_range[1]):
            paths.append(path)

    # Collate data.
    bins_list, coords_list, counts_list = [], [], []
    measurement_raw_list, measurement_shotnoise_list = [], []

    for path in paths:
        data_ = np.loadtxt(path, usecols=USECOLS[stat], unpack=True)

        bins_list.append(data_[0])
        coords_list.append(data_[1])
        counts_list.append(data_[2])

        measurement_raw_list.append(data_[3] + 1j * data_[4])
        try:
            measurement_shotnoise_list.append(data_[5] + 1j * data_[6])
        except IndexError:
            measurement_shotnoise_list.append(np.zeros(data_.shape[-1]))

    # Check data agreement.
    for name, array in zip(
        ['centres', 'effective coordinates', 'contributing counts'],
        [bins_list, coords_list, counts_list],
    ):
        if not np.allclose(np.diff(array, axis=0), 0.):
            warnings.warn(
                "Data bin {} do not match for measurements from '{}'."
                .format(name, fnpattern)
            )

    # Organise data arrays.
    bins, coords, counts = bins_list[-1], coords_list[-1], counts_list[-1]
    measurements_raw = np.vstack(measurement_raw_list)
    measurements_shotnoise = np.vstack(measurement_shotnoise_list)

    return bins, coords, counts, measurements_raw, measurements_shotnoise


def insert_component(component, directory,
                     components_bins, components_coords, components_counts,
                     components_vectors, components_vectors_shotnoise,
                     tag='', sn_range=None):
    """Insert data vector component to a list to be combined.

    Parameters
    ----------
    component : str
        Statistic component.
    directory : path-like
        Source directory of the component.
    components_bins : list of :class:`numpy.ndarray`
        Data bin centres for all components to append to.
    components_coords : list of :class:`numpy.ndarray`
        Data bin effective coordinates for all components to append to.
    components_counts : list of :class:`numpy.ndarray`
        Data bin contribution counts for all components to append to.
    components_vectors : list of :class:`numpy.ndarray`
        Data vectors of raw measurements for all components to append to.
    components_vectors_shotnoise : list of :class:`numpy.ndarray`
        Data vectors of shot noise measurements for all components
        to append to.
    tag : str, optional
        Data file common tag (default is '').
    sn_range : sequence of int of length 2, optional
        Data file serial number range (default is `None`).

    """
    # Extract component.
    stat, _ = split_component(component)
    fnpattern = component + tag

    print(f"Reading {component} measurements.")

    bins, coords, counts, measurements_raw, measurements_shotnoise = \
        read_measurements(stat, directory, fnpattern, sn_range=sn_range)

    print(f"Read {len(measurements_raw)} {component} measurement files.")

    # Append component.
    components_bins.append(bins)
    components_coords.append(coords)
    components_counts.append(counts)
    components_vectors.append(measurements_raw)
    components_vectors_shotnoise.append(measurements_shotnoise)


def combine_components(components_bins, components_coords, components_counts,
                       components_vectors, components_vectors_shotnoise,
                       save=None):
    """Combine statistic components.

    Parameters
    ----------
    components_bins : list of :class:`numpy.ndarray`
        Data bin centres for all components.
    components_coords : list of :class:`numpy.ndarray`
        Data bin effective coordinates for all components.
    components_counts : list of :class:`numpy.ndarray`
        Data bin contribution counts for all components.
    components_vectors : list of :class:`numpy.ndarray`
        Data vectors of raw measurements for all components.
    components_vectors_shotnoise : list of :class:`numpy.ndarray`
        Data vectors of shot noise measurements for all components.

    Returns
    -------
    tuple of :class:`numpy.ndarray` of length 5
        Combined flattened bin centres,
        flattened bin effective coordinates,
        flattened bin contribution counts,
        stacked raw measurement data vectors
        and stacked shot noise measurement data vectors.

    """
    # Check numbers of components agree.
    if not (
        len(components_bins)
        == len(components_coords)
        == len(components_counts)
        == len(components_vectors)
        == len(components_vectors_shotnoise)
    ):
        raise ValueError("Number of components do not agree.")

    # Combine components.
    data_bin_centres = np.concatenate(components_bins)
    data_bin_coords = np.concatenate(components_coords)
    data_bin_counts = np.concatenate(components_counts)
    data_vectors = np.hstack(components_vectors)
    data_vectors_shotnoise = np.hstack(components_vectors_shotnoise)

    # Save combined data vector.
    if save:
        np.savez(
            save,
            data_bin_centres=data_bin_centres,
            data_bin_coords=data_bin_coords,
            data_bin_counts=data_bin_counts,
            data_vectors=data_vectors,
            data_vectors_shotnoise=data_vectors_shotnoise
        )

    return data_bin_centres, data_bin_coords, data_bin_counts, \
        data_vectors, data_vectors_shotnoise


if __name__ == '__main__':

    # Configure.
    config = configure()

    # Extract components.
    components_bins, components_coords = [], []
    components_counts = []
    components_vectors, components_vectors_shotnoise = [], []
    for component in config.components:
        insert_component(
            component, config.input_dir,
            components_bins, components_coords,
            components_counts,
            components_vectors, components_vectors_shotnoise,
            tag=config.tag, sn_range=config.range
        )

    # Combine components.
    component_names = ''.join(config.components)
    components_tag = config.tag.replace('(.*)', '').replace('__', '_')

    combine_components(
        components_bins, components_coords, components_counts,
        components_vectors, components_vectors_shotnoise,
        save=f"{config.output_dir}/{component_names}{components_tag}.npz"
    )
