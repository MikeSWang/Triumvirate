"""Reshape full-form three-point clustering measurement data table into
2-d arrays.

"""
import numpy as np


def reshape_threept_datatab(paired_coords, flat_multipole):
    """Reshape flattended full-form three-point clustering measurements
    into 2-d arrays.

    Parameters
    ----------
    paired_coords : 2-d array of float
        Paired coordinate points.
    flat_multipole : 1-d array of float
        Flattened multipole values.

    Returns
    -------
    common_coords : 1-d array of float
        Common coordinate points.
    multipole : 2-d array of float
        Multipole values.

    Raises
    ------
    ValueError
        If `paired_coords` is not a 2-d array.
    ValueError
        If `paired_coords` and `flat_multipole` have different lengths.
    ValueError
        If `paired_coords` and `flat_multipole` do not correspond to an
        upper triangular matrix.

    """
    # Check dimensions and reshape.
    if np.array(paired_coords).ndim != 2:
        raise ValueError("Paired coordinate points must be a 2-d array.")
    if (flat_len := len(paired_coords)) != len(flat_multipole):
        raise ValueError(
            "Paired coordinate points and flattened multipole "
            "array must have the same length."
        )

    nbins = (np.sqrt(8*flat_len + 1) - 1) / 2.
    if not np.isclose(nbins, int(nbins)):
        raise ValueError(
            "Input arrays do not correspond to an upper triangular matrix."
        )
    nbins = int(nbins)

    triu_indices = np.triu_indices(nbins)

    common_coords = paired_coords[:nbins, 1]

    multipole = np.zeros((nbins, nbins))
    multipole[triu_indices] = flat_multipole
    multipole.T[triu_indices] = flat_multipole

    return common_coords, multipole
