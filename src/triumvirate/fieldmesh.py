"""
Field and Mesh (:mod:`~triumvirate.fieldmesh`)
==========================================================================

.. versionadded:: 0.3.0


Handle fields and mesh grids.

.. autosummary::
    record_binned_vectors

"""
import numpy as np

from ._field import _record_binned_vectors
from .parameters import (
    _modify_sampling_parameters,
    fetch_paramset_template,
    ParameterSet,
)


def record_binned_vectors(binning, paramset=None, boxsize=None, ngrid=None):
    """Record binned vectors given a binning scheme and mesh grid
    parameters.

    Parameters
    ----------
    binning : :class:`~triumvirate.binning.Binning`
        Binning.
    paramset : :class:`~triumvirate.parameters.ParameterSet`, optional
        Parameters set including the mesh grid parameters.  If `None`
        (default), `boxsize` and `ngrid` must be provided.
    boxsize : float or sequence of [float, float, float], optional
        Box size in each dimension (default is `None`).  Must be provided
        if `paramset` is `None`.  This will override the corresponding
        entries in `paramset`.
    ngrid : int or sequence of [int, int, int], optional
        Grid number in each dimension (default is `None`).  Must be
        provided if `paramset` is set.  This will override the
        corresponding entries in `paramset`.

    Returns
    -------
    binned_vectors : :class:`numpy.ndarray`
        Binned vectors as a structured array with the following fields---

        - ``'index'``: bin index;
        - ``'lower_edge'``: lower edge of the bin;
        - ``'upper_edge'``: upper edge of the bin;
        - ``'vecx'``: x-component of the vector;
        - ``'vecy'``: y-component of the vector;
        - ``'vecz'``: z-component of the vector,

    """
    if paramset is None:
        paramset = fetch_paramset_template('dict')

        paramset.update(range=[binning.bin_min, binning.bin_max])  # filler
        paramset.update(num_bins=binning.num_bins)  # filler

        paramset.update(
            degrees={'ell1': 0, 'ell2': 0, 'ELL': 0}
        )  # placeholder

        if boxsize is None or ngrid is None:
            raise ValueError(
                "Either `paramset` or both `boxsize` and `ngrid` "
                "must be provided."
            )

    params_sampling = {}
    if boxsize is not None:
        params_sampling['boxsize'] = [boxsize] * 3 if np.isscalar(boxsize) \
            else list(boxsize)
    if ngrid is not None:
        params_sampling['ngrid'] = [ngrid] * 3 if np.isscalar(ngrid) \
            else list(ngrid)

    if params_sampling:
        paramset = _modify_sampling_parameters(paramset, params_sampling)

    if isinstance(paramset, dict):  # likely redundant but safe
        paramset = ParameterSet(param_dict=paramset)

    binned_vectors = _record_binned_vectors(binning, paramset)

    return binned_vectors
