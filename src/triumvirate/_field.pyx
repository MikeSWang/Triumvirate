"""
Mesh field and statistics (:mod:`~triumvirate._field`)
==========================================================================

Parse mesh fields and their statistics.

"""
# STYLE: Un-Pythonic import order to prioritise Cython imports.
from cython.operator cimport dereference as deref

import numpy as np
cimport numpy as np

from ._field cimport CppFieldStats
from .dataobjs cimport Binning
from .parameters cimport ParameterSet

np.import_array()


def _record_binned_vectors(Binning binning not None,
                           ParameterSet paramset not None):
    """Record binned vectors.

    Parameters
    ----------
    binning : :class:`~triumvirate.dataobjs.Binning`
        Binning for the measurements.
    paramset : :class:`~triumvirate.parameters.ParameterSet`
        Parameter set for the sampling mesh grid.

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
    fieldstats_ptr = new CppFieldStats(deref(paramset.thisptr), False)

    binned_vectors_struct = fieldstats_ptr.record_binned_vectors(
        deref(binning.thisptr), ''.encode('utf-8')
    )

    binned_vectors = np.empty(
        binned_vectors_struct.count,
        dtype=[
            ('index', 'i4'),
            ('lower_edge', 'f8'),
            ('upper_edge', 'f8'),
            ('vecx', 'f8'),
            ('vecy', 'f8'),
            ('vecz', 'f8'),
        ]
    )

    binned_vectors['index'] = np.asarray(binned_vectors_struct.indices)
    binned_vectors['lower_edge'] = np.asarray(
        binned_vectors_struct.lower_edges
    )
    binned_vectors['upper_edge'] = np.asarray(
        binned_vectors_struct.upper_edges
    )
    binned_vectors['vecx'] = np.asarray(binned_vectors_struct.vecx)
    binned_vectors['vecy'] = np.asarray(binned_vectors_struct.vecy)
    binned_vectors['vecz'] = np.asarray(binned_vectors_struct.vecz)

    return binned_vectors
