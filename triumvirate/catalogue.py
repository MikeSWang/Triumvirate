"""
Catalogue I/O (:mod:`~triumvirate.measurements.catalogue`)
===========================================================================

Read and process catalogues.

"""
import numpy as np
from astropy.table import Table


class ParticleCatalogue:
    """Particle catalogue holding coordinates and weights.

    Parameters
    ----------
    filepath : str or :class:`pathlib.Path`
        Catalogue file path.
    names : sequence of str
        Catalogue field names.
    format : str
        File format specifier.
    alt_names : dict of {str: str}, optional
        Mapping between the default column names 'x', 'y', 'z', 'ws' and
        'wc' (keys) and the corresponding names (values) in `names`.

    Attributes
    ----------
    nparticles : int
        Number of particles.
    bbox : list of tuple
        Bounding coordinates for each dimension.

    See Also
    --------
    :class:`~astropy.table.Table`

    """

    def __init__(self, filepath, names, format='ascii', alt_names={}):

        # Construct data array.
        self._data_array = Table.read(filepath, names=names, format=format)

        for def_name, alt_name in alt_names.items():
            self._data_array.rename_column(alt_name, def_name)

        # Compute catalogue properties.
        self.nparticles = len(self._data_array)

        bbox = []
        for axis in ['x', 'y', 'z']:
            bbox.append((
                np.min(self._data_array[axis]), np.max(self._data_array[axis])
            ))
        self.bbox = bbox

    def offset_coords(self, origin):
        """Offset coordinates for a given origin or purpose.

        Parameters
        ----------
        origin : array of float
            Coordinates of the new origin.

        """
        for axis, centre in zip(['x', 'y', 'z'], origin):
            self._data_array[axis] -= centre
