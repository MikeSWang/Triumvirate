"""
Catalogue (:mod:`~triumvirate.catalogue`)
===========================================================================

Handle catalogue I/O and processing.

"""
import numpy as np
from astropy.table import Table


class ParticleCatalogue:
    """Particle catalogue holding coordinates and weights.

    Notes
    -----
    There are two types of weights: systematic weights 'ws' and
    clustering weights 'wc'.

    Parameters
    ----------
    x, y, z, ws, wc : array of float
        The Cartesian coordinates, systematic weights and clustering
        weights of particles.

    Attributes
    ----------
    num_particles : int
        Number of particles.
    wgt_particles : float
        Systematic-weight sum of particles.
    bounds : dict of tuple
        Particle coordinate bounds.

    """

    def __init__(self, x, y, z, ws, wc):

        # Initialise particle data.
        self._data = Table()
        self._data['x'] = x
        self._data['y'] = y
        self._data['z'] = z
        self._data['ws'] = ws
        self._data['wc'] = wc

        # Compute catalogue properties.
        self.num_particles = len(self._data)
        self.wgt_particles = np.sum(self._data['ws'])

        self._calc_bounds()

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        if isinstance(key, (int, slice)):
            return self._data[key]
        return [self._data[i] for i in key]

    def __setitem__(self, key, value):
        if isinstance(key, (int, slice)):
            self._data[key] = value
        else:
            for i in key:
                if isinstance(i, int):
                    self._data[i] = value
                else:
                    raise ValueError('Cannot interpret item key.')

    def __iter__(self):
        return self._data.__iter__()

    def __next__(self):
        return self._data.__next__()

    @classmethod
    def read_from_file(cls, filepath, names, name_mapping={},
                       format='ascii', tab_kwargs={}):
        """Read particle data from file.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`
            Catalogue file path.
        names : sequence of str
            Catalogue file field names.
        name_mapping : dict of {str: str}, optional
            Mapping between the default column names 'x', 'y', 'z', 'ws' and
            'wc' (keys) and the corresponding names (values) in `names`.
        format : str, optional
            File format specifier (default is 'ascii').  See also
            :class:`astropy.table.Table`.
        tab_kwargs : dict, optional
            Keyword arguments to be passed to
            :meth:`astropy.table.Table.read`.

        """
        self = object.__new__(cls)

        if not {'x', 'y', 'z', 'ws', 'wc'} < set(name_mapping.keys()):
            raise ValueError(
                "`name_mapping` must contain all default column names."
            )

        # Initialise particle data.
        self._data = Table.read(
            filepath, names=names, format=format, **tab_kwargs
        )

        for name_, name_alt in name_mapping.items():
            self._data.rename_column(name_alt, name_)

        # Compute catalogue properties.
        self.num_particles = len(self._data)
        self.wgt_particles = np.sum(self._data['ws'])

        self._calc_bounds()

    def offset_coords(self, origin):
        """Offset particle coordinates for a given origin.

        Parameters
        ----------
        origin : (3,) array of float
            Coordinates of the new origin.

        """
        for axis, coord in zip(['x', 'y', 'z'], origin):
            self._data[axis] -= coord

        self._calc_bounds()

    def centre(self, boxsize):
        """Place particles in a periodic box of given box size.

        Parameters
        ----------
        boxsize : (3,) array of float
            Box size in each dimension.

        """
        origin = [
            np.mean(self.bounds[axis]) - boxsize[idx_axis] / 2.
            for idx_axis, axis in enumerate(['x', 'y', 'z'])
        ]

        self.offset_coords(origin)

    def periodise(self, boxsize):
        """Place particles in a periodic box of given box size.

        Parameters
        ----------
        boxsize : (3,) array of float
            Box size in each dimension.

        """
        self._data['x'] %= boxsize[0]
        self._data['y'] %= boxsize[1]
        self._data['z'] %= boxsize[2]

        self._calc_bounds()

    @staticmethod
    def align_catalogues_for_fft(catalogue_data, catalogue_rand,
                                 boxsize, ngrid, shift_factor=3.):
        """Align a pair of catalogues by offsetting particle positions
        for FFT (though mesh grid shift).

        Parameters
        ----------
        catalogue_data : :class:`triumvirate.catalogue.ParticleCatalogue`
            (Data-source) particle catalogue.
        catalogue_rand : :class:`triumvirate.catalogue.ParticleCatalogue`
            (Random-source) particle catalogue.
        boxsize : (3,) array of float
            Box size in each dimension.
        ngrid : (3,) array of int
            Grid number in each dimension.
        shift_factor : float, optional
            Offset grid shift factor (default is 3.).

        """
        dpos = np.subtract(
            [catalogue_rand.bounds[axis][0] for axis in ['x', 'y', 'z']],
            shift_factor * np.divide(boxsize, ngrid)
        )

        catalogue_data.offset_coords(dpos)
        catalogue_rand.offset_coords(dpos)

    def _calc_bounds(self):
        """Calculate coordinate bounds in each dimension.

        """
        bounds = {}
        for axis in ['x', 'y', 'z']:
            bounds[axis] = (
                np.min(self._data[axis]), np.max(self._data[axis])
            )

        self.bounds = bounds
