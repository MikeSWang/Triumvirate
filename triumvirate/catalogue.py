"""
Catalogue (:mod:`~triumvirate.catalogue`)
==========================================================================

Handle catalogue I/O and processing.

"""
import warnings

import numpy as np
from astropy.table import Table

from triumvirate._particles import _ParticleCatalogue

try:
    from nbodykit.base.catalog import CatalogSource
    from nbodykit.source.catalog import CSVCatalog
    _nbkt_imported = True
except Exception:
    _nbkt_imported = False


class MissingField(ValueError):
    """Value error raised when a mandatory field is missing/empty in
    a catalogue.

    """
    pass


class ParticleCatalogue:
    """Catalogue holding particle coordinates, weights and
    redshift-dependent mean number densities.

    Notes
    -----
    There are two types of weights: sample weights ``ws``
    (e.g. completeness weights) and clustering weights ``wc``
    (e.g. Feldman--Kaiser--Peacock weights).

    .. attention::

        Note the naming convention above: in particular, ``wc`` is not
        the completeness weight (which is a component of ``ws`` instead).

    Parameters
    ----------
    x, y, z : 1-d array of float
        Cartesian coordinates of particles.  `x`, `y` and `z` must have
        the same length.
    nz : (1-d array of) float, optional
        Redshift-dependent mean number density of
        (defaults is `None`).  If an array, it must be of the same length
        as `x`, `y` and `z`.
    ws, wc : (1-d array of) float, optional
        Sample weights and clustering weights of particles (defaults
        are 1.).  If an array, it must be of the same length as
        `x`, `y` and `z`.
    logger : :class:`logging.Logger`, optional
        Program logger (default is `None`).

    Attributes
    ----------
    bounds : dict of {str: tuple}
        Particle coordinate bounds.
    ntotal : int
        Total particle number.
    wtotal : float
        Total sample weight.

    """

    def __init__(self, x, y, z, nz=None, ws=1., wc=1., logger=None):

        self._logger = logger
        self._source = 'extdata'

        # Initialise particle data.
        if _nbkt_imported:
            self._pdata = CatalogSource(comm=None)
            self._arch = 'nbodykit'
        else:
            self._pdata = Table()
            self._arch = 'astropy'

        if nz is None:
            warnings.warn(
                "Catalogue 'nz' field is None and thus set to zero, "
                "which may raise errors in some computations."
            )
            nz = 0.

        self._pdata['x'] = x
        self._pdata['y'] = y
        self._pdata['z'] = z
        self._pdata['nz'] = nz
        self._pdata['ws'] = ws
        self._pdata['wc'] = wc

        # Compute catalogue properties.
        self._calc_bounds(init=True)

        self.ntotal = len(self._pdata)
        self.wtotal = self._compute(self._pdata['ws'].sum())

        if self._logger:
            self._logger.info(
                "Catalogue initialised: %d particles with "
                "total sample weights %.3f (%s).",
                self.ntotal, self.wtotal, self
            )

    @classmethod
    def read_from_file(cls, filepath, reader='astropy', names=None,
                       format='ascii.no_header', table_kwargs=None,
                       name_mapping=None, logger=None):
        """Read particle data from file.

        Notes
        -----
        If any of 'nz', 'ws' and 'wc' columns are not provided,
        these columns are initialised with the default values in
        :meth:`~triumvirate.catalogue.__init__`.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`
            Catalogue file path.
        reader : {'astropy', 'nbodykit'}, optional
            If 'astropy' (default), :class:`astropy.table.Table`
            is used for reading in the catalogue file; else if 'nbodykit',
            :mod:`nbodykit.source.catalog` is used (if ``nbodykit``
            is installed).
        names : sequence of str, optional
            Catalogue file field names.  Cannot be `None` (default)
            if `reader` is 'nbodykit'.  If `None`, the header in the file
            is used to infer the names.
        format : str, optional
            File format specifier (default is 'ascii.no_header', with any
            header included as comment lines).  Used only when `reader`
            is 'astropy'.  See
            `<https://docs.astropy.org/en/stable/io/ascii/
            index.html#supported-formats>`_ for supported file formats.
        name_mapping : dict of {str: str}, optional
            Mapping between any of the default column names 'x', 'y', 'z',
            'nz', 'ws' and 'wc' (keys) and the corresponding field names
            (values) in `names` (default is `None`). Used only
            when `reader` is 'astropy'.
        table_kwargs : dict, optional
            Keyword arguments to be passed to
            :attr:`astropy.table.Table.read` (default is `None`).
            Used only when `reader` is 'astropy'.
        logger : :class:`logging.Logger`, optional
            Program logger (default is `None`).

        """
        self = object.__new__(cls)

        self._logger = logger
        self._source = f'extfile:{filepath}'

        if reader.lower() == 'nbodykit':
            if _nbkt_imported:
                self._arch = 'nbodykit'
            else:
                warnings.warn(
                    "'astropy' is used for catalogue I/O as "
                    "'nbodykit' is unavailable",
                    category=warnings.RuntimeWarning
                )
                self._arch = 'astropy'
        elif reader.lower() == 'astropy':
            self._arch = 'astropy'
        else:
            raise ValueError(
                f"Unsupported catalogue file reader: {reader}. "
                "Possible options: {'astropy', 'nbodykit'}."
            )

        # Initialise particle data.
        if self._arch == 'nbodykit':
            self._pdata = CSVCatalog(str(filepath), names)

        if self._arch == 'astropy':
            if table_kwargs is None:
                table_kwargs = {}

            self._pdata = Table.read(
                filepath, format=format, names=names, **table_kwargs
            )

            if name_mapping:
                for name_, name_alt_ in name_mapping.items():
                    self._pdata.rename_column(name_alt_, name_)

        # Validate particle data columns.
        if self._arch == 'nbodykit':
            colnames = self._pdata.columns
        if self._arch == 'astropy':
            colnames = self._pdata.colnames

        for axis_name in ['x', 'y', 'z']:
            if axis_name not in colnames:
                raise MissingField(f"Mandatory field {axis_name} is missing.")

        if 'nz' not in colnames:
            self._pdata['nz'] = 0.
            warnings.warn(
                "Catalogue 'nz' field is not provided and thus set to zero, "
                "which may raise errors in some computations."
            )

        for name_wgt in ['ws', 'wc']:
            if name_wgt not in colnames:
                warnings.warn(
                    f"Catalogue '{name_wgt}' field is not provided, "
                    "so is set to unity."
                )
                self._pdata[name_wgt] = 1.

        # Compute catalogue properties.
        self._calc_bounds(init=True)

        self.ntotal = len(self._pdata)
        self.wtotal = self._compute(self._pdata['ws'].sum())

        if self._logger:
            self._logger.info(
                "Catalogue loaded: %d particles with "
                "total sample weights %.3f (%s).",
                self.ntotal, self.wtotal, self
            )

        return self

    def __str__(self):
        try:
            return "ParticleCatalogue(source={})".format(self._source)
        except NameError:
            return "ParticleCatalogue(address={})".format(
                object.__str__(self).split()[-1].lstrip().rstrip('>')
            )

    def __getitem__(self, key):
        """Return data entry.

        Parameters
        ----------
        key : (list of) str, int, slice
            Data entry key.

        Returns
        -------
        :class:`numpy.ndarray`
            Data entry.

        """
        return self._compute(self._pdata[key])

    def __setitem__(self, key, value):
        if isinstance(key, str):
            self._pdata[key] = value
        else:
            raise ValueError('Cannot set values for non-string-type key.')

    def __len__(self):
        return len(self._pdata)

    def __iter__(self):
        return self._pdata.__iter__()

    def __next__(self):
        return self._pdata.__next__()

    def compute_los(self):
        """Compute the line of sight to each particle.

        Returns
        -------
        los : (N, 3) :class:`numpy.ndarray`
            Normalised line-of-sight vectors.

        """
        los_norm = np.sqrt(
            self._pdata['x']**2 + self._pdata['y']**2 + self._pdata['z']**2
        )

        los_norm[los_norm == 0.] = 1.

        los = np.transpose([
            self._pdata['x'] / los_norm,
            self._pdata['y'] / los_norm,
            self._pdata['z'] / los_norm
        ])

        return self._compute(los)

    def compute_mean_density(self, volume=None, boxsize=None):
        """Compute the mean density over a volume.

        This sets the 'nz' column to :attr:`ntotal` divided by
        `volume` or cubic (product of) `boxsize`, and is typically used
        for calculating the homogeneous background number density
        in a simulation box.

        Parameters
        ----------
        volume : double, optional
            Volume over which the mean density is calculated.
        boxsize : ((3,) array of) double, optional
            Box size (in each dimension) over which the mean density
            is calculated.  Used only when `volume` is `None`.

        """
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        volume = volume if volume else np.product(boxsize)

        self._pdata['nz'] = self.ntotal / volume

    def centre(self, boxsize, catalogue_ref=None):
        """Centre a (pair of) catalogue(s) in a box.

        Parameters
        ----------
        boxsize : ((3,) array of) float
            Box size (in each dimension).
        catalogue_ref : :class:`~triumvirate.catalogue.ParticleCatalogue`, \
                        optional
            Reference catalogue used for box alignment, also to be centred
            in the same box.  If `None` (default), the current catalogue
            itself is used as the reference catalogue.

        Notes
        -----
        The reference catalogue is typically the random-source catalogue
        (if provided).  The box corner closest to the minimal particle
        extents is used as the new coordinate origin.

        """
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        if catalogue_ref is None:
            catalogue_ref_ = self
        else:
            catalogue_ref_ = catalogue_ref

        origin = np.array([
            np.mean(catalogue_ref_.bounds[axis]) - boxsize[iaxis]/2.
            for iaxis, axis in enumerate(['x', 'y', 'z'])
        ])

        self.offset_coords(origin)
        if catalogue_ref is not None:
            catalogue_ref.offset_coords(origin)

    def pad(self, boxsize, ngrid=None, boxsize_pad=None, ngrid_pad=None,
            catalogue_ref=None):
        """Pad a (pair of) catalogue(s) in a box.

        The minimum particle extents in each dimension are shifted away
        from the origin of the box depending on the amount of padding,
        which can be set as a multiple of either the grid sizes or
        the box sizes (i.e. number of grids or fraction of boxsizes).

        Parameters
        ----------
        boxsize : (3,) array of float
            Box size in each dimension.
        ngrid : (3,) array of int, optional
            Grid number in each dimension (default is `None`).
            Must be provided if `ngrid_pad` is set.
        boxsize_pad : ((3,) array of) float, optional
            Box size padding factor.  If not `None` (default), then
            `ngrid_pad` must be `None`.
        ngrid_pad : ((3,) array of) float, optional
            Grid padding factor.  If not `None` (default), then
            `boxsize_pad` must be `None`.
        catalogue_ref : :class:`~triumvirate.catalogue.ParticleCatalogue`, \
                        optional
            Reference catalogue used for box alignment, also to be put in
            the same box.  If `None` (default), the current catalogue
            itself is used as the reference catalogue.

        Raises
        ------
        ValueError
            If `boxsize_pad` and `ngrid_pad` are both set to not `None`.

        Notes
        -----
        The reference catalogue is typically the random-source catalogue
        (if provided).  Padding is applied at the box corner closest to
        the minimal particle extents, which is used as the new origin.

        """
        if boxsize_pad is None and ngrid_pad is None:
            warnings.warn(
                "`boxsize_pad` and `ngrid_pad` are both None. "
                "No padding is applied to the catalogue."
            )
            return
        if boxsize_pad is not None and ngrid_pad is not None:
            raise ValueError(
                "Conflicting padding as `boxsize_pad` and `ngrid_pad` "
                "are both set (not None)."
            )

        if catalogue_ref is None:
            catalogue_ref_ = self
        else:
            catalogue_ref_ = catalogue_ref

        origin = np.array([
            catalogue_ref_.bounds[axis][0] for axis in ['x', 'y', 'z']
        ])

        if boxsize_pad:
            origin -= np.multiply(boxsize_pad, boxsize)
        if ngrid_pad:
            origin -= np.multiply(ngrid_pad, np.divide(boxsize, ngrid))

        self.offset_coords(origin)
        if catalogue_ref is not None:
            catalogue_ref.offset_coords(origin)

    def periodise(self, boxsize):
        """Place particles in a periodic box of given box size.

        Parameters
        ----------
        boxsize : ((3,) array of) float
            Box size (in each dimension).

        """
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        for iaxis, axis in enumerate(['x', 'y', 'z']):
            self._pdata[axis] %= boxsize[iaxis]

        self._calc_bounds()

    def offset_coords(self, origin):
        """Offset particle coordinates for a given origin.

        Parameters
        ----------
        origin : (3,) array of float
            Coordinates of the new origin.

        """
        for axis, coord in zip(['x', 'y', 'z'], origin):
            self._pdata[axis] -= coord

        self._calc_bounds()

    def _calc_bounds(self, init=False):
        """Calculate coordinate bounds in each dimension.

        Parameters
        ----------
        init : bool, optional
            If `True` (default is `False`), particle positions are
            original and have not been offset previously.

        """
        self.bounds = {}
        for axis in ['x', 'y', 'z']:
            self.bounds[axis] = (
                self._compute(self._pdata[axis].min()),
                self._compute(self._pdata[axis].max())
            )

        if self._logger:
            if init:
                self._logger.info(
                    "Original extents of particle coordinates: "
                    "{'x': (%.3f, %.3f), 'y': (%.3f, %.3f), 'z': (%.3f, %.3f)} "
                    "(%s).",
                    *self.bounds['x'], *self.bounds['y'], *self.bounds['z'],
                    self
                )
            else:
                self._logger.info(
                    "Offset extents of particle coordinates: "
                    "{'x': (%.3f, %.3f), 'y': (%.3f, %.3f), 'z': (%.3f, %.3f)} "
                    "(%s).",
                    *self.bounds['x'], *self.bounds['y'], *self.bounds['z'],
                    self
                )

    def _convert_to_cpp_catalogue(self, verbose=-1):
        """Convert to a C++-wrapped catalogue.

        Parameters
        ----------
        verbose : int, optional
            Verbosity level to be passed to the backend C++ logger
            (default is -1, i.e. unused).

        Returns
        -------
        :class:`~triumvirate._particles._ParticleCatalogue`
            C++-wrapped catalogue.

        """
        x = self._compute(self._pdata['x'])
        y = self._compute(self._pdata['y'])
        z = self._compute(self._pdata['z'])
        nz = self._compute(self._pdata['nz'])
        ws = self._compute(self._pdata['ws'])
        wc = self._compute(self._pdata['wc'])

        return _ParticleCatalogue(x, y, z, nz, ws, wc, verbose=verbose)

    def _compute(self, quant):
        """Return a quantity in standard form (i.e. apply
        :meth:`dask.array.Array.compute` to any Dask array).

        Parameters
        ----------
        quant
            Quantity.

        Returns
        -------
        The same quantity in a non-Dask-array form.

        """
        try:
            return quant.compute()
        except AttributeError:
            return quant

    def write_attrs_as_header(self, catalogue_ref=None):
        """Write catalogue attributes as a header.

        Parameters
        ----------
        catalogue_ref : :class:`~triumvirate.catalogue.ParticleCatalogue`, \
                        optional
            Reference catalogue (default is `None`) whose attributes are
            also written out.  This is typically the
            random-source catalogue.

        Returns
        -------
        text_header : str
            Catalogue attributes as a header string.

        """
        if catalogue_ref is None:
            text_lines = [
                "Catalogue source: {}"
                    .format(self._source),
                "Catalogue size: {:d} particles of total weight {:.3f}"
                    .format(self.ntotal, self.wtotal),
                "Catalogue particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                    .format(
                        *self.bounds['x'], *self.bounds['y'], *self.bounds['z']
                    ),
            ]
        else:
            text_lines = [
                "Data catalogue source: {}"
                    .format(self._source),
                "Data catalogue size: {:d} particles of total weight {:.3f}"
                    .format(self.ntotal, self.wtotal),
                "Data-source particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                    .format(
                        *self.bounds['x'], *self.bounds['y'], *self.bounds['z']
                    ),
                "Random catalogue source: {}"
                    .format(catalogue_ref._source),
                "Random catalogue size: {:d} particles of total weight {:.3f}"
                    .format(catalogue_ref.ntotal, catalogue_ref.wtotal),
                "Random-source particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                    .format(
                        *catalogue_ref.bounds['x'],
                        *catalogue_ref.bounds['y'],
                        *catalogue_ref.bounds['z']
                    ),
            ]

        text_header = "\n".join(text_lines)

        return text_header
