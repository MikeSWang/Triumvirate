"""
Catalogue (:mod:`~triumvirate.catalogue`)
==========================================================================

Handle catalogue I/O and processing.

.. autosummary::
    MissingValueError
    DefaultValueWarning
    ParticleCatalogue

"""
import warnings

import numpy as np
from astropy.table import Table

from ._particles import _ParticleCatalogue

try:
    from nbodykit.source.catalog import (
        ArrayCatalog,
        BinaryCatalog,
        CSVCatalog,
        # FITSCatalog,
        HDFCatalog,
    )
    _nbkt_imported = True
except Exception:
    _nbkt_imported = False


class MissingValueError(ValueError):
    """Value error raised when a mandatory field is missing/empty in
    a catalogue.

    """
    pass


class DefaultValueWarning(UserWarning):
    """Warning issued when values of a field are not provided and set
    to the default.

    """
    pass


class ParticleCatalogue:
    """Catalogue holding particle coordinates, weights and
    redshift-dependent mean number density.

    Parameters
    ----------
    x, y, z : 1-d array of float
        Cartesian coordinates of particles.  `x`, `y` and `z` must have
        the same length.
    nz : (1-d array of) float, optional
        Redshift-dependent mean number density (defaults is `None`).
        If an array, it must be of the same length as `x`, `y` and `z`.
        This quantity should include sample weighting if sample weights
        are used; see the :ref:`note <sample_and_clustering_weights>`
        below for more details.
    ws, wc : (1-d array of) float, optional
        Sample weight and clustering weight of particles (defaults
        are 1.).  If an array, it must be of the same length as
        `x`, `y` and `z`.  See the
        :ref:`note <sample_and_clustering_weights>` below
        for more details.
    logger : :class:`logging.Logger`, optional
        Program logger (default is `None`).

    Attributes
    ----------
    bounds : dict of {str: tuple of (float, float)}
        Particle coordinate bounds.
    ntotal : int
        Total particle number.
    wtotal : float
        Total particle overall weight.
    wstotal : float
        Total particle sample weight.


    .. versionchanged:: 0.2.0
       Separate `wtotal` and `wstotal` attributes.

    .. _sample_and_clustering_weights:

    .. admonition:: Sample and clustering weights

        There are two types of weights: sample weight ``ws``
        (e.g. completeness weights) and clustering weight ``wc``
        (e.g. Feldman--Kaiser--Peacock weights).  The overall weight
        is the product of the two for each particle.

        Note the naming convention above: in particular, ``wc`` is not
        the completeness weight (which is a component of ``ws`` instead).

        When sample weights are applied (e.g. the catalogue has been
        corrected for incompleteness), the redshift-dependent mean
        number density ``nz`` is the weighted mean number density
        (e.g. including the completeness weights).

    """

    def __init__(self, x, y, z, nz=None, ws=1., wc=1., logger=None):

        self._logger = logger

        # Initialise particle data (with 'astropy').
        if nz is None:
            warnings.warn(
                "Catalogue 'nz' field is None and thus set to zero, "
                "which may raise errors in some computations.",
                category=DefaultValueWarning
            )
            nz = 0.

        self._pdata = Table()
        self._pdata['x'] = x
        self._pdata['y'] = y
        self._pdata['z'] = z
        self._pdata['nz'] = nz
        self._pdata['ws'] = ws
        self._pdata['wc'] = wc

        if _nbkt_imported:
            self._pdata = ArrayCatalog(self._pdata)
            self._backend = 'nbodykit'
        else:
            self._backend = 'astropy'

        self._source = 'extdata:{}'.format(id(self._pdata))

        # Compute catalogue properties.
        self._calc_bounds(init=True)

        self.ntotal = len(self)
        self.wtotal = self._compute(
            (self._pdata['ws'] * self._pdata['wc']).sum()
        )
        self.wstotal = self._compute(self._pdata['ws'].sum())

        if self._logger:
            self._logger.info(
                "Catalogue initialised: "
                "ntotal = %d, wtotal = %.3f, wstotal = %.3f (%s).",
                self.ntotal, self.wtotal, self.wstotal, self
            )

    @classmethod
    def read_from_file(cls, filepath, reader='astropy', names=None,
                       format='ascii.no_header', table_kwargs=None,
                       file_kwargs=None, name_mapping=None, logger=None):
        """Read particle data from file.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`
            Catalogue file path.
        reader : {'astropy', 'nbodykit'}, optional
            If 'astropy' (default), |Table| is used for reading in the
            catalogue file; else if 'nbodykit', child classes of
            |FileCatalogBase| is used (if :mod:`nbodykit` is available).
        names : sequence of str, list of tuple or :class:`numpy.dtype`, \
                optional
            Catalogue file field names or data types.  If `None`
            (default), the header in the file may be used to infer the
            field names or data types.  See the
            :ref:`note <format_and_names_arguments>` below
            for more details.
        format : str, optional
            File format specifier (default is 'ascii.no_header' for
            default `reader` value 'astropy').  See the
            :ref:`note <format_and_names_arguments>` below
            for more details.
        name_mapping : dict of {str: str}, optional
            Mapping between any of the default column names 'x', 'y', 'z',
            'nz', 'ws' and 'wc' (keys) and the corresponding field names
            (values) in `names` (default is `None`).
        table_kwargs : dict, optional
            Keyword arguments to be passed to
            :attr:`astropy.table.Table.read` (default is `None`).
            Used only when ``reader='astropy'``.
        file_kwargs : dict, optional
            Keyword arguments to be passed to child classes of
            :class:`nbodykit.source.catalog.file.FileCatalogBase` (default
            is `None`). Used only when ``reader='nbodykit'``.  See also
            the hint below.
        logger : :class:`logging.Logger`, optional
            Program logger (default is `None`).


        .. _format_and_names_arguments:

        .. admonition:: `format` and `names` arguments

            For ``reader='astropy'``, supported `format` can be found in
            `'Built-In Table Readers/Writers'`_, and `names` correspond to
            the ``names`` keyword argument (sequence of str) in
            :attr:`astropy.table.Table.read` for a subset of formats.
            See |Table| for more details including appropriate
            `table_kwargs`.

            For ``reader='nbodykit'``, supported `format` and the
            argument corresponding to `names` in the reader are---

            - ``'text'`` and ``names`` (sequence of str): plain-text files
              read in by |CSVCatalog|;
            - ``'binary'`` and ``dtype`` (list of tuple or
              :class:`numpy.dtype`): binary files read in
              by |BinaryCatalog|;
            - ``'hdf5'`` and not applicable: HDF files read in
              by |HDFCatalog| (pass `file_kwargs` as appropriate).

            See `'Reading Catalogs from Disk'`_ for more details.


        .. hint::

            If any of the 'nz', 'ws' and 'wc' columns are not provided,
            these columns are initialised with the default values in
            :class:`~triumvirate.catalogue.ParticleCatalogue`.


        Examples
        --------
        >>> ParticleCatalogue.read_from_file(
        ...     filepath,
        ...     reader='astropy',
        ...     format='fits',
        ...     # Match original data column names to the default names.
        ...     name_mapping={
        ...         'x': 'X', 'y': 'Y', 'z': 'Z',
        ...         'nz': 'NZ', 'ws': 'WEIGHT_SYS', 'wc': 'WEIGHT_FKP',
        ...     }
        ... )
        >>> ParticleCatalogue.read_from_file(
        ...     filepath,
        ...     reader='nbodykit',
        ...     format='hdf5',
        ...     # `root` keyword argument is passed to `HDFCatalog`.
        ...     file_kwargs={root='particles'}
        ... )


        .. |Table| replace:: :class:`astropy.table.Table`

        .. |FileCatalogBase| replace:: \
            :class:`nbodykit.source.catalog.file.FileCatalogBase`

        .. |CSVCatalog| replace:: \
            :class:`nbodykit.source.catalog.file.CSVCatalog`

        .. |FITSCatalog| replace:: \
            :class:`nbodykit.source.catalog.file.FITSCatalog`

        .. |BinaryCatalog| replace::  \
            :class:`nbodykit.source.catalog.file.BinaryCatalog`

        .. |HDFCatalog| replace:: \
            :class:`nbodykit.source.catalog.file.HDFCatalog`

        .. _'Built-In Table Readers/Writers': \
            https://docs.astropy.org/en/latest/io/unified.html
            #built-in-table-readers-writers

        .. _'Reading Catalogs from Disk': \
            https://nbodykit.readthedocs.io/en/latest/catalogs/reading.html

        """
        self = object.__new__(cls)

        self._logger = logger
        self._source = f'extfile:{filepath}'

        if reader.lower() == 'nbodykit':
            if _nbkt_imported:
                self._backend = 'nbodykit'
            else:
                warnings.warn(
                    "'astropy' is used for catalogue I/O as "
                    "'nbodykit' is unavailable",
                    category=RuntimeWarning
                )
                self._backend = 'astropy'
        elif reader.lower() == 'astropy':
            self._backend = 'astropy'
        else:
            raise ValueError(
                f"Unsupported catalogue file reader: {reader}. "
                "Possible options: {'astropy', 'nbodykit'}."
            )

        # Initialise particle data.
        if self._backend == 'nbodykit':
            if file_kwargs is None:
                file_kwargs = {}

            filepath = str(filepath)
            if format == 'text':
                self._pdata = CSVCatalog(filepath, names, **file_kwargs)
            elif format == 'fits':
                raise ValueError(
                    "'nbodykit' FITS reader is broken. "
                    "Use 'astropy' FITS reader instead."
                )
            elif format == 'binary':
                self._pdata = \
                    BinaryCatalog(filepath, dtype=names, **file_kwargs)
            elif format == 'hdf5':
                self._pdata = HDFCatalog(filepath, **file_kwargs)
            else:
                raise ValueError(
                    "Unsupported `format` for ``reader='nbodykit'``: {}."
                    .format(filepath)
                )

            if name_mapping:
                for name_, name_alt_ in name_mapping.items():
                    self._pdata[name_] = self._pdata[name_alt_]
                    try:
                        del self._pdata[name_alt_]
                    except ValueError:
                        pass

        if self._backend == 'astropy':
            if table_kwargs is None:
                table_kwargs = {}
            if format.startswith('ascii'):
                table_kwargs.update(names=names)

            self._pdata = Table.read(
                filepath, format=format, **table_kwargs
            )

            if name_mapping:
                for name_, name_alt_ in name_mapping.items():
                    self._pdata.rename_column(name_alt_, name_)

        # Validate particle data columns.
        if self._backend == 'nbodykit':
            colnames = self._pdata.columns
        if self._backend == 'astropy':
            colnames = self._pdata.colnames

        for axis_name in ['x', 'y', 'z']:
            if axis_name not in colnames:
                raise MissingValueError(
                    f"Mandatory field {axis_name} is missing."
                )

        if 'nz' not in colnames:
            self._pdata['nz'] = 0.
            warnings.warn(
                "Catalogue 'nz' field is not provided and thus set to zero, "
                "which may raise errors in some computations.",
                category=DefaultValueWarning
            )

        for name_wgt in ['ws', 'wc']:
            if name_wgt not in colnames:
                warnings.warn(
                    f"Catalogue '{name_wgt}' field is not provided, "
                    "so is set to unity.",
                    category=DefaultValueWarning
                )
                self._pdata[name_wgt] = 1.

        # Compute catalogue properties.
        self._calc_bounds(init=True)

        self.ntotal = len(self._pdata)
        self.wtotal = self._compute(
            (self._pdata['ws'] * self._pdata['wc']).sum()
        )
        self.wstotal = self._compute(self._pdata['ws'].sum())

        if self._logger:
            self._logger.info(
                "Catalogue loaded: "
                "ntotal = %d, wtotal = %.3f, wstotal = %.3f (%s).",
                self.ntotal, self.wtotal, self.wstotal, self
            )

        return self

    def __str__(self):
        try:
            return "ParticleCatalogue(source={})".format(self._source)
        except NameError:
            return "ParticleCatalogue(address={})".format(
                object.__str__(self).split()[-1].lstrip().rstrip('>')
            )

    def __len__(self):
        return len(self._pdata)

    def __getitem__(self, key):
        """Return one or more data entries.

        Parameters
        ----------
        key : (list of) str, (list of) int or slice
            Data entry key(s) or selectors (indices or slices).

        Returns
        -------
        :class:`numpy.ndarray`
            Data entry/entries.

        """
        return self._compute(self._pdata[key])

    def __setitem__(self, key, val):
        """Set data column.

        Parameters
        ----------
        key : str
            Data column name.
        val : 1-d array_like
            Data column array.

        """
        if isinstance(key, str):
            self._pdata[key] = val
        else:
            raise ValueError('Cannot set values for non-string-type key.')

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
        boxsize : float or sequence of [float, float, float], optional
            Box size (in each dimension) over which the mean density
            is calculated.  Used only when `volume` is `None`.


        .. attention::

            The invocation of this method resets the particle data
            column ``'nz'``.

        """
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        volume = volume or np.prod(boxsize)

        self._pdata['nz'] = self.ntotal / volume

    def centre(self, boxsize, catalogue_ref=None):
        """Centre a (pair of) catalogue(s) in a box.

        Parameters
        ----------
        boxsize : float or sequence of [float, float, float]
            Box size (in each dimension).
        catalogue_ref : :class:`~triumvirate.catalogue.ParticleCatalogue`, \
                        optional
            Reference catalogue used for box alignment, also to be centred
            in the same box.  If `None` (default), the current catalogue
            itself is used as the reference catalogue.


        .. note::

            The reference catalogue is typically the random-source
            catalogue (if provided).  Particle coordinates in both
            catalogues are shifted by the same displacement vector such
            that the mid-point of particle coordinate extents in the
            reference catalogue is at the centre of the box.

        """
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        if catalogue_ref is None:
            catalogue_ref_ = self
        else:
            catalogue_ref_ = catalogue_ref

        _axes_overflow = self._check_bounds_in_boxsize(boxsize)
        if _axes_overflow:
            warnings.warn(
                "Catalogue extent exceeds the box size along axis {} ({})."
                "Some partcles may lie outside the box after centring."
                .format(set(_axes_overflow), self)
            )
        if catalogue_ref is not None:
            _axes_overflow_ref = \
                catalogue_ref_._check_bounds_in_boxsize(boxsize)
            if _axes_overflow_ref:
                warnings.warn(
                    "Catalogue extent exceeds the box size along axis {} ({})."
                    "Some partcles may lie outside the box after centring."
                    .format(set(_axes_overflow_ref), catalogue_ref_)
                )

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

        The particle coordinates are shifted away from the box corner
        at the origin such that in each dimension the minimum particle
        coordinate is given by the amount of padding, which can be set as
        a fraction of the box size or a multiple of the grid cell size.

        Parameters
        ----------
        boxsize : float or sequence of [float, float, float]
            Box size in each dimension.
        ngrid : int or sequence of [int, int, int], optional
            Grid number in each dimension (default is `None`).
            Must be provided if `ngrid_pad` is set.
        boxsize_pad : float or sequence of [float, float, float], optional
            Box size padding factor.  If not `None` (default), then
            `ngrid_pad` must be `None`.
        ngrid_pad : float or sequence of [float, float, float], optional
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


        .. note::

            The reference catalogue is typically the random-source
            catalogue (if provided).  Particle coordinates in both
            catalogues are shifted by the same displacement vector such
            that the minimum particle coordinates in the reference
            catalogue are the amounts of padding specified.

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
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        if catalogue_ref is None:
            catalogue_ref_ = self
        else:
            catalogue_ref_ = catalogue_ref

        _axes_overflow = self._check_bounds_in_boxsize(boxsize)
        if _axes_overflow:
            warnings.warn(
                "Catalogue extent exceeds the box size along axis {} ({})."
                "Some partcles may lie outside the box after padding."
                .format(set(_axes_overflow), self)
            )
        if catalogue_ref is not None:
            _axes_overflow_ref = \
                catalogue_ref_._check_bounds_in_boxsize(boxsize)
            if _axes_overflow_ref:
                warnings.warn(
                    "Catalogue extent exceeds the box size along axis {} ({})."
                    "Some partcles may lie outside the box after padding."
                    .format(set(_axes_overflow_ref), catalogue_ref_)
                )

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

        for iaxis, axis in enumerate(['x', 'y', 'z']):
            if (max(self.bounds[axis]) > boxsize[iaxis]
                    or max(catalogue_ref_.bounds[axis]) > boxsize[iaxis]):
                warnings.warn(
                    "`boxsize` is smaller than the particle extents "
                    f"in {axis} axis. Some particles now lie outside "
                    "the box after padding."
                )

    def periodise(self, boxsize):
        """Place particles in a periodic box of given box size.

        Parameters
        ----------
        boxsize : float or sequence of [float, float, float]
            Box size (in each dimension).

        """
        if np.isscalar(boxsize):
            boxsize = [boxsize, boxsize, boxsize]

        _axes_overflow = self._check_bounds_in_boxsize(boxsize)
        if _axes_overflow:
            warnings.warn(
                "Box size is smaller than particle coordinate extents "
                "along axis: {}."
                .format(_axes_overflow)
            )

        for iaxis, axis in enumerate(['x', 'y', 'z']):
            # Also centre.
            self._pdata[axis] = (
                self._pdata[axis] + boxsize[iaxis] / 2. - (
                    self.bounds[axis][1] + self.bounds[axis][0]
                ) / 2.
            ) % boxsize[iaxis]

        self._calc_bounds()

    def offset_coords(self, origin):
        """Offset particle coordinates for a given origin.

        Parameters
        ----------
        origin : float or sequence of [float, float, float]
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
        self.spans = {}
        for axis in ['x', 'y', 'z']:
            self.bounds[axis] = (
                self._compute(self._pdata[axis].min()),
                self._compute(self._pdata[axis].max())
            )
            self.spans[axis] = self.bounds[axis][1] - self.bounds[axis][0]

        if self._logger:
            if init:
                self._logger.info(
                    "Original extents of particle coordinates: "
                    "{'x': (%.3f, %.3f | %.3f),"
                    " 'y': (%.3f, %.3f | %.3f),"
                    " 'z': (%.3f, %.3f | %.3f)}"
                    " (%s).",
                    *self.bounds['x'], self.spans['x'],
                    *self.bounds['y'], self.spans['y'],
                    *self.bounds['z'], self.spans['z'],
                    self
                )
            else:
                self._logger.info(
                    "Offset extents of particle coordinates: "
                    "{'x': (%.3f, %.3f | %.3f),"
                    " 'y': (%.3f, %.3f | %.3f),"
                    " 'z': (%.3f, %.3f | %.3f)}"
                    " (%s).",
                    *self.bounds['x'], self.spans['x'],
                    *self.bounds['y'], self.spans['y'],
                    *self.bounds['z'], self.spans['z'],
                    self
                )

    def _check_bounds_in_boxsize(self, boxsize):
        """Check if partice coordinate extents are covered by the box size
        in all dimensions, and return the axis name(s) if not.

        Parameters
        ----------
        boxsize : sequence of float
            Box size in all dimensions.

        Returns
        -------
        list of str
            Axis if and where box size is smaller than the particle
            coordinate extent.

        """
        _axes = []
        for iaxis, axis in enumerate(['x', 'y', 'z']):
            if np.abs(np.diff(self.bounds[axis])) > boxsize[iaxis]:
                _axes.append(axis)

        return _axes

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
        x = np.ascontiguousarray(self._compute(self._pdata['x']))
        y = np.ascontiguousarray(self._compute(self._pdata['y']))
        z = np.ascontiguousarray(self._compute(self._pdata['z']))
        nz = np.ascontiguousarray(self._compute(self._pdata['nz']))
        ws = np.ascontiguousarray(self._compute(self._pdata['ws']))
        wc = np.ascontiguousarray(self._compute(self._pdata['wc']))

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
                    .format(self._source),  # noqa: E131
                "Catalogue size: "
                "ntotal = {:d}, wtotal = {:.3f}, wstotal = {:.3f}"
                    .format(self.ntotal, self.wtotal, self.wstotal),
                "Catalogue particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                    .format(
                        *self.bounds['x'], *self.bounds['y'], *self.bounds['z']
                    ),
            ]
        else:
            text_lines = [
                "Data catalogue source: {}"
                    .format(self._source),  # noqa: E131
                "Data catalogue size: "
                "ntotal = {:d}, wtotal = {:.3f}, wstotal = {:.3f}"
                    .format(self.ntotal, self.wtotal, self.wstotal),
                "Data-source particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                    .format(
                        *self.bounds['x'], *self.bounds['y'], *self.bounds['z']
                    ),
                "Random catalogue source: {}"
                    .format(catalogue_ref._source),  # noqa: E131
                "Random catalogue size: "
                "ntotal = {:d}, wtotal = {:.3f}, wstotal = {:.3f}"
                    .format(  # noqa: E131
                        catalogue_ref.ntotal,
                        catalogue_ref.wtotal,
                        catalogue_ref.wstotal,
                    ),
                "Random-source particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                    .format(  # noqa: E131
                        *catalogue_ref.bounds['x'],
                        *catalogue_ref.bounds['y'],
                        *catalogue_ref.bounds['z']
                    ),
            ]

        text_header = "\n".join(text_lines)

        return text_header
