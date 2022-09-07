"""
Catalogue (:mod:`~triumvirate.catalogue`)
===========================================================================

Handle catalogue I/O and processing.

"""
import warnings

import numpy as np
from astropy.table import Table

from triumvirate._catalogue import _ParticleCatalogue


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
    There are two types of weights: systematic weights 'ws' and
    clustering weights 'wc'.

    Parameters
    ----------
    x, y, z : 1-d array of float
        Cartesian coordinates of particles.  Lengths of `x`, `y` and `z`
        must agree.
    nz : (1-d array of) float, optional
        Redshift-dependent mean number density (defaults is `None`).
        If an array, its length must be the same as those of `x`, `y`
        and `z`.
    ws, wc : (1-d array of) float, optional
        Systematic weights and clustering weights of particles (defaults
        are 1).  If an array, its length must be the same as those of
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
        Total systematic weight.

    """

    def __init__(self, x, y, z, nz=None, ws=1., wc=1., logger=None):

        self._logger = logger

        # Initialise particle data.
        if nz is None:
            warnings.warn(
                "Catalogue 'nz' field is set to `None`, "
                "which may raise errors in some computations."
            )

        self._pdata = Table()
        self._pdata['x'] = x
        self._pdata['y'] = y
        self._pdata['z'] = z
        self._pdata['nz'] = nz
        self._pdata['ws'] = ws
        self._pdata['wc'] = wc

        # Compute catalogue properties.
        self._calc_bounds(init=True)

        self.ntotal = len(self._pdata)
        self.wtotal = np.sum(self._pdata['ws'])

        if self._logger:
            self._logger.info(
                "Catalogue initialised: %d particles with "
                "total systematic weights %.3f (%s).",
                self.ntotal, self.wtotal, self
            )

    @classmethod
    def read_from_file(cls, filepath, format='ascii.no_header', names=None,
                       name_mapping=None, table_kwargs=None, logger=None):
        """Read particle data from file.

        Notes
        -----
        If any of 'nz', 'ws' and 'wc' columns are not provided,
        these columns are initialised with the default values in
        :meth:`~triumvirate.catalogue.__init__`.

        For supported file data formats, see
        `<https://docs.astropy.org/en/stable/io/ascii/index.html#supported-formats>`_.

        Parameters
        ----------
        filepath : str or :class:`pathlib.Path`
            Catalogue file path.
        format : str, optional
            File format specifier (default is 'ascii.no_header', with any
            header included as comment lines).
        names : sequence of str, optional
            Catalogue file field names.  If `None` (default), the header
            row in the file is used to provide the names.
        name_mapping : dict of {str: str}, optional
            Mapping between any of the default column names 'x', 'y', 'z',
            'nz', 'ws' and 'wc' (keys) and the corresponding field names
            (values) in `names` (default is `None`).
        table_kwargs : dict, optional
            Keyword arguments to be passed to
            :meth:`astropy.table.Table.read` (default is `None`).
        logger : :class:`logging.Logger`, optional
            Program logger (default is `None`).

        """
        self = object.__new__(cls)

        self._logger = logger
        self._source = str(filepath)  # non-existent in `__init__`

        # Initialise particle data.
        if table_kwargs is None:
            table_kwargs = {}

        self._pdata = Table.read(
            filepath, format=format, names=names, **table_kwargs
        )

        # Validate particle data columns.
        if name_mapping is not None:
            for name_, name_alt in name_mapping.items():
                self._pdata.rename_column(name_alt, name_)

        for name_axis in ['x', 'y', 'z']:
            if name_axis not in self._pdata.colnames:
                raise MissingField(f"Mandatory field {name_axis} is missing.")

        if 'nz' not in self._pdata.colnames:
            self._pdata['nz'] = None
            warnings.warn(
                "Catalogue 'nz' field is set to `None`, "
                "which may raise errors in some computations."
            )

        for name_wgt in ['ws', 'wc']:
            if name_wgt not in self._pdata.colnames:
                self._pdata[name_wgt] = 1.

        # Compute catalogue properties.
        self._calc_bounds(init=True)

        self.ntotal = len(self._pdata)
        self.wtotal = np.sum(self._pdata['ws'])

        if self._logger:
            self._logger.info(
                "Catalogue loaded: %d particles with "
                "total systematic weights %.3f (%s).",
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
        if isinstance(key, (int, slice)):
            return self._pdata[key]
        return [self._pdata[key_] for key_ in key]

    def __setitem__(self, key, value):
        if isinstance(key, (int, str, slice)):
            self._pdata[key] = value
        else:
            for key_ in key:
                if isinstance(key_, (int, str)):
                    self._pdata[key_] = value
                else:
                    raise ValueError('Cannot interpret item key.')

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
            Normalised lines of sight.

        """
        los_norm = np.sqrt(
            self._pdata['x']**2 + self._pdata['y']**2 + self._pdata['z']**2
        )

        los = np.transpose([
            self._pdata['x'] / los_norm,
            self._pdata['y'] / los_norm,
            self._pdata['z'] / los_norm
        ])

        return los

    def centre(self, boxsize, catalogue_ref=None):
        """Centre a (pair of) catalogue(s) in a box.

        Parameters
        ----------
        boxsize : (3,) array of float
            Box size in each dimension.
        catalogue_ref : :class:`~triumvirate.catalogue.ParticleCatalogue`, optional
            Reference catalogue used for box alignment, to be put in the
            same box.  If `None` (default), the current catalogue itself
            is used as the reference catalogue.

        Notes
        -----
        The reference catalogue is typically the random-source catalogue
        (if provided).  Padding is applied at the box corner used as the
        new coordinate origin.

        """
        if catalogue_ref is None:
            origin = [
                np.mean(self.bounds[axis]) - boxsize[iaxis] / 2.
                for iaxis, axis in enumerate(['x', 'y', 'z'])
            ]

            self.offset_coords(origin)
        else:
            origin = np.array([
                np.mean(catalogue_ref.bounds[axis]) - boxsize[iaxis]/2.
                for iaxis, axis in enumerate(['x', 'y', 'z'])
            ])

            self.offset_coords(origin)
            catalogue_ref.offset_coords(origin)

    def pad(self, boxsize, ngrid=None, boxsize_pad=None, ngrid_pad=None,
            catalogue_ref=None):
        """Pad a (pair of) catalogue(s) in a box.

        The minima of the particle coordinates in each dimension are
        shifted away from the origin of the box depending on the amount
        of padding, which can be set either as a multiple of the grid
        sizes (i.e. number of grids) or a multipole of the box sizes
        (i.e. a percentage).

        Parameters
        ----------
        boxsize : (3,) array of float
            Box size in each dimension.
        ngrid : (3,) array of int
            Grid number in each dimension (default is `None`).
        boxsize_pad : (sequence of) float, optional
            Box size padding factor.  If `None` (default), then
            `ngrid_pad` should be set; otherwise, `ngrid_pad` must
            be `None`.
        ngrid_pad : (sequence of) float, optional
            Grid padding factor.  If `None` (default), then `boxsize_pad`
            should be set; otherwise, `boxsize_pad` must be `None`.
        catalogue_ref : :class:`~triumvirate.catalogue.ParticleCatalogue`, optional
            Reference catalogue used for box alignment, to be put in the
            same box.  If `None` (default), the current catalogue itself
            is used as the reference catalogue.

        Raises
        ------
        ValueError
            If `boxsize_pad` and `ngrid_pad` are both set to not `None`.

        Notes
        -----
        The reference catalogue is typically the random-source catalogue
        (if provided).  Padding is applied at the box corner used as the
        new coordinate origin.

        """
        if boxsize_pad is None and ngrid_pad is None:
            warnings.warn(
                "`boxsize_pad` and `ngrid_pad` are both `None`. "
                "No padding is applied to the catalogue."
            )
            return
        if boxsize_pad is not None and ngrid_pad is not None:
            raise ValueError(
                "Conflicting padding as `boxsize_pad` and `ngrid_pad` "
                "are both set (not `None`)."
            )

        if catalogue_ref is None:
            origin = np.array([
                self.bounds[axis][0] for axis in ['x', 'y', 'z']
            ])
        else:
            origin = np.array([
                catalogue_ref.bounds[axis][0] for axis in ['x', 'y', 'z']
            ])

        if ngrid_pad:
            origin -= np.multiply(ngrid_pad, np.divide(boxsize, ngrid))
        if boxsize_pad:
            origin -= np.multiply(boxsize_pad, boxsize)

        self.offset_coords(origin)
        if catalogue_ref is not None:
            catalogue_ref.offset_coords(origin)

    def periodise(self, boxsize):
        """Place particles in a periodic box of given box size.

        Parameters
        ----------
        boxsize : (3,) array of float
            Box size in each dimension.

        """
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
            If `True` (default is `False`), particles positions are
            original and have not been offset previously.

        """
        self.bounds = {}
        for axis in ['x', 'y', 'z']:
            self.bounds[axis] = (
                np.min(self._pdata[axis]), np.max(self._pdata[axis])
            )

        if self._logger:
            if init:
                self._logger.info(
                    "Original extents of particle coordinates: %s (%s).",
                    self.bounds, self
                )
            else:
                self._logger.info(
                    "Offset extents of particle coordinates: %s (%s).",
                    self.bounds, self
                )

    def _calc_powspec_normalisation_from_particles(self, alpha=1.):
        """Calculate particle-based power spectrum normalisation.

        Parameters
        ----------
        alpha : float, optional
            Alpha ratio (default is 1.).

        Returns
        -------
        float
            Power spectrum normalisation factor.

        """
        if None in self._pdata['nz']:
            raise MissingField(
                "Cannot calculate power spectrum normalisation "
                "because of missing 'nz' value(s) in catalogue ({})."
                .format(self)
            )

        return 1. / alpha / np.sum(
            self._pdata['nz'] * self._pdata['ws'] * self._pdata['wc']**2
        )

    def _calc_powspec_shotnoise_from_particles(self, alpha=1.):
        """Calculate (unnormalised) particle-based power spectrum
        shot noise.

        Parameters
        ----------
        alpha : float, optional
            Alpha ratio (default is 1.).

        Returns
        -------
        float
            Shot noise power.

        """
        if None in self._pdata['nz']:
            raise MissingField(
                "Cannot calculate power spectrum shot noise "
                "because of missing 'nz' value(s) in catalogue ({})."
                .format(self)
            )

        return alpha**2 * np.sum(self._pdata['ws']**2 * self._pdata['wc']**2)


def _prepare_catalogue(catalogue):
    return _ParticleCatalogue(
        catalogue._pdata['x'],
        catalogue._pdata['y'],
        catalogue._pdata['z'],
        np.nan_to_num(np.array(catalogue._pdata['nz'], dtype=float)),
        catalogue._pdata['ws'],
        catalogue._pdata['wc']
    )
