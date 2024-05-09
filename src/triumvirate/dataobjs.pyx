"""
Data Objects (:mod:`~triumvirate.dataobjs`)
==========================================================================

Define data objects.

.. autosummary::
    Binning

"""
import numpy as np

from .dataobjs cimport CppBinning


cdef class Binning:
    """Binning.

    Parameters
    ----------
    space : {'config', 'fourier'}
        Coordinate space.
    scheme : {'lin', 'log', 'linpad', 'logpad'}
        Binning scheme.  See the :ref:`note <binning_scheme>` below
        for more details.
    bin_min, bin_max : float, optional
        Minimum and maximum of the bin range (defaults are `None`).
    num_bins : int, optional
        Number of bins (default is `None`).

    Attributes
    ----------
    space : {'config', 'fourier'}
        Coordinate space.
    scheme : {'lin', 'log', 'linpad', 'logpad'}
        Binning scheme.  See the :ref:`note <binning_scheme>` below
        for more details.
    bin_min : float or None
        Minimum of the bin range.
    bin_max : float or None
        Maximum of the bin range.
    num_bins : int or None
        Number of bins.
    bin_edges : list of float
        Bin edges of length (:attr:`num_bins` + 1).
    bin_centres : list of float
        Bin centres of length :attr:`num_bins`.
    bin_widths : list of float
        Bin widths of length :attr:`num_bins`.


    .. _binning_scheme:

    .. admonition:: Binning `scheme`

        The bin setting method supports linear ('lin') and
        log-linear/exponential ('log') binning, with possible linear
        padding from zero for the first 5 bins ('linpad' or 'logpad').
        The padding is either 1.e-3 ('fourier' space) and 10. ('config'
        space), or determined by the mesh grid resolution if set using
        mesh grid sizes.

    """

    def __cinit__(self, space, scheme,
                  bin_min=None, bin_max=None, num_bins=None):

        self.thisptr = new CppBinning(
            space.encode('utf-8'), scheme.encode('utf-8')
        )

        self.scheme = scheme
        self.space = space

        _settable = True
        if bin_min is not None:
            self.bin_min = bin_min
        else:
            _settable = False
        if bin_max is not None:
            self.bin_max = bin_max
        else:
            _settable = False
        if num_bins is not None:
            self.num_bins = num_bins
        else:
            _settable = False

        if _settable:
            self.set_bins(self.bin_min, self.bin_max, self.num_bins)

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def from_parameter_set(cls, paramset):
        """Create binning from a parameter set.

        This sets :attr:`scheme` and :attr:`space`.

        Parameters
        ----------
        paramset : :class:`~triumvirate.parameters.ParameterSet`
            Parameter set.

        Raises
        ------
        ValueError
            When the 'space' parameter is unset or unrecognised in
            `paramset`.  See the :ref:`note <binning_space>` below
            for more details.


        .. _binning_space:

        .. attention::

            If `paramset` has been initialised without an appropriate
            'statistic_type' parameter value, the derived 'space'
            parameter is unset and must be specified to be the either of
            {'fourier', 'config'} before `paramset` is to this method.
            See also :class:`~triumvirate.parameters.ParameterSet`.

        """
        if paramset['space'] not in ['fourier', 'config']:
            raise ValueError(
                "'space' parameter is unset or unrecognised in `paramset`."
            )

        self = cls(
            paramset['space'], paramset['binning'],
            bin_min=paramset['range'][0],
            bin_max=paramset['range'][-1],
            num_bins=paramset['num_bins']
        )

        return self

    def set_bins(self, bin_min, bin_max, num_bins):
        """Set binning based on the bin range and number using the
        current binning scheme.

        Bin edges, centres and widths are recalculated.

        Parameters
        ----------
        bin_min, bin_max : float
            Minimum and maximum of the bin range.
        num_bins : int
            Number of bins.

        """
        self.bin_min, self.bin_max = bin_min, bin_max
        self.num_bins = num_bins

        self.thisptr.set_bins(self.bin_min, self.bin_max, self.num_bins)

        self.bin_edges = self.thisptr.bin_edges
        self.bin_centres = self.thisptr.bin_centres
        self.bin_widths = self.thisptr.bin_widths

    def set_grid_based_bins(self, boxsize, ngrid):
        """Set linear binning based on a mesh grid.

        The binning :attr:`scheme` is overriden to 'lin'.  The bin width
        is given by the grid resolution in configuration space or the
        fundamental wavenumber in Fourier space.  The bin minimum is zero
        and the bin maximum is half the box size in configuration space or
        the Nyquist wavenumber in Fourier space.

        Parameters
        ----------
        boxsize : float or sequence of [float, float, float]
            Mesh box size (maximum).
        ngrid : int or sequence of [int, int, int]
            Mesh grid number (minimum).

        """
        boxsize = np.asarray(boxsize).max()
        ngrid = np.asarray(ngrid).min()

        self.thisptr.set_bins(boxsize, ngrid)

        self.bin_min, self.bin_max = self.thisptr.bin_min, self.thisptr.bin_max
        self.num_bins = self.thisptr.num_bins
        self.bin_edges = self.thisptr.bin_edges
        self.bin_centres = self.thisptr.bin_centres
        self.bin_widths = self.thisptr.bin_widths
        self.scheme = self.thisptr.scheme

    def set_custom_bins(self, bin_edges):
        """Set customised binning using bin edges.

        The binning :attr:`scheme` is set to 'custom'.
        :attr:`bin_centres` and :attr:`bin_widths` are automatically
        calculated.

        Parameters
        ----------
        bin_edges : 1-d array of float
            Bin edges of length (:attr:`num_bins` + 1).

        """
        self.bin_min, self.bin_max = bin_edges[0], bin_edges[-1]
        self.num_bins = len(bin_edges) - 1

        self.bin_edges = np.array(bin_edges)
        self.bin_centres = np.add(bin_edges[:-1], bin_edges[1:]) / 2.
        self.bin_widths = np.subtract(bin_edges[1:], bin_edges[:-1])

        self.scheme = 'custom'

        self.thisptr.bin_min = self.bin_min
        self.thisptr.bin_max = self.bin_max
        self.thisptr.num_bins = self.num_bins
        self.thisptr.bin_edges = self.bin_edges
        self.thisptr.bin_centres = self.bin_centres
        self.thisptr.bin_widths = self.bin_widths
        self.thisptr.scheme = self.scheme.encode('utf-8')


# NOTE: Remove docstrings from Cython pseudo special functions.
vars()['__reduce_cython__'].__doc__ = None
vars()['__setstate_cython__'].__doc__ = None
