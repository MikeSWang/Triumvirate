"""
Data Objects (:mod:`~triumvirate.dataobjs`)
===========================================================================

Define data objects.

"""
# import numpy as np

from .dataobjs cimport CppBinning


cdef class Binning:
    """Binning.

    Parameters
    ----------
    bin_min, bin_max : float
        Minimum and maximum of the bin range.
    num_bins : int
        Number of bins.
    scheme : str, optional
        Binning scheme, {'lin', 'log', 'linpad', 'logpad'}.
    space : str, optional
        Coordinate space, {'config', 'fourier'}.

    Attributes
    ----------
    bin_min, bin_max : float
        Minimum and maximum of the bin range.
    num_bins : int
        Number of bins.
    scheme : str or None
        Binning scheme, {'lin', 'log', 'linpad', 'logpad'}.
    space : str or None
        Coordinate space, {'config', 'fourier'}.
    bin_edges : (:attr:`num_bins` + 1,) list of float
        Bin edges.
    bin_centres : (:attr:`num_bins`,) list of float
        Bin centres.
    bin_widths : (:attr:`num_bins`,) list of float
        Bin widths.

    Notes
    -----
    The bin setting method supports linear ('lin') and
    log-linear/exponential ('log') binning, with possible linear padding
    from zero using 5 bins of width 1.e-3 ('fourier' space) or 10.
    ('config' space) ('linpad' or 'logpad').

    """

    def __cinit__(self, bin_min, bin_max, num_bins, scheme=None, space=None):

        self.thisptr = new CppBinning(bin_min, bin_max, num_bins)

        self.bin_min = bin_min
        self.bin_max = bin_max
        self.num_bins = num_bins

        self.scheme = scheme
        self.space = space

    @classmethod
    def from_parameter_set(cls, parameter_set):
        """Create binning scheme from a parameter set.

        This sets :attr:`scheme` and :attr:`space`.

        Parameters
        ----------
        parameter_set : :class:`~triumvirate.parameters.ParameterSet`
            Parameter set.

        """
        self = cls(
            *parameter_set['range'], parameter_set['num_bins'],
            scheme=parameter_set['binning'], space=parameter_set['space']
        )

        return self

    def set_bins(self, scheme=None, space=None):
        """Set bin quantities including bin edges, centres and widths.

        Parameters
        ----------
        scheme : {'lin', 'log', 'linpad', 'logpad'}, optional
            Binning scheme.  If `None` (default), :attr:`scheme` is used
            which must not also be `None`.
        space : {'config', 'fourier'}, optional
            Coordinate space.  If `None` (default), :attr:`space` is used
            which must not also be `None`.

        Raises
        ------
        ValueError
            If `scheme` (or `space`) is `None` when :attr:`scheme`
            (or :attr:`space`) is also `None`.

        """
        scheme = scheme if scheme else self.scheme
        space = space if space else self.space

        if scheme is None:
            raise ValueError(
              "`scheme` cannot be None as the corresponding attribute is unset."
            )
        if space is None:
            raise ValueError(
              "`space` cannot be None as the corresponding attribute is unset."
            )

        self.thisptr.set_bins(scheme.encode('utf-8'), space.encode('utf-8'))

        self.bin_edges = self.thisptr.bin_edges
        self.bin_centres = self.thisptr.bin_centres
        self.bin_widths = self.thisptr.bin_widths

    # @property
    # def bin_edges(self):
    #     """Return :attr:`bin_edges`.
    #
    #     """
    #     return self.bin_edges

    # @bin_edges.setter
    # def bin_edges(self, bin_edges):
    #     """Set customised :attr:`bin_edges`.
    #
    #     :attr:`bin_centres` and :attr:`bin_widths` are also set
    #     for internal consistency.
    #
    #     Parameters
    #     ----------
    #     bin_edges : array of float
    #         Bin edges of length (:attr:`num_bins` + 1).
    #
    #     """
    #     bin_centres = np.add(bin_edges[:-1], bin_edges[1:]) / 2.
    #     bin_widths = np.subtract(bin_edges[1:], bin_edges[:-1])

    #     self.bin_edges = bin_edges
    #     self.bin_centres = bin_centres
    #     self.bin_widths = bin_widths

    #     self.thisptr.bin_edges = self.bin_edges
    #     self.thisptr.bin_centres = self.bin_centres
    #     self.thisptr.bin_widths = self.bin_widths
