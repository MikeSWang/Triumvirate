"""Interface with data objects.

"""
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport numpy as np

np.import_array()


cdef extern from "include/dataobjs.hpp":
    # --------------------------------------------------------------------
    # Binning schemes
    # --------------------------------------------------------------------

    cdef cppclass CppBinning "trv::Binning":
        string space
        string scheme
        double bin_min
        double bin_max
        int num_bins
        vector[double] bin_edges
        vector[double] bin_centres
        vector[double] bin_widths

        CppBinning(string space, string scheme) except +

        void set_bins(double bin_min, double bin_max, int nbins) except +
        void set_bins(double boxsize_max, int ngrid_min) except +


    # --------------------------------------------------------------------
    # Mesh grids
    # --------------------------------------------------------------------

    struct BinnedVectors "trv::BinnedVectors":
        int count
        int num_bins
        vector[int] indices
        vector[double] lower_edges
        vector[double] upper_edges
        vector[double] vecx
        vector[double] vecy
        vector[double] vecz


    # --------------------------------------------------------------------
    # Line of sight
    # --------------------------------------------------------------------

    struct LineOfSight "trv::LineOfSight":
        double pos[3]


    # --------------------------------------------------------------------
    # Clustering statistics
    # --------------------------------------------------------------------

    # -- Two-point statistics --------------------------------------------

    struct PowspecMeasurements "trv::PowspecMeasurements":
        int dim
        vector[double] kbin
        vector[double] keff
        vector[int] nmodes
        vector[double complex] pk_raw
        vector[double complex] pk_shot

    struct TwoPCFMeasurements "trv::TwoPCFMeasurements":
        int dim
        vector[double] rbin
        vector[double] reff
        vector[int] npairs
        vector[double complex] xi

    struct TwoPCFWindowMeasurements "trv::TwoPCFWindowMeasurements":
        int dim
        vector[double] rbin
        vector[double] reff
        vector[int] npairs
        vector[double complex] xi

    # -- Three-point statistics ------------------------------------------

    struct BispecMeasurements "trv::BispecMeasurements":
        int dim
        vector[double] k1_bin
        vector[double] k2_bin
        vector[double] k1_eff
        vector[double] k2_eff
        vector[int] nmodes_1
        vector[int] nmodes_2
        vector[double complex] bk_raw
        vector[double complex] bk_shot

    struct ThreePCFMeasurements "trv::ThreePCFMeasurements":
        int dim
        vector[double] r1_bin
        vector[double] r2_bin
        vector[double] r1_eff
        vector[double] r2_eff
        vector[int] npairs_1
        vector[int] npairs_2
        vector[double complex] zeta_raw
        vector[double complex] zeta_shot

    struct ThreePCFWindowMeasurements "trv::ThreePCFWindowMeasurements":
        int dim
        vector[double] r1_bin
        vector[double] r2_bin
        vector[double] r1_eff
        vector[double] r2_eff
        vector[int] npairs_1
        vector[int] npairs_2
        vector[double complex] zeta_raw
        vector[double complex] zeta_shot


cdef class Binning:
    cdef CppBinning* thisptr
    cdef public str space
    cdef public str scheme
    cdef public double bin_min
    cdef public double bin_max
    cdef public int num_bins
    cdef public vector[double] bin_edges
    cdef public vector[double] bin_centres
    cdef public vector[double] bin_widths
