"""Declaration of data objects.

"""
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport numpy as np

# from parameters cimport CppParameterSet


cdef extern from "include/dataobjs.hpp":
    # --------------------------------------------------------------------
    # Binning schemes
    # --------------------------------------------------------------------

    cdef cppclass CppBinning "trv::Binning":
        string scheme
        string space
        double bin_min
        double bin_max
        int num_bins
        vector[double] bin_edges
        vector[double] bin_centres
        vector[double] bin_widths

        CppBinning(double coord_min, double coord_max, int nbin) except +
        # CppBinning(CppParameterSet& params)

        void set_bins(string scheme, string space) except +
        void set_bins() except +


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
        vector[double] kbin
        vector[double] keff
        vector[int] nmodes
        vector[np.complex128_t] pk_raw
        vector[np.complex128_t] pk_shot

    struct TwoPCFMeasurements "trv::TwoPCFMeasurements":
        vector[double] rbin
        vector[double] reff
        vector[int] npairs;
        vector[np.complex128_t] xi

    struct TwoPCFWindowMeasurements "trv::TwoPCFWindowMeasurements":
        vector[double] rbin
        vector[double] reff
        vector[int] npairs;
        vector[np.complex128_t] xi

    # -- Three-point statistics ------------------------------------------

    struct BispecMeasurements "trv::BispecMeasurements":
        vector[double] k1bin
        vector[double] k2bin
        vector[double] k1eff
        vector[double] k2eff
        vector[int] nmodes
        vector[np.complex128_t] bk_raw
        vector[np.complex128_t] bk_shot

    struct ThreePCFMeasurements "trv::ThreePCFMeasurements":
        vector[double] r1bin
        vector[double] r2bin
        vector[double] r1eff
        vector[double] r2eff
        vector[int] npairs
        vector[np.complex128_t] zeta_raw
        vector[np.complex128_t] zeta_shot

    struct ThreePCFWindowMeasurements "trv::ThreePCFWindowMeasurements":
        vector[double] r1bin
        vector[double] r2bin
        vector[double] r1eff
        vector[double] r2eff
        vector[int] npairs
        vector[np.complex128_t] zeta_raw
        vector[np.complex128_t] zeta_shot


cdef class Binning:
    cdef CppBinning* thisptr
    cdef public str scheme
    cdef public str space
    cdef public double bin_min
    cdef public double bin_max
    cdef public int num_bins
    cdef public vector[double] bin_edges
    cdef public vector[double] bin_centres
    cdef public vector[double] bin_widths
