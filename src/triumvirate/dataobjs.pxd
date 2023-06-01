"""Declare data objects.

"""
from libcpp.string cimport string
from libcpp.vector cimport vector

cimport numpy as np


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
        vector[int] npairs
        vector[np.complex128_t] xi

    struct TwoPCFWindowMeasurements "trv::TwoPCFWindowMeasurements":
        vector[double] rbin
        vector[double] reff
        vector[int] npairs
        vector[np.complex128_t] xi

    # -- Three-point statistics ------------------------------------------

    struct BispecMeasurements "trv::BispecMeasurements":
        vector[double] k1_bin
        vector[double] k2_bin
        vector[double] k1_eff
        vector[double] k2_eff
        vector[int] nmodes
        vector[np.complex128_t] bk_raw
        vector[np.complex128_t] bk_shot

    struct ThreePCFMeasurements "trv::ThreePCFMeasurements":
        vector[double] r1_bin
        vector[double] r2_bin
        vector[double] r1_eff
        vector[double] r2_eff
        vector[int] npairs
        vector[np.complex128_t] zeta_raw
        vector[np.complex128_t] zeta_shot

    struct ThreePCFWindowMeasurements "trv::ThreePCFWindowMeasurements":
        vector[double] r1_bin
        vector[double] r2_bin
        vector[double] r1_eff
        vector[double] r2_eff
        vector[int] npairs
        vector[np.complex128_t] zeta_raw
        vector[np.complex128_t] zeta_shot


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
