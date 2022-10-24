"""Declaration of :cpp:class:`trv::ParameterSet` and
its members and methods.

"""
from libcpp cimport bool as bool_t
from libcpp.string cimport string


cdef extern from "include/parameters.hpp":
    cdef cppclass CppParameterSet "trv::ParameterSet":
        # ----------------------------------------------------------------
        # Members
        # ----------------------------------------------------------------

        # -- I/O ---------------------------------------------------------

        string catalogue_dir
        string measurement_dir
        string data_catalogue_file
        string rand_catalogue_file
        # string catalogue_columns
        string output_tag

        # -- Mesh sampling -----------------------------------------------

        double boxsize[3]
        int ngrid[3]

        double volume
        int nmesh

        string alignment
        string padscale
        double padfactor

        string assignment
        string interlace

        # -- Measurement -------------------------------------------------

        string catalogue_type
        string statistic_type

        string norm_convention

        string binning
        string form

        string npoint
        string space

        int ell1
        int ell2
        int ELL

        int i_wa
        int j_wa

        double bin_min
        double bin_max

        int num_bins
        int idx_bin

        # -- Misc --------------------------------------------------------

        int verbose

        # ----------------------------------------------------------------
        # Methods
        # ----------------------------------------------------------------

        int validate()

        # int print_to_file(char* out_parameter_filepath)
        # int print_to_file()


cdef class ParameterSet:
    cdef CppParameterSet* thisptr
    cdef object _logger
    cdef string _source
    cdef bool_t _original
    cdef bool_t _validity
    cdef dict _params
