"""Interface with the FFTLog algoritms for the Hankel and
associated transforms.

"""
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector


cdef extern from "include/fftlog.hpp":
    cdef cppclass CppHankelTransform "trv::maths::HankelTransform":
        double order
        double bias
        int nsamp
        int nsamp_trans
        double logres
        double pivot

        vector[double] pre_sampts
        vector[double] post_sampts

        CppHankelTransform(double mu, double q, bool_t threaded)

        void initialise(
            vector[double] sample_pts, double kr_c, bool_t lowring,
            int extrap, double extrap_exp
        ) except +

        void biased_transform(
            double complex* a, double complex* b
        ) except +


cdef class HankelTransform:
    cdef CppHankelTransform* thisptr
    cdef public double order
    cdef public double bias
    cdef public int _nsamp
    cdef public int _nsamp_trans
    cdef public double _logres
    cdef public double _pivot
    cdef public vector[double] _pre_sampts
    cdef public vector[double] _post_sampts
