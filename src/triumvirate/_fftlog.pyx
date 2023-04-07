"""
FFTLog algorithm (:mod:`~triumvirate._fftlog`)
==========================================================================

Implement the FFTLog algorithm for Hankel-related transforms.

"""
import numpy as np
cimport numpy as np


cdef extern from "include/fftlog.hpp":
    int sj_transform_cpp "trv::maths::sj_transform" (
        int ell, int m, int N, double* r, double* a, double* k, double* b
    )

    int sj_transform_symm_biased_cpp "trv::maths::sj_transform_symm_biased" (
        int ell, int i, int N, double* r, double* a, double* k, double* b
    )

    int transform_powspec_to_corrfunc_multipole_cpp \
        "trv::transform_powspec_to_corrfunc_multipole" (
            int ell, int N, double* k, double* pk, double* r, double* xi
        )

    int transform_corrfunc_to_powspec_multipole_cpp \
        "trv::transform_corrfunc_to_powspec_multipole" (
            int ell, int N, double* r, double* xi, double* k, double* pk
        )


def sj_transform(
        int ell, int m, int N,
        np.ndarray[double, ndim=1, mode='c'] x not None,
        np.ndarray[double, ndim=1, mode='c'] fx not None
    ):
    cdef np.ndarray[double, ndim=1, mode='c'] y = np.zeros(N)
    cdef np.ndarray[double, ndim=1, mode='c'] gy = np.zeros(N)

    sj_transform_cpp(ell, m, N, &x[0], &fx[0], &y[0], &gy[0])

    return y, gy


def sj_transform_symm_biased(
        int ell, int i, int N,
        np.ndarray[double, ndim=1, mode='c'] x not None,
        np.ndarray[double, ndim=1, mode='c'] fx not None
    ):
    cdef np.ndarray[double, ndim=1, mode='c'] y = np.zeros(N)
    cdef np.ndarray[double, ndim=1, mode='c'] gy = np.zeros(N)

    sj_transform_symm_biased_cpp(ell, i, N, &x[0], &fx[0], &y[0], &gy[0])

    return y, gy


def trans_powspec_to_corrfunc_multipole(
        int ell, int N,
        np.ndarray[double, ndim=1, mode='c'] k not None,
        np.ndarray[double, ndim=1, mode='c'] pk not None
    ):
    cdef np.ndarray[double, ndim=1, mode='c'] r = np.zeros(N)
    cdef np.ndarray[double, ndim=1, mode='c'] xi = np.zeros(N)

    transform_powspec_to_corrfunc_multipole_cpp(
        ell, N, &k[0], &pk[0], &r[0], &xi[0]
    )

    return r, xi


def trans_corrfunc_to_powspec_multipole(
        int ell, int N,
        np.ndarray[double, ndim=1, mode='c'] r not None,
        np.ndarray[double, ndim=1, mode='c'] xi not None
    ):
    cdef np.ndarray[double, ndim=1, mode='c'] k = np.zeros(N)
    cdef np.ndarray[double, ndim=1, mode='c'] pk = np.zeros(N)

    transform_corrfunc_to_powspec_multipole_cpp(
        ell, N, &r[0], &xi[0], &k[0], &pk[0]
    )

    return k, pk
