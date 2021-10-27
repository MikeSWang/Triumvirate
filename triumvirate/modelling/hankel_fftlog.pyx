# cython: c_string_type=unicode, c_string_encoding=utf8
import numpy as np
cimport numpy as np

# distutils: sources=["include/hankel_fftlog.hpp"], language=['c++']
cdef extern from "include/hankel_fftlog.hpp":
    int cosmo_forward_hankel_transform(
        int ell, int m, int N, double* k, double* pk, double* r, double* xi
    )
    int wide_angle_forward_hankel_transform(
        int ell, int i, int N, double* k, double* pk, double* r, double* xi
    )
    int transform_powspec_to_corrfunc(
        int N, double* k, double* pk, double* r, double* xi
    )
    int transform_corrfunc_to_powspec(
        int N, double* r, double* xi, double* k, double* pk
    )


def hankel_transform(int ell, int m, int N,
                     np.ndarray[double, ndim=1, mode='c'] x not None,
                     np.ndarray[double, ndim=1, mode='c'] fx not None,
                     np.ndarray[double, ndim=1, mode='c'] y not None,
                     np.ndarray[double, ndim=1, mode='c'] gy not None):
    return cosmo_forward_hankel_transform(
        ell, m, N, &x[0], &fx[0], &y[0], &gy[0]
    )


def wide_angle_hankel_transform(int ell, int m, int N,
                                np.ndarray[double, ndim=1, mode='c'] x
                                    not None,
                                np.ndarray[double, ndim=1, mode='c'] fx
                                    not None,
                                np.ndarray[double, ndim=1, mode='c'] y
                                    not None,
                                np.ndarray[double, ndim=1, mode='c'] gy
                                    not None):
    return wide_angle_forward_hankel_transform(
        ell, i, N, &x[0], &fx[0], &y[0], &gy[0]
    )


def trans_powspec_to_corrfunc(int N,
                              np.ndarray[double, ndim=1, mode='c'] k not None,
                              np.ndarray[double, ndim=1, mode='c'] pk not None,
                              np.ndarray[double, ndim=1, mode='c'] r not None,
                              np.ndarray[double, ndim=1, mode='c'] xi not None):
    return transform_powspec_to_corrfunc(N, &k[0], &pk[0], &r[0], &xi[0])


def trans_corrfunc_to_powspec(int N,
                              np.ndarray[double, ndim=1, mode='c'] r not None,
                              np.ndarray[double, ndim=1, mode='c'] xi not None,
                              np.ndarray[double, ndim=1, mode='c'] k not None,
                              np.ndarray[double, ndim=1, mode='c'] pk not None):
    return transform_corrfunc_to_powspec(N, &r[0], &xi[0], &k[0], &pk[0])
