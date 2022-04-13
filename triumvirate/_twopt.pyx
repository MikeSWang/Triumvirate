"""
Two-Point Correlators (:mod:`~triumvirate._twopt`)
===========================================================================

Compute two-point correlator statistics.

"""
from cython.operator cimport dereference as deref
from libc.stdlib cimport malloc
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

from parameters cimport CppParameterSet, ParameterSet
from _catalogue cimport CppParticleCatalogue, _ParticleCatalogue


cdef extern from "include/particles.hpp":
    struct LineOfSight:
        double pos[3]


cdef extern from "include/twopt.hpp":
    struct PowspecMeasurements:
        vector[double] kbin
        vector[double] keff
        vector[int] nmode
        vector[np.complex128_t] pk_raw
        vector[np.complex128_t] pk_shot
    struct CorrfuncMeasurements:
        vector[double] rbin
        vector[double] reff
        vector[int] npair;
        vector[np.complex128_t] xi
    struct PowspecWindowMeasurements:
        vector[double] kbin
        vector[double] keff
        vector[int] nmode;
        vector[np.complex128_t] pk
    struct CorrfuncWindowMeasurements:
        vector[double] rbin
        vector[double] reff
        vector[int] npair;
        vector[np.complex128_t] xi

    double calc_powspec_normalisation_from_mesh_cpp \
        "calc_powspec_normalisation_from_mesh" (
            CppParticleCatalogue& catalogue,
            CppParameterSet& params,
            double alpha
        )
    double calc_powspec_normalisation_from_particles_cpp \
        "calc_powspec_normalisation_from_particles" (
            CppParticleCatalogue& catalogue,
            double alpha
        )

    PowspecMeasurements compute_powspec_cpp "compute_powspec" (
        CppParticleCatalogue& particles_data,
        CppParticleCatalogue& particles_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        double* kbin,
        double alpha,
        double norm,
        bool_t save
    )
    CorrfuncMeasurements compute_corrfunc_cpp "compute_corrfunc" (
        CppParticleCatalogue& particles_data,
        CppParticleCatalogue& particles_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        double* rbin,
        double alpha,
        double norm,
        bool_t save
    )
    # PowspecWindowMeasurements compute_powspec_window_cpp \
    #     "compute_powspec_window" (
    #         CppParticleCatalogue& particles_rand,
    #         LineOfSight* los_rand,
    #         CppParameterSet& params,
    #         double* kbin,
    #         double alpha,
    #         double norm,
    #         bool_t save
    #     )
    CorrfuncWindowMeasurements compute_corrfunc_window_cpp \
        "compute_corrfunc_window" (
            CppParticleCatalogue& particles_rand,
            LineOfSight* los_rand,
            CppParameterSet& params,
            double* rbin,
            double alpha,
            double norm,
            bool_t save
        )
    PowspecMeasurements compute_powspec_in_box_cpp "compute_powspec_in_box" (
        CppParticleCatalogue& particles_data,
        CppParameterSet& params,
        double* kbin,
        bool_t save
    )
    CorrfuncMeasurements compute_corrfunc_in_box_cpp "compute_corrfunc_in_box" (
        CppParticleCatalogue& particles_data,
        CppParameterSet& params,
        double* rbin,
        bool_t save
    )


def _calc_powspec_normalisation_from_mesh(
        _ParticleCatalogue particles not None,
        ParameterSet params not None, double alpha
    ):
    return calc_powspec_normalisation_from_mesh_cpp(
        deref(particles.thisptr), deref(params.thisptr), alpha
    )

def _calc_powspec_normalisation_from_particles(
        _ParticleCatalogue particles not None, double alpha
    ):
    return calc_powspec_normalisation_from_particles_cpp(
        deref(particles.thisptr), alpha
    )

def _compute_powspec(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_data not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] kbin not None,
        double alpha,
        double norm,
        bool_t save
    ):

    # Parse lines of sight per particle.
    cdef LineOfSight* los_data_cpp = <LineOfSight*>malloc(
        len(los_data) * sizeof(LineOfSight)
    )
    for pid, (los_x, los_y, los_z) in enumerate(los_data):
        los_data_cpp[pid].pos[0] = los_x
        los_data_cpp[pid].pos[1] = los_y
        los_data_cpp[pid].pos[2] = los_z

    cdef LineOfSight* los_rand_cpp = <LineOfSight*>malloc(
        len(los_rand) * sizeof(LineOfSight)
    )
    for pid, (los_x, los_y, los_z) in enumerate(los_rand):
        los_rand_cpp[pid].pos[0] = los_x
        los_rand_cpp[pid].pos[1] = los_y
        los_rand_cpp[pid].pos[2] = los_z

    # Run algorithm.
    cdef PowspecMeasurements meas
    meas = compute_powspec_cpp(
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr),
        &kbin[0], alpha, norm,
        save
    )

    return {
        'kbin': np.asarray(meas.kbin),
        'keff': np.asarray(meas.keff),
        'nmode': np.asarray(meas.nmode),
        'pk_raw': np.asarray(meas.pk_raw),
        'pk_shot': np.asarray(meas.pk_shot),
    }

def _compute_corrfunc(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_data not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] rbin not None,
        double alpha,
        double norm,
        bool_t save
    ):

    # Parse lines of sight per particle.
    cdef LineOfSight* los_data_cpp = <LineOfSight*>malloc(
        len(los_data) * sizeof(LineOfSight)
    )
    for pid, (los_x, los_y, los_z) in enumerate(los_data):
        los_data_cpp[pid].pos[0] = los_x
        los_data_cpp[pid].pos[1] = los_y
        los_data_cpp[pid].pos[2] = los_z

    cdef LineOfSight* los_rand_cpp = <LineOfSight*>malloc(
        len(los_rand) * sizeof(LineOfSight)
    )
    for pid, (los_x, los_y, los_z) in enumerate(los_rand):
        los_rand_cpp[pid].pos[0] = los_x
        los_rand_cpp[pid].pos[1] = los_y
        los_rand_cpp[pid].pos[2] = los_z

    # Run algorithm.
    cdef CorrfuncMeasurements meas
    meas = compute_corrfunc_cpp(
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr),
        &rbin[0], alpha, norm,
        save
    )

    return {
        'rbin': np.asarray(meas.rbin),
        'reff': np.asarray(meas.reff),
        'npair': np.asarray(meas.npair),
        'xi': np.asarray(meas.xi),
    }

# def _compute_powspec_window(
#         _ParticleCatalogue particles_rand not None,
#         np.ndarray[double, ndim=2, mode='c'] los_rand not None,
#         ParameterSet params not None,
#         np.ndarray[double, ndim=1, mode='c'] kbin not None,
#         double alpha,
#         double norm,
#         bool_t save
#     ):
#
#     # Parse lines of sight per particle.
#     cdef LineOfSight* los_rand_cpp = <LineOfSight*>malloc(
#         len(los_rand) * sizeof(LineOfSight)
#     )
#     for pid, (los_x, los_y, los_z) in enumerate(los_rand):
#         los_rand_cpp[pid].pos[0] = los_x
#         los_rand_cpp[pid].pos[1] = los_y
#         los_rand_cpp[pid].pos[2] = los_z
#
#     # Run algorithm.
#     cdef PowspecWindowMeasurements meas
#     meas = compute_powspec_window_cpp(
#         deref(particles_rand.thisptr), los_rand_cpp, deref(params.thisptr),
#         &kbin[0], alpha, norm,
#         save
#     )
#
#     return {
#         'kbin': np.asarray(meas.kbin),
#         'keff': np.asarray(meas.keff),
#         'nmode': np.asarray(meas.nmode),
#         'pk': np.asarray(meas.pk_raw),
#     }

def _compute_corrfunc_window(
        _ParticleCatalogue particles_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] rbin not None,
        double alpha,
        double norm,
        bool_t save
    ):

    # Parse lines of sight per particle.
    cdef LineOfSight* los_rand_cpp = <LineOfSight*>malloc(
        len(los_rand) * sizeof(LineOfSight)
    )
    for pid, (los_x, los_y, los_z) in enumerate(los_rand):
        los_rand_cpp[pid].pos[0] = los_x
        los_rand_cpp[pid].pos[1] = los_y
        los_rand_cpp[pid].pos[2] = los_z

    # Run algorithm.
    cdef CorrfuncWindowMeasurements meas
    meas = compute_corrfunc_window_cpp(
        deref(particles_rand.thisptr), los_rand_cpp, deref(params.thisptr),
        &rbin[0], alpha, norm,
        save
    )

    return {
        'rbin': np.asarray(meas.rbin),
        'reff': np.asarray(meas.reff),
        'npair': np.asarray(meas.npair),
        'xi': np.asarray(meas.xi),
    }

def _compute_powspec_in_box(
        _ParticleCatalogue particles_data not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] kbin not None,
        double norm,
        bool_t save
    ):

    cdef PowspecMeasurements meas
    meas = compute_powspec_in_box_cpp(
        deref(particles_data.thisptr), deref(params.thisptr),
        &kbin[0], norm,
        save
    )

    return {
        'kbin': np.asarray(meas.kbin),
        'keff': np.asarray(meas.keff),
        'nmode': np.asarray(meas.nmode),
        'pk_raw': np.asarray(meas.pk_raw),
        'pk_shot': np.asarray(meas.pk_shot),
    }

def _compute_corrfunc_in_box(
        _ParticleCatalogue particles_data not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] rbin not None,
        double norm,
        bool_t save
    ):

    cdef CorrfuncMeasurements meas
    meas = compute_corrfunc_in_box_cpp(
        deref(particles_data.thisptr), deref(params.thisptr),
        &rbin[0], norm,
        save
    )

    return {
        'rbin': np.asarray(meas.rbin),
        'reff': np.asarray(meas.reff),
        'npair': np.asarray(meas.npair),
        'xi': np.asarray(meas.xi),
    }
