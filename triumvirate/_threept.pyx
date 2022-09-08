"""
Three-Point Correlators (:mod:`~triumvirate._threept`)
===========================================================================

Compute three-point correlator statistics.

"""
from cython.operator cimport dereference as deref
from libc.stdlib cimport malloc
from libcpp cimport bool as bool_t
from libcpp.vector cimport vector

import numpy as np
cimport numpy as np

from triumvirate.parameters cimport CppParameterSet, ParameterSet
from triumvirate._catalogue cimport CppParticleCatalogue, _ParticleCatalogue


cdef extern from "include/dataobjs.hpp":
    struct LineOfSight "trv::LineOfSight":
        double pos[3]


cdef extern from "include/threept.hpp":
    struct BispecMeasurements "trv::algo::BispecMeasurements":
        vector[double] kbin1
        vector[double] kbin2
        vector[double] keff1
        vector[double] keff2
        vector[int] nmode
        vector[np.complex128_t] bk_raw
        vector[np.complex128_t] bk_shot

    struct ThreePCFMeasurements "trv::algo::ThreePCFMeasurements":
        vector[double] rbin1
        vector[double] rbin2
        vector[double] reff1
        vector[double] reff2
        vector[int] npair
        vector[np.complex128_t] zeta_raw
        vector[np.complex128_t] zeta_shot

    struct ThreePCFWindowMeasurements "trv::algo::ThreePCFWindowMeasurements":
        vector[double] rbin1
        vector[double] rbin2
        vector[double] reff1
        vector[double] reff2
        vector[int] npair
        vector[np.complex128_t] zeta_raw
        vector[np.complex128_t] zeta_shot

    double calc_bispec_normalisation_from_mesh_cpp \
        "trv::algo::calc_bispec_normalisation_from_mesh" (
            CppParticleCatalogue& catalogue,
            CppParameterSet& params,
            double alpha
        )

    double calc_bispec_normalisation_from_particles_cpp \
        "trv::algo::calc_bispec_normalisation_from_particles" (
            CppParticleCatalogue& catalogue,
            double alpha
        )

    BispecMeasurements compute_bispec_cpp "trv::algo::compute_bispec" (
        CppParticleCatalogue& particles_data,
        CppParticleCatalogue& particles_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        vector[double] kbin,
        double alpha,
        double norm,
        double norm_alt,
        bool_t save
    )

    ThreePCFMeasurements compute_3pcf_cpp "trv::algo::compute_3pcf" (
        CppParticleCatalogue& particles_data,
        CppParticleCatalogue& particles_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        vector[double] rbin,
        double alpha,
        double norm,
        double norm_alt,
        bool_t save
    )

    BispecMeasurements compute_bispec_in_box_cpp \
        "trv::algo::compute_bispec_in_box" (
            CppParticleCatalogue& particles_data,
            CppParameterSet& params,
            vector[double] kbin,
            double norm,
            double norm_alt,
            bool_t save
        )

    ThreePCFMeasurements compute_3pcf_in_box_cpp \
        "trv::algo::compute_3pcf_in_box" (
            CppParticleCatalogue& particles_data,
            CppParameterSet& params,
            vector[double] rbin,
            double norm,
            double norm_alt,
            bool_t save
        )

    ThreePCFWindowMeasurements compute_3pcf_window_cpp \
        "trv::algo::compute_3pcf_window" (
            CppParticleCatalogue& particles_rand,
            LineOfSight* los_rand,
            CppParameterSet& params,
            vector[double] rbin,
            double alpha,
            double norm,
            double norm_alt,
            bool_t wide_angle,
            bool_t save
        )

    # BispecMeasurements compute_bispec_for_los_choice_cpp \
    #     "trv::algo::compute_bispec_for_los_choice" (
    #         CppParticleCatalogue& particles_data,
    #         CppParticleCatalogue& particles_rand,
    #         LineOfSight* los_data,
    #         LineOfSight* los_rand,
    #         int los_choice,
    #         CppParameterSet& params,
    #         vector[double] kbin,
    #         double alpha,
    #         double norm,
    #         double norm_alt,
    #         bool_t save
    #     )


def _calc_bispec_normalisation_from_mesh(
        _ParticleCatalogue catalogue not None,
        ParameterSet params not None,
        double alpha
    ):
    return calc_bispec_normalisation_from_mesh_cpp(
        deref(catalogue.thisptr),
        deref(params.thisptr),
        alpha
    )


def _calc_bispec_normalisation_from_particles(
        _ParticleCatalogue catalogue not None,
        double alpha
    ):
    return calc_bispec_normalisation_from_particles_cpp(
        deref(catalogue.thisptr),
        alpha
    )


def _compute_bispec(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_data not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] kbin not None,
        double alpha,
        double norm,
        double norm_alt,
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
    cdef BispecMeasurements meas
    meas = compute_bispec_cpp(
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr),
        kbin,
        alpha, norm, norm_alt,
        save
    )

    return {
        'kbin1': np.asarray(meas.kbin1),
        'kbin2': np.asarray(meas.kbin2),
        'keff1': np.asarray(meas.keff1),
        'keff2': np.asarray(meas.keff2),
        'nmode': np.asarray(meas.nmode),
        'bk_raw': np.asarray(meas.bk_raw),
        'bk_shot': np.asarray(meas.bk_shot),
    }


def _compute_3pcf(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_data not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] rbin not None,
        double alpha,
        double norm,
        double norm_alt,
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
    cdef ThreePCFMeasurements meas
    meas = compute_3pcf_cpp(
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr),
        rbin,
        alpha, norm, norm_alt,
        save
    )

    return {
        'rbin1': np.asarray(meas.rbin1),
        'rbin2': np.asarray(meas.rbin2),
        'reff1': np.asarray(meas.reff1),
        'reff2': np.asarray(meas.reff2),
        'nmode': np.asarray(meas.nmode),
        'zeta_raw': np.asarray(meas.zeta_raw),
        'zeta_shot': np.asarray(meas.zeta_shot),
    }


def _compute_bispec_in_box(
        _ParticleCatalogue particles_data not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] kbin not None,
        double norm,
        double norm_alt,
        bool_t save
    ):

    cdef BispecMeasurements meas
    meas = compute_bispec_in_box_cpp(
        deref(particles_data.thisptr), deref(params.thisptr),
        kbin,
        norm, norm_alt,
        save
    )

    return {
        'kbin1': np.asarray(meas.kbin1),
        'kbin2': np.asarray(meas.kbin2),
        'keff1': np.asarray(meas.keff1),
        'keff2': np.asarray(meas.keff2),
        'nmode': np.asarray(meas.nmode),
        'bk_raw': np.asarray(meas.bk_raw),
        'bk_shot': np.asarray(meas.bk_shot),
    }


def _compute_3pcf_in_box(
        _ParticleCatalogue particles_data not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] rbin not None,
        double norm,
        double norm_alt,
        bool_t save
    ):

    cdef ThreePCFMeasurements meas
    meas = compute_3pcf_in_box_cpp(
        deref(particles_data.thisptr), deref(params.thisptr),
        rbin,
        norm, norm_alt,
        save
    )

    return {
        'rbin1': np.asarray(meas.rbin1),
        'rbin2': np.asarray(meas.rbin2),
        'reff1': np.asarray(meas.reff1),
        'reff2': np.asarray(meas.reff2),
        'nmode': np.asarray(meas.nmode),
        'zeta_raw': np.asarray(meas.zeta_raw),
        'zeta_shot': np.asarray(meas.zeta_shot),
    }


def _compute_3pcf_window(
        _ParticleCatalogue particles_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        np.ndarray[double, ndim=1, mode='c'] rbin not None,
        double alpha,
        double norm,
        double norm_alt,
        bool_t wide_angle,
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
    cdef ThreePCFWindowMeasurements meas
    meas = compute_3pcf_window_cpp(
        deref(particles_rand.thisptr),
        los_rand_cpp,
        deref(params.thisptr),
        rbin,
        alpha, norm, norm_alt,
        wide_angle,
        save
    )

    return {
        'rbin1': np.asarray(meas.rbin1),
        'rbin2': np.asarray(meas.rbin2),
        'reff1': np.asarray(meas.reff1),
        'reff2': np.asarray(meas.reff2),
        'nmode': np.asarray(meas.nmode),
        'zeta_raw': np.asarray(meas.zeta_raw),
        'zeta_shot': np.asarray(meas.zeta_shot),
    }


# def _compute_bispec_for_los_choice(
#         _ParticleCatalogue particles_data not None,
#         _ParticleCatalogue particles_rand not None,
#         np.ndarray[double, ndim=2, mode='c'] los_data not None,
#         np.ndarray[double, ndim=2, mode='c'] los_rand not None,
#         int los_choice,
#         ParameterSet params not None,
#         np.ndarray[double, ndim=1, mode='c'] kbin not None,
#         double alpha,
#         double norm,
#         double norm_alt,
#         bool_t save
#     ):
#
#     # Parse lines of sight per particle.
#     cdef LineOfSight* los_data_cpp = <LineOfSight*>malloc(
#         len(los_data) * sizeof(LineOfSight)
#     )
#     for pid, (los_x, los_y, los_z) in enumerate(los_data):
#         los_data_cpp[pid].pos[0] = los_x
#         los_data_cpp[pid].pos[1] = los_y
#         los_data_cpp[pid].pos[2] = los_z
#
#     cdef LineOfSight* los_rand_cpp = <LineOfSight*>malloc(
#         len(los_rand) * sizeof(LineOfSight)
#     )
#     for pid, (los_x, los_y, los_z) in enumerate(los_rand):
#         los_rand_cpp[pid].pos[0] = los_x
#         los_rand_cpp[pid].pos[1] = los_y
#         los_rand_cpp[pid].pos[2] = los_z
#
#     # Run algorithm.
#     cdef BispecMeasurements meas
#     meas = compute_bispec_for_los_choice_cpp(
#         deref(particles_data.thisptr), deref(particles_rand.thisptr),
#         los_data_cpp, los_rand_cpp,
#         los_choice,
#         deref(params.thisptr),
#         kbin,
#         alpha, norm, norm_alt,
#         save
#     )
#
#     return {
#         'kbin1': np.asarray(meas.kbin1),
#         'kbin2': np.asarray(meas.kbin2),
#         'keff1': np.asarray(meas.keff1),
#         'keff2': np.asarray(meas.keff2),
#         'nmode': np.asarray(meas.nmode),
#         'bk_raw': np.asarray(meas.bk_raw),
#         'bk_shot': np.asarray(meas.bk_shot),
#     }
