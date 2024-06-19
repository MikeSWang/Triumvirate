"""
Three-Point Correlator Algorithms (:mod:`~triumvirate._threept`)
==========================================================================

Declare and interface with three-point correlator algorithms.

"""
# STYLE: Un-Pythonic import order to prioritise Cython imports.
from cython.operator cimport dereference as deref
from libc.stdlib cimport free, malloc
from libcpp cimport bool as bool_t

import numpy as np
cimport numpy as np

from ._particles cimport CppParticleCatalogue, _ParticleCatalogue
from .dataobjs cimport (
    Binning, CppBinning,
    LineOfSight,
    BispecMeasurements, ThreePCFMeasurements, ThreePCFWindowMeasurements
)
from .parameters cimport CppParameterSet, ParameterSet

np.import_array()


cdef extern from "include/threept.hpp":
    # --------------------------------------------------------------------
    # Normalisation
    # --------------------------------------------------------------------

    double calc_bispec_normalisation_from_particles_cpp \
        "trv::calc_bispec_normalisation_from_particles" (
            CppParticleCatalogue& particles,
            double alpha
        ) except +

    double calc_bispec_normalisation_from_mesh_cpp \
        "trv::calc_bispec_normalisation_from_mesh" (
            CppParticleCatalogue& particles,
            CppParameterSet& params,
            double alpha
        )


    # --------------------------------------------------------------------
    # Full correlator
    # --------------------------------------------------------------------

    BispecMeasurements compute_bispec_cpp "trv::compute_bispec" (
        CppParticleCatalogue& catalogue_data,
        CppParticleCatalogue& catalogue_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        CppBinning& kbinning,
        double norm_factor
    )

    ThreePCFMeasurements compute_3pcf_cpp "trv::compute_3pcf" (
        CppParticleCatalogue& catalogue_data,
        CppParticleCatalogue& catalogue_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        CppBinning& rbinning,
        double norm_factor
    )

    BispecMeasurements compute_bispec_in_gpp_box_cpp \
        "trv::compute_bispec_in_gpp_box" (
            CppParticleCatalogue& catalogue_data,
            CppParameterSet& params,
            CppBinning& kbinning,
            double norm_factor
        )

    ThreePCFMeasurements compute_3pcf_in_gpp_box_cpp \
        "trv::compute_3pcf_in_gpp_box" (
            CppParticleCatalogue& catalogue_data,
            CppParameterSet& params,
            CppBinning& rbinning,
            double norm_factor
        )

    ThreePCFWindowMeasurements compute_3pcf_window_cpp \
        "trv::compute_3pcf_window" (
            CppParticleCatalogue& catalogue_rand,
            LineOfSight* los_rand,
            CppParameterSet& params,
            CppBinning& rbinning,
            double alpha,
            double norm_factor,
            bool_t wide_angle
        )

    # BispecMeasurements compute_bispec_for_los_choice_cpp \
    #     "trv::compute_bispec_for_los_choice" (
    #         CppParticleCatalogue& catalogue_data,
    #         CppParticleCatalogue& catalogue_rand,
    #         LineOfSight* los_data,
    #         LineOfSight* los_rand,
    #         int los_choice,
    #         CppParameterSet& params,
    #         CppBinning& kbinning,
    #         double norm_factor
    #     )


def _calc_bispec_normalisation_from_particles(
        _ParticleCatalogue particles not None, double alpha
    ):
    return calc_bispec_normalisation_from_particles_cpp(
        deref(particles.thisptr), alpha
    )


def _calc_bispec_normalisation_from_mesh(
        _ParticleCatalogue particles not None,
        ParameterSet params not None,
        double alpha
    ):
    return calc_bispec_normalisation_from_mesh_cpp(
        deref(particles.thisptr), deref(params.thisptr), alpha
    )


def _compute_bispec(
        _ParticleCatalogue catalogue_data not None,
        _ParticleCatalogue catalogue_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_data not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        Binning kbinning not None,
        double norm_factor
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
    cdef BispecMeasurements results
    results = compute_bispec_cpp(
        deref(catalogue_data.thisptr), deref(catalogue_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr), deref(kbinning.thisptr),
        norm_factor
    )

    free(los_data_cpp); free(los_rand_cpp)

    return {
        'k1_bin': np.asarray(results.k1_bin),
        'k2_bin': np.asarray(results.k2_bin),
        'k1_eff': np.asarray(results.k1_eff),
        'k2_eff': np.asarray(results.k2_eff),
        'nmodes_1': np.asarray(results.nmodes_1),
        'nmodes_2': np.asarray(results.nmodes_2),
        'bk_raw': np.asarray(results.bk_raw),
        'bk_shot': np.asarray(results.bk_shot),
    }


def _compute_3pcf(
        _ParticleCatalogue catalogue_data not None,
        _ParticleCatalogue catalogue_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_data not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        Binning rbinning not None,
        double norm_factor
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
    cdef ThreePCFMeasurements results
    results = compute_3pcf_cpp(
        deref(catalogue_data.thisptr), deref(catalogue_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr), deref(rbinning.thisptr),
        norm_factor
    )

    free(los_data_cpp); free(los_rand_cpp)

    return {
        'r1_bin': np.asarray(results.r1_bin),
        'r2_bin': np.asarray(results.r2_bin),
        'r1_eff': np.asarray(results.r1_eff),
        'r2_eff': np.asarray(results.r2_eff),
        'npairs_1': np.asarray(results.npairs_1),
        'npairs_2': np.asarray(results.npairs_2),
        'zeta_raw': np.asarray(results.zeta_raw),
        'zeta_shot': np.asarray(results.zeta_shot),
    }


def _compute_bispec_in_gpp_box(
        _ParticleCatalogue catalogue_data not None,
        ParameterSet params not None,
        Binning kbinning not None,
        double norm_factor
    ):
    cdef BispecMeasurements results
    results = compute_bispec_in_gpp_box_cpp(
        deref(catalogue_data.thisptr),
        deref(params.thisptr), deref(kbinning.thisptr),
        norm_factor
    )

    return {
        'k1_bin': np.asarray(results.k1_bin),
        'k2_bin': np.asarray(results.k2_bin),
        'k1_eff': np.asarray(results.k1_eff),
        'k2_eff': np.asarray(results.k2_eff),
        'nmodes_1': np.asarray(results.nmodes_1),
        'nmodes_2': np.asarray(results.nmodes_2),
        'bk_raw': np.asarray(results.bk_raw),
        'bk_shot': np.asarray(results.bk_shot),
    }


def _compute_3pcf_in_gpp_box(
        _ParticleCatalogue catalogue_data not None,
        ParameterSet params not None,
        Binning rbinning not None,
        double norm_factor
    ):
    cdef ThreePCFMeasurements results
    results = compute_3pcf_in_gpp_box_cpp(
        deref(catalogue_data.thisptr),
        deref(params.thisptr), deref(rbinning.thisptr),
        norm_factor
    )

    return {
        'r1_bin': np.asarray(results.r1_bin),
        'r2_bin': np.asarray(results.r2_bin),
        'r1_eff': np.asarray(results.r1_eff),
        'r2_eff': np.asarray(results.r2_eff),
        'npairs_1': np.asarray(results.npairs_1),
        'npairs_2': np.asarray(results.npairs_2),
        'zeta_raw': np.asarray(results.zeta_raw),
        'zeta_shot': np.asarray(results.zeta_shot),
    }


def _compute_3pcf_window(
        _ParticleCatalogue catalogue_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        Binning rbinning not None,
        double alpha,
        double norm_factor,
        bool_t wide_angle
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
    cdef ThreePCFWindowMeasurements results
    results = compute_3pcf_window_cpp(
        deref(catalogue_rand.thisptr), los_rand_cpp,
        deref(params.thisptr), deref(rbinning.thisptr),
        alpha, norm_factor,
        wide_angle
    )

    free(los_rand_cpp)

    return {
        'r1_bin': np.asarray(results.r1_bin),
        'r2_bin': np.asarray(results.r2_bin),
        'r1_eff': np.asarray(results.r1_eff),
        'r2_eff': np.asarray(results.r2_eff),
        'npairs_1': np.asarray(results.npairs_1),
        'npairs_2': np.asarray(results.npairs_2),
        'zeta_raw': np.asarray(results.zeta_raw),
        'zeta_shot': np.asarray(results.zeta_shot),
    }


# def _compute_bispec_for_los_choice(
#         _ParticleCatalogue catalogue_data not None,
#         _ParticleCatalogue catalogue_rand not None,
#         np.ndarray[double, ndim=2, mode='c'] los_data not None,
#         np.ndarray[double, ndim=2, mode='c'] los_rand not None,
#         int los_choice,
#         ParameterSet params not None,
#         Binning kbinning not None,
#         double norm_factor
#     ):
#     # Parse lines of sight per particle.
#     cdef LineOfSight* los_data_cpp = <LineOfSight*>malloc(
#         len(los_data) * sizeof(LineOfSight)
#     )
#     for pid, (los_x, los_y, los_z) in enumerate(los_data):
#         los_data_cpp[pid].pos[0] = los_x
#         los_data_cpp[pid].pos[1] = los_y
#         los_data_cpp[pid].pos[2] = los_z

#     cdef LineOfSight* los_rand_cpp = <LineOfSight*>malloc(
#         len(los_rand) * sizeof(LineOfSight)
#     )
#     for pid, (los_x, los_y, los_z) in enumerate(los_rand):
#         los_rand_cpp[pid].pos[0] = los_x
#         los_rand_cpp[pid].pos[1] = los_y
#         los_rand_cpp[pid].pos[2] = los_z

#     # Run algorithm.
#     cdef BispecMeasurements results
#     results = compute_bispec_for_los_choice_cpp(
#         deref(catalogue_data.thisptr), deref(catalogue_rand.thisptr),
#         los_data_cpp, los_rand_cpp, los_choice,
#         deref(params.thisptr), deref(kbinning.thisptr),
#         norm_factor
#     )

#     free(los_data_cpp); free(los_rand_cpp)

#     return {
#         'k1_bin': np.asarray(results.k1_bin),
#         'k2_bin': np.asarray(results.k2_bin),
#         'k1_eff': np.asarray(results.k1_eff),
#         'k2_eff': np.asarray(results.k2_eff),
#         'nmodes_1': np.asarray(results.nmodes_1),
#         'nmodes_2': np.asarray(results.nmodes_2),
#         'bk_raw': np.asarray(results.bk_raw),
#         'bk_shot': np.asarray(results.bk_shot),
#     }
