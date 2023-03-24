"""
Three-Point Correlator Algorithms (:mod:`~triumvirate._threept`)
==========================================================================

Declaration and wrapping of three-point correlator algorithms.

"""
from cython.operator cimport dereference as deref
from libc.stdlib cimport free, malloc
from libcpp cimport bool as bool_t

import numpy as np
cimport numpy as np

from triumvirate._particles cimport CppParticleCatalogue, _ParticleCatalogue
from triumvirate.dataobjs cimport (
    Binning, CppBinning,
    LineOfSight,
    BispecMeasurements, ThreePCFMeasurements, ThreePCFWindowMeasurements
)
from triumvirate.parameters cimport CppParameterSet, ParameterSet


cdef extern from "include/threept.hpp":
    # --------------------------------------------------------------------
    # Normalisation
    # --------------------------------------------------------------------

    double calc_bispec_normalisation_from_mesh_cpp \
        "trv::calc_bispec_normalisation_from_mesh" (
            CppParticleCatalogue& catalogue,
            CppParameterSet& params,
            double alpha
        )

    double calc_bispec_normalisation_from_particles_cpp \
        "trv::calc_bispec_normalisation_from_particles" (
            CppParticleCatalogue& catalogue,
            double alpha
        ) except +


    # --------------------------------------------------------------------
    # Full correlator
    # --------------------------------------------------------------------

    BispecMeasurements compute_bispec_cpp "trv::compute_bispec" (
        CppParticleCatalogue& particles_data,
        CppParticleCatalogue& particles_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        CppBinning& kbinning,
        double norm_factor
    )

    ThreePCFMeasurements compute_3pcf_cpp "trv::compute_3pcf" (
        CppParticleCatalogue& particles_data,
        CppParticleCatalogue& particles_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        CppBinning& rbinning,
        double norm_factor
    )

    BispecMeasurements compute_bispec_in_gpp_box_cpp \
        "trv::compute_bispec_in_gpp_box" (
            CppParticleCatalogue& particles_data,
            CppParameterSet& params,
            CppBinning& kbinning,
            double norm_factor
        )

    ThreePCFMeasurements compute_3pcf_in_gpp_box_cpp \
        "trv::compute_3pcf_in_gpp_box" (
            CppParticleCatalogue& particles_data,
            CppParameterSet& params,
            CppBinning& rbinning,
            double norm_factor
        )

    ThreePCFWindowMeasurements compute_3pcf_window_cpp \
        "trv::compute_3pcf_window" (
            CppParticleCatalogue& particles_rand,
            LineOfSight* los_rand,
            CppParameterSet& params,
            CppBinning& rbinning,
            double alpha,
            double norm_factor,
            bool_t wide_angle
        )

    # BispecMeasurements compute_bispec_for_los_choice_cpp \
    #     "trv::compute_bispec_for_los_choice" (
    #         CppParticleCatalogue& particles_data,
    #         CppParticleCatalogue& particles_rand,
    #         LineOfSight* los_data,
    #         LineOfSight* los_rand,
    #         int los_choice,
    #         CppParameterSet& params,
    #         CppBinning& kbinning,
    #         double norm_factor
    #     )


def _calc_bispec_normalisation_from_mesh(
        _ParticleCatalogue catalogue not None,
        ParameterSet params not None,
        double alpha
    ):
    return calc_bispec_normalisation_from_mesh_cpp(
        deref(catalogue.thisptr), deref(params.thisptr), alpha
    )


def _calc_bispec_normalisation_from_particles(
        _ParticleCatalogue catalogue not None, double alpha
    ):
    return calc_bispec_normalisation_from_particles_cpp(
        deref(catalogue.thisptr), alpha
    )


def _compute_bispec(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
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
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr), deref(kbinning.thisptr),
        norm_factor
    )

    free(los_data_cpp); free(los_rand_cpp)

    return {
        'k1_bin': np.array(results.k1_bin),
        'k2_bin': np.array(results.k2_bin),
        'k1_eff': np.array(results.k1_eff),
        'k2_eff': np.array(results.k2_eff),
        'nmodes': np.array(results.nmodes),
        'bk_raw': np.array(results.bk_raw),
        'bk_shot': np.array(results.bk_shot),
    }


def _compute_3pcf(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
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
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr), deref(rbinning.thisptr),
        norm_factor
    )

    free(los_data_cpp); free(los_rand_cpp)

    return {
        'r1_bin': np.array(results.r1_bin),
        'r2_bin': np.array(results.r2_bin),
        'r1_eff': np.array(results.r1_eff),
        'r2_eff': np.array(results.r2_eff),
        'npairs': np.array(results.npairs),
        'zeta_raw': np.array(results.zeta_raw),
        'zeta_shot': np.array(results.zeta_shot),
    }


def _compute_bispec_in_gpp_box(
        _ParticleCatalogue particles_data not None,
        ParameterSet params not None,
        Binning kbinning not None,
        double norm_factor
    ):
    cdef BispecMeasurements results
    results = compute_bispec_in_gpp_box_cpp(
        deref(particles_data.thisptr),
        deref(params.thisptr), deref(kbinning.thisptr),
        norm_factor
    )

    return {
        'k1_bin': np.array(results.k1_bin),
        'k2_bin': np.array(results.k2_bin),
        'k1_eff': np.array(results.k1_eff),
        'k2_eff': np.array(results.k2_eff),
        'nmodes': np.array(results.nmodes),
        'bk_raw': np.array(results.bk_raw),
        'bk_shot': np.array(results.bk_shot),
    }


def _compute_3pcf_in_gpp_box(
        _ParticleCatalogue particles_data not None,
        ParameterSet params not None,
        Binning rbinning not None,
        double norm_factor
    ):
    cdef ThreePCFMeasurements results
    results = compute_3pcf_in_gpp_box_cpp(
        deref(particles_data.thisptr),
        deref(params.thisptr), deref(rbinning.thisptr),
        norm_factor
    )

    return {
        'r1_bin': np.array(results.r1_bin),
        'r2_bin': np.array(results.r2_bin),
        'r1_eff': np.array(results.r1_eff),
        'r2_eff': np.array(results.r2_eff),
        'npairs': np.array(results.npairs),
        'zeta_raw': np.array(results.zeta_raw),
        'zeta_shot': np.array(results.zeta_shot),
    }


def _compute_3pcf_window(
        _ParticleCatalogue particles_rand not None,
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
        deref(particles_rand.thisptr), los_rand_cpp,
        deref(params.thisptr), deref(rbinning.thisptr),
        alpha, norm_factor,
        wide_angle
    )

    free(los_rand_cpp)

    return {
        'r1_bin': np.array(results.r1_bin),
        'r2_bin': np.array(results.r2_bin),
        'r1_eff': np.array(results.r1_eff),
        'r2_eff': np.array(results.r2_eff),
        'npairs': np.array(results.npairs),
        'zeta_raw': np.array(results.zeta_raw),
        'zeta_shot': np.array(results.zeta_shot),
    }


# def _compute_bispec_for_los_choice(
#         _ParticleCatalogue particles_data not None,
#         _ParticleCatalogue particles_rand not None,
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
#         deref(particles_data.thisptr), deref(particles_rand.thisptr),
#         los_data_cpp, los_rand_cpp, los_choice,
#         deref(params.thisptr), deref(kbinning.thisptr),
#         norm_factor
#     )

#     free(los_data_cpp); free(los_rand_cpp)

#     return {
#         'k1_bin': np.array(results.k1_bin),
#         'k2_bin': np.array(results.k2_bin),
#         'k1_eff': np.array(results.k1_eff),
#         'k2_eff': np.array(results.k2_eff),
#         'nmodes': np.array(results.nmodes),
#         'bk_raw': np.array(results.bk_raw),
#         'bk_shot': np.array(results.bk_shot),
#     }
