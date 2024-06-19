"""
Two-Point Correlator Algorithms (:mod:`~triumvirate._twopt`)
==========================================================================

Declare and interface with two-point correlator algorithms.

"""
# STYLE: Un-Pythonic import order to prioritise Cython imports.
from cython.operator cimport dereference as deref
from libc.stdlib cimport free, malloc
from libcpp.string cimport string

import numpy as np
cimport numpy as np

from ._particles cimport CppParticleCatalogue, _ParticleCatalogue
from .dataobjs cimport (
    Binning, CppBinning,
    LineOfSight,
    PowspecMeasurements, TwoPCFMeasurements, TwoPCFWindowMeasurements
)
from .parameters cimport CppParameterSet, ParameterSet

np.import_array()


cdef extern from "include/twopt.hpp":
    # --------------------------------------------------------------------
    # Normalisation
    # --------------------------------------------------------------------

    double calc_powspec_normalisation_from_particles_cpp \
        "trv::calc_powspec_normalisation_from_particles" (
            CppParticleCatalogue& particles,
            double alpha
        ) except +

    double calc_powspec_normalisation_from_mesh_cpp \
        "trv::calc_powspec_normalisation_from_mesh" (
            CppParticleCatalogue& particles,
            CppParameterSet& params,
            double alpha
        )

    double calc_powspec_normalisation_from_meshes_cpp \
        "trv::calc_powspec_normalisation_from_meshes" (
            CppParticleCatalogue& particles_data,
            CppParticleCatalogue& particles_rand,
            CppParameterSet& params,
            double alpha
        )

    double calc_powspec_normalisation_from_meshes_cpp \
        "trv::calc_powspec_normalisation_from_meshes" (
            CppParticleCatalogue& particles_data,
            CppParticleCatalogue& particles_rand,
            CppParameterSet& params,
            double alpha,
            double padding,
            double cellsize,
            string assignment
        )


    # --------------------------------------------------------------------
    # Full correlator
    # --------------------------------------------------------------------

    PowspecMeasurements compute_powspec_cpp "trv::compute_powspec" (
        CppParticleCatalogue& catalogue_data,
        CppParticleCatalogue& catalogue_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        CppBinning& kbinning,
        double norm_factor
    )

    TwoPCFMeasurements compute_corrfunc_cpp "trv::compute_corrfunc" (
        CppParticleCatalogue& catalogue_data,
        CppParticleCatalogue& catalogue_rand,
        LineOfSight* los_data,
        LineOfSight* los_rand,
        CppParameterSet& params,
        CppBinning& rbinning,
        double norm_factor
    )

    PowspecMeasurements compute_powspec_in_gpp_box_cpp \
        "trv::compute_powspec_in_gpp_box" (
            CppParticleCatalogue& catalogue_data,
            CppParameterSet& params,
            CppBinning& kbinning,
            double norm_factor
        )

    TwoPCFMeasurements compute_corrfunc_in_gpp_box_cpp \
        "trv::compute_corrfunc_in_gpp_box" (
            CppParticleCatalogue& catalogue_data,
            CppParameterSet& params,
            CppBinning& rbinning,
            double norm_factor
        )

    TwoPCFWindowMeasurements compute_corrfunc_window_cpp \
        "trv::compute_corrfunc_window" (
            CppParticleCatalogue& catalogue_rand,
            LineOfSight* los_rand,
            CppParameterSet& params,
            CppBinning& rbinning,
            double alpha,
            double norm_factor
        )


def _calc_powspec_normalisation_from_particles(
        _ParticleCatalogue particles not None, double alpha
    ):
    return calc_powspec_normalisation_from_particles_cpp(
        deref(particles.thisptr), alpha
    )


def _calc_powspec_normalisation_from_mesh(
        _ParticleCatalogue particles not None,
        ParameterSet params not None,
        double alpha
    ):
    return calc_powspec_normalisation_from_mesh_cpp(
        deref(particles.thisptr), deref(params.thisptr), alpha
    )


def _calc_powspec_normalisation_from_meshes(
        _ParticleCatalogue particles_data not None,
        _ParticleCatalogue particles_rand not None,
        ParameterSet params not None,
        double alpha,
        padding=None, cellsize=None, assignment=None
    ):
    if None in [padding, cellsize, assignment]:
        return calc_powspec_normalisation_from_meshes_cpp(
            deref(particles_data.thisptr), deref(particles_rand.thisptr),
            deref(params.thisptr), alpha
        )
    else:  # STYLE: non-Pythonic use of `else`
        if not (
            isinstance(padding, float)
            and isinstance(cellsize, float)
            and isinstance(assignment, str)
        ):
            raise TypeError(
                "`padding`, `cellsize` and `assignment` must be "
                "of type float, float and str: received {}, {} and {}."
                .format(type(padding), type(cellsize), type(assignment))
            )
        return calc_powspec_normalisation_from_meshes_cpp(
            deref(particles_data.thisptr), deref(particles_rand.thisptr),
            deref(params.thisptr), alpha,
            padding, cellsize, assignment.encode('utf-8')
        )


def _compute_powspec(
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
    cdef PowspecMeasurements results
    results = compute_powspec_cpp(
        deref(catalogue_data.thisptr), deref(catalogue_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr), deref(kbinning.thisptr),
        norm_factor
    )

    free(los_data_cpp); free(los_rand_cpp)

    return {
        'kbin': np.asarray(results.kbin),
        'keff': np.asarray(results.keff),
        'nmodes': np.asarray(results.nmodes),
        'pk_raw': np.asarray(results.pk_raw),
        'pk_shot': np.asarray(results.pk_shot),
    }


def _compute_corrfunc(
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
    cdef TwoPCFMeasurements results
    results = compute_corrfunc_cpp(
        deref(catalogue_data.thisptr), deref(catalogue_rand.thisptr),
        los_data_cpp, los_rand_cpp,
        deref(params.thisptr), deref(rbinning.thisptr),
        norm_factor
    )

    free(los_data_cpp); free(los_rand_cpp)

    return {
        'rbin': np.asarray(results.rbin),
        'reff': np.asarray(results.reff),
        'npairs': np.asarray(results.npairs),
        'xi': np.asarray(results.xi),
    }


def _compute_powspec_in_gpp_box(
        _ParticleCatalogue catalogue_data not None,
        ParameterSet params not None,
        Binning kbinning not None,
        double norm_factor
    ):
    cdef PowspecMeasurements results
    results = compute_powspec_in_gpp_box_cpp(
        deref(catalogue_data.thisptr),
        deref(params.thisptr), deref(kbinning.thisptr),
        norm_factor
    )

    return {
        'kbin': np.asarray(results.kbin),
        'keff': np.asarray(results.keff),
        'nmodes': np.asarray(results.nmodes),
        'pk_raw': np.asarray(results.pk_raw),
        'pk_shot': np.asarray(results.pk_shot),
    }


def _compute_corrfunc_in_gpp_box(
        _ParticleCatalogue catalogue_data not None,
        ParameterSet params not None,
        Binning rbinning not None,
        double norm_factor
    ):
    cdef TwoPCFMeasurements results
    results = compute_corrfunc_in_gpp_box_cpp(
        deref(catalogue_data.thisptr),
        deref(params.thisptr), deref(rbinning.thisptr),
        norm_factor
    )

    return {
        'rbin': np.asarray(results.rbin),
        'reff': np.asarray(results.reff),
        'npairs': np.asarray(results.npairs),
        'xi': np.asarray(results.xi),
    }


def _compute_corrfunc_window(
        _ParticleCatalogue catalogue_rand not None,
        np.ndarray[double, ndim=2, mode='c'] los_rand not None,
        ParameterSet params not None,
        Binning rbinning not None,
        double alpha,
        double norm_factor
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
    cdef TwoPCFWindowMeasurements results
    results = compute_corrfunc_window_cpp(
        deref(catalogue_rand.thisptr), los_rand_cpp,
        deref(params.thisptr), deref(rbinning.thisptr),
        alpha, norm_factor
    )

    free(los_rand_cpp)

    return {
        'rbin': np.asarray(results.rbin),
        'reff': np.asarray(results.reff),
        'npairs': np.asarray(results.npairs),
        'xi': np.asarray(results.xi),
    }
