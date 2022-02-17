"""

"""
import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from cython.operator cimport dereference as deref

from params cimport MeasurementParameters, Parameters
from _morph cimport ParticleCatalogue, CppParticleCatalogue

cdef extern from "include/threept.hpp":
    int calc_bispec_(
        CppParticleCatalogue, CppParticleCatalogue,
        vector[vector[double]], vector[vector[double]],
        Parameters&,
        double, double*
    )

def calc_bispec_py(
        ParticleCatalogue particles_data, ParticleCatalogue particles_rand,
        np.ndarray[double, ndim=2, mode='c'] los_data,
        np.ndarray[double, ndim=2, mode='c'] los_rand,
        MeasurementParameters params,
        double alpha,
        np.ndarray[double, ndim=1, mode='c'] kbin
    ):
    return calc_bispec_(
        deref(particles_data.thisptr), deref(particles_rand.thisptr),
        los_data, los_rand,
        deref(params.thisptr),
        alpha,
        &kbin[0]
    )