"""
Object Morphism (:mod:`~triumvirate.measurements._morph`)
===========================================================================

Convert relevant Python objects into C++ objects.

"""
import numpy as np
cimport numpy as np
from cython.operator cimport dereference as deref


cdef class ParticleCatalogue:
    def __cinit__(self,
                  np.ndarray[double, ndim=1, mode='c'] x not None,
                  np.ndarray[double, ndim=1, mode='c'] y not None,
                  np.ndarray[double, ndim=1, mode='c'] z not None,
                  np.ndarray[double, ndim=1, mode='c'] w not None):

        self.thisptr = new CppParticleCatalogue()
        self.thisptr.read_particles_catalogue(x, y, z, w)
