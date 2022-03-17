"""
Catalogue Parser (:mod:`~triumvirate._catalogue`)
===========================================================================

Parse Python catalogue objects into C++ catalogue objects.

"""
import numpy as np
cimport numpy as np


cdef class _ParticleCatalogue:

    def __cinit__(
        self,
        np.ndarray[double, ndim=1, mode='c'] x not None,
        np.ndarray[double, ndim=1, mode='c'] y not None,
        np.ndarray[double, ndim=1, mode='c'] z not None,
        np.ndarray[double, ndim=1, mode='c'] nz not None,
        np.ndarray[double, ndim=1, mode='c'] ws not None,
        np.ndarray[double, ndim=1, mode='c'] wc not None
    ):

        self.thisptr = new CppParticleCatalogue()

        self.thisptr.read_particle_data(x, y, z, nz, ws, wc)

    def __dealloc__(self):
        del self.thisptr
