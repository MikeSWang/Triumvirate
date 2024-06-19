"""
Catalogue Parser (:mod:`~triumvirate._particles`)
==========================================================================

Parse Python catalogue objects into C++ particle catalogues.

"""
cimport numpy as np

from ._particles cimport CppParticleCatalogue

np.import_array()


cdef class _ParticleCatalogue:

    def __cinit__(
        self,
        np.ndarray[double, ndim=1, mode='c'] x not None,
        np.ndarray[double, ndim=1, mode='c'] y not None,
        np.ndarray[double, ndim=1, mode='c'] z not None,
        np.ndarray[double, ndim=1, mode='c'] nz not None,
        np.ndarray[double, ndim=1, mode='c'] ws not None,
        np.ndarray[double, ndim=1, mode='c'] wc not None,
        verbose=-1
    ):

        self.thisptr = new CppParticleCatalogue(verbose)

        self.thisptr.load_particle_data(x, y, z, nz, ws, wc)

    def __dealloc__(self):
        del self.thisptr
