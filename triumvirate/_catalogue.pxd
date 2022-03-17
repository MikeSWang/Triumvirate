from libcpp.vector cimport vector


cdef extern from "include/particles.hpp":
    cppclass CppParticleCatalogue "ParticleCatalogue":
        CppParticleCatalogue()
        int read_particle_data(
            vector[double], vector[double], vector[double],
            vector[double], vector[double], vector[double]
        )


cdef class _ParticleCatalogue:
    cdef CppParticleCatalogue* thisptr
