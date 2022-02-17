from libcpp.vector cimport vector


cdef extern from "include/particles.hpp":
    cppclass CppParticleCatalogue "ParticleCatalogue":
        CppParticleCatalogue()
        int read_particles_catalogue(
            vector[double], vector[double], vector[double], vector[double]
        )


cdef class ParticleCatalogue:
    cdef CppParticleCatalogue* thisptr
