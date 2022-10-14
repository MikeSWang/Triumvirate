from libcpp.vector cimport vector


cdef extern from "include/particles.hpp":
    cppclass CppParticleCatalogue "trv::ParticleCatalogue":
        CppParticleCatalogue()
        int load_particle_data(
            vector[double] x,
            vector[double] y,
            vector[double] z,
            vector[double] nz,
            vector[double] ws,
            vector[double] wc
        )


cdef class _ParticleCatalogue:
    cdef CppParticleCatalogue* thisptr
