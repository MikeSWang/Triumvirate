from libcpp.string cimport string

cdef class ParameterSet:
    cdef string _source
    cdef string _status
    cdef dict _params
