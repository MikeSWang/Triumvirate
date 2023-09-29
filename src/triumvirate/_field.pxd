"""Interface with field and mesh objects and statistics.

"""
from libcpp.string cimport string

from .dataobjs cimport BinnedVectors, CppBinning
from .parameters cimport CppParameterSet


cdef extern from "include/field.hpp":
    cdef cppclass CppFieldStats "trv::FieldStats":
        CppFieldStats(CppParameterSet& params)

        BinnedVectors record_binned_vectors(
            CppBinning& binning,
            string save_file
        ) except +
