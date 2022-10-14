# Copyright (C) [GPLv3 Licence]
#
# This file is part of the Triumvirate program. See the COPYRIGHT
# and LICENCE files at the top-level directory of this distribution
# for details of copyright and licensing.
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program. If not, see <https:#www.gnu.org/licenses/>.

from libcpp.string cimport string


cdef extern from "include/parameters.hpp":
    cppclass CppParameterSet "trv::ParameterSet":
        # ----------------------------------------------------------------
        # Members
        # ----------------------------------------------------------------

        # -- I/O ---------------------------------------------------------

        string catalogue_dir
        string measurement_dir
        string data_catalogue_file
        string rand_catalogue_file
        # string catalogue_columns
        string output_tag

        # -- Mesh sampling -----------------------------------------------

        double boxsize[3]
        int ngrid[3]

        double volume
        int nmesh

        string alignment
        string padscale
        double padfactor

        string assignment
        string interlace

        # -- Measurement -------------------------------------------------

        string catalogue_type
        string measurement_type

        string norm_convention

        string binning
        string form

        string npoint
        string space

        int ell1
        int ell2
        int ELL

        int i_wa
        int j_wa

        double bin_min
        double bin_max
        int num_bins
        int idx_bin

        # -- Misc --------------------------------------------------------

        int verbose

        # ----------------------------------------------------------------
        # Methods
        # ----------------------------------------------------------------

        int validate()


cdef class ParameterSet:
    cdef CppParameterSet* thisptr
    cdef string _source
    cdef string _status
    cdef dict _params
    cdef object _logger
