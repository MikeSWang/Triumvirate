// Copyright (C) [GPLv3 Licence]
//
// This file is part of the Triumvirate program. See the COPYRIGHT
// and LICENCE files at the top-level directory of this distribution
// for details of copyright and licensing.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

/**
 * @file io.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "io.hpp"

namespace trv {

/// **********************************************************************
/// System
/// **********************************************************************

namespace sys {

bool if_filepath_is_set(std::string pathstr) {
  /// Check if the string is empty.
  if (pathstr.empty()) {return false;}

  /// Check if the path is a directory not file.
  std::string endchar = "/";
  int strcomp = pathstr.compare(
    pathstr.length() - endchar.length(), endchar.length(), endchar
  );
  if (strcomp == 0) {return false;}  // `pathstr` ends in "/"

  /// Check if the string contains non-whitespace characters.  If so,
  /// the path is set, otherwise not.
  for (int ichar = 0; ichar < pathstr.length(); ichar++){
    if (!std::isspace(pathstr[ichar])) {return true;}
  }

  return false;
}

}  // namespace trv::sys


/// **********************************************************************
/// Program
/// **********************************************************************

void print_premeasurement_info(
  std::FILE* fileptr, trv::ParameterSet& params,
  trv::ParticleCatalogue& catalogue_data, trv::ParticleCatalogue& catalogue_rand
) {
  std::fprintf(
    fileptr,
    "# Data catalogue: %d particles of total weight %.3f\n",
    catalogue_data.ntotal, catalogue_data.wtotal
  );
  std::fprintf(
    fileptr,
    "# Random catalogue: %d particles of total weight %.3f\n",
    catalogue_rand.ntotal, catalogue_rand.wtotal
  );
  std::fprintf(
    fileptr,
    "# Box size: (%.3f, %.3f, %.3f)\n",
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    fileptr,
    "# Mesh number: (%d, %d, %d)\n",
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    fileptr,
    "# Data-source particle extents: "
    "([%.3f, %.3f], [%.3f, %.3f], [%.3f, %.3f])\n",
    catalogue_data.pos_min[0], catalogue_data.pos_max[0],
    catalogue_data.pos_min[1], catalogue_data.pos_max[1],
    catalogue_data.pos_min[2], catalogue_data.pos_max[2]
  );
  std::fprintf(
    fileptr,
    "# Random-source particle extents: "
    "([%.3f, %.3f], [%.3f, %.3f], [%.3f, %.3f])\n",
    catalogue_rand.pos_min[0], catalogue_rand.pos_max[0],
    catalogue_rand.pos_min[1], catalogue_rand.pos_max[1],
    catalogue_rand.pos_min[2], catalogue_rand.pos_max[2]
  );
  std::fprintf(
    fileptr,
    "# Box alignment: %s,\n",
    params.alignment.c_str()
  );
  std::fprintf(
    fileptr,
    "# Mesh assignment and interlacing: %s, %s\n",
    params.assignment.c_str(), params.interlace.c_str()
  );
}

void print_premeasurement_info(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ParticleCatalogue& catalogue
) {
  std::fprintf(
    fileptr,
    "# Catalogue: %d particles of total weight %.3f\n",
    catalogue.ntotal, catalogue.wtotal
  );
  std::fprintf(
    fileptr,
    "# Box size: (%.3f, %.3f, %.3f)\n",
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    fileptr,
    "# Mesh number: (%d, %d, %d)\n",
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    fileptr,
    "# Particle extents: "
    "([%.3f, %.3f], [%.3f, %.3f], [%.3f, %.3f])\n",
    catalogue.pos_min[0], catalogue.pos_max[0],
    catalogue.pos_min[1], catalogue.pos_max[1],
    catalogue.pos_min[2], catalogue.pos_max[2]
  );
  std::fprintf(
    fileptr,
    "# Box alignment: %s,\n",
    params.alignment.c_str()
  );
  std::fprintf(
    fileptr,
    "# Mesh assignment and interlacing: %s, %s\n",
    params.assignment.c_str(), params.interlace.c_str()
  );
}



}  // namespace trv
