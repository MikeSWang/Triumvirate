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
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "io.hpp"

namespace trv {

// ***********************************************************************
// Paths
// ***********************************************************************

namespace sys {

bool if_filepath_is_set(std::string pathstr) {
  // Check if the string is empty.
  if (pathstr.empty()) {return false;}

  // Check if the path is a directory not file.
  std::string endchar = "/";
  int strcomp = pathstr.compare(
    pathstr.length() - endchar.length(), endchar.length(), endchar
  );
  if (strcomp == 0) {return false;}  // `pathstr` ends in "/"

  // Check if the string contains non-whitespace characters.  If so,
  // the path is set, otherwise not.
  for (std::string::size_type ichar = 0; ichar < pathstr.length(); ichar++){
    if (!std::isspace(pathstr[ichar])) {return true;}
  }

  return false;
}

void make_write_dir(std::string dirstr) {
  // Skip special directories.
  if (dirstr.empty() || dirstr == "." || dirstr == "./" || dirstr == "/") {
    return;
  }
  if (mkdir(dirstr.c_str(), 0777) && errno != EEXIST) {
    trv::sys::logger.error(
      "Failed to create output measurement directory: %s.", dirstr.c_str()
    );
    throw trv::sys::IOError(
      "Failed to create output measurement directory: %s.\n", dirstr.c_str()
    );
  }
}

}  // namespace trv::sys


// ***********************************************************************
// Files
// ***********************************************************************

namespace io {

// -----------------------------------------------------------------------
// Measurement header
// -----------------------------------------------------------------------

void print_measurement_header_to_file(
  std::FILE* fileptr, trv::ParameterSet& params,
  trv::ParticleCatalogue& catalogue_data,
  trv::ParticleCatalogue& catalogue_rand,
  double norm_factor_part, double norm_factor_mesh, double norm_factor_meshes
) {
  std::fprintf(
    fileptr,
    "%s Data catalogue source: %s\n",
    comment_delimiter, catalogue_data.source.c_str()
  );
  std::fprintf(
    fileptr,
    "%s Data catalogue size: ntotal = %d, wtotal = %.3f, wstotal = %.3f\n",
    comment_delimiter,
    catalogue_data.ntotal, catalogue_data.wtotal, catalogue_data.wstotal
  );
  std::fprintf(
    fileptr,
    "%s Data-source particle extents: "
    "[(%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f)]\n",
    comment_delimiter,
    catalogue_data.pos_min[0], catalogue_data.pos_max[0],
    catalogue_data.pos_min[1], catalogue_data.pos_max[1],
    catalogue_data.pos_min[2], catalogue_data.pos_max[2]
  );
  std::fprintf(
    fileptr,
    "%s Random catalogue source: %s\n",
    comment_delimiter, catalogue_rand.source.c_str()
  );
  std::fprintf(
    fileptr,
    "%s Random catalogue size: ntotal = %d, wtotal = %.3f, wstotal = %.3f\n",
    comment_delimiter,
    catalogue_rand.ntotal, catalogue_rand.wtotal, catalogue_rand.wstotal
  );
  std::fprintf(
    fileptr,
    "%s Random-source particle extents: "
    "[(%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f)]\n",
    comment_delimiter,
    catalogue_rand.pos_min[0], catalogue_rand.pos_max[0],
    catalogue_rand.pos_min[1], catalogue_rand.pos_max[1],
    catalogue_rand.pos_min[2], catalogue_rand.pos_max[2]
  );

  std::fprintf(
    fileptr,
    "%s Box size: [%.3f, %.3f, %.3f]\n",
    comment_delimiter,
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    fileptr,
    "%s Box alignment: %s\n",
    comment_delimiter,
    params.alignment.c_str()
  );
  std::fprintf(
    fileptr,
    "%s Mesh number: [%d, %d, %d]\n",
    comment_delimiter,
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    fileptr,
    "%s Mesh assignment and interlacing: %s, %s\n",
    comment_delimiter,
    params.assignment.c_str(), params.interlace.c_str()
  );

  if (params.norm_convention == "none") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, 1., params.norm_convention.c_str()
    );
  } else
  if (params.norm_convention == "particle") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, norm_factor_part, params.norm_convention.c_str()
    );
  } else
  if (params.norm_convention == "mesh") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, norm_factor_mesh, params.norm_convention.c_str()
    );
  } else
  if (params.norm_convention == "mesh-mixed") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, norm_factor_meshes, params.norm_convention.c_str()
    );
  }

  std::fprintf(
    fileptr,
    "%s Normalisation factor alternatives: "
    "%.9e (particle), %.9e (mesh), %.9e (mesh-mixed)\n",
    comment_delimiter,
    norm_factor_part, norm_factor_mesh, norm_factor_meshes
  );
}

void print_measurement_header_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ParticleCatalogue& catalogue,
  double norm_factor_part, double norm_factor_mesh, double norm_factor_meshes
) {
  std::fprintf(
    fileptr,
    "%s Catalogue source: %s\n",
    comment_delimiter, catalogue.source.c_str()
  );
  std::fprintf(
    fileptr,
    "%s Catalogue size: ntotal = %d, wtotal = %.3f, wstotal = %.3f\n",
    comment_delimiter, catalogue.ntotal, catalogue.wtotal, catalogue.wstotal
  );
  std::fprintf(
    fileptr,
    "%s Catalogue particle extents: "
    "[(%.3f, %.3f), (%.3f, %.3f), (%.3f, %.3f)]\n",
    comment_delimiter,
    catalogue.pos_min[0], catalogue.pos_max[0],
    catalogue.pos_min[1], catalogue.pos_max[1],
    catalogue.pos_min[2], catalogue.pos_max[2]
  );

  std::fprintf(
    fileptr,
    "%s Box size: [%.3f, %.3f, %.3f]\n",
    comment_delimiter,
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    fileptr,
    "%s Box alignment: %s\n",
    comment_delimiter,
    params.alignment.c_str()
  );
  std::fprintf(
    fileptr,
    "%s Mesh number: [%d, %d, %d]\n",
    comment_delimiter,
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    fileptr,
    "%s Mesh assignment and interlacing: %s, %s\n",
    comment_delimiter,
    params.assignment.c_str(), params.interlace.c_str()
  );

  if (params.norm_convention == "none") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, 1., params.norm_convention.c_str()
    );
  } else
  if (params.norm_convention == "particle") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, norm_factor_part, params.norm_convention.c_str()
    );
  } else
  if (params.norm_convention == "mesh") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, norm_factor_mesh, params.norm_convention.c_str()
    );
  } else
  if (params.norm_convention == "mesh-mixed") {
    std::fprintf(
      fileptr,
      "%s Normalisation factor: %.9e (%s)\n",
      comment_delimiter, norm_factor_meshes, params.norm_convention.c_str()
    );
  }

  std::fprintf(
    fileptr,
    "%s Normalisation factor alternatives: "
    "%.9e (particle), %.9e (mesh), %.9e (mesh-mixed)\n",
    comment_delimiter,
    norm_factor_part, norm_factor_mesh, norm_factor_meshes
  );
}


// -----------------------------------------------------------------------
// Two-point measurement data table
// -----------------------------------------------------------------------

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::PowspecMeasurements& meas_powspec
) {
  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] k_cen, [1] k_eff, [2] nmodes, "
    "[3] Re{pk%d_raw}, [4] Im{pk%d_raw}, "
    "[5] Re{pk%d_shot}, [6] Im{pk%d_shot}\n",
    comment_delimiter,
    params.ELL, params.ELL, params.ELL, params.ELL
  );

  // Print data table.
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_powspec.kbin[ibin],
      meas_powspec.keff[ibin],
      meas_powspec.nmodes[ibin],
      meas_powspec.pk_raw[ibin].real(), meas_powspec.pk_raw[ibin].imag(),
      meas_powspec.pk_shot[ibin].real(), meas_powspec.pk_shot[ibin].imag()
    );
  }
}

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::TwoPCFMeasurements& meas_2pcf
) {
  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s [0] r_cen, [1] r_eff, [2] npairs, [3] Re{xi%d}, [4] Im{xi%d}\n",
    comment_delimiter,
    params.ELL, params.ELL
  );

  // Print data table.
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t% .9e\t% .9e\n",
      meas_2pcf.rbin[ibin],
      meas_2pcf.reff[ibin],
      meas_2pcf.npairs[ibin],
      meas_2pcf.xi[ibin].real(), meas_2pcf.xi[ibin].imag()
    );
  }
}

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::TwoPCFWindowMeasurements& meas_2pcf_win
) {
  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s [0] r_cen, [1] r_eff, [2] npairs, [3] Re{xi%d}, [4] Im{xi%d}\n",
    comment_delimiter,
    params.ELL, params.ELL
  );

  // Print data table.
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t% .9e\t% .9e\n",
      meas_2pcf_win.rbin[ibin],
      meas_2pcf_win.reff[ibin],
      meas_2pcf_win.npairs[ibin],
      meas_2pcf_win.xi[ibin].real(), meas_2pcf_win.xi[ibin].imag()
    );
  }
};


// -----------------------------------------------------------------------
// Three-point measurement data table
// -----------------------------------------------------------------------

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::BispecMeasurements& meas_bispec
) {
  char multipole_str[4];
  std::sprintf(multipole_str, "%d%d%d", params.ell1, params.ell2, params.ELL);

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] k1_cen, [1] k1_eff, [2] k2_cen, [3] k2_eff, [4] nmodes, "
    "[5] Re{bk%s_raw}, [6] Im{bk%s_raw}, "
    "[7] Re{bk%s_shot}, [8] Im{bk%s_shot}\n",
    comment_delimiter,
    multipole_str, multipole_str, multipole_str, multipole_str
  );

  // Print data table.
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_bispec.k1_bin[ibin], meas_bispec.k1_eff[ibin],
      meas_bispec.k2_bin[ibin], meas_bispec.k2_eff[ibin],
      meas_bispec.nmodes[ibin],
      meas_bispec.bk_raw[ibin].real(), meas_bispec.bk_raw[ibin].imag(),
      meas_bispec.bk_shot[ibin].real(), meas_bispec.bk_shot[ibin].imag()
    );
  }
}

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ThreePCFMeasurements& meas_3pcf
) {
  char multipole_str[4];
  std::sprintf(multipole_str, "%d%d%d", params.ell1, params.ell2, params.ELL);

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] r1_cen, [1] r1_eff, [2] r2_cen, [3] r2_eff, [4] npairs, "
    "[5] Re{zeta%s_raw}, [6] Im{zeta%s_raw}, "
    "[7] Re{zeta%s_shot}, [8] Im{zeta%s_shot}\n",
    comment_delimiter,
    multipole_str, multipole_str, multipole_str, multipole_str
  );

  // Print data table.
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_3pcf.r1_bin[ibin], meas_3pcf.r1_eff[ibin],
      meas_3pcf.r2_bin[ibin], meas_3pcf.r2_eff[ibin],
      meas_3pcf.npairs[ibin],
      meas_3pcf.zeta_raw[ibin].real(), meas_3pcf.zeta_raw[ibin].imag(),
      meas_3pcf.zeta_shot[ibin].real(), meas_3pcf.zeta_shot[ibin].imag()
    );
  }
}

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ThreePCFWindowMeasurements& meas_3pcf_win
) {
  char multipole_str[4];
  std::sprintf(multipole_str, "%d%d%d", params.ell1, params.ell2, params.ELL);

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] r1_cen, [1] r1_eff, [2] r2_cen, [3] r2_eff, [4] npairs, "
    "[5] Re{zeta%s_raw}, [6] Im{zeta%s_raw}, "
    "[7] Re{zeta%s_shot}, [8] Im{zeta%s_shot}\n",
    comment_delimiter,
    multipole_str, multipole_str, multipole_str, multipole_str
  );

  // Print data table.
  for (int ibin = 0; ibin < params.num_bins; ibin++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_3pcf_win.r1_bin[ibin], meas_3pcf_win.r1_eff[ibin],
      meas_3pcf_win.r2_bin[ibin], meas_3pcf_win.r2_eff[ibin],
      meas_3pcf_win.npairs[ibin],
      meas_3pcf_win.zeta_raw[ibin].real(),
      meas_3pcf_win.zeta_raw[ibin].imag(),
      meas_3pcf_win.zeta_shot[ibin].real(),
      meas_3pcf_win.zeta_shot[ibin].imag()
    );
  }
}

}  // namespace trv::io
}  // namespace trv
