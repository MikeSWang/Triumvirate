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

bool if_filepath_is_set(const std::string& pathstr) {
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

  while (!dirstr.empty() && (dirstr.back() == '/')) {
    dirstr.pop_back();
  }

  std::filesystem::path dir(dirstr);
  std::error_code ec;
  if (std::filesystem::create_directories(dir, ec)) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.info("Directory created: %s", dirstr.c_str());
    }
  } else {
    if (ec) {
      if (trv::sys::currTask == 0) {
        trv::sys::logger.error(
          "Failed to create directory: %s", dirstr.c_str()
        );
      }
      throw trv::sys::IOError(
        "Failed to create directory: %s\n", dirstr.c_str()
      );
    }
  }
}

}  // namespace trv::sys


// ***********************************************************************
// Files
// ***********************************************************************

namespace io {

// -----------------------------------------------------------------------
// Pre-measurement header
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
// Binning details
// -----------------------------------------------------------------------

void print_binned_vectors_to_file(
  std::FILE* fileptr, trv::ParameterSet& params,
  trv::BinnedVectors& binned_vectors
) {
  // Print header.
  std::fprintf(
    fileptr,
    "%s Box size: [%.3f, %.3f, %.3f]\n",
    comment_delimiter,
    params.boxsize[0], params.boxsize[1], params.boxsize[2]
  );
  std::fprintf(
    fileptr,
    "%s Mesh number: [%d, %d, %d]\n",
    comment_delimiter,
    params.ngrid[0], params.ngrid[1], params.ngrid[2]
  );
  std::fprintf(
    fileptr,
    "%s Vector count: %d\n", comment_delimiter, binned_vectors.count
  );
  std::fprintf(
    fileptr,
    "%s Bin number: %d\n", comment_delimiter, binned_vectors.num_bins
  );

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] bin_index, [1] bin_lower, [2] bin_upper, "
    "[3] vec_x, [4] vec_y, [5] vec_z\n",
    comment_delimiter
  );

  // Print data table.
  for (int ivec = 0; ivec < binned_vectors.count; ivec++) {
    std::fprintf(
      fileptr,
      "%d\t%.9e\t%.9e\t% .9e\t% .9e\t% .9e\n",
      binned_vectors.indices[ivec],
      binned_vectors.lower_edges[ivec],
      binned_vectors.upper_edges[ivec],
      binned_vectors.vecx[ivec],
      binned_vectors.vecy[ivec],
      binned_vectors.vecz[ivec]
    );
  }
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
  for (int idx_dv = 0; idx_dv < meas_powspec.dim; idx_dv++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_powspec.kbin[idx_dv],
      meas_powspec.keff[idx_dv],
      meas_powspec.nmodes[idx_dv],
      meas_powspec.pk_raw[idx_dv].real(), meas_powspec.pk_raw[idx_dv].imag(),
      meas_powspec.pk_shot[idx_dv].real(), meas_powspec.pk_shot[idx_dv].imag()
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
  for (int idx_dv = 0; idx_dv < meas_2pcf.dim; idx_dv++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t% .9e\t% .9e\n",
      meas_2pcf.rbin[idx_dv],
      meas_2pcf.reff[idx_dv],
      meas_2pcf.npairs[idx_dv],
      meas_2pcf.xi[idx_dv].real(), meas_2pcf.xi[idx_dv].imag()
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
  for (int idx_dv = 0; idx_dv < meas_2pcf_win.dim; idx_dv++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t% .9e\t% .9e\n",
      meas_2pcf_win.rbin[idx_dv],
      meas_2pcf_win.reff[idx_dv],
      meas_2pcf_win.npairs[idx_dv],
      meas_2pcf_win.xi[idx_dv].real(), meas_2pcf_win.xi[idx_dv].imag()
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
  char multipole_str[8];
  std::snprintf(
    multipole_str, sizeof(multipole_str), "%d%d%d",
    params.ell1, params.ell2, params.ELL
  );

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] k1_cen, [1] k1_eff, [2] nmodes_1, "
    "[3] k2_cen, [4] k2_eff, [5] nmodes_2, "
    "[6] Re{bk%s_raw}, [7] Im{bk%s_raw}, "
    "[8] Re{bk%s_shot}, [9] Im{bk%s_shot}\n",
    comment_delimiter,
    multipole_str, multipole_str, multipole_str, multipole_str
  );

  // Print data table.
  for (int idx_dv = 0; idx_dv < meas_bispec.dim; idx_dv++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_bispec.k1_bin[idx_dv], meas_bispec.k1_eff[idx_dv],
      meas_bispec.nmodes_1[idx_dv],
      meas_bispec.k2_bin[idx_dv], meas_bispec.k2_eff[idx_dv],
      meas_bispec.nmodes_2[idx_dv],
      meas_bispec.bk_raw[idx_dv].real(), meas_bispec.bk_raw[idx_dv].imag(),
      meas_bispec.bk_shot[idx_dv].real(), meas_bispec.bk_shot[idx_dv].imag()
    );
  }
}

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ThreePCFMeasurements& meas_3pcf
) {
  char multipole_str[8];
  std::snprintf(
    multipole_str, sizeof(multipole_str), "%d%d%d",
    params.ell1, params.ell2, params.ELL
  );

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] r1_cen, [1] r1_eff, [2] npairs_1, "
    "[3] r2_cen, [4] r2_eff, [5] npairs_2, "
    "[6] Re{zeta%s_raw}, [7] Im{zeta%s_raw}, "
    "[8] Re{zeta%s_shot}, [9] Im{zeta%s_shot}\n",
    comment_delimiter,
    multipole_str, multipole_str, multipole_str, multipole_str
  );

  // Print data table.
  for (int idx_dv = 0; idx_dv < meas_3pcf.dim; idx_dv++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_3pcf.r1_bin[idx_dv], meas_3pcf.r1_eff[idx_dv],
      meas_3pcf.npairs_1[idx_dv],
      meas_3pcf.r2_bin[idx_dv], meas_3pcf.r2_eff[idx_dv],
      meas_3pcf.npairs_2[idx_dv],
      meas_3pcf.zeta_raw[idx_dv].real(), meas_3pcf.zeta_raw[idx_dv].imag(),
      meas_3pcf.zeta_shot[idx_dv].real(), meas_3pcf.zeta_shot[idx_dv].imag()
    );
  }
}

void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ThreePCFWindowMeasurements& meas_3pcf_win
) {
  char multipole_str[8];
  std::snprintf(
    multipole_str, sizeof(multipole_str), "%d%d%d",
    params.ell1, params.ell2, params.ELL
  );

  // Print data table columns.
  std::fprintf(
    fileptr,
    "%s "
    "[0] r1_cen, [1] r1_eff, [2] npairs_1, "
    "[3] r2_cen, [4] r2_eff, [5] npairs_2, "
    "[6] Re{zeta%s_raw}, [7] Im{zeta%s_raw}, "
    "[8] Re{zeta%s_shot}, [9] Im{zeta%s_shot}\n",
    comment_delimiter,
    multipole_str, multipole_str, multipole_str, multipole_str
  );

  // Print data table.
  for (int idx_dv = 0; idx_dv < meas_3pcf_win.dim; idx_dv++) {
    std::fprintf(
      fileptr,
      "%.9e\t%.9e\t%10d\t%.9e\t%.9e\t%10d\t% .9e\t% .9e\t% .9e\t% .9e\n",
      meas_3pcf_win.r1_bin[idx_dv], meas_3pcf_win.r1_eff[idx_dv],
      meas_3pcf_win.npairs_1[idx_dv],
      meas_3pcf_win.r2_bin[idx_dv], meas_3pcf_win.r2_eff[idx_dv],
      meas_3pcf_win.npairs_2[idx_dv],
      meas_3pcf_win.zeta_raw[idx_dv].real(),
      meas_3pcf_win.zeta_raw[idx_dv].imag(),
      meas_3pcf_win.zeta_shot[idx_dv].real(),
      meas_3pcf_win.zeta_shot[idx_dv].imag()
    );
  }
}

}  // namespace trv::io
}  // namespace trv
