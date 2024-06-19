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
 * @file parameters.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Program parameter configuration.
 *
 * This module provides a parameter set object with methods for reading,
 * validating and printing parameters (from/to a file).
 */

#ifndef TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_

#include <sys/stat.h>

#include <fftw3.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>

#include "monitor.hpp"

namespace trv {

/**
 * @brief Parameter set.
 *
 * This reads parameters from a file, stores and prints out the extracted
 * parameters, and validates the parameters.
 *
 */
class ParameterSet {
 public:
  // ---------------------------------------------------------------------
  // I/O
  // ---------------------------------------------------------------------

  /// catalogue directory
  std::string catalogue_dir;
  /// measurement/output directory
  std::string measurement_dir;
  /// data catalogue file
  std::string data_catalogue_file;
  /// random catalogue file
  std::string rand_catalogue_file;
  /// catalogue data columns (comma-separated without space)
  std::string catalogue_columns;
  /// output tag
  std::string output_tag;

  // ---------------------------------------------------------------------
  // Mesh sampling
  // ---------------------------------------------------------------------

  // Mesh properties.
  /// box size (in Mpc/h) in each dimension
  double boxsize[3] = {0., 0., 0.};
  /// grid number in each dimension
  int ngrid[3] = {0, 0, 0};

  // Mesh alignment.
  /// box alignment: {"centre" (default), "pad"}
  std::string alignment = "centre";
  /// padding scale (if @c alignment is "pad"): {"box" (default), "grid"}
  std::string padscale = "box";
  /// padding factor
  double padfactor = 0.;

  // Mesh assignment.
  /// mesh assignment scheme: {"ngp", "cic", "tsc" (default), "pcs"}
  std::string assignment = "tsc";
  /// interlacing switch: {"true"/"on", "false"/"off" (default)}
  std::string interlace = "false";

  // Derived mesh quantities.
  double volume;    ///< box volume (in Mpc^3/h^3)
  long long nmesh;  ///< number of mesh grid cells

  int assignment_order = 0;  ///< order of the assignment scheme

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

  // Measurement type.
  /// catalogue type: {"survey", "random", "sim", "none"}
  std::string catalogue_type;
  /// statistic type: {"powspec", "2pcf", "2pcf-win", "bispec", "3pcf",
  ///                  "3pcf-win", "3pcf-win-wa", "modes", "pairs"}
  std::string statistic_type;

  // Derived measurement type.
  std::string npoint;  ///< <i>N</i>-point case: {"2pt", "3pt", "none"}
  std::string space;   ///< coordinte space: {"fourier", "config"}

  // Measurement indexing.
  /// spherical degree associated with the first wavevector
  int ell1 = 0;
  /// spherical degree associated with the second wavevector
  int ell2 = 0;
  /// spherical degree associated with the line of sight
  int ELL = 0;

  int i_wa = 0;  ///< first order of the wide-angle correction term
  int j_wa = 0;  ///< second order of the wide-angle correction term

  // Measurement choices.
  /// form of the bispectrum measurement: {"full",
  ///                                      "diag" (default), "off-diag",
  ///                                      "row"}
  std::string form = "diag";

  /// normalisation convention: {"none", "particle" (default), "mesh",
  ///                            "mesh-mixed"}
  std::string norm_convention = "particle";

  // Derived measurement choices.
  /// shape of the 3PCF measurement: {"full", "diag" (default), "off-diag",
  ///                                 "row", "triu"}
  std::string shape = "diag";

  // Measurement parameters.
  /// binning scheme: {"lin" (default), "log",
  ///                  "linpad", "logpad", "custom"}
  std::string binning = "lin";

  double bin_min = 0.;  ///< measurement range minimum (in Mpc/h or h/Mpc)
  double bin_max = 0.;  ///< measurement range maximum (in Mpc/h or h/Mpc)

  /// number of measurement bins
  int num_bins = 0;
  /// fixed bin index in "off-fiag"/"row" @c form three-point measurements
  int idx_bin = 0;

  // ---------------------------------------------------------------------
  // Misc
  // ---------------------------------------------------------------------

  /// FFTW scheme: {"estimate", "measure" (default), "patient"}
  std::string fftw_scheme = "measure";

  /// derived FFTW planner flag
  unsigned fftw_planner_flag = FFTW_MEASURE;

  /// use FFTW wisdom: {"false" (default), <path-to-dir>}
  std::string use_fftw_wisdom = "false";

  /// derived FFTW wisdom file paths
  std::string fftw_wisdom_file_f;  ///< forward-transform wisdom file path
  std::string fftw_wisdom_file_b;  ///< backward-transform wisdom file path

  /// save flag/path for detailed binning of vectors: {"true",
  ///                                                  "false" (default),
  ///                                                  <relpath-to-file>}
  std::string save_binned_vectors = "false";

  /// logging verbosity level: {0  (NSET), 10 (DBUG), 20 (STAT) (default),
  ///                           30 (INFO), 40 (WARN), 50 (ERRO)}
  int verbose = 20;

  /**
   * @brief Construct a parameter set.
   */
  ParameterSet() = default;

  /**
   * @brief Destroy the parameter set.
   *
   */
  ~ParameterSet() = default;

  /**
   * @brief Construct a copied parameter set.
   *
   * @param other Parameter set to be copied.
   */
  ParameterSet(const ParameterSet& other);

  /**
   * @brief Read parameters from a file.
   *
   * @param parameter_filepath Parameter file path.
   * @returns Validation exit status.
   */
  int read_from_file(char* parameter_filepath);

  /**
   * @brief Validate parameters.
   *
   * @returns Exit status.
   * @throws trv::sys::InvalidParameterError When a parameter is invalid.
   *
   * @note This method is called by
   *       @ref trv::ParameterSet::read_from_file().
   */
  int validate();

  /**
   * @brief Print out extracted parameters to a file in the
   *        output measurement directory.
   *
   * @param out_parameter_filepath Printout parameter file path.
   * @returns Exit status.
   */
  int print_to_file(char* out_parameter_filepath);

  /**
   * @brief Print out extracted parameters to the default file path in the
   *        output measurement directory.
   *
   * @returns Exit status.
   *
   * @overload
   */
  int print_to_file();
};

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
