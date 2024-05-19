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
 * @file io.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief I/O support including custom exceptions and utility functions.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_IO_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_IO_HPP_INCLUDED_

#include <sys/stat.h>
#include <filesystem>
#include <string>
#include <system_error>

#include "parameters.hpp"
#include "particles.hpp"
#include "dataobjs.hpp"

namespace trv {

// ***********************************************************************
// Paths
// ***********************************************************************

namespace sys {

/**
 * @brief Check if a file path is set.
 *
 * @param pathstr File path string.
 * @returns { @c true , @c false }
 */
bool if_filepath_is_set(const std::string& pathstr);

/**
 * @brief Make write directory.
 *
 * @param dirstr Directory path string.
 * @throws trv::sys::IOError When write directory cannot be created.
 */
void make_write_dir(std::string dirstr);

}  // namespace trv::sys


// ***********************************************************************
// Files
// ***********************************************************************

namespace io {

/// @cond DOXYGEN_DOC_MISC
const char comment_delimiter[] = "#";  ///< header comment delimiter
/// @endcond

// -----------------------------------------------------------------------
// Pre-measurement header
// -----------------------------------------------------------------------

/**
 * @brief Print the pre-measurement header to a file including information
 *        about the catalogue(s) and mesh grid assignment.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param catalogue_data (Data-source) particle catalogue.
 * @param catalogue_rand (Random-source) particle catalogue.
 * @param norm_factor_part Particle-based normalisation factor.
 * @param norm_factor_mesh Mesh-based normalisation factor.
 * @param norm_factor_meshes Normalisation factor based on mixed meshes.
 */
void print_measurement_header_to_file(
  std::FILE* fileptr, trv::ParameterSet& params,
  trv::ParticleCatalogue& catalogue_data,
  trv::ParticleCatalogue& catalogue_rand,
  double norm_factor_part, double norm_factor_mesh, double norm_factor_meshes
);

/**
 * @brief Print the pre-measurement header to a file including information
 *        about the catalogue(s) and mesh grid assignment.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param catalogue Particle catalogue.
 * @param norm_factor_part Particle-based normalisation factor.
 * @param norm_factor_mesh Mesh-based normalisation factor.
 * @param norm_factor_meshes Normalisation factor based on mixed meshes.
 *
 * @overload
 */
void print_measurement_header_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ParticleCatalogue& catalogue,
  double norm_factor_part, double norm_factor_mesh, double norm_factor_meshes
);


// -----------------------------------------------------------------------
// Binning details
// -----------------------------------------------------------------------

/**
 * @brief Print binned vectors to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param binned_vectors Binned vectors.
 */
void print_binned_vectors_to_file(
  std::FILE* fileptr, trv::ParameterSet& params,
  trv::BinnedVectors& binned_vectors
);


// -----------------------------------------------------------------------
// Two-point measurement data table
// -----------------------------------------------------------------------

/**
 * @brief Print measurements as a data table to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param meas_powspec Power spectrum measurements.
 */
void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::PowspecMeasurements& meas_powspec
);

/**
 * @brief Print measurements as a data table to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param meas_2pcf Two-point correlation function measurements.
 *
 * @overload
 */
void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::TwoPCFMeasurements& meas_2pcf
);

/**
 * @brief Print measurements as a data table to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param meas_2pcf_win Two-point correlation function window measurements.
 *
 * @overload
 */
void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::TwoPCFWindowMeasurements& meas_2pcf_win
);


// -----------------------------------------------------------------------
// Three-point measurement data table
// -----------------------------------------------------------------------

/**
 * @brief Print measurements as a data table to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param meas_bispec Bispectrum measurements.
 *
 * @overload
 */
void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::BispecMeasurements& meas_bispec
);

/**
 * @brief Print measurements as a data table to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param meas_3pcf Three-point correlation function measurements.
 *
 * @overload
 */
void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ThreePCFMeasurements& meas_3pcf
);

/**
 * @brief Print measurements as a data table to a file.
 *
 * @param fileptr File to print to.
 * @param params Parameter set.
 * @param meas_3pcf_win Three-point correlation function window
 *                      measurements.
 *
 * @overload
 */
void print_measurement_datatab_to_file(
  std::FILE* fileptr,
  trv::ParameterSet& params, trv::ThreePCFWindowMeasurements& meas_3pcf_win
);

}  // namespace trv::io
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_IO_HPP_INCLUDED_
