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
 * @author Mike S Wang (https://github.com/MikeSWang)
 * @brief Program parameter configuration.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_

#include <cstdio>
#include <cstring>
#include <fstream>
#include <string>
#include <sys/stat.h>

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
  /// --------------------------------------------------------------------
  /// I/O
  /// --------------------------------------------------------------------

  std::string catalogue_dir;        ///< catalogue directory
  std::string measurement_dir;      ///< measurement/output directory
  std::string data_catalogue_file;  ///< data catalogue file
  std::string rand_catalogue_file;  ///< random catalogue file
  std::string catalogue_columns;    ///< catalogue data columns
                                    ///< (comma-separated)
  std::string output_tag;           ///< output tag

  /// --------------------------------------------------------------------
  /// Mesh sampling
  /// --------------------------------------------------------------------

  /// Mesh properties.
  double boxsize[3];      ///< box size (in Mpc/h) in each dimension
  unsigned int ngrid[3];  ///< grid number in each dimension

  /// Derived mesh quantities.
  double volume;       ///< box volume (in Mpc^3/h^3)
  unsigned int nmesh;  ///< number of mesh grids

  /// Mesh alignment.
  std::string alignment = "centre";  ///< box alignment:
                                     ///< {'centre' (default), 'pad'}
  std::string padscale = "box";      ///< padding scale
                                     ///< (if `alignment` is 'pad'):
                                     ///< {'box' (default), 'grid'}
  double padfactor;                  ///< padding factor

  /// Mesh assignment.
  std::string assignment = "tsc";   ///< mesh assignment scheme:
                                    ///< {'ngp', 'cic',
                                    ///<  'tsc' (default), 'pcs'}
  std::string interlace = "false";  ///< interlacing switch:
                                    ///< {'true'/'on',
                                    ///<  'false'/'off' (default)}

  /// --------------------------------------------------------------------
  /// Measurement
  /// --------------------------------------------------------------------

  /// Measurement specification.
  std::string catalogue_type;    ///< catalogue type:
                                 ///< {'survey', 'random', 'sim'}
  std::string measurement_type;  ///< measurement_type:
                                 ///< {'powspec', '2pcf', '2pcf-win',
                                 ///<  'bispec', '3pcf', '3pcf-win',
                                 ///<  '3pcf-win-wa'}

  std::string norm_convention = "particle";   ///< normalisation
                                              ///< convention:
                                              ///< {'mesh',
                                              ///<  'particle' (default)}
  std::string shotnoise_convention = "mesh";  ///< shot noise convention:
                                              ///< {'mesh' (default),
                                              ///<  'particle'}

  std::string binning = "lin";  ///< binning scheme:
                                ///< {'lin' (default), 'log',
                                ///<  'linpad', 'logpad', 'custom'}
  std::string form = "diag";    ///< form of the bispectrum measurement:
                                ///< {'diag' (default), 'full'}

  /// Derived measurement specification.
  std::string npoint;  ///< npoint case: {'2pt', '3pt'}
  std::string space;   ///< coordinte space: {'fourier', 'config'}

  /// Measurement parameters.
  int ell1;  ///< spherical degree associated with the first wavevector
  int ell2;  ///< spherical degree associated with the second wavevector
  int ELL;   ///< spherical degree associated with the line of sight

  int i_wa;  ///< first order of the wide-angle correction term
  int j_wa;  ///< second order of the wide-angle correction term

  double bin_min;  ///< measurement range minimum (in Mpc/h or h/Mpc)
  double bin_max;  ///< measurement range maximum (in Mpc/h or h/Mpc)

  unsigned int num_bins;  ///< number of measurement bins
  unsigned int idx_bin;   ///< fixed bin index in 'full' `form`
                          ///< bispectrum measurements

  /**
   * @brief Read parameters from a file.
   *
   * @param parameter_file Command-line argument list.
   * @returns Validation exit status.
   */
  int read_from_file(char* parameter_file);

  /**
   * @brief Validate parameters.
   *
   * @returns Exit status.
   */
  int validate();

  /**
   * @brief Print out extracted parameters to a file in the
   *        output measurement directory.
   *
   * @returns Exit status.
   */
  int printout();
};

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
