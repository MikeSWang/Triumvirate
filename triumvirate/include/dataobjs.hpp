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
 * @file dataobjs.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang)
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Clustering measurement data objects.
 *
 * Clustering measurement data objects provided include:
 * - binning schemes;
 * - line of sight; and
 * - clustering statistics.
 */

#ifndef TRIUMVIRATE_INCLUDE_DATAOBJS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_DATAOBJS_HPP_INCLUDED_

#include <cmath>
#include <vector>

#include "monitor.hpp"
#include "parameters.hpp"

namespace trv {

/// **********************************************************************
/// Binning schemes
/// **********************************************************************

/**
 * @brief Isotropic coordinate binning.
 *
 * This sets up isotropic wavenumber or separation bins
 * in configuration or Fourier space.
 *
 */
class Binning {
 public:
  std::string scheme;               ///< binning scheme
  std::string space;                ///< coordinate space
  double bin_min;                   ///< lowest bin edge
  double bin_max;                   ///< highest bin edge
  int num_bins;                     ///< number of bins
  std::vector<double> bin_centres;  ///< bin centres
  std::vector<double> bin_edges;    ///< bin edges

  /**
   * @brief Construct binnng from bin specification.
   *
   * @param coord_min Minimum coordinate in binning range.
   * @param coord_max Maximum coordinate in binning range.
   * @param nbin Number of bins.
   */
  Binning(double coord_min, double coord_max, int nbin);

  /**
   * @brief Construct binning from a parameter set.
   *
   * @param params Paramater set.
   *
   * @overload
   */
  Binning(trv::ParameterSet& params);

  /**
   * @brief Set bins.
   *
   * @param scheme Binning scheme, one of
   *               {"lin", "log", "linpad", "logpad", "custom"}.
   * @param space Coordinate space, one of {"fourier", "config"}.
   * @throws trv::sys::UnimplementedError When @c binning is set to
   *                                      "custom" but not implemented.
   * @throws trv::sys::InvalidParameter When @c scheme is not one of the
   *                                    options above.
   */
  void set_bins(std::string scheme, std::string space);

  /**
   * @brief Set bins.
   *
   * The binning scheme is inferred from @c params when initialised
   * with @ref trv::Binning::Binning(trv::ParameterSet&).
   *
   * @overload
   */
  void set_bins();
};


/// **********************************************************************
/// Line of sight
/// **********************************************************************

/**
 * @brief Line-of-sight vector.
 *
 */
struct LineOfSight {
  double pos[3];  ///< 3-d position vector
};

/// **********************************************************************
/// Clustering statistics
/// **********************************************************************

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_DATAOBJS_HPP_INCLUDED_
