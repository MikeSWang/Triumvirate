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
 * @file tools.hpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 * @brief Miscellaneous numerical tools.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_

#include <cmath>
#include <complex>
#include <cstdarg>
#include <stdexcept>
#include <vector>

#include "monitor.hpp"
#include "parameters.hpp"

namespace trv {

/// **********************************************************************
/// Arrays
/// **********************************************************************

namespace mon {

/**
 * @brief Exception raised when an extrapolation error occurs.
 *
 */
class ExtrapError: public std::runtime_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an ExtrapError exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  ExtrapError(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @return String representation of the exception.
   */
  virtual const char* what() const noexcept;
};

}  // namespace trv::mon

namespace ops {

/**
 * @brief Extrapolate sample series exponentially (i.e. log-linearly).
 *
 * @param[in] a Samples series.
 * @param[in] N Sample number.
 * @param[in] N_ext Extrapolation number on either end.
 * @param[out] a_ext Extrapolated sample series.
 * @throws trv::mon::ExtrapError Extrapolation error when sign change
 *                               or zero occurs.
 */
void extrap_loglin(double* a, int N, int N_ext, double* a_ext);

/**
 * @brief Extrapolate sample bi-series exponentially
 *        (i.e. log-bilinearly).
 *
 * @param[in] a Samples bi-series.
 * @param[in] N Sample number (in both dimensions).
 * @param[in] N_ext Extrapolation number on either end
 *                  (in both dimensions).
 * @param[out] a_ext Extrapolated sample bi-series.
 * @throws trv::mon::ExtrapError Extrapolation error when sign change
 *                               or zero occurs.
 */
void extrap2d_logbilin(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
);

/**
 * @brief Extrapolate sample bi-series linearly (i.e. bilinearly).
 *
 * @param[in] a Samples bi-series.
 * @param[in] N Sample number (in both dimensions).
 * @param[in] N_ext Extrapolation number on either end
 *                  (in both dimensions).
 * @param[out] a_ext Extrapolated sample bi-series.
 */
void extrap2d_bilin(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
);

/**
 * @brief Extrapolate sample bi-series by zero padding.
 *
 * @param[in] a Samples bi-series.
 * @param[in] N Sample number (in both dimensions).
 * @param[in] N_ext Extrapolation number on either end
 *                  (in both dimensions).
 * @param[out] a_ext Extrapolated sample bi-series.
 */
void extrap2d_bizeros(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
);

}  // namespace trv::ops

/// **********************************************************************
/// Algorithms
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
   * @param nbin Nunber of bins.
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
   * @throws trv::mon::UnimplementedError When @c "custom" binning is not
   *                                      implemented by the user.
   * @throws trv::mon::InvalidParameter When @c scheme is not one of the
   *                                    options above.
   */
  void set_bins(std::string scheme, std::string space);

  /**
   * @brief Set bins.
   *
   * The binning scheme is inferred from @c params when initialised
   * with @ref Binning(trv::ParameterSet&).
   *
   * @overload
   */
  void set_bins();
};

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_TOOLS_HPP_INCLUDED_
