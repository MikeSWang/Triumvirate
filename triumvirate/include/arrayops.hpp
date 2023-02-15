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
 * @file arrayops.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Array operations.
 *
 * Array operations provided include:
 * - data extrapolations, both linearly and logarithmically, in 1- or 2-d.
 */

#ifndef TRIUMVIRATE_INCLUDE_ARRAYOPS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_ARRAYOPS_HPP_INCLUDED_

#include <cmath>
#include <cstdarg>
#include <stdexcept>
#include <vector>

#include "monitor.hpp"

namespace trv {

// ***********************************************************************
// Extrapolation
// ***********************************************************************

namespace sys {

/**
 * @brief Exception raised when an extrapolation error occurs.
 *
 */
class ExtrapError: public std::runtime_error {
 public:
  std::string err_mesg;  ///< error message

  /**
   * @brief Construct an @ref trv::sys::ExtrapError exception.
   *
   * @param fmt_string Error message format string.
   * @param ... An arbitrary number of substitution arguments.
   */
  ExtrapError(const char* fmt_string, ...);

  /**
   * @brief Exception string representation.
   *
   * @returns String representation of the exception.
   */
  virtual const char* what() const noexcept;
};

}  // namespace trv::sys

namespace array {

/**
 * @brief Extrapolate sample series exponentially (i.e. log-linearly).
 *
 * @param[in] a Samples series.
 * @param[in] N Sample number.
 * @param[in] N_ext Extrapolation number on either end.
 * @param[out] a_ext Extrapolated sample series.
 * @throws trv::sys::ExtrapError When sign change or zero occurs.
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
 * @throws trv::sys::ExtrapError When sign change or zero occurs.
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

}  // namespace trv::array
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_ARRAYOPS_HPP_INCLUDED_
