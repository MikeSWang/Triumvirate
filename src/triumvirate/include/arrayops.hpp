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
 * 1- or 2-d array operations provided include:
 * - data extrapolation;
 * - data sorting.
 */

#ifndef TRIUMVIRATE_INCLUDE_ARRAYOPS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_ARRAYOPS_HPP_INCLUDED_

#include <algorithm>
#include <cfenv>
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
 * @brief Extrapolation scheme.
 *
 */
enum ExtrapOption {
  NONE = 0,    ///< 0: none
  ZERO = 1,    ///< 1: zero
  PAD = 2,     ///< 2: pad
  LIN = 3,     ///< 3: linear
  LOGLIN = 4,  ///< 4: log-linear
};

/**
 * @brief Check if all elements of a 1-d array are close to a given value.
 *
 * @tparam T A @c typename.
 * @param arr 1-d array as a vector.
 * @param val Value to check against.
 * @param atol Absolute tolerance (default is 1.e-8).
 * @param rtol Relative tolerance (default is 1.e-5).
 * @returns @c true if all elements are close to @p val,
 *          @c false otherwise.
 */
template <typename T>
bool check_isclose(
  const std::vector<T>& arr, const T val,
  T atol = 1.e-8, T rtol = 1.e-5
) {
  for (T elem : arr) {
      if (std::abs(elem - val) > atol + rtol * std::abs(val)) {return false;}
  }
  return true;
}

/**
 * @brief Check the linearity or log-linearity of a 1-d array.
 *
 * @param a 1-d array as a vector.
 * @param check_lin If true, check linearity.
 * @param check_loglin If true, check log-linearity.  This also checks
 *                     for any sign change or zeros.
 * @param check_sign If true, check for any sign change or zeros.
 * @returns 0 if no check fails, 1 if linearity check fails,
 *          2 if log-linearity check fails, or 3 if sign check fails.
 */
int check_1d_array(
  std::vector<double>& a, bool check_lin, bool check_loglin, bool check_sign
);

/**
 * @brief Extrapolate a 1-d array linearly.
 *
 * @param[in] a 1-d array as a vector.
 * @param[in] N_ext Number of extra elements on either side.
 * @param[out] a_ext Extrapolated 1-d array.
 */
void extrap_lin(std::vector<double>& a, int N_ext, std::vector<double>& a_ext);

/**
 * @brief Extrapolate a 1-d array exponentially (i.e. log-linearly).
 *
 * @param[in] a 1-d array as a vector.
 * @param[in] N_ext Number of extra elements on either side.
 * @param[out] a_ext Extrapolated 1-d array.
 */
void extrap_loglin(
  std::vector<double>& a, int N_ext, std::vector<double>& a_ext
);

/**
 * @brief Extrapolate a 1-d array by constant padding.
 *
 * @param[in] a 1-d array as a vector.
 * @param[in] N_ext Number of extra elements on either side.
 * @param[in] c_lower Lower-end padding constant.
 * @param[in] c_upper Upper-end padding constant.
 * @param[out] a_ext Extrapolated 1-d array.
 */
void extrap_pad(
  std::vector<double>& a,
  int N_ext, double c_lower, double c_upper,
  std::vector<double>& a_ext
);

/**
 * @brief Extrapolate a 2-d array bi-linearly.
 *
 * @param[in] a 2-d array.
 * @param[in] N_row_ext Number of extra elements on either side
 *                      of each row.
 * @param[in] N_col_ext Number of extra elements on either side
 *                      of each column.
 * @param[out] a_ext Extrapolated 2-d array.
 */
void extrap2d_lin(
  std::vector< std::vector<double> >& a,
  int N_row_ext, int N_col_ext,
  std::vector< std::vector<double> >& a_ext
);

/**
 * @brief Extrapolate a 2-d array bi-exponentially (i.e. log-bilinearly).
 *
 * @param[in] a 2-d array.
 * @param[in] N_row_ext Number of extra elements on either side
 *                      of each row.
 * @param[in] N_col_ext Number of extra elements on either side
 *                      of each column.
 * @param[out] a_ext Extrapolated 2-d array.
 */
void extrap2d_loglin(
  std::vector< std::vector<double> >& a,
  int N_row_ext, int N_col_ext,
  std::vector< std::vector<double> >& a_ext
);

/**
 * @brief Extrapolate a 2-d array by constant padding.
 *
 * @param[in] a 2-d array.
 * @param[in] N_row_ext Number of extra elements on either side
 *                      of each row.
 * @param[in] N_col_ext Number of extra elements on either side
 *                      of each column.
 * @param[in] c_row_lower Lower-end padding constant for each row.
 * @param[in] c_row_upper Upper-end padding constant for each row.
 * @param[in] c_col_lower Lower-end padding constant for each column.
 * @param[in] c_col_upper Upper-end padding constant for each column.
 * @param[out] a_ext Extrapolated 2-d array.
 */
void extrap2d_pad(
  std::vector< std::vector<double> >& a,
  int N_row_ext, int N_col_ext,
  double c_row_lower, double c_row_upper,
  double c_col_lower, double c_col_upper,
  std::vector< std::vector<double> >& a_ext
);


// ***********************************************************************
// Sorting
// ***********************************************************************

/**
 * @brief Get the sorted indices.
 *
 * @param sorting_vector Sorting vector.
 * @return Sorted indices.
 */
std::vector<int> get_sorted_indices(std::vector<int> sorting_vector);

}  // namespace trv::array

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_ARRAYOPS_HPP_INCLUDED_
