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
 * @file arrayops.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "arrayops.hpp"

namespace trvs = trv::sys;

namespace trv {

// ***********************************************************************
// Extrapolation
// ***********************************************************************

namespace sys {

ExtrapError::ExtrapError(const char* fmt_string, ...) : std::runtime_error(
    "Extrapolation error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsnprintf(err_mesg_buf, sizeof(err_mesg_buf), fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* ExtrapError::what() const noexcept {return this->err_mesg.c_str();}

}  // namespace trv::sys

namespace array {

int check_1d_array(
  std::vector<double>& a, bool check_lin, bool check_loglin, bool check_sign
) {
  std::size_t N = a.size();

  if (check_lin) {
    std::vector<double> diff_a(N - 1);
    for (std::size_t i = 0; i < N - 1; i++) {
      diff_a[i] = a[i + 1] - a[i];
    }
    if (!check_isclose(diff_a, diff_a[0])) {return 1;}
  }

  if (check_loglin) {
    std::vector<double> diff_loga(N - 1);
    std::feclearexcept(FE_ALL_EXCEPT);
    for (std::size_t i = 0; i < N - 1; i++) {
      diff_loga[i] = std::log(a[i + 1] / a[i]);
      if (std::fetestexcept(FE_DIVBYZERO) || std::fetestexcept(FE_INVALID)) {
        return 3;
      }
    }
    if (!check_isclose(diff_loga, diff_loga[0])) {return 2;}
  }

  if (check_sign) {
    for (std::size_t i = 0; i < N - 1; i++) {
      if (a[i] * a[i + 1] <= 0.0) {return 3;}
    }
  }

  return 0;
}

void extrap_lin(
  std::vector<double>& a, int N_ext, std::vector<double>& a_ext
) {
  // Prepare arrays.
  int N = a.size();
  int N_notch_left = N_ext;
  int N_notch_right = N_ext + N;
  int N_extrap = N_ext + N + N_ext;

  a_ext.resize(N_extrap);

  // Extrapolate the lower end.
  double da_left = a[1] - a[0];
  for (int i = 0; i < N_notch_left; i++) {
    a_ext[i] = a.front() + (i - N_notch_left) * da_left;
  }

  // Fill in the middle part.
  for (int i = N_notch_left; i < N_notch_right; i++) {
    a_ext[i] = a[i - N_notch_left];
  }

  // Extrapolate the upper end.
  double da_right = a[N - 1] - a[N - 2];
  for (int i = N_notch_right; i < N_extrap; i++) {
    a_ext[i] = a.back() + (i - (N_notch_right - 1)) * da_right;
  }
}

void extrap_loglin(
  std::vector<double>& a, int N_ext, std::vector<double>& a_ext
) {
  std::feclearexcept(FE_ALL_EXCEPT);

  // Prepare arrays.
  int N = a.size();
  int N_notch_left = N_ext;
  int N_notch_right = N_ext + N;
  int N_extrap = N_ext + N + N_ext;

  if (check_1d_array(a, false, false, true) != 0) {
    throw std::invalid_argument(
      "Sign change or zeros detected in the input array."
    );
  }

  a_ext.resize(N_extrap);

  // Extrapolate the lower end.
  double dlna_left = std::log(a[1] / a[0]);
  if (std::isnan(dlna_left)) {
    throw std::invalid_argument(
      "NaN detected at lower-end log-linear spacing."
    );
  }

  for (int i = 0; i < N_notch_left; i++) {
    a_ext[i] = a.front() * std::exp((i - N_notch_left) * dlna_left);
    if (std::fetestexcept(FE_DIVBYZERO) || std::fetestexcept(FE_INVALID)) {
      throw std::invalid_argument(
        "NaN detected at lower-end log-linear extrapolation."
      );
    }
  }

  // Fill in the middle part.
  for (int i = N_notch_left; i < N_notch_right; i++) {
    a_ext[i] = a[i - N_ext];
  }

  // Extrapolate the upper end.
  double dlna_right = std::log(a[N - 1] / a[N - 2]);
  if (std::isnan(dlna_right)) {
    throw std::invalid_argument(
      "NaN detected at upper-end log-linear spacing."
    );
  }

  for (int i = N_notch_right; i < N_extrap; i++) {
    a_ext[i] = a.back() * std::exp((i - (N_notch_right - 1)) * dlna_right);
    if (std::fetestexcept(FE_DIVBYZERO) || std::fetestexcept(FE_INVALID)) {
      throw std::invalid_argument(
        "NaN detected at upper-end log-linear extrapolation."
      );
    }
  }
}

void extrap_pad(
  std::vector<double>& a,
  int N_ext, double c_lower, double c_upper,
  std::vector<double>& a_ext
) {
  int N = a.size();
  int N_notch_left = N_ext;
  int N_notch_right = N_ext + N;
  int N_extrap = N_ext + N + N_ext;

  a_ext.resize(N_extrap);

  for (int i = 0; i < N_notch_left; i++) {
    a_ext[i] = c_lower;
  }
  for (int i = N_notch_left; i < N_notch_right; i++) {
    a_ext[i] = a[i - N_notch_left];
  }
  for (int i = N_notch_right; i < N_extrap; i++) {
    a_ext[i] = c_upper;
  }
}

void extrap2d_lin(
  std::vector< std::vector<double> >& a,
  int N_row_ext, int N_col_ext,
  std::vector< std::vector<double> >& a_ext
) {
  // FUTURE: Implement this function.
  throw trvs::UnimplementedError(
    "Function `extrap2d_lin` is not yet implemented in the C++ backend."
  );
}

void extrap2d_loglin(
  std::vector< std::vector<double> >& a,
  int N_row_ext, int N_col_ext,
  std::vector< std::vector<double> >& a_ext
) {
  // FUTURE: Implement this function.
  throw trvs::UnimplementedError(
    "Function `extrap2d_loglin` is not yet implemented in the C++ backend."
  );
}

void extrap2d_pad(
  std::vector< std::vector<double> >& a,
  int N_row_ext, int N_col_ext,
  double c_row_lower, double c_row_upper,
  double c_col_lower, double c_col_upper,
  std::vector< std::vector<double> >& a_ext
) {
  // FUTURE: Implement this function.
  throw trvs::UnimplementedError(
    "Function `extrap2d_pad` is not yet implemented in the C++ backend."
  );
}


// ***********************************************************************
// Sorting
// ***********************************************************************

std::vector<int> get_sorted_indices(std::vector<int> sorting_vector) {
  // Create an index vector to store the indices of the elements.
  std::vector<int> indices(sorting_vector.size());

#ifdef TRV_USE_OMP
#pragma omp parallel for simd
#endif  // TRV_USE_OMP
  for (int i = 0; i < int(sorting_vector.size()); i++) {
      indices[i] = i;
  }

  // Sort the index vector based on the sorting vector.
  std::sort(indices.begin(), indices.end(), [&](int a, int b) {
      return sorting_vector[a] < sorting_vector[b];;
  });

  return indices;
}

}  // namespace trv::array

}  // namespace trv
