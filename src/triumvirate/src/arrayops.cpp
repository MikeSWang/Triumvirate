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

namespace trv {

// ***********************************************************************
// Extrapolation
// ***********************************************************************

namespace sys {

ExtrapError::ExtrapError(const char* fmt_string, ...): std::runtime_error(
    "Extrapolation error."  // mandatory default error message
) {
  std::va_list args;

  char err_mesg_buf[4096];
  va_start(args, fmt_string);
  std::vsprintf(err_mesg_buf, fmt_string, args);
  va_end(args);

  this->err_mesg = std::string(err_mesg_buf);
}

const char* ExtrapError::what() const noexcept {return this->err_mesg.c_str();}

}  // namespace trv::sys

namespace array {

void extrap_loglin(double* a, int N, int N_ext, double* a_ext) {
  // Check for sign change or zero.
  double dlna_left = std::log(a[1] / a[0]);
  if (std::isnan(dlna_left)) {
    throw trv::sys::ExtrapError(
      "[ERRO] Sign change or zero detected in log-linear extrapolation "
      "at the lower end."
    );
  }

  double dlna_right = std::log(a[N - 1] / a[N - 2]);
  if (std::isnan(dlna_left)) {
    throw trv::sys::ExtrapError(
      "[ERRO] Sign change or zero detected in log-linear extrapolation "
      "at the upper end."
    );
  }

  // Extrapolate the lower end.
  if (a[0] > 0) {
    for (int i = 0; i < N_ext; i++) {
      a_ext[i] = a[0] * std::exp((i - N_ext) * dlna_left);
    }
  } else {
    for (int i = 0; i < N_ext; i++) {
      a_ext[i] = 0.;
    }
  }

  // Fill in the middle part.
  for (int i = N_ext; i < N_ext + N; i++) {
    a_ext[i] = a[i - N_ext];
  }

  // Extrapolate the upper end.
  if (a[N - 1] > 0) {
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      a_ext[i] = a[N - 1] * std::exp((i - N_ext - (N - 1)) * dlna_right);
    }
  } else {
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      a_ext[i] = 0.;
    }
  }
}

void extrap2d_logbilin(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
) {
  double dlna_left, dlna_right;
  for (int i = N_ext; i < N_ext + N; i++) {
    // Check for sign change or zero.
    dlna_left = std::log(a[i - N_ext][1] / a[i - N_ext][0]);
    if (std::isnan(dlna_left)) {
      throw trv::sys::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the left end."
      );
    }

    dlna_right = std::log(a[i - N_ext][N - 1] / a[i - N_ext][N - 2]);
    if (std::isnan(dlna_right)) {
      throw trv::sys::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the right end."
      );
    }

    // Extrapolate the middle left edge.
    for (int j = 0; j < N_ext; j++) {
      if (a[i - N_ext][0] > 0) {
        a_ext[i][j] = a[i - N_ext][0] * std::exp((j - N_ext) * dlna_left);
      } else {
        a_ext[i][j] = 0.;
      }
    }

    // Fill in the central part.
    for (int j = N_ext; j < N_ext + N; j++) {
      a_ext[i][j] = a[i - N_ext][j - N_ext];
    }

    // Extrapolate the middle right edge.
    for (int j = N_ext + N; j < N + 2*N_ext; j++) {
      if (a[i - N_ext][N - 1] > 0) {
        a_ext[i][j] = a[i - N_ext][N - 1] * std::exp(
          (j - N_ext - (N - 1)) * dlna_right
        );
      } else {
        a_ext[i][j] = 0.;
      }
    }
  }

  double dlna_up, dlna_down;
  for (int j = 0; j < N + 2*N_ext; j++) {
    // Check for sign change or zero.
    dlna_up = std::log(a_ext[N_ext + 1][j] / a_ext[N_ext][j]);
    if (std::isnan(dlna_up)) {
      throw trv::sys::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the top end."
      );
    }

    dlna_down = std::log(a_ext[N_ext + N - 1][j] / a_ext[N_ext + N - 2][j]);
    if (std::isnan(dlna_down)) {
      throw trv::sys::ExtrapError(
        "[ERRO] Sign change or zero detected in log-linear extrapolation "
        "at the bottom end."
      );
    }

    // Extrapolate the upper edge.
    for (int i = 0; i < N_ext; i++) {
      if (a_ext[N_ext][j] > 0) {
        a_ext[i][j] = a_ext[N_ext][j] * std::exp((i - N_ext) * dlna_up);
      } else {
        a_ext[i][j] = 0.;
      }
    }

    // Extrapolate the lower edge.
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      if (a_ext[N_ext + N - 1][j] > 0) {
        a_ext[i][j] = a_ext[N_ext + N - 1][j] * std::exp(
          (i - N_ext - (N - 1)) * dlna_down
        );
      } else {
        a_ext[i][j] = 0.;
      }
    }
  }
}

void extrap2d_bilin(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
) {
  double da_left, da_right;
  for (int i = N_ext; i < N + N_ext; i++) {
    da_left = a[i - N_ext][1] -  a[i - N_ext][0];
    da_right = a[i - N_ext][N - 1] - a[i - N_ext][N - 2];

    // Extrapolate the middle left edge.
    for (int j = 0; j < N_ext; j++) {
      a_ext[i][j] = a[i - N_ext][0] + (j - N_ext) * da_left;
    }

    // Fill in the central part.
    for (int j = N_ext; j < N_ext + N; j++) {
      a_ext[i][j] = a[i - N_ext][j - N_ext];
    }

    // Extrapolate the middle right edge.
    for (int j = N_ext + N; j < N + 2*N_ext; j++) {
      a_ext[i][j] = a[i - N_ext][N-1] + (j - N_ext - (N - 1)) * da_right;
    }
  }

  double da_up, da_down;
  for (int j = 0; j < N + 2*N_ext; j++) {
    da_up = a_ext[N_ext + 1][j] - a_ext[N_ext][j];
    da_down = a_ext[N_ext + N - 1][j] - a_ext[N_ext + N - 2][j];

    // Extrapolate the upper edge.
    for (int i =  0; i < N_ext; i++) {
      a_ext[i][j] = a_ext[N_ext][j] + (i - N_ext) * da_up;
    }

    // Extrapolate the lower edge.
    for (int i =  N_ext + N; i < 2*N_ext + N; i++) {
      a_ext[i][j] = a_ext[N_ext + N - 1][j] + (i - N_ext - (N - 1)) * da_down;
    }
  }
}

void extrap2d_bizeros(
  std::vector< std::vector<double> > a,
  int N, int N_ext,
  std::vector< std::vector<double> >& a_ext
) {
  for (int i = N_ext; i < N + N_ext; i++) {
    for (int j = 0; j < N_ext; j++) {
      a_ext[i][j] = 0.;
    }
    for (int j = N_ext; j < N_ext + N; j++) {
      a_ext[i][j] = a[i - N_ext][j - N_ext];
    }
    for (int j = N_ext + N; j < N + 2*N_ext; j++) {
      a_ext[i][j] = 0.;
    }
  }

  for (int j = 0; j < N + 2*N_ext; j++) {
    for (int i = 0; i < N_ext; i++) {
      a_ext[i][j] = 0.;
    }
    for (int i = N_ext + N; i < N + 2*N_ext; i++) {
      a_ext[i][j] = 0.;
    }
  }
}

}  // namespace trv::array
}  // namespace trv
