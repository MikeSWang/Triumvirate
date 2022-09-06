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
 * @file bessel.cpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 *
 */

#include "bessel.hpp"

namespace trv {
namespace maths {

SphericalBesselCalculator::SphericalBesselCalculator(const int ell) {
  /// Set up sampling range and number.
  /// CAVEAT: Discretionary choices.
  const double xmin = 0.;       ///< minimum of interpolation range
  const double xmax = 10000.;   ///< maximum of interpolation range
  const int nsample = 1000000;  ///< interpolation sample number

  /// Initialise and evaluate at sample points.
  double dx = (xmax - xmin) / (nsample - 1);

  double* x = new double[nsample];
  double* j_ell = new double[nsample];
  for (int i = 0; i < nsample; i++) {
    x[i] = xmin + dx * i;
    j_ell[i] = gsl_sf_bessel_jl(ell, x[i]);
  }

  /// Initialise the interpolator using cubic spline and the accelerator.
  this->accel = gsl_interp_accel_alloc();
  this->spline = gsl_spline_alloc(gsl_interp_cspline, nsample);

  gsl_spline_init(this->spline, x, j_ell, nsample);

  /// Delete sample points.
  delete[] x; x = nullptr;
  delete[] j_ell; j_ell = nullptr;
}

SphericalBesselCalculator::~SphericalBesselCalculator() {
  if (this->accel != nullptr) {
    gsl_interp_accel_free(this->accel); this->accel = nullptr;
  }

  if (this->spline != nullptr) {
    gsl_spline_free(this->spline); this->spline = nullptr;
  }
}

double SphericalBesselCalculator::eval(double x) {
  return gsl_spline_eval(this->spline, x, this->accel);
}

}  // namespace trv::maths
}  // namespace trv
