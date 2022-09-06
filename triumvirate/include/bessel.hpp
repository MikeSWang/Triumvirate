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
 * @file bessel.hpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 * @brief Spherical Bessel function calculations.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_BESSEL_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_BESSEL_HPP_INCLUDED_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>

namespace trv {
namespace maths {

/**
 * @brief Interpolated spherical Bessel function @f$ j_\ell(x) @f$
 *        of the first kind.
 *
 */
class SphericalBesselCalculator {
 public:
  /**
   * @brief Construct the interpolated function.
   *
   * @param ell Order @f$ \ell @f$.
   */
  SphericalBesselCalculator(const int ell);

  /**
   * @brief Destruct the interpolated function.
   */
  ~SphericalBesselCalculator();

  /**
   * @brief Evaluate the interpolated function.
   *
   * @param x Argument f@$ x f@$.
   * @returns Value of @f$ j_\ell @f$.
   */
  double eval(double x);

 private:
  gsl_interp_accel* accel;  ///< interpolation accelerator
  gsl_spline* spline;       ///< interpolation scheme
};

}  // namespace trv::maths
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_BESSEL_HPP_INCLUDED_
