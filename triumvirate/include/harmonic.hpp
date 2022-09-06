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
 * @file harmonic.hpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 * @brief (Reduced) spherical harmonic calculations.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_HARMONIC_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_HARMONIC_HPP_INCLUDED_

#include <cmath>
#include <complex>

#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>

#include "monitor.hpp"

namespace trv {
namespace maths {

/**
 * @brief Calculate Wigner 3-j symbol.
 *
 * @param j1, j2, j3, m1, m2, m3 Wigner 3-j symbol components.
 * @returns Value of the Wigner 3-j symbol.
 */
double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3);

/**
 * @brief Reduced spherical harmonics.
 *
 * The 'reduced' (conjugated and unit normalised) spherical harmonics are
 * given by
 *
 * @f[
 *   y_\ell^m = \sqrt(\frac{4*\pi}{2\ell + 1}) {Y_\ell^m}^\ast
 * @f]
 *
 * with f@$ y_0^0 = 1 f@$.
 */
class SphericalHarmonicCalculator {
 public:
  /**
   * @brief Calculate the reduced spherical harmonic.
   *
   * @param ell Degree @f$ ell @f$.
   * @param m Order @f$ m @f$.
   * @param pos 3-d position vector.
   * @returns Value of @f$ y_\ell^m f@$.
   */
  static std::complex<double> calc_reduced_spherical_harmonic(
    const int ell, const int m, double pos[3]
  );

  /**
   * @brief Store reduced spherical harmonics computed in Fourier space.
   *
   * @param[in] ell Degree @f$ ell @f$.
   * @param[in] m Order @f$ m @f$.
   * @param[in] boxsize Box size in each dimension.
   * @param[in] ngrid Grid number in each dimension.
   * @param[out] ylm_out Stored @f$ y_\ell^m f@$ values.
   * @throws trv::mon::InvalidData Exception thrown when the output variable
   *                          is not provided.
   */
  static void store_reduced_spherical_harmonic_in_fourier_space(
    const int ell, const int m,
    const double boxsize[3], const unsigned int ngrid[3],
    std::complex<double>* ylm_out
  );

  /**
   * @brief Store reduced spherical harmonics computed in
   *        configuration space.
   *
   * @param[in] ell Degree @f$ ell @f$.
   * @param[in] m Order @f$ m @f$.
   * @param[in] boxsize Box size in each dimension.
   * @param[in] ngrid Grid number in each dimension.
   * @param[out] ylm_out Stored @f$ y_\ell^m f@$ values.
   * @throws trv::mon::InvalidData Exception thrown when the output variable
   *                          is not provided.
   */
  static void store_reduced_spherical_harmonic_in_config_space(
    const int ell, const int m,
    const double boxsize[3], const unsigned int ngrid[3],
    std::complex<double>* ylm_out
  );
};

}  // namespace trv::maths
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_HARMONIC_HPP_INCLUDED_
