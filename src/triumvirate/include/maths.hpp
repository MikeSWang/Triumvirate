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
 * @file maths.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief Mathematical calculations.
 *
 * Mathematical calculations provided include:
 * - spherical Bessel functions of the first kind with interpolation;
 * - (reduced) spherical harmonics include 3-d mesh grid storage;
 * - Wigner 3-j symbols;
 * - the gamma function and related quantities with Lanzcos approximation.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_MATHS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_MATHS_HPP_INCLUDED_

#include <gsl/gsl_interp.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_spline.h>

#include <cmath>
#include <complex>
#include <vector>

#include "monitor.hpp"

namespace trv {
namespace maths {

// ***********************************************************************
// Complex numbers
// ***********************************************************************

extern const std::complex<double> M_I;  ///< imaginary unit

/**
 * @brief Evaluate a complex number @f$ r e^{i \theta} @f$
 *        in the polar form.
 *
 * @param r Modulus @f$ r @f$.
 * @param theta Argument @f$ \theta @f$.
 * @returns Value of the complex number.
 */
std::complex<double> eval_complex_in_polar(double r, double theta);


// ***********************************************************************
// Vectors
// ***********************************************************************

/**
 * @brief Return the magnitude of a 3-d vector.
 *
 * @param vec A 3-d vector.
 * @returns Vector magnitude.
 */
double get_vec3d_magnitude(std::vector<double> vec);

/**
 * @brief Return the magnitude of a 3-d vector.
 *
 * @param vec A 3-d vector.
 * @returns Vector magnitude.
 *
 * @overload
 */
double get_vec3d_magnitude(double* vec);


// ***********************************************************************
// Gamma function
// ***********************************************************************

/**
 * @brief Evaluate the Lanczos approximation series @f$ A_g(z) @f$
 *        for the gamma function.
 *
 * The series approximate the gamma function by
 * @f[
 *   \Gamma(z + 1) = \sqrt{2\pi}
 *     (z + g + 1/2)^{z + 1/2} \mathrm{e}^{- z - g - 1/2} A_g(z)
 * @f]
 * where @f$ g @f$ is the Lanczos constant.
 *
 * @param z Complex argument.
 * @returns Value of @f$ A_g(z) @f$.
 */
std::complex<double> eval_lanczos_approx_series(std::complex<double> z);

/**
 * @brief Evaluate the gamma function @f$ \Gamma(z) @f$ on
 *        the complex plane using the Lanczos approximation.
 *
 * @param z Complex argument.
 * @returns Value of @f$ \Gamma(z) @f$.
 *
 * @see trv::eval_lanczos_approx_series()
 */
std::complex<double> eval_gamma(std::complex<double> z);

/**
 * @brief Evaluate the log-gamma function @f$ \ln\Gamma(z) @f$
 *        on the complex plane using the Lanczos approximation.
 *
 * @param z Complex argument.
 * @returns Value of @f$ \ln\Gamma(z) @f$.
 *
 * @see trv::eval_lanczos_approx_series()
 */
std::complex<double> eval_lngamma(std::complex<double> z);

/**
 * @brief Evaluate the logarithm of the ratio of two gamma functions.
 *
 * The ratio is of the form
 * @f$
 *   \Gamma\big(\frac{\nu + \mu + 1}{2}\big)
 *     \big/ \Gamma\big(\frac{\nu - \mu + 1}{2}\big)
 * @f$
 * in the asymptotic limit @f$ \mu, \nu \to \infty @f$.
 *
 * @param mu Real variable.
 * @param nu Complex variable.
 * @returns Approximate logarithmic ratio.
 */
std::complex<double> eval_gamma_ratio_asymp(
  double mu, std::complex<double> nu
);

/**
 * @brief Return the real and imaginary parts of the log-gamma function.
 *
 * @param[in] x Real part of the complex argument.
 * @param[in] y Imaginary part of the complex argument.
 * @param[out] lnr Real part of the log-gamma function value.
 * @param[out] theta Imaginary part of the log-gamma function value.
 */
void get_lngamma_components(double x, double y, double& lnr, double& theta);


// ***********************************************************************
// Spherical harmonics
// ***********************************************************************

/// zero-tolerance for Wigner 3-j coupling coefficients
extern const double eps_coupling;

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
 * @f[
 *   y_\ell^m = \sqrt{\frac{4\pi}{2\ell + 1}} {Y_\ell^m}^\ast
 * @f]
 * with @f$ y_0^0 = 1 @f$.
 */
class SphericalHarmonicCalculator {
 public:
  /**
   * @brief Calculate the reduced spherical harmonic.
   *
   * @param ell Degree @f$ \ell @f$.
   * @param m Order @f$ m @f$.
   * @param pos 3-d position vector.
   * @returns Value of @f$ y_\ell^m @f$.
   */
  static std::complex<double> calc_reduced_spherical_harmonic(
    const int ell, const int m, double pos[3]
  );

  /**
   * @brief Store reduced spherical harmonics computed in Fourier space.
   *
   * @param[in] ell Degree @f$ \ell @f$.
   * @param[in] m Order @f$ m @f$.
   * @param[in] boxsize Box size in each dimension.
   * @param[in] ngrid Grid number in each dimension.
   * @param[out] ylm_out Stored @f$ y_\ell^m @f$ values.
   * @throws trv::sys::InvalidDataError When the output variable is
   *                                    not provided.
   */
  static void store_reduced_spherical_harmonic_in_fourier_space(
    const int ell, const int m,
    const double boxsize[3], const int ngrid[3],
    std::vector< std::complex<double> >& ylm_out
  );

  /**
   * @brief Store reduced spherical harmonics computed in
   *        configuration space.
   *
   * @param[in] ell Degree @f$ \ell @f$.
   * @param[in] m Order @f$ m @f$.
   * @param[in] boxsize Box size in each dimension.
   * @param[in] ngrid Grid number in each dimension.
   * @param[out] ylm_out Stored @f$ y_\ell^m @f$ values.
   * @throws trv::sys::InvalidDataError When the output variable
   *                                    is not provided.
   */
  static void store_reduced_spherical_harmonic_in_config_space(
    const int ell, const int m,
    const double boxsize[3], const int ngrid[3],
    std::vector< std::complex<double> >& ylm_out
  );
};


// ***********************************************************************
// Spherical Bessel function
// ***********************************************************************

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
   * @param x Argument @f$ x @f$.
   * @returns Value of @f$ j_\ell @f$.
   */
  double eval(double x);

 private:
  gsl_interp_accel* accel;  ///< interpolation accelerator
  gsl_spline* spline;       ///< interpolation scheme
};

}  // namespace trv::maths
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_MATHS_HPP_INCLUDED_
