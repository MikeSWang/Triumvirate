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
 * @file fftlog.hpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 * @brief FFTLog algorithm for Hankel and related transforms and its
 *        application to cosmological functions.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_

#include <fftw3.h>

#include <cmath>
#include <complex>

#include "maths.hpp"
#include "arrayops.hpp"

namespace trv {

namespace maths {

/**
 * @brief Calculate a low-ringing pivot value @f$ k r = k_c r_c @f$
 *        for the FFTLog transform.
 *
 * @param[in] mu Order of the Hankel transform.
 * @param[in] q FFTLog power-law bias index.
 * @param[in] L Logarithmic interval.
 * @param[in] N Sample number.
 * @param[out] kr_c Initial pivot value.
 * @returns Low-ringing pivot value.
 */
double calc_kr_pivot_lowring(
  double mu, double q, double L, int N, double kr_c=1.
);

/**
 * @brief Compute the FFTLog transform kernel coefficients @f$ u @f$.
 *
 * @param[in] mu Order of the Hankel transform.
 * @param[in] q FFTLog power-law bias index.
 * @param[in] L Logarithmic interval.
 * @param[in] N Sample number.
 * @param[in] kr_c Pivot value.
 * @param[out] u Kernel coefficients.
 */
void compute_u_kernel_coeff(
  double mu, double q, double L, int N, double kr_c, std::complex<double>* u
);

/**
 * @brief Perform the (forward) Hankel transform.
 *
 * The transform is defined here as
 * @f[
 *   b(k) = \int_0^\infty k \mathrm{d}r \, (k r)^q J_\mu(k r) a(r) \,.
 * @f]
 *
 * @param[in] mu Order of the Hankel transform.
 * @param[in] q FFTLog power-law bias index.
 * @param[in] N Sample number.
 * @param[in] kr_c (Initial) pivot value.
 * @param[in] lowring Boolean low-ringing condition.
 * @param[in] r Pre-transform sample points.
 * @param[in] a Pre-transform sample values.
 * @param[out] k Post-transform sample points.
 * @param[out] b Post-trasform sample values.
 * @param[out] u FFTLog transform kernel coefficients.
 */
void hankel_transform(
  double mu, double q, double kr_c, int N, bool lowring,
  double* r, const std::complex<double>* a,
  double* k, std::complex<double>* b,
  std::complex<double>* u
);

/**
 * @brief Perform the (forward) spherical Fourier--Bessel transform.
 *
 * The transform is defined here as
 * @f[
 *   b_\ell(k) = (2\pi k)^{-3/2} \int_0^\infty k \mathrm{d}r \,
 *     r^{m - 1/2} J_\ell(k r) a_\ell(r) \,.
 * @f]
 *
 * @note When ``m = 2``, this is the Hankel transform definition used in
 *       cosmological correlators,
 *       @f[
 *         \xi_\ell(r) = 4\pi \mathrm{i}^\ell \int_0^\infty
 *           \frac{\mathrm{d}k \, k^2}{(2\pi)^3}
 *           j_\ell(kr) P_\ell(k) \,,
 *       @f]
 *       i.e. @f$ a(r) @f$ corresponds to @f$ P_\ell(k) @f$ and
 *       @f$ b(k) @f$ corresponds to @f$ \mathrm{i}^{-\ell} \xi_\ell(r)@f$.
 *
 * @param[in] ell Order of the transform.
 * @param[in] m Dimensional power-law index.
 * @param[in] N Sample number.
 * @param[in] r Pre-transform sample points.
 * @param[in] a Pre-transform sample values.
 * @param[out] k Post-transform sample points.
 * @param[out] b Post-transform sample values.
 */
void sj_transform(
  int ell, int m, int N, double* r, double* a, double* k, double* b
);

/**
 * @brief Perform the (forward) biased symmetric spherical
 *        Bessel transform.
 *
 * The transform is defined here as
 * @f[
 *   b_\ell(k) = (2\pi k)^{-3/2} \int_0^\infty k \mathrm{d}r \,
 *     r^{3/2 + i} J_\ell(k r) a_\ell(r) \,.
 * @f]
 *
 * @note This is the Hankel transform definition used in cosmological
 *       wide-angle expansions,
 *       @f[
 *         \xi_\ell^{(n)}(r) = 4\pi \mathrm{i}^\ell \int_0^\infty
 *           \frac{\mathrm{d}k \, k^2}{(2\pi)^3}
 *           (k r)^{-n} j_\ell(kr) P_\ell^{(n)}(k) \,,
 *       @f]
 *       i.e. @f$ a(r) @f$ corresponds to @f$ P_\ell^{(n)}(k) @f$ and
 *       @f$ b(k) @f$ corresponds to
 *       @f$ \mathrm{i}^{-\ell} \xi_\ell^{(n)}(r)@f$ with @f$ n @f$
 *       identified with `i`.
 *
 * @param[in] ell Order of the transform.
 * @param[in] i Power-law bias index.
 * @param[in] N Sample number.
 * @param[in] r Pre-transform sample points.
 * @param[in] a Pre-transform sample values.
 * @param[out] k Post-transform sample points.
 * @param[out] b Post-transform sample values.
 */
void sj_transform_symm_biased(
  int ell, int i, int N, double* r, double* a, double* k, double* b
);

}  // namespace trv::maths

/**
 * @brief Transform power spectrum to correlation function
 *        multipole samples.
 *
 * @param[in] ell Multipole degree.
 * @param[in] N Sample number.
 * @param[in] k Wavenumber sample points.
 * @param[in] pk Power spectrum samples.
 * @param[out] r Separation sample points.
 * @param[out] xi Correlation function samples.
 */
void transform_powspec_to_corrfunc_multipole(
  int ell, int N, double* k, double* pk, double* r, double* xi
);

/**
 * @brief Transform correlation function to power spectrum
 *        multipole samples.
 *
 * @param[in] ell Multipole degree.
 * @param[in] N Sample number.
 * @param[in] r (Pointer) Separation sample points.
 * @param[in] xi (Pointer) Correlation function samples.
 * @param[out] k (Pointer) Wavenumber sample points.
 * @param[out] pk (Pointer) Power spectrum samples.
 */
void transform_corrfunc_to_powspec_multipole(
  int ell, int N, double* r, double* xi, double* k, double* pk
);

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_
