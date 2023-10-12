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
#include <vector>

#include "maths.hpp"
#include "arrayops.hpp"

namespace trv {

namespace maths {

/**
 * @brief Perform the (forward biased) Hankel and associated transforms
 *        using the FFTLog algorithm.
 *
 */
class HankelTransform {
 public:
  double order;        ///< order of the Hankel transform
  double bias;         ///< power-law bias index
  int nsamp = 0;       ///< number of samples
  double logres = 0.;  ///< logarithmic interval sample spacing
  double pivot = 1.;   ///< pivot value

  /// logarithmically linearly-spaced sample points pre-transform
  std::vector<double> pre_sampts;

  /// logarithmically linearly-spaced sample points post-transform
  std::vector<double> post_sampts;

  /**
   * @brief Construct the Hankel transform.
   *
   * @param mu Order of the Hankel transform.
   * @param q Power-law bias index.
   */
  HankelTransform(double mu, double q);

  /**
   * @brief Initialise the Hankel transform.
   *
   * This sets up the pre- and post-transform samples as well as the
   * transform kernel and the pivot value.
   *
   * @param sample_pts Logarithmically linearly-spaced sample points.
   * @param kr_c Pivot value.
   * @param lowring If true (default), set the pivot value by the
   *                low-ringing condition.
   * @throws trv::sys::InvalidParameterError When the size of `sample_pts`
   *                                         is less than 2.
   * @throws trv::sys::InvalidParameterError When `kr_c` is non-positive
   *                                         while `lowring` is false.
   */
  void initialise(
    std::vector<double> sample_pts, double kr_c, bool lowring = true
  );

  /**
   * @brief Perform the (forward biased) Hankel transform.
   *
   * The transform is defined here as
   * @f[
   *   b(k) = \int_0^\infty k \mathrm{d}r \, (k r)^q J_\mu(k r) a(r) \,.
   * @f]
   *
   * @param[in] a Pre-transform sample values.
   * @param[out] b Post-trasform sample values.
   */
  void biased_transform(std::complex<double>* a, std::complex<double>* b);

  /**
   * @brief Compute the FFTLog transform kernel coefficients @f$ u @f$.
   *
   * @returns Kernel coefficients.
   */
  std::vector< std::complex<double> > compute_kernel_coeff();

  /**
   * @brief Calculate a low-ringing FFTLog transform pivot value
   *        @f$ k r = k_c r_c @f$.
   *
   * @param delta Logarithmic-interval sample spacing.
   * @param kr_c Initial pivot value (default is 1.).  If zero, it is
   *             calculated directly; otherwise, it is adjusted from
   *             the initial value.
   * @returns Low-ringing pivot value.
   *
   * @note This is not bound to the logarithmic-interval sample spacing
   *       initialised.
   */
  double calc_lowring_pivot(double delta, double kr_c = 1.);

 private:
  /// FFTLog transform kernel coefficients
  std::vector< std::complex<double> > kernel;
};

class SphericalBesselTransform: public HankelTransform {
 public:
  int degree;  ///< degree of the spherical Bessel transform

  /**
   * @brief Construct the spherical Bessel transform.
   *
   * @param ell Degree of the spherical Bessel transform.
   * @param n Power-law bias index.
   */
  SphericalBesselTransform(int ell, int n);

  /**
   * @brief Initialise the spherical Bessel transform.
   *
   * This sets up the underlying Hankel transform with order
   * @f$ \mu = \ell + 1/2 @f$.  The low-ringing pivot value is enforced
   * from an initial value of 1.
   *
   * @param sample_pts Logarithmically linearly-spaced sample points.
   * @param kr_c Pivot value.
   * @param lowring If true (default), set the pivot value by the
   *                low-ringing condition.
   */
  void initialise(
    std::vector<double> sample_pts, double kr_c, bool lowring = true
  );

  /**
   * @brief Perform the (forward biased) spherical Bessel transform.
   *
   * The transform is defined here as
   * @f[
   *   b(k) = 4\pi \int_0^\infty r^2 \mathrm{d}r \,
   *     (k r)^q j_\ell(k r) a(r) \,.
   * @f]
   *
   * @note This is equivalent to the (forward biased) Hankel transform for
   *       @f$ A(r) = r^{3/2} a(r) @f$ and $B(k) = (2\pi / k)^{3/2} b(k)$
   *       with @f$ \mu = \ell + 1/2 @f$ and the same @f$ q @f$.
   *
   * @param[in] a Pre-transform sample values.
   * @param[out] b Post-transform sample values.
   */
  void biased_transform(
    std::vector< std::complex<double> >& a,
    std::vector< std::complex<double> >& b
  );

  /**
   * @brief Transform csomological multipole samples.
   *
   * @param[in] dir Transform direction: +1 (forward) for configuration to
   *                Fourier space, -1 (backward) for Fourier to
   *                configuration space.
   * @param[in] pre_samples Pre-transform multipole samples.
   * @param[out] post_samples Post-transform multipoles samples.
   */
  void transform_cosmological_multipole(
    int dir,
    std::vector< std::complex<double> >& pre_samples,
    std::vector< std::complex<double> >& post_samples
  );
};

}  // namespace trv::maths

}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_
