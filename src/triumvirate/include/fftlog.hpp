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
#include <cstring>
#include <vector>

#include "maths.hpp"
#include "arrayops.hpp"

namespace trva = trv::array;

namespace trv {

namespace maths {

/**
 * @brief Perform the (forward biased) Hankel and associated transforms
 *        using the FFTLog algorithm.
 *
 */
class HankelTransform {
 public:
  double order;         ///< order of the Hankel transform
  double bias;          ///< power-law bias index
  int nsamp = 0;        ///< number of samples provided
  int nsamp_trans = 0;  ///< number of samples transformed
  double logres = 0.;   ///< logarithmic interval sample spacing
  double pivot = 1.;    ///< pivot value

  /// logarithmically linearly-spaced sample points pre-transform
  std::vector<double> pre_sampts;

  /// logarithmically linearly-spaced sample points post-transform
  std::vector<double> post_sampts;

  /**
   * @brief Construct the Hankel transform.
   *
   * @param mu Order of the Hankel transform.
   * @param q Power-law bias index.
   * @param threaded If `true` (default), use multi-threads FFT.
   */
  HankelTransform(double mu, double q, bool threaded = true);

  /**
   * @brief Destruct the Hankel transform.
   */
  ~HankelTransform();

  /**
   * @brief Initialise the Hankel transform.
   *
   * This sets up the pre- and post-transform samples as well as the
   * transform kernel and the pivot value.
   *
   * @param sample_pts Logarithmically linearly-spaced sample points.
   *                   Must be even in length if extrapolation is enabled.
   * @param kr_c Pivot value.
   * @param lowring If true (default), set the pivot value by the
   *                low-ringing condition.
   * @param extrap Extrapolation option.  If not
   *               @ref trv::array::ExtrapOption::NONE (default),
   *               the sample size for the transform is the smallest
   *               power of 2 that is greater than or equal to
   *               `extrap_exp` times the original number of
   *               sample points; the pre-transform samples are assumed
   *               to be real and must be even in length.
   * @param extrap_exp Sample size expansion factor (default is 2.) for
   *                   extrapolation.  The smallest power of 2 greater
   *                   than or equal to this times the original number
   *                   of sample points is used as the sample size for
   *                   the transform.
   * @throws trv::sys::InvalidParameterError When the size of `sample_pts`
   *                                         is less than 2.
   * @throws trv::sys::InvalidParameterError When `sample_pts` is not
   *                                         logarithmically spaced.
   * @throws trv::sys::InvalidParameterError When the size of `sample_pts`
   *                                         is not even with
   *                                         extrapolation enabled.
   * @throws trv::sys::InvalidParameterError When `extrap_exp` results in
   *                                         a shrunken sample size.
   * @throws trv::sys::InvalidParameterError When `kr_c` is non-positive
   *                                         while `lowring` is false.
   */
  void initialise(
    std::vector<double> sample_pts, double kr_c,
    bool lowring = true,
    trva::ExtrapOption extrap = trva::ExtrapOption::NONE,
    double extrap_exp = 2.
  );

  /**
   * @brief Initialise the Hankel transform.
   *
   * This sets up the pre- and post-transform samples as well as the
   * transform kernel and the pivot value.
   *
   * @param sample_pts Logarithmically linearly-spaced sample points.
   *                   Must be even in length if extrapolation is enabled.
   * @param kr_c Pivot value.
   * @param lowring If true (default), set the pivot value by the
   *                low-ringing condition.
   * @param extrap Extrapolation option.  If not 0 (default),
   *               the sample size for the transform is the smallest
   *               power of 2 that is greater than or equal to
   *               `extrap_exp` times the original number of
   *               sample points; the pre-transform samples are assumed
   *               to be real and must be even in length.
   * @param extrap_exp Sample size expansion factor (default is 2.) for
   *                   extrapolation.  The smallest power of 2 greater
   *                   than or equal to this times the original number
   *                   of sample points is used as the sample size for
   *                   the transform.
   *
   * @overload
   */
  void initialise(
    std::vector<double> sample_pts, double kr_c,
    bool lowring = true,
    int extrap = 0,
    double extrap_exp = 2.
  );

  /**
   * @brief Perform the (forward biased) Hankel transform.
   *
   * The transform is defined here as
   * @f[
   *   b(k) = \int_0^\infty k \mathrm{d}r \, (k r)^q J_\mu(k r) a(r) \,.
   * @f]
   *
   * @param[in] a Pre-transform sample values. Must be even in length
   *              if extrapolation is enabled.
   * @param[out] b Post-trasform sample values.
   *
   * @attention If extrapolation is enabled (`extrap` not
   *            @ref trv::array::ExtrapOption::NONE, the pre-transform
   *            samples are assumed to be real.
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
  /// extrapolation option (default is none)
  trva::ExtrapOption extrap = trva::ExtrapOption::NONE;

  /// number of extra sample points on either side
  int n_ext = 0;

  /// pre-transform sample points with extrapolation option
  std::vector<double> pre_sampts_extrap;

  /// post-transform sample points with extrapolation option
  std::vector<double> post_sampts_extrap;

  /// FFTLog transform kernel coefficients
  std::vector< std::complex<double> > kernel;

  /// pre-kernel FFTW plan and array
  fftw_plan pre_plan;
  fftw_complex* pre_buffer = nullptr;

  /// post-kernel FFTW plan and array
  fftw_plan post_plan;
  fftw_complex* post_buffer = nullptr;

  /// FFTW plan initialisation flag
  bool plan_init = false;

  /// FFTW multi-threading flag
  bool threaded = true;

  /**
   * @brief Reset FFTW plans and buffers.
   *
   */
  void reset();
};

class SphericalBesselTransform: public HankelTransform {
 public:
  int degree;  ///< degree of the spherical Bessel transform

  /**
   * @brief Construct the spherical Bessel transform.
   *
   * @param ell Degree of the spherical Bessel transform.
   * @param n Power-law bias index.
   * @param threaded If `true` (default), use multi-threads FFT.
   */
  SphericalBesselTransform(int ell, int n, bool threaded = true);

  /**
   * @brief Initialise the spherical Bessel transform.
   *
   * This sets up the underlying Hankel transform with order
   * @f$ \mu = \ell + 1/2 @f$.  The low-ringing pivot value is enforced
   * from an initial value of 1.
   *
   * @param sample_pts Logarithmically linearly-spaced sample points.
   *                   Must be even in length if extrapolation is enabled.
   * @param kr_c Pivot value.
   * @param lowring If true (default), set the pivot value by the
   *                low-ringing condition.
   * @param extrap Extrapolation option.  If not
   *               @ref trv::array::ExtrapOption::NONE (default),
   *               the sample size for the transform is the smallest
   *               power of 2 that is greater than or equal to
   *               `extrap_exp` times the original number of
   *               sample points; the pre-transform samples are assumed
   *               to be real and must be even in length.
   * @param extrap_exp Sample size expansion factor (default is 2) for
   *                   extrapolation.  The smallest power of 2 greater
   *                   than or equal to this times the original number
   *                   of sample points is used as the sample size for
   *                   the transform.
   */
  void initialise(
    std::vector<double> sample_pts, double kr_c,
    bool lowring = true,
    trva::ExtrapOption extrap = trva::ExtrapOption::NONE,
    double extrap_exp = 2.
  );

  /**
   * @brief Initialise the spherical Bessel transform.
   *
   * This sets up the underlying Hankel transform with order
   * @f$ \mu = \ell + 1/2 @f$.  The low-ringing pivot value is enforced
   * from an initial value of 1.
   *
   * @param sample_pts Logarithmically linearly-spaced sample points.
   *                   Must even in length if extrapolation is enabled.
   * @param kr_c Pivot value.
   * @param lowring If true (default), set the pivot value by the
   *                low-ringing condition.
   * @param extrap Extrapolation option.  If not 0 (default),
   *               the sample size for the transform is the smallest
   *               power of 2 that is greater than or equal to
   *               `extrap_exp` times the original number of
   *               sample points; the pre-transform samples are assumed
   *               to be real and must be even in length.
   * @param extrap_exp Sample size expansion factor (default is 2) for
   *                   extrapolation.  The smallest power of 2 greater
   *                   than or equal to this times the original number
   *                   of sample points is used as the sample size for
   *                   the transform.
   *
   * @overload
   */
  void initialise(
    std::vector<double> sample_pts, double kr_c,
    bool lowring = true,
    int extrap = 0,
    double extrap_exp = 2.
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
   *       @f$ A(r) = r^{3/2} a(r) @f$ and
   *       @f$ B(k) = (2\pi / k)^{3/2} b(k) @f$
   *       with @f$ \mu = \ell + 1/2 @f$ and the same @f$ q @f$.
   *
   * @param[in] a Pre-transform sample values.  Must even in length
   *              if extrapolation is enabled.
   * @param[out] b Post-transform sample values.
   *
   * @attention If extrapolation is enabled by `extrap`, the pre-transform
   *            samples are assumed to be real.
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
   * @param[in] pre_samples Pre-transform multipole samples.  Must be even
   *                        in length if extrapolation is enabled.
   * @param[out] post_samples Post-transform multipoles samples.
   *
   * @attention If extrapolation is enabled (`extrap` not
   *            @ref trv::array::ExtrapOption::NONE, the pre-transform
   *            samples are assumed to be real.
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
