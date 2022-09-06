/**
 * @file fftlog.hpp
 * @brief FFTLog algorithm for Hankel and related transforms and its
 *        application to cosmological functions.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_

#include <cmath>
#include <complex>

#include <fftw3.h>

#include "tools.hpp"
#include "gamma.hpp"

using namespace trv::maths;

namespace trv {
namespace ops {

/**
 * Calculate a low-ringing pivot value f@$ k r = k_c r_c f@$ for the
 * FFTLog Hankel transform.
 *
 * @param[in] mu Order of the Hankel transform.
 * @param[in] q FFTLog power-law bias index.
 * @param[in] L Logarithmic interval.
 * @param[in] N Sample number.
 * @param[out] kr_c Initial pivot value.
 */
double calc_kr_pivot_lowring(
  double mu, double q, double L, int N, double kr_c=1.
) {
  double x_p = (mu + 1 + q)/2;
  double x_m = (mu + 1 - q)/2;
  double y = M_PI * N / (2*L);

  /// Note that no minus is involved by complex conjugation of
  /// the gamma function.
  double _lnr, phi_p, phi_m;  // placeholder `_lnr`
  get_lngamma_components(x_p, y, _lnr, phi_p);
  get_lngamma_components(x_m, y, _lnr, phi_m);

  double arg_ = std::log(2/kr_c) * N / L + (phi_p + phi_m) / M_PI;
  double arg_int = std::round(arg_);
  if (arg_ != arg_int) {
    kr_c *= std::exp(L / N * (arg_ - arg_int));
  }

  return kr_c;
}

/**
 * Compute the FFTLog Hankel transform kernel coefficients f@$ u f@$.
 *
 * @param[in] mu Order of the Hankel transform.
 * @param[in] q FFTLog power-law bias index.
 * @param[in] L Logarithmic interval.
 * @param[in] N Sample number.
 * @param[in] kr_c Pivot value.
 * @param[out] u Kernel coefficients.
 */
void compute_u_kernel_coeff(
  double mu, double q, double L, int N, double kr_c, std::complex<double> u[]
) {
  double y = M_PI / L;
  double kr_0 = kr_c * std::exp(-L);

  double t = -2 * y * std::log(kr_0/2);

  if (q == 0) {
    double x = (mu + 1)/2;

    double lnr, phi;
    for (int m = 0; m <= N/2; m++) {
      get_lngamma_components(x, m*y, lnr, phi);
      u[m] = eval_complex_in_polar(1., m*t + 2*phi);
    }
  } else {
    double x_p = (mu + 1 + q)/2;
    double x_m = (mu + 1 - q)/2;

    double lnr_p, phi_p;
    double lnr_m, phi_m;
    for (int m = 0; m <= N/2; m++) {
      get_lngamma_components(x_p, m*y, lnr_p, phi_p);
      get_lngamma_components(x_m, m*y, lnr_m, phi_m);

      u[m] = eval_complex_in_polar(
        std::exp(std::log(2) * q + lnr_p - lnr_m), m*t + phi_p - phi_m
      );
    }
  }

  for (int m = N/2 + 1; m < N; m++) {
    u[m] = conj(u[N - m]);
  }

  if (N % 2 == 0) {
    u[N/2] = u[N/2].real() + M_I * 0.;  // make real by log-periodicity
  }
}

/**
 * Perform the (forward) Hankel transform, defined here as
 * f@$
 *   b(k) = \int_0^\infty k \mathrm{d}r \, (k r)^q J_\mu(k r) a(r) \,.
 * f@$
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
  const double r[], std::complex<double> a[],
  double k[], std::complex<double> b[],
  std::complex<double>* u
) {
  /// Calculate the logarithmic interval.
  double L = N * std::log(r[N - 1] / r[0]) / (N - 1.);

  /// Compute the forward transform kernel.
  std::complex<double>* u_ = NULL;
  if (u == NULL) {
    u_ = new std::complex<double>[N];

    if (lowring) {
      kr_c = calc_kr_pivot_lowring(mu, q, L, N, kr_c);
    }

    compute_u_kernel_coeff(mu, q, L, N, kr_c, u_);
    u = u_;
  }

  /// Compute output sample points corresponding to the input sample points.
  double kr_0 = kr_c * std::exp(-L);

  k[0] = kr_0 / r[0];
  for (int j = 1; j < N; j++) {
    k[j] = k[0] * std::exp(j * L / N);
  }

  /// Compute the convolution b = a * u using FFT.
  fftw_plan forward_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) a, (fftw_complex*) b, -1, FFTW_ESTIMATE
  );
  fftw_plan reverse_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) b, (fftw_complex*) b, +1, FFTW_ESTIMATE
  );
    /** ``(`` and ``)`` necessary
        (see https://www.fftw.org/doc/Complex-numbers.html) */

  fftw_execute(forward_plan);
  for (int m = 0; m < N; m++) {
    b[m] *= u[m] / double(N);  // divide by `N` to normalise the inverse DFT
  }
  fftw_execute(reverse_plan);

  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(reverse_plan);

  /// Reverse the array `b`.
  std::complex<double> b_;
  for (int n = 0; n < N/2; n++) {
    b_ = b[n];
    b[n] = b[N - n - 1];
    b[N - n - 1] = b_;
  }

  delete[] u_;
}

/**
 * Perform the (forward) spherical Fourier--Bessel transform,
 * defined here as
 * f@$
 *   b_\ell(k) = (2\pi k)^{-3/2} \int_0^\infty k \mathrm{d}r \,
 *     r^{m - 1/2} J_\ell(k r) a_\ell(r) \,.
 * f@$
 *
 * @note When ``m = 2``, this is the Hankel transform definition used in
 *       cosmological correlators,
 *       f@$
 *         \xi_\ell(r) = 4\pi \mathrm{i}^\ell \int_0^\infty
 *           \frac{\mathrm{d}k \, k^2}{(2\pi)^3}
 *           j_\ell(kr) P_\ell(k) \,,
 *       f@$ i.e. f@$ a(r) f@$ corresponds to f@$ P_\ell(k) f@$ and
 *       f@$ b(k) f@$ corresponds to f@$ \mathrm{i}^{-\ell} \xi_\ell(r)f@$.
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
) {
  double mu = ell + 1./2;
  double q = 0.;
  double kr_c = 1.;
  bool lowring = true;

  std::complex<double>* A = new std::complex<double>[N];
  std::complex<double>* B = new std::complex<double>[N];

  for (int j = 0; j < N; j++) {
    A[j] = std::pow(r[j], m - 1./2) * a[j];
  }

  hankel_transform(mu, q, kr_c, N, lowring, r, A, k, B, NULL);

  for (int j = 0; j < N; j++) {
    std::complex<double> bj = std::pow(2*M_PI * k[j], -3./2) * B[j];
    b[j] = bj.real();
  }

  delete[] A; delete[] B;
}

/**
 * Perform the (forward) biased symmetric spherical Bessel transform,
 * defined here as
 * f@$
 *   b_\ell(k) = (2\pi k)^{-3/2} \int_0^\infty k \mathrm{d}r \,
 *     r^{3/2 + i} J_\ell(k r) a_\ell(r) \,.
 * f@$
 *
 * @note This is the Hankel transform definition used in cosmological
 *       wide-angle expansions,
 *       f@$
 *         \xi_\ell^{(n)}(r) = 4\pi \mathrm{i}^\ell \int_0^\infty
 *           \frac{\mathrm{d}k \, k^2}{(2\pi)^3}
 *           (k r)^{-n} j_\ell(kr) P_\ell^{(n)}(k) \,,
 *       f@$ i.e. f@$ a(r) f@$ corresponds to f@$ P_\ell^{(n)}(k) f@$ and
 *       f@$ b(k) f@$ corresponds to
 *       f@$ \mathrm{i}^{-\ell} \xi_\ell^{(n)}(r)f@$ with f@$ n f@$
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
) {
  double mu = ell + 1./2;
  double m = 2.;
  double q = 0.;
  double kr_c = 1.;
  bool lowring = true;

  std::complex<double>* A = new std::complex<double>[N];
  std::complex<double>* B = new std::complex<double>[N];

  for (int j = 0; j < N; j++) {
    A[j] = std::pow(r[j], m + i - 1./2) * a[j];
  }

  hankel_transform(mu, q, kr_c, N, lowring, r, A, k, B, NULL);

  for (int j = 0; j < N; j++) {
    std::complex<double> bj = std::pow(2*M_PI * k[j], -3./2) * B[j];
    b[j] = bj.real();
  }

  delete[] A; delete[] B;
}

/**
 * Transform power spectrum to correlation function multipole samples.
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
) {
  int m = 2;

  sj_transform(ell, m, N, k, pk, r, xi);
}

/**
 * Transform correlation function to power spectrum multipole samples.
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
) {
  int m = 2;

  sj_transform(ell, m, N, r, xi, k, pk);

  /// Factors of pi are needed here the forward transform is used
  /// in the backend as the backward transform.
  for (int j = 0; j < N; j++) {
    pk[j] *= 8 * M_PI * M_PI * M_PI;
  }
}

}  // trv::ops::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_FFTLOG_HPP_INCLUDED_
