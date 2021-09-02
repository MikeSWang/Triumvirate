#ifndef HANKEL_FFTLOG_H_INCLUDED_
#define HANKEL_FFTLOG_H_INCLUDED_

#include <complex>

#include <fftw3.h>

const double PI = 3.14159265358979323846;
const std::complex<double> I(0., 1.);

/**
 * Evaluate a complex number f@$ r e^{i \phi} f@$ in the polar form.
 *
 * @param r Modulus.
 * @param phi Argument.
 * @returns Value.
 */
std::complex<double> eval_polar_complex (double r, double phi) {
  return r * (cos(phi) + I * sin(phi));
}

/**
 * Evaluate the gamma function f@$ \Gamma(z) f@$ with complex arguments
 * using the Lanczos approximation.
 *
 * @param z Complex argument of the function.
 * @returns Value of the function.
 */
std::complex<double> eval_gamma(std::complex<double> z) {
  /// Set up N = 9 Lanczos coefficients with Lanzcos parameter g = 7.
  static double g = 7.;
  static int N = 9;
  static double lanczos_coeff[] = {
    0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
    -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
    12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7
  };

  /// Explot Euler's reflection formula as the Lanczos approximation
  /// is only valid for Re{z} > 1/2.
  if (z.real() < 1./2) {
    return PI / (sin(PI * z) * eval_gamma(1. - z));
  }

  /// Evaluate the Lanczos approximation series.
  z -= 1;
  std::complex<double> t = z + g + 1./2;
  std::complex<double> x = lanczos_coeff[0];
  for (int n = 1; n < N; n++) {
	  x += lanczos_coeff[n] / (z + double(n));  // double-conversion essential
  }

  return sqrt(2*PI) * pow(t, z + 1./2) * exp(-t) * x;
}

/**
 * Evaluate the log-gamma function f@$ \ln\Gamma(z) f@$.
 *
 * @param z Complex argument of the function.
 * @returns Value of the function.
 */
std::complex<double> eval_lngamma(std::complex<double> z) {
  return log(eval_gamma(z));
}

/**
 * Extract the polar components of the gamma function value (i.e. the real
 * and imaginary parts of the log-gamma function value).
 *
 * @param[in] x Real part of the complex argument.
 * @param[in] y Imaginary part of the complex argument.
 * @param[out] lnr (Pointer to) log-modulus of the gamma function value.
 * @param[out] phi (Pointer to) argument of the gamma function value.
 */
static void get_lngamma_components(
  double x, double y, double* lnr, double* phi
) {
  std::complex<double> lngamma = eval_lngamma(x + I * y);
  if (lnr) {*lnr = lngamma.real();}
  if (phi) {*phi = lngamma.imag();}
}

/**
 * Calculate a low-ringing pivot value f@$ k_c r_c f@$.
 *
 * @param N Sampling number.
 * @param mu Orderof the Hankel transform.
 * @param q FFTLog power-law bias index.
 * @param L Log-interval period.
 * @param kr_c Initial pivot value.
 */
static double calc_lowring_kr_pivot(
  int N, double mu, double q, double L, double kr_c
) {
  double x_plus = (mu + 1 + q)/2;
  double x_minus = (mu + 1 - q)/2;
  double y = PI * N / (2*L);

  /// Note no minus is involved by complex conjugation of the gamma
  /// function.
  double lnr, phi_minus, phi_plus;
  get_lngamma_components(x_plus, y, &lnr, &phi_plus);
  get_lngamma_components(x_minus, y, &lnr, &phi_minus);
  double arg = log(2 / kr_c) * N/L + (phi_plus + phi_minus) / PI;

  double arg_int = round(arg);
  if (arg != arg_int) {
    kr_c *= exp(L/N * (arg - arg_int));
  }

  return kr_c;
}

/**
 * Calculate the FFTLog Hankel transform kernel f@$ U f@$.
 *
 * @param N Sampling number.
 * @param mu Orderof the Hankel transform.
 * @param q FFTLog power-law bias index.
 * @param L Log-interval period.
 * @param kr_c Pivot f@$ kr f@$ value.
 * @param[out] u Kernel coefficients.
 */
void calc_hankel_kernel_u_coeff(
  int N, double mu, double q, double L, double kr_c, std::complex<double> u[]
) {
  double y = PI / L;
  double kr_0 = kr_c * exp(-L);
  double t = -2 * y * log(kr_0 / 2);

  if (q == 0) {
    double x = (mu + 1)/2;

    double lnr, phi;
    for (int m = 0; m <= N/2; m++) {
      get_lngamma_components(x, m*y, &lnr, &phi);

      u[m] = eval_polar_complex(1., m*t + 2*phi);
    }
  } else {
    double x_plus = (mu + 1 + q)/2;
    double x_minus = (mu + 1 - q)/2;

    double lnr_plus, phi_plus, lnr_minus, phi_minus;
    for (int m = 0; m <= N/2; m++) {
      get_lngamma_components(x_plus, m*y, &lnr_plus, &phi_plus);
      get_lngamma_components(x_minus, m*y, &lnr_minus, &phi_minus);

      u[m] = eval_polar_complex(
        exp(log(2) * q + lnr_plus - lnr_minus),
        m*t + phi_plus - phi_minus
      );
    }
  }

  for (int m = N/2 + 1; m < N; m++) {
    u[m] = conj(u[N - m]);
  }

  if (N % 2 == 0) {
    u[N/2] = u[N/2].real() + I * 0.;  // make real by log-periodicity
  }
}


/**
 * Perform the forward Hankel transform, defined as
 *
 * f@$
 *   a(k) = \int_0^\infty k dr (k r)^q J_\mu(k r) a(r) \,.
 * f@$
 *
 * @param[in] N Sampling number.
 * @param[in] mu Orderof the Hankel transform.
 * @param[in] q FFTLog power-law bias index.
 * @param[in] kr_c Pivot f@$ kr f@$ value.
 * @param[in] lowring Low-ringing condition (binary, 0 or 1).
 * @param[in] r Pre-transform sample points.
 * @param[out] k Post-transform sample points.
 * @param[out] a Pre-transform samples.
 * @param[out] b Post-trasform samples.
 * @param[out] u FFTLog transform kernel coefficients.
 */
void forward_hankel_transform(
  int N, double mu, double q, double kr_c, int lowring,
  const double r[], double k[],
  std::complex<double> a[], std::complex<double> b[],
  std::complex<double>* u
) {
  double L = N * log(r[N - 1] / r[0]) / (N - 1.);

  /// Compute the forward transform kernel.
  std::complex<double>* u_ = NULL;
  if (u == NULL) {
    if (lowring) {
      kr_c = calc_lowring_kr_pivot(N, mu, q, L, kr_c);
    }

    u_ = new std::complex<double>[N];
    calc_hankel_kernel_u_coeff(N, mu, q, L, kr_c, u_);
    u = u_;
  }

  /// Compute the convolution ``b = a * u`` using FFT.
  fftw_plan forward_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) a, (fftw_complex*) b, -1, FFTW_ESTIMATE
  );
  fftw_plan reverse_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) b, (fftw_complex*) b, +1, FFTW_ESTIMATE
  );
    // `(` and `)` necessary (see www.fftw.org/doc/Complex-numbers.html)

  fftw_execute(forward_plan);
  for (int m = 0; m < N; m++) {
    b[m] *= u[m] / double(N);
      // division by N as FFTW doesn't normalise the inverse FFT
  }
  fftw_execute(reverse_plan);

  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(reverse_plan);

  /// Reverse the array `b`.
  std::complex<double> b_tmp;
  for (int n = 0; n < N/2; n++) {
    b_tmp = b[n];
    b[n] = b[N - n - 1];
    b[N - n - 1] = b_tmp;
  }

  /// Compute sample points `k` corresponding to input `r`.
  double kr_0 = kr_c * exp(-L);

  k[0] = kr_0 / r[0];
  for (int j = 1; j < N; j++) {
    k[j] = k[0] * exp(j*L/N);
  }

  delete[] u_;
}

/**
 * Perform the cosmological forward Hankel transform (spherical
 * Bessel transform), defined as
 *
 * f@$
 *   \xi_\ell(r) = (2\pi r)^{-3/2} \int_0^\infty
 *     r dk k^{m - 1/2} J_\ell(k r) P_\ell(k) \,.
 * f@$
 *
 * Note the comparison with the usual (but generalised) definition
 *
 * f@$
 *   \xi_\ell^{(n)}(r) = 4\pi i^\ell \int_0^\infty dk k^2 / (2\pi)^3
 *     (k r)^{-n} j_\ell(kr) P_\ell^{(n)}(k) \,,
 * f@$
 *
 * i.e. f@$ \xi_\ell(r; m = 2) = i^{-\ell} \xi_\ell^{(n = 0)}(r) f@$,
 * where 'n' is identified with `q` in the standard Hankel transform.
 *
 * @param[in] ell Order of the transform.
 * @param[in] m Dimensional power-law index.
 * @param[in] N Sampling number.
 * @param[in] k Wavenumber sample points.
 * @param[in] pk Power spectrum samples.
 * @param[out] r Separation sample points.
 * @param[out] xi Correlation function samples.
 */
int cosmo_forward_hankel_transform(
  int ell, int m, int N, double* k, double* pk, double* r, double* xi
) {
  static double q = 0.;

  std::complex<double>* a = new std::complex<double>[N];
  std::complex<double>* b = new std::complex<double>[N];

  for (int j = 0; j < N; j++) {
    a[j] = pow(k[j], m - 1./2) * pk[j];
  }

  forward_hankel_transform(N, ell + 1./2, q, 1., 1, k, r, a, b, NULL);

  for (int j = 0; j < N; j++) {
    std::complex<double> xi_j = pow(2 * PI * r[j], -3./2) * b[j];
    xi[j] = xi_j.real();
  }

  delete[] a;
  delete[] b;

  return 0;
}

/**
 * Perform the Hankel transform (spherical Bessel transform) for
 * wide-angle correction terms.
 *
 * @param[in] ell Order of the transform.
 * @param[in] i Wide-angle correction order.
 * @param[in] N Sampling number.
 * @param[in] k Wavenumber sample points.
 * @param[in] pk Power spectrum samples.
 * @param[out] r Separation sample points.
 * @param[out] xi Correlation function samples.
 */
int wide_angle_hankel_transform(
  int ell, int i, int N, double* k, double* pk, double* r, double* xi
) {
  static double m = 2.;
  static double q = 0.;

  std::complex<double>* a = new std::complex<double>[N];
  std::complex<double>* b = new std::complex<double>[N];

  for (int j = 0; j < N; j++) {
    a[j] = pow(k[j], m + i - 1./2) * pk[j];
  }

  forward_hankel_transform(N, ell + 1./2, q, 1., 1, k, r, a, b, NULL);

  for (int j = 0; j < N; j++) {
    std::complex<double> xi_j = pow(2 * PI * r[j], -3./2) * b[j];
    xi[j] = xi_j.real();
  }

  delete[] a;
  delete[] b;

  return 0;
}

/**
 * Transform power spectrum samples to correlation function samples
 * (isotropic or monopole only).
 *
 * @param[in] N Sample number.
 * @param[in] k Wavenumber sample points.
 * @param[in] pk Power spectrum samples.
 * @param[out] r Separation sample points.
 * @param[out] xi Correlation function samples.
 */
int transform_power_spec_to_corr_func(
  int N, double* k, double* pk, double* r, double* xi
) {
  cosmo_forward_hankel_transform(0, 2, N, k, pk, r, xi);

  return 0;
}

/**
 * Transform correlation function samples to power spectrum samples
 * (isotropic or monopole only).
 *
 * @param[in] N Sample number.
 * @param[in] r Separation sample points.
 * @param[in] xi Correlation function samples.
 * @param[out] k Wavenumber sample points.
 * @param[out] pk Power spectrum samples.
 */
int transform_corr_func_to_power_spec(
  int N, double* r, double* xi, double* k, double* pk
) {
  cosmo_forward_hankel_transform(0, 2, N, r, xi, k, pk);

  /// Factors of pi are needed here the forward transform is used
  /// in the backend as the backward transform.
  static const double cubic_2pi = 8 * PI * PI * PI;
  for (int j = 0; j < N; j++) {
    pk[j] *= cubic_2pi;
  }

  return 0;
}

#endif
