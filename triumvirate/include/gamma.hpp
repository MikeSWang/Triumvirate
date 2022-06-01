/**
 * @file gamma.hpp
 * @brief Gamma function calculations.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_GAMMA_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_GAMMA_HPP_INCLUDED_

#include <cmath>
#include <complex>

#include "tools.hpp"

namespace trv {
namespace maths {

/**
 * Evaluate the gamma function f@$ \Gamma(z) f@$ on the complex plane
 * using the Lanczos approximation (with nine terms and f@$ g = 7 f@$).
 *
 * @param z Complex argument.
 * @returns Function value.
 */
std::complex<double> eval_gamma(std::complex<double> z) {
  /// Set up ``N = 9`` Lanczos coefficients with Lanzcos parameter ``g = 7``.
  int N = 9;
  double g = 7.;

  double coeff_lanczos[] = {
    0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
    -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
    12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7
  };  // varying number of significant figures

  /// Exploit Euler's reflection formula as the Lanczos approximation
  /// is only valid for Re{z} > 1/2.
  if (z.real() < 1./2) {
    return M_PI / (std::sin(M_PI * z) * eval_gamma(1. - z));
  }

  /// Evaluate the Lanczos approximation series.
  z -= 1;
  std::complex<double> t = z + g + 1./2;

  std::complex<double> term = coeff_lanczos[0];
  for (int n = 1; n < N; n++) {
    term += coeff_lanczos[n] / (z + double(n));  // double-conversion essential
  }

  return std::sqrt(2*M_PI) * std::pow(t, z + 1./2) * std::exp(-t) * term;
}

/**
 * Evaluate the log-gamma function f@$ \ln\Gamma(z) f@$
 * on the complex plane using the Lanczos approximation (with nine terms
 * and f@$ g = 7 f@$).
 *
 * @param z Complex argument.
 * @returns Function value.
 * @see eval_gamma(std::complex<double>)
 */
std::complex<double> eval_lngamma(std::complex<double> z) {
  int N = 9;
  double g = 7.;

  double coeff_lanczos[] = {
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

  if (z.real() < 1./2) {
    return std::log(M_PI) - std::log(std::sin(M_PI * z)) - eval_lngamma(1. - z);
  }

  z -= 1;
  std::complex<double> t = z + g + 1./2;

  std::complex<double> term = coeff_lanczos[0];
  for (int n = 1; n < N; n++) {
    term += coeff_lanczos[n] / (z + double(n));
  }
  std::complex<double> res = std::log(2*M_PI) / 2.
    + (z + 1./2) * std::log(t) - t + std::log(term);

  return std::log(2*M_PI) / 2. + (z + 1./2) * std::log(t) - t + std::log(term);
}

/**
 * Evaluate logarithm of the ratio of two gamma functions, in the form
 * f@$
 *   \Gamma(\frac{\nu + \mu + 1}{2}) / \Gamma(\frac{\nu - \mu + 1}{2})
 * f@$
 * in the asymptotic limit f@$ \mu, \nu \to \infty f@$.
 *
 * @param mu Argument in the nominator.
 * @param nu Argument in the .
 * @return ln_ratio Approximate logarithmic ratio.
 */
std::complex<double> eval_gamma_ratio_asymp(
  double mu, std::complex<double> nu
) {
  std::complex<double> x_p = (mu + 1 + nu)/2.;
  std::complex<double> x_m = (mu + 1 - nu)/2.;

  std::complex<double> ln_ratio = - nu
    + (x_p - 1./2) * std::log(x_p) - (x_m - 1./2) * std::log(x_m)
    + 1./12 * (1./x_p - 1./x_m)
    - 1./360 * (1./std::pow(x_p, 3) - 1./std::pow(x_m, 3))
    + 1./1260 * (1./std::pow(x_p, 5) - 1./std::pow(x_m, 5));

  return ln_ratio;
}

/**
 * Return the real and imaginary parts of the log-gamma function.
 *
 * @param[in] x Real part of the complex argument.
 * @param[in] y Imaginary part of the complex argument.
 * @param[out] lnr Real part of the log-gamma function value.
 * @param[out] theta Imaginary part of the log-gamma function value.
 */
void get_lngamma_components(
  double x, double y, double* lnr, double* theta
) {
  std::complex<double> lngamma = eval_lngamma(x + M_I * y);
  if (lnr) {
    *lnr = lngamma.real();
  }
  if (theta) {
    *theta = lngamma.imag();
  }
}

}  // trv::maths::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_GAMMA_HPP_INCLUDED_
