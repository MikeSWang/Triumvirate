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
 * @file gamma.cpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 *
 */

#include "gamma.hpp"

namespace trv {
namespace maths {

/// Lanzcos approximation parameters.
/// CAVEAT: Discretionary choices.
const int nterm_lanczos = 9;  ///< number of terms in approximation series
const double gconst_lanczos = 7.;  ///< Lanczos approximation constant
const double pcoeff_lanczos[] = {
      0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
  -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
   -176.61502916214059906584551354,
     12.507343278686904814458936853,
     -0.13857109526572011689554707,
      9.984369578019570859563e-6,
      1.50563273514931155834e-7
};  ///< Lanczos approximation coefficients for `gconst` & `nterm_lanczos`
    // varying number of significant figures

std::complex<double> eval_lanczos_approx_series(std::complex<double> z) {
  std::complex<double> series = pcoeff_lanczos[0];
  for (int i = 1; i < nterm_lanczos; i++) {
    series += pcoeff_lanczos[i] / (z + double(i));
  }

  return series;
}

std::complex<double> eval_gamma(std::complex<double> z) {
  /// Exploit Euler's reflection formula as the Lanczos approximation
  /// is only valid for Re{z} > 1/2.
  if (z.real() < 1./2) {
    return M_PI / (std::sin(M_PI * z) * eval_gamma(1. - z));
  }

  /// Substitute variables into the approximation formula.
  z -= 1.;

  std::complex<double> t = z + gconst_lanczos + 1./2;
  std::complex<double> series = eval_lanczos_approx_series(z);

  std::complex<double> gamma = std::sqrt(2*M_PI)
    * std::pow(t, z + 1./2) * std::exp(-t) * series;

  return gamma;
}

std::complex<double> eval_lngamma(std::complex<double> z) {
  /// Exploit Euler's reflection formula as the Lanczos approximation
  /// is only valid for Re{z} > 1/2.
  if (z.real() < 1./2) {
    return std::log(M_PI) - std::log(std::sin(M_PI * z)) - eval_lngamma(1. - z);
  }

  /// Substitute variables into the approximation formula.
  z -= 1;

  std::complex<double> t = z + gconst_lanczos + 1./2;
  std::complex<double> series = eval_lanczos_approx_series(z);

  std::complex<double> lngamma = std::log(2*M_PI) / 2.
    + (z + 1./2) * std::log(t) - t + std::log(series);

  return lngamma;
}

std::complex<double> eval_gamma_ratio_asymp(
  double mu, std::complex<double> nu
) {
  std::complex<double> x_p = (mu + 1 + nu)/2.;
  std::complex<double> x_m = (mu + 1 - nu)/2.;

  std::complex<double> lnratio =
    - nu
    + (x_p - 1./2) * std::log(x_p) - (x_m - 1./2) * std::log(x_m)
    + 1./12 * (1./x_p - 1./x_m)
    - 1./360 * (1./std::pow(x_p, 3) - 1./std::pow(x_m, 3))
    + 1./1260 * (1./std::pow(x_p, 5) - 1./std::pow(x_m, 5));

  return lnratio;
}

void get_lngamma_components(double x, double y, double& lnr, double& theta) {
  /// Define complex unit.
  const std::complex<double> M_I_(0., 1.);

  std::complex<double> lngamma = eval_lngamma(x + M_I_ * y);

  if (lnr) {lnr = lngamma.real();}
  if (theta) {theta = lngamma.imag();}
}

}  // namespace trv::maths
}  // namespace trv
