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
 * @file gamma.hpp
 * @author Mike S Wang (https://github.com/MikeSWang)
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
 * @brief Evaluate the Lanczos approximation series f@$ A_g(z) f@$
 *        for the gamma function.
 *
 * The series approximate the gamma function by
 *
 * f@[
 *   \Gamma(z + 1) = \sqrt{2\pi}
 *     (z + g + 1/2)^{z + 1/2} \mathrm{e}^{- z - g - 1/2} A_g(z)
 * f@]
 *
 * where f@$ g f@$ is the Lanczos constant.
 *
 * @param z Complex argument.
 * @returns Value of f@$ A_g(z) f@$.
 */
std::complex<double> eval_lanczos_approx_series(std::complex<double> z);

/**
 * @brief Evaluate the gamma function f@$ \Gamma(z) f@$ on
 *        the complex plane using the Lanczos approximation.
 *
 * @param z Complex argument.
 * @returns Value of f@$ \Gamma(z) f@$.
 *
 * @see trv::eval_lanczos_approx_series()
 */
std::complex<double> eval_gamma(std::complex<double> z);

/**
 * @brief Evaluate the log-gamma function f@$ \ln\Gamma(z) f@$
 *        on the complex plane using the Lanczos approximation.
 *
 * @param z Complex argument.
 * @returns Value of f@$ \ln\Gamma(z) f@$.
 *
 * @see trv::eval_lanczos_approx_series()
 */
std::complex<double> eval_lngamma(std::complex<double> z);

/**
 * @brief Evaluate the logarithm of the ratio of two gamma functions.
 *
 * The ratio is of the form
 * f@$
 *   \Gamma(\frac{\nu + \mu + 1}{2}) / \Gamma(\frac{\nu - \mu + 1}{2})
 * f@$
 * in the asymptotic limit f@$ \mu, \nu \to \infty f@$.
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

}  // namespace trv::maths
}  // namespace trv

#endif  // !TRIUMVIRATE_INCLUDE_GAMMA_HPP_INCLUDED_
