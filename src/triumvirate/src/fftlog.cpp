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
 * @file fftlog.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "fftlog.hpp"

namespace trvm = trv::maths;
namespace trvs = trv::sys;

namespace trv {

namespace maths {

HankelTransform::HankelTransform(double mu, double q) {
  this->order = mu;
  this->bias = q;
}

void HankelTransform::initialise(
  std::vector<double> sample_pts, double kr_c, bool lowring
) {
  // Initialise pre-sample points.
  if (sample_pts.size() < 2) {
    throw trvs::InvalidParameterError(
      "The number of sample points must be at least 2."
    );
  }
  if (
    trv::array::check_1d_array(
      sample_pts.data(), sample_pts.size(), false, true
    ) != 0
  ) {
    throw trvs::InvalidParameterError(
      "The sample points are not log-linearly spaced."
    );
  }

  this->pre_sampts = sample_pts;
  this->nsamp = sample_pts.size();
  this->logres =
    std::log(sample_pts.back() / sample_pts.front()) / (this->nsamp - 1);

  // Initialise the transform.
  if (lowring) {
    this->pivot = this->calc_lowring_pivot(this->logres, kr_c);
  } else {
    if (kr_c <= 0.) {
      throw trvs::InvalidParameterError(
        "Pivot value must be positive."
      );
    }
    this->pivot = kr_c;
  }

  this->kernel = this->compute_kernel_coeff();

  // Initialise post-sample points.
  this->post_sampts.resize(this->nsamp);
  for (int j = 0; j < this->nsamp; j++) {
    this->post_sampts[j] = this->pivot / this->pre_sampts[this->nsamp - j - 1];
  }

  // // Alternative to the for-loop block above, note that k_c and r_c are
  // // both exp(L/2) times `k0` and `r0`.
  // // STYLE: Standard naming convention is not followed below.
  // int N = this->nsamp;
  // double dL = this->logres;
  // double kr_c_used = this->pivot;
  // double r0 = this->pre_sampts[0];

  // double kr_0 = kr_c_used * std::exp(- N * dL);
  // double k0 = kr_0 / r0;
}

double HankelTransform::calc_lowring_pivot(double delta, double kr_c) {
  // STYLE: Standard naming convention is not followed below.
  double mu = this->order;
  double q = this->bias;
  double dL = delta;

  double x_p = (mu + 1. + q)/2.;
  double x_m = (mu + 1. - q)/2.;
  double y = M_PI / (2.*dL);

  if (kr_c > 0.) {
    // Note that no minus sign is involved in the phase by
    // complex conjugation of the gamma function.
    double lnr_p, lnr_m, phi_p, phi_m;
    trvm::get_lngamma_parts(x_p, y, lnr_p, phi_p);
    trvm::get_lngamma_parts(x_m, y, lnr_m, phi_m);

    double argphase = std::log(2./kr_c) / dL + (phi_p + phi_m) / M_PI;

    kr_c *= std::exp(dL * (argphase - std::round(argphase)));
  } else {
    // Alternatively, simply set the value as follows.
    std::complex<double> gamma_lnratio = eval_gamma_lnratio(mu, q + 2*y*M_I);

    double lnkr_c = std::log(2) + gamma_lnratio.imag() / (2*y);

    kr_c = std::exp(lnkr_c - std::floor(lnkr_c / dL) * dL);
  }

  return kr_c;
}

std::vector< std::complex<double> > HankelTransform::compute_kernel_coeff() {
  // STYLE: Standard naming convention is not followed below.
  double mu = this->order;
  double q = this->bias;
  int N = this->nsamp;
  double dL = this->logres;
  double kr_c = this->pivot;

  if (N <= 0 || dL <= 0. || kr_c <= 0.) {
    throw std::runtime_error(
      "This instance of trv::maths::HankelTransform has not been "
      "initialised with `initialise`."
    );
  }

  double x_p = (mu + 1. + q)/2.;
  double x_m = (mu + 1. - q)/2.;
  double x_eq = (mu + 1.)/2.;

  double y = M_PI / (N * dL);

  double t = -2. * y * std::log(kr_c/2.);

  // Note that no minus sign is involved in the phase by
  // complex conjugation of the gamma function.
  std::vector< std::complex<double> > u(N);
  if (q == 0.) {
    for (int m = 0; m <= N/2; m++) {
      double lnr_, phi_eq;
      trvm::get_lngamma_parts(x_eq, m*y, lnr_, phi_eq);

      u[m] = trvm::eval_complex_in_polar(1., m*t + 2*phi_eq);
    }
  } else {
    double qln2 = q * std::log(2.);
    for (int m = 0; m <= N/2; m++) {
      double lnr_p, phi_p;
      double lnr_m, phi_m;
      trvm::get_lngamma_parts(x_p, m*y, lnr_p, phi_p);
      trvm::get_lngamma_parts(x_m, m*y, lnr_m, phi_m);

      u[m] = trvm::eval_complex_in_polar(
        std::exp(qln2 + lnr_p - lnr_m), m*t + phi_p + phi_m
      );
    }
  }

  for (int m = N/2 + 1; m < N; m++) {
    u[m] = conj(u[N - m]);
  }

  if (N % 2 == 0) {
    // Make the mid-point real by log-periodicity.
    u[N/2] = u[N/2].real() + trvm::M_I * 0.;
  }

  return u;
}

void HankelTransform::biased_transform(
  std::complex<double>* a, std::complex<double>* b
) {
  // STYLE: Standard naming convention is not followed below.
  int N = this->nsamp;

  if (this->kernel.empty() || N <= 2) {
    throw std::runtime_error(
      "This instance of trv::maths::HankelTransform has not been "
      "initialised with `initialise`."
    );
  }

  // Compute the convolution b = a * u using FFT.
  // NOTE: ``(`` and ``)`` are necessary
  // (see https://www.fftw.org/doc/Complex-numbers.html).
  fftw_plan forward_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) a, (fftw_complex*) b, FFTW_FORWARD, FFTW_ESTIMATE
  );
  fftw_execute(forward_plan);
  fftw_destroy_plan(forward_plan);

  for (int m = 0; m < N; m++) {
    // Divide by `N` to normalise the inverse DFT.
    b[m] *= this->kernel[m] / double(N);
  }

  fftw_plan reverse_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) b, (fftw_complex*) b, FFTW_BACKWARD, FFTW_ESTIMATE
  );
  fftw_execute(reverse_plan);
  fftw_destroy_plan(reverse_plan);

  // Reverse the array `b` as inverse FFT is used above instead of FFT.
  std::complex<double> b_;
  for (int n = 0; n < N/2; n++) {
    b_ = b[n];
    b[n] = b[N - n - 1];
    b[N - n - 1] = b_;
  }
}

SphericalBesselTransform::SphericalBesselTransform(int ell, int n): \
HankelTransform(ell + 1./2, double(n)) {
  this->degree = ell;
}

void SphericalBesselTransform::initialise(
  std::vector<double> sample_pts, double kr_c, bool lowring
) {
  HankelTransform::initialise(sample_pts, kr_c, lowring);
}

void SphericalBesselTransform::biased_transform(
  std::vector< std::complex<double> >& a,
  std::vector< std::complex<double> >& b
) {
  // STYLE: Standard naming convention is not followed below.
  int N = this->nsamp;

  if (int(a.size()) != N) {
    throw trvs::InvalidParameterError(
      "The size of array `a` must be equal to the number of samples."
    );
  }
  if (int(a.size()) != N) {
    b.resize(N);
  }

  std::complex<double> A[N];
  std::complex<double> B[N];

  for (int j = 0; j < N; j++) {
    A[j] = std::pow(this->pre_sampts[j], 3./2) * a[j];
  }

  HankelTransform::biased_transform(A, B);

  for (int j = 0; j < N; j++) {
    b[j] = std::pow(2*M_PI / this->post_sampts[j], 3./2) * B[j];
  }
}

void SphericalBesselTransform::transform_cosmological_multipole(
    int dir,
    std::vector< std::complex<double> >& pre_samples,
    std::vector< std::complex<double> >& post_samples
) {
  if (abs(dir) != 1) {
    throw trvs::InvalidParameterError(
      "The transform direction must be either +1 (forward) or -1 (backward)."
    );
  }

  double pifactor = (dir == -1) ? std::pow(2*M_PI, -3) : 1.;
  std::complex<double> parity = std::pow(trvm::M_I, - dir * this->degree);
  std::complex<double> prefactor = pifactor * parity;

  this->biased_transform(pre_samples, post_samples);

  for (int j = 0; j <this->nsamp; j++) {
    post_samples[j] *= prefactor;
  }
}

}  // namespace trv::maths

}  // namespace trv
