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

namespace trva = trv::array;
namespace trvm = trv::maths;
namespace trvs = trv::sys;

namespace trv {

namespace maths {

HankelTransform::HankelTransform(double mu, double q, bool threaded) {
  // Check transform parameters.
  if (
    ((mu + 1. + q) <= 0. && std::floor(mu + 1. + q) == (mu + 1. + q)) ||
    ((mu + 1. - q) <= 0. && std::floor(mu + 1. - q) == (mu + 1. - q))
  ) {
    throw trvs::InvalidParameterError(
      "Invalid values of order `mu` and bias `q`: "
      "negative poles in the gamma function."
    );
  }

  this->order = mu;
  this->bias = q;

  this->threaded = threaded;

  // Initialise FFTW plans.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  if (this->threaded) {
    fftw_init_threads();
    fftw_plan_with_nthreads(omp_get_max_threads());
  }
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
}

HankelTransform::~HankelTransform() {
  this->reset();

  // Do NOT call `fftw_cleanup` in reusable code here
  // as it may cause memory leaks.
// #if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
//   if (this->threaded) {
//     fftw_cleanup_threads();
//   } else {
//     fftw_cleanup();
//   }
// #else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
//   fftw_cleanup();
// #endif  // TRV_USE_OMP && TRV_USE_FFTWOMP
}

void HankelTransform::reset() {
  if (this->plan_init) {
    fftw_destroy_plan(this->pre_plan);
    fftw_destroy_plan(this->post_plan);
    this->plan_init = false;
  }
  if (this->pre_buffer != nullptr) {
    fftw_free(this->pre_buffer); this->pre_buffer = nullptr;
  }
  if (this->post_buffer != nullptr) {
    fftw_free(this->post_buffer); this->post_buffer = nullptr;
  }
}

void HankelTransform::initialise(
  std::vector<double> sample_pts, double kr_c, bool lowring,
  trva::ExtrapOption extrap, double extrap_exp
) {
  // Initialise pre-sample points.
  if (sample_pts.size() < 2) {
    throw trvs::InvalidParameterError(
      "The number of sample points must be at least 2."
    );
  }
  if (trva::check_1d_array(sample_pts, false, true, false) != 0) {
    throw trvs::InvalidParameterError(
      "The sample points are not logarithmically spaced."
    );
  }

  this->pre_sampts = sample_pts;
  this->nsamp = sample_pts.size();
  this->logres =
    std::log(sample_pts.back() / sample_pts.front()) / (this->nsamp - 1);

  this->extrap = extrap;
  if (this->extrap != trva::ExtrapOption::NONE) {
    if (this->nsamp % 2 != 0) {
      throw trvs::InvalidParameterError(
        "The number of sample points must be even for extrapolation."
      );
    }
    this->nsamp_trans = std::pow(
      2, std::ceil(std::log2(extrap_exp * this->nsamp))
    );
    if (this->nsamp_trans < this->nsamp) {
      throw trvs::InvalidParameterError(
        "The sample size expansion factor results in a shrunken sample size."
      );
    }

    this->n_ext = (this->nsamp_trans - this->nsamp) / 2;
    trva::extrap_loglin(
      this->pre_sampts, this->n_ext, this->pre_sampts_extrap
    );
  } else {
    this->nsamp_trans = this->nsamp;
    this->n_ext = 0;
  }

  // Initialise the transform.
  if (lowring) {
    this->pivot = this->calc_lowring_pivot(this->logres, kr_c);
  } else {
    if (kr_c <= 0.) {
      throw trvs::InvalidParameterError("Pivot value must be positive.");
    }
    this->pivot = kr_c;
  }

  this->kernel = this->compute_kernel_coeff();

  // Initialise post-sample points.
  // Note the end-point of the periodic interval is not included in
  // the sample points, thus the shift in the product from the pivot
  // by one sample point.
  double kr_aprod = this->pivot * std::exp(- this->logres);

  this->post_sampts.resize(this->nsamp);
  for (int j = 0; j < this->nsamp; j++) {
    this->post_sampts[j] = kr_aprod / this->pre_sampts[this->nsamp - j - 1];
  }

  if (this->extrap != trva::ExtrapOption::NONE) {
    this->post_sampts_extrap.resize(this->nsamp_trans);
    for (int j = 0; j < this->nsamp_trans; j++) {
      this->post_sampts_extrap[j] =
        kr_aprod / this->pre_sampts[this->nsamp_trans - j - 1];
    }
  }

  // ---->
  // // Alternative to the for-loop block above, note that k_c and r_c are
  // // both exp(L/2) times `k0` and `r0`.
  // double kr_0 = kr_aprod * std::exp(- this->nsamp * this->logres);
  // double r0 = this->pre_sampts[0];
  // double k0 = kr_0 / r0;
  // ...
  // ----<

  // Initialise FFTW plans.
  this->reset();

  this->pre_buffer = fftw_alloc_complex(this->nsamp_trans);
  this->pre_plan = fftw_plan_dft_1d(
    this->nsamp_trans, this->pre_buffer, this->pre_buffer,
    FFTW_FORWARD, FFTW_ESTIMATE
  );

  this->post_buffer = fftw_alloc_complex(this->nsamp_trans);
  this->post_plan = fftw_plan_dft_1d(
    this->nsamp_trans, this->post_buffer, this->post_buffer,
    FFTW_FORWARD, FFTW_ESTIMATE
  );

  this->plan_init = true;
}

void HankelTransform::initialise(
  std::vector<double> sample_pts, double kr_c, bool lowring,
  int extrap, double extrap_exp
) {
  this->initialise(
    sample_pts, kr_c, lowring, trva::ExtrapOption(extrap), extrap_exp
  );
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
  int N_trans = this->nsamp_trans;
  double dL = this->logres;
  double kr_c = this->pivot;

  if (N_trans <= 0 || dL <= 0. || kr_c <= 0.) {
    throw std::runtime_error(
      "This instance of trv::maths::HankelTransform has not been "
      "initialised with `initialise`."
    );
  }

  double x_p = (mu + 1. + q)/2.;
  double x_m = (mu + 1. - q)/2.;
  double x_eq = (mu + 1.)/2.;

  double y = M_PI / (N_trans * dL);

  double t = -2. * y * std::log(kr_c/2.);

  // Note that no minus sign is involved in the phase by
  // complex conjugation of the gamma function.
  std::vector< std::complex<double> > u(N_trans);
  if (q == 0.) {
    for (int m = 0; m <= N_trans/2; m++) {
      double lnr_, phi_eq;
      trvm::get_lngamma_parts(x_eq, m*y, lnr_, phi_eq);

      u[m] = trvm::eval_complex_in_polar(1., m*t + 2*phi_eq);
    }
  } else {
    double qln2 = q * std::log(2.);
    for (int m = 0; m <= N_trans/2; m++) {
      double lnr_p, phi_p;
      double lnr_m, phi_m;
      trvm::get_lngamma_parts(x_p, m*y, lnr_p, phi_p);
      trvm::get_lngamma_parts(x_m, m*y, lnr_m, phi_m);

      u[m] = trvm::eval_complex_in_polar(
        std::exp(qln2 + lnr_p - lnr_m), m*t + phi_p + phi_m
      );
    }
  }

  for (int m = N_trans/2 + 1; m < N_trans; m++) {
    u[m] = conj(u[N_trans - m]);
  }

  if (N_trans % 2 == 0) {
    // Make the mid-point real by log-periodicity.
    u[N_trans/2] = u[N_trans/2].real() + trvm::M_I * 0.;
  }

  return u;
}

void HankelTransform::biased_transform(
  std::complex<double>* a, std::complex<double>* b
) {
  // STYLE: Standard naming convention is not followed below.
  int N = this->nsamp;
  int N_trans = this->nsamp_trans;

  if (this->kernel.empty() || N <= 2 || N_trans <= 2) {
    throw std::runtime_error(
      "This instance of trv::maths::HankelTransform has not been "
      "initialised with `initialise`."
    );
  }

  // Perform any extrapolation required.
  if (this->extrap == trva::ExtrapOption::NONE) {
    for (int j = 0; j < N_trans; j++) {
      this->pre_buffer[j][0] = (a[j]).real();
      this->pre_buffer[j][1] = (a[j]).imag();
      this->post_buffer[j][0] = 0.;
      this->post_buffer[j][1] = 0.;
    }
  } else {
    // Assume reality with extrapolation.
    std::vector<double> a_vec(N);
    for (int j = 0; j < N; j++) {
      a_vec[j] = a[j].real();
    }

    std::vector<double> a_trans_vec(N_trans);
    switch (this->extrap) {
      case trva::ExtrapOption::NONE:
        break;
      case trva::ExtrapOption::ZERO:
        trva::extrap_pad(
          a_vec, this->n_ext, 0., 0., a_trans_vec
        );
        break;
      case trva::ExtrapOption::PAD:
        trva::extrap_pad(
          a_vec, this->n_ext, a_vec.front(), a_vec.back(), a_trans_vec
        );
        break;
      case trva::ExtrapOption::LIN:
        trva::extrap_lin(a_vec, this->n_ext, a_trans_vec);
        break;
      case trva::ExtrapOption::LOGLIN:
        trva::extrap_loglin(a_vec, this->n_ext, a_trans_vec);
        break;
      default:
        throw trvs::InvalidParameterError("Unsupported extrapolation option.");
    }

    for (int j = 0; j < N_trans; j++) {
      this->pre_buffer[j][0] = a_trans_vec[j];
      this->pre_buffer[j][1] = 0.;
      this->post_buffer[j][0] = 0.;
      this->post_buffer[j][1] = 0.;
    }
  }

  // Compute the convolution b = a * u using FFT.
  fftw_execute(this->pre_plan);

  for (int m = 0; m < N_trans; m++) {
    // Divide by `N` to normalise the inverse DFT.
    std::complex<double> a_(this->pre_buffer[m][0], this->pre_buffer[m][1]);
    std::complex<double> b_ = a_ * this->kernel[m] / double(N_trans);
    this->post_buffer[m][0] = b_.real();
    this->post_buffer[m][1] = b_.imag();
  }

  fftw_execute(this->post_plan);

  // Trim any extrapolation.
  for (int j = 0; j < N; j++) {
    b[j] = std::complex<double>(
      this->post_buffer[j + this->n_ext][0],
      this->post_buffer[j + this->n_ext][1]
    );
  }

  // ---->
  // // Reverse and shift the array `b_trans` when inverse FFT is used above
  // // instead of FFT, and trim any extrapolation.
  // for (int m = 0; m < N_trans/2; m++) {
  //   double b_real_ = this->post_buffer[m][0];
  //   double b_imag_ = this->post_buffer[m][1];

  //   this->post_buffer[m][0] = this->post_buffer[N_trans - m - 1][0];
  //   this->post_buffer[m][1] = this->post_buffer[N_trans - m - 1][1];

  //   this->post_buffer[N_trans - m - 1][0] = b_real_;
  //   this->post_buffer[N_trans - m - 1][1] = b_imag_;
  // }
  // for (int j = 0; j < N; j++) {
  //   int j_ = (j + N_trans - 1) % N_trans + this->n_ext;
  //   b[j] = std::complex<double>(
  //     this->post_buffer[j_][0],
  //     this->post_buffer[j_][1]
  //   );
  // }
  // ----<
}

SphericalBesselTransform::SphericalBesselTransform(
  int ell, int n, bool threaded
) : HankelTransform(ell + 1./2, double(n), threaded) {
  this->degree = ell;
}

void SphericalBesselTransform::initialise(
  std::vector<double> sample_pts, double kr_c, bool lowring,
  trva::ExtrapOption extrap, double extrap_exp
) {
  HankelTransform::initialise(sample_pts, kr_c, lowring, extrap, extrap_exp);
}

void SphericalBesselTransform::initialise(
  std::vector<double> sample_pts, double kr_c, bool lowring,
  int extrap, double extrap_exp
) {
  HankelTransform::initialise(sample_pts, kr_c, lowring, extrap, extrap_exp);
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
      "The transform direction must be either +1 (forward) or -1 (backward).\n"
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
