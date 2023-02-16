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

namespace trv {

namespace maths {

double calc_kr_pivot_lowring(
  double mu, double q, double L, int N, double kr_c
) {
  double x_p = (mu + 1 + q)/2;
  double x_m = (mu + 1 - q)/2;
  double y = M_PI * N / (2*L);

  // Note that no minus is involved by complex conjugation of
  // the gamma function.
  double lnr_, phi_p, phi_m;
  trvm::get_lngamma_components(x_p, y, lnr_, phi_p);
  trvm::get_lngamma_components(x_m, y, lnr_, phi_m);

  double arg_ = std::log(2/kr_c) * N / L + (phi_p + phi_m) / M_PI;
  double arg_int = std::round(arg_);
  if (arg_ != arg_int) {
    kr_c *= std::exp(L / N * (arg_ - arg_int));
  }

  return kr_c;
}

void compute_u_kernel_coeff(
  double mu, double q, double L, int N, double kr_c, std::complex<double>* u
) {
  double y = M_PI / L;
  double kr_0 = kr_c * std::exp(-L);

  double t = -2 * y * std::log(kr_0/2);

  if (q == 0) {
    double x = (mu + 1)/2;

    double lnr, phi;
    for (int m = 0; m <= N/2; m++) {
      trvm::get_lngamma_components(x, m*y, lnr, phi);
      u[m] = trvm::eval_complex_in_polar(1., m*t + 2*phi);
    }
  } else {
    double x_p = (mu + 1 + q)/2;
    double x_m = (mu + 1 - q)/2;

    double lnr_p, phi_p;
    double lnr_m, phi_m;
    for (int m = 0; m <= N/2; m++) {
      trvm::get_lngamma_components(x_p, m*y, lnr_p, phi_p);
      trvm::get_lngamma_components(x_m, m*y, lnr_m, phi_m);

      u[m] = trvm::eval_complex_in_polar(
        std::exp(std::log(2) * q + lnr_p - lnr_m), m*t + phi_p - phi_m
      );
    }
  }

  for (int m = N/2 + 1; m < N; m++) {
    u[m] = conj(u[N - m]);
  }

  if (N % 2 == 0) {
    u[N/2] = u[N/2].real() + trvm::M_I * 0.;  // make real by log-periodicity
  }
}

void hankel_transform(
  double mu, double q, double kr_c, int N, bool lowring,
  double* r, const std::complex<double>* a,
  double* k, std::complex<double>* b,
  std::complex<double>* u
) {
  // Calculate the logarithmic interval.
  double L = N * std::log(r[N - 1] / r[0]) / (N - 1.);

  // Compute the forward transform kernel.
  std::complex<double>* u_ = nullptr;
  if (u == nullptr) {
    u_ = new std::complex<double>[N];

    if (lowring) {
      kr_c = calc_kr_pivot_lowring(mu, q, L, N, kr_c);
    }

    compute_u_kernel_coeff(mu, q, L, N, kr_c, u_);
    u = u_;
  }

  // Compute output sample points corresponding to the input sample points.
  double kr_0 = kr_c * std::exp(-L);

  k[0] = kr_0 / r[0];
  for (int j = 1; j < N; j++) {
    k[j] = k[0] * std::exp(j * L / N);
  }

  // Compute the convolution b = a * u using FFT.
  // ``(`` and ``)`` are necessary
  // (see https://www.fftw.org/doc/Complex-numbers.html).
  fftw_plan forward_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) a, (fftw_complex*) b, -1, FFTW_ESTIMATE
  );
  fftw_plan reverse_plan = fftw_plan_dft_1d(
    N, (fftw_complex*) b, (fftw_complex*) b, +1, FFTW_ESTIMATE
  );

  fftw_execute(forward_plan);
  for (int m = 0; m < N; m++) {
    b[m] *= u[m] / double(N);  // divide by `N` to normalise the inverse DFT
  }
  fftw_execute(reverse_plan);

  fftw_destroy_plan(forward_plan);
  fftw_destroy_plan(reverse_plan);

  // Reverse the array `b`.
  std::complex<double> b_;
  for (int n = 0; n < N/2; n++) {
    b_ = b[n];
    b[n] = b[N - n - 1];
    b[N - n - 1] = b_;
  }

  delete[] u_;
}

void sj_transform(
  int ell, int m, int N, double* r, double* a, double* k, double* b
) {
  double mu = ell + 1./2;
  double q = 0.;
  double kr_c = 1.;
  bool lowring = true;

  std::complex<double> A[N];
  std::complex<double> B[N];

  for (int j = 0; j < N; j++) {
    A[j] = std::pow(r[j], m - 1./2) * a[j];
  }

  hankel_transform(mu, q, kr_c, N, lowring, r, A, k, B, nullptr);

  for (int j = 0; j < N; j++) {
    std::complex<double> bj = std::pow(2*M_PI * k[j], -3./2) * B[j];
    b[j] = bj.real();
  }
}

void sj_transform_symm_biased(
  int ell, int i, int N, double* r, double* a, double* k, double* b
) {
  double mu = ell + 1./2;
  double m = 2.;
  double q = 0.;
  double kr_c = 1.;
  bool lowring = true;

  std::complex<double> A[N];
  std::complex<double> B[N];

  for (int j = 0; j < N; j++) {
    A[j] = std::pow(r[j], m + i - 1./2) * a[j];
  }

  hankel_transform(mu, q, kr_c, N, lowring, r, A, k, B, nullptr);

  for (int j = 0; j < N; j++) {
    std::complex<double> bj = std::pow(2*M_PI * k[j], -3./2) * B[j];
    b[j] = bj.real();
  }
}

}  // namespace trv::maths

void transform_powspec_to_corrfunc_multipole(
  int ell, int N, double* k, double* pk, double* r, double* xi
) {
  int m = 2;

  trvm::sj_transform(ell, m, N, k, pk, r, xi);
}

void transform_corrfunc_to_powspec_multipole(
  int ell, int N, double* r, double* xi, double* k, double* pk
) {
  int m = 2;

  trvm::sj_transform(ell, m, N, r, xi, k, pk);

  // Factors of π are needed here the forward transform is used
  // in the backend as the backward transform.
  for (int j = 0; j < N; j++) {
    pk[j] *= 8 * M_PI * M_PI * M_PI;
  }
}

}  // namespace trv
