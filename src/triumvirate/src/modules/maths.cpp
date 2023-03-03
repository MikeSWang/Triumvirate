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
 * @file maths.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "maths.hpp"

namespace trvs = trv::sys;

namespace trv {
namespace maths {

// ***********************************************************************
// Complex numbers
// ***********************************************************************

/// @cond DOXYGEN_DOC_MISC
const std::complex<double> M_I(0., 1.);  ///< imaginary unit
/// @endcond

std::complex<double> eval_complex_in_polar(double r, double theta) {
  return r * (std::cos(theta) + M_I * std::sin(theta));
}


// ***********************************************************************
// Vectors
// ***********************************************************************

double get_vec3d_magnitude(std::vector<double> vec) {
  return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

double get_vec3d_magnitude(double* vec) {
  return std::sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}


// ***********************************************************************
// Gamma function
// ***********************************************************************

// Lanzcos approximation parameters.
// CAVEAT: Discretionary choices (for the Lanczos approximation series).
const int nterm_lanczos = 9;  ///< number of terms in approximation series
const double gconst_lanczos = 7.;  ///< Lanczos approximation constant
/// Lanczos approximation coefficients for @ref trv::maths::gconst_lanczos
/// and @ref trv::maths::nterm_lanczos
const double pcoeff_lanczos[] = {
      0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
  -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
   -176.61502916214059906584551354,
     12.507343278686904814458936853,
     -0.13857109526572011689554707,
      9.984369578019570859563e-6,
      1.50563273514931155834e-7,
};  // varying number of significant figures

std::complex<double> eval_lanczos_approx_series(std::complex<double> z) {
  std::complex<double> series = pcoeff_lanczos[0];
  for (int i = 1; i < nterm_lanczos; i++) {
    series += pcoeff_lanczos[i] / (z + double(i));
  }

  return series;
}

std::complex<double> eval_gamma(std::complex<double> z) {
  // Exploit Euler's reflection formula as the Lanczos approximation
  // is only valid for Re{z} > 1/2.
  if (z.real() < 1./2) {
    return M_PI / (std::sin(M_PI * z) * eval_gamma(1. - z));
  }

  // Substitute variables into the approximation formula.
  z -= 1.;

  std::complex<double> t = z + gconst_lanczos + 1./2;
  std::complex<double> series = eval_lanczos_approx_series(z);

  std::complex<double> gamma = std::sqrt(2*M_PI)
    * std::pow(t, z + 1./2) * std::exp(-t) * series;

  return gamma;
}

std::complex<double> eval_lngamma(std::complex<double> z) {
  // Exploit Euler's reflection formula as the Lanczos approximation
  // is only valid for Re{z} > 1/2.
  if (z.real() < 1./2) {
    return
      std::log(M_PI) - std::log(std::sin(M_PI * z)) - eval_lngamma(1. - z);
  }

  // Substitute variables into the approximation formula.
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
  std::complex<double> lngamma = eval_lngamma(x + M_I * y);

  if (lnr) {lnr = lngamma.real();}
  if (theta) {theta = lngamma.imag();}
}


// ***********************************************************************
// Spherical harmonics
// ***********************************************************************

// CAVEAT: Discretionary choice such that eps = 1.e-9.
const double eps_coupling = 1.e-9;

double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  return gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3);
}

std::complex<double> \
SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
  const int ell, const int m, double pos[3]
) {
  // CAVEAT: Discretionary choice such that eps = 1.e-9.
  const double eps = 1.e-9;

  // Return unity in the trivial case.
  if (ell == 0 && m == 0) {return 1.;}

  // Calculate modulus.
  double xyz_mod_sq = 0.;
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    xyz_mod_sq += pos[iaxis] * pos[iaxis];
  }  // r² = x² + y² + z²

  double xyz_mod = std::sqrt(xyz_mod_sq);  // r = √(x² + y² + z²)

  // Return zero in the trivial case.
  if (std::fabs(xyz_mod) < eps) {return 0.;}

  // Calculate the angular variable μ = cos(θ).
  double mu = pos[2] / xyz_mod;  // μ = z / r

  // Calculate the angular variable ϕ.
  double xy_mod = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]); // r_xy =
                                                                // √(x² + y²)

  double phi = 0.;
  if (std::fabs(xy_mod) >= eps) {
    phi = std::acos(pos[0] / xy_mod);  // ϕ = arccos(x / r_xy)
    if (pos[1] < 0.) {
      phi = - phi + 2.*M_PI;  // ϕ ∈ [π, 2π] if y < 0
    }
  }

  // Calculate spherical harmonics with m >= 0 via the normalised
  // associated Legendre polynomial, i.e. Y_lm = √((2l + 1)/(4π))
  // √((l - |m|)!/(l + |m|)!) P_l^|m|(μ) * exp(imϕ).
  std::complex<double> ylm = std::exp(M_I * double(m) * phi)
    * gsl_sf_legendre_sphPlm(ell, std::abs(m), mu);

  // Impose parity and conjugation.
  ylm = std::pow(-1, (m - std::abs(m))/2) * std::conj(ylm);

  // Normalise to the reduced form.
  ylm *= std::sqrt(4.*M_PI / (2.*ell + 1.));

  return ylm;
}

void SphericalHarmonicCalculator::\
store_reduced_spherical_harmonic_in_fourier_space(
  const int ell, const int m,
  const double boxsize[3], const int ngrid[3],
  std::vector< std::complex<double> >& ylm_out
) {
  // Determine the fundamental wavenumber in each dimension.
  double dk[3] = {
    2.*M_PI / boxsize[0], 2.*M_PI / boxsize[1], 2.*M_PI / boxsize[2]
  };

  // Assign a wavevector to each grid cell.
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < ngrid[0]; i++) {
    for (int j = 0; j < ngrid[1]; j++) {
      for (int k = 0; k < ngrid[2]; k++) {
        // Lay the 'bricks' vertically, then inwards, then to
        // the right, i.e. along z-axis, y-axis and then x-axis.
        // The assigned flattened-grid array index is
        // (i * ngrid_y * ngrid_z + j * ngrid_z + k)
        // where ngrid is the grid number along each axis.
        long long idx_grid = (i * ngrid[1] + j) * ngrid[2] + k;

        // This conforms to the (absurd) FFT array-ordering convention
        // that negative wavenumbers/frequencies come after zero and
        // positive wavenumbers/frequencies.
        double kvec[3];
        kvec[0] = (i < ngrid[0]/2) ? i * dk[0] : (i - ngrid[0]) * dk[0];
        kvec[1] = (j < ngrid[1]/2) ? j * dk[1] : (j - ngrid[1]) * dk[1];
        kvec[2] = (k < ngrid[2]/2) ? k * dk[2] : (k - ngrid[2]) * dk[2];

        ylm_out[idx_grid] = calc_reduced_spherical_harmonic(ell, m, kvec);
      }
    }
  }
}

void SphericalHarmonicCalculator::\
store_reduced_spherical_harmonic_in_config_space(
  const int ell, const int m,
  const double boxsize[3], const int ngrid[3],
  std::vector< std::complex<double> >& ylm_out
) {
  // Determine the grid cell size in each dimension.
  double dr[3] = {
    boxsize[0] / double(ngrid[0]),
    boxsize[1] / double(ngrid[1]),
    boxsize[2] / double(ngrid[2])
  };

  // Assign a position vector to each grid cell.
#ifdef TRV_USE_OMP
#pragma omp parallel for collapse(3)
#endif  // TRV_USE_OMP
  for (int i = 0; i < ngrid[0]; i++) {
    for (int j = 0; j < ngrid[1]; j++) {
      for (int k = 0; k < ngrid[2]; k++) {
        // Lay the 'bricks' vertically, then inwards, then to
        // the right, i.e. along z-axis, y-axis and then x-axis.
        // The assigned flattened-grid array index is
        // (i * ngrid_y * ngrid_z + j * ngrid_z + k)
        // where ngrid is the grid number along each axis.
        long long idx_grid = (i * ngrid[1] + j) * ngrid[2] + k;

        // This conforms to the (absurd) FFT array-ordering convention
        // that negative wavenumbers/frequencies come after zero and
        // positive wavenumbers/frequencies.
        double rvec[3];
        rvec[0] = (i < ngrid[0]/2) ? i * dr[0] : (i - ngrid[0]) * dr[0];
        rvec[1] = (j < ngrid[1]/2) ? j * dr[1] : (j - ngrid[1]) * dr[1];
        rvec[2] = (k < ngrid[2]/2) ? k * dr[2] : (k - ngrid[2]) * dr[2];

        ylm_out[idx_grid] = calc_reduced_spherical_harmonic(ell, m, rvec);
      }
    }
  }
}


// ***********************************************************************
// Spherical Bessel function
// ***********************************************************************

SphericalBesselCalculator::SphericalBesselCalculator(const int ell) {
  // Set up sampling range and number.
  // CAVEAT: Discretionary choices such that max(kr) > 4096π, Δ(kr) = 0.01.
  const double xmin = 0.;       // minimum of interpolation range
  const double xmax = 15000.;   // maximum of interpolation range
  const int nsample = 1500000;  // interpolation sample number

  // Initialise and evaluate at sample points.
  double dx = (xmax - xmin) / (nsample - 1);

  double* x = new double[nsample];
  double* j_ell = new double[nsample];

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
  for (int i = 0; i < nsample; i++) {
    x[i] = xmin + dx * i;
    j_ell[i] = gsl_sf_bessel_jl(ell, x[i]);
  }

  // Initialise the interpolator using cubic spline and the accelerator.
  this->accel = gsl_interp_accel_alloc();
  this->spline = gsl_spline_alloc(gsl_interp_cspline, nsample);

  gsl_spline_init(this->spline, x, j_ell, nsample);

  delete[] x; delete[] j_ell;
}

SphericalBesselCalculator::~SphericalBesselCalculator() {
  if (this->accel != nullptr) {
    gsl_interp_accel_free(this->accel); this->accel = nullptr;
  }

  if (this->spline != nullptr) {
    gsl_spline_free(this->spline); this->spline = nullptr;
  }
}

double SphericalBesselCalculator::eval(double x) {
  return gsl_spline_eval(this->spline, x, this->accel);
}

}  // namespace trv::maths
}  // namespace trv
