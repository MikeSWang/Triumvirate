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
 * @file harmonic.cpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 *
 */

#include "harmonic.hpp"

namespace trv {
namespace maths {

double wigner_3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  return gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3);
}

std::complex<double> SphericalHarmonicCalculator::calc_reduced_spherical_harmonic(
  const int ell, const int m, double pos[3]
) {
  /// Define zero precision.
  /// CAVEAT: Discretionary choice.
  const double eps = 1.e-15;

  /// Define complex unit.
  const std::complex<double> M_I_(0., 1.);

  /// Return unity in the trivial case.
  if (ell == 0 && m == 0) {return 1.;}

  /// Calculate modulus.
  double xyz_mod_sq = 0.;
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    xyz_mod_sq += pos[iaxis] * pos[iaxis];
  }  // r^2 = x^2 + y^2 + z^2

  double xyz_mod = std::sqrt(xyz_mod_sq);  // r = √(x^2 + y^2 + z^2)

  /// Return zero in the trivial case.
  if (std::fabs(xyz_mod) < eps) {return 0.;}

  /// Calculate the angular variable μ = cos(θ).
  double mu = pos[2] / xyz_mod;  // μ = z / r

  /// Calculate the angular variable ϕ.
  double xy_mod = std::sqrt(pos[0] * pos[0] + pos[1] * pos[1]); // r_xy =
                                                                // √(x^2 + y^2)

  double phi = 0.;
  if (std::fabs(xy_mod) >= eps) {
    phi = std::acos(pos[0] / xy_mod);  // ϕ = arccos(x / r_xy)
    if (pos[1] < 0.) {
      phi = - phi + 2.*M_PI;  // ϕ ∈ [π, 2π] if y < 0
    }
  }

  /// Calculate spherical harmonics with m >= 0 via the normalised
  /// associated Legendre polynomial, i.e. Y_lm = √((2l + 1)/(4π))
  /// √((l - |m|)!/(l + |m|)!) P_l^|m|(μ) * exp(imϕ).
  std::complex<double> ylm = std::exp(M_I_ * double(m) * phi)
    * gsl_sf_legendre_sphPlm(ell, std::abs(m), mu);

  /// Impose parity and conjugation.
  ylm = std::pow(-1, (m - std::abs(m))/2) * std::conj(ylm);

  /// Normalise to the reduced form.
  ylm *= std::sqrt(4.*M_PI / (2.*ell + 1.));

  return ylm;
}

void SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_fourier_space(
  const int ell, const int m,
  const double boxsize[3], const unsigned int ngrid[3],
  std::complex<double>* ylm_out
) {
  /// Throw an exception when the output variable is provided.
  if (ylm_out == nullptr) {
    if (trv::mon::currTask == 0) {
      throw trv::mon::InvalidData(
        "[%s ERRO] Cannot store computed spherical harmonics. "
        "Output array is not provided.\n",
        trv::mon::show_timestamp().c_str()
      );
    }
  }

  /// Determine the fundamental wavenumber in each dimension.
  double dk[3] = {
    2.*M_PI / boxsize[0], 2.*M_PI / boxsize[1], 2.*M_PI / boxsize[2]
  };

  /// Assign a wavevector to each grid cell.
  double kvec[3];
  for (int i = 0; i < ngrid[0]; i++) {
    for (int j = 0; j < ngrid[1]; j++) {
      for (int k = 0; k < ngrid[2]; k++) {
        /// Lay the 'bricks' vertically, then inwards, then to
        /// the right, i.e. along z-axis, y-axis and then x-axis.
        /// The assigned flattened-grid array index is
        /// (i * ngrid_y * ngrid_z + j * ngrid_z + k)
        /// where ngrid is the grid number along each axis.
        long long idx_grid = (i * ngrid[1] + j) * ngrid[2] + k;

        /// This conforms to the (absurd) FFT array-ordering convention
        /// that negative wavenumbers/frequencies come after zero and
        /// positive wavenumbers/frequencies.
        kvec[0] = (i < ngrid[0]/2) ? i * dk[0] : (i - ngrid[0]) * dk[0];
        kvec[1] = (j < ngrid[1]/2) ? j * dk[1] : (j - ngrid[1]) * dk[1];
        kvec[2] = (k < ngrid[2]/2) ? k * dk[2] : (k - ngrid[2]) * dk[2];

        ylm_out[idx_grid] = calc_reduced_spherical_harmonic(ell, m, kvec);
      }
    }
  }
}

void SphericalHarmonicCalculator::store_reduced_spherical_harmonic_in_config_space(
  const int ell, const int m,
  const double boxsize[3], const unsigned int ngrid[3],
  std::complex<double>* ylm_out
) {
  /// Throw an exception when the output variable is provided.
  if (ylm_out == nullptr) {
    if (trv::mon::currTask == 0) {
      throw trv::mon::InvalidData(
        "[%s ERRO] Cannot store computed spherical harmonics. "
        "Output array is not provided.\n",
        trv::mon::show_timestamp().c_str()
      );
    }
  }

  /// Determine the grid cell size in each dimension.
  double dr[3] = {
    boxsize[0] / double(ngrid[0]),
    boxsize[1] / double(ngrid[1]),
    boxsize[2] / double(ngrid[2]),
  };

  /// Assign a position vector to each grid cell.
  double rvec[3];
  for (int i = 0; i < ngrid[0]; i++) {
    for (int j = 0; j < ngrid[1]; j++) {
      for (int k = 0; k < ngrid[2]; k++) {
        /// Lay the 'bricks' vertically, then inwards, then to
        /// the right, i.e. along z-axis, y-axis and then x-axis.
        /// The assigned flattened-grid array index is
        /// (i * ngrid_y * ngrid_z + j * ngrid_z + k)
        /// where ngrid is the grid number along each axis.
        long long idx_grid = (i * ngrid[1] + j) * ngrid[2] + k;

        /// This conforms to the (absurd) FFT array-ordering convention
        /// that negative wavenumbers/frequencies come after zero and
        /// positive wavenumbers/frequencies.
        rvec[0] = (i < ngrid[0]/2) ? (i * dr[0]) : (i - ngrid[0]) * dr[0];
        rvec[1] = (j < ngrid[1]/2) ? (j * dr[1]) : (j - ngrid[1]) * dr[1];
        rvec[2] = (k < ngrid[2]/2) ? (k * dr[2]) : (k - ngrid[2]) * dr[2];

        ylm_out[idx_grid] = calc_reduced_spherical_harmonic(ell, m, rvec);
      }
    }
  }
}

}  // namespace trv::maths
}  // namespace trv
