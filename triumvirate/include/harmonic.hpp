/**
 * @file harmonic.hpp
 * @brief (Reduced) spherical harmonic calculations.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_HARMONIC_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_HARMONIC_HPP_INCLUDED_

#include <cmath>
#include <complex>
#include <cstdio>

#include <gsl/gsl_sf_legendre.h>

#include "common.hpp"

/**
 * Calculation and storage of spherical harmonics.
 *
 * This includes the calculation of the 'reduced' (conjugated and
 * suitably normalised) spherical harmonics
 * @f$ \sqrt(\frac{4*\pi}{2\ell + 1}) Y_{\ell m}^\ast @f$.
 *
 * This also includes the storage of computed values on a regular mesh
 * grid in configuration or Fourier space.
 *
 */
class SphericalHarmonicCalculator {
 public:
  /**
   * Calculate reduced spherical harmonics.
   *
   * @param ell Degree of the spherical harmonic.
   * @param m Order of the spherical harmonic.
   * @param pos 3-d position vector.
   * @returns ylm Normalised conjugated value.
   */
  static std::complex<double> calc_reduced_spherical_harmonic(
    const int ell, const int m, double pos[3]
  ) {
    /// Return value in the trivial case.
    if (ell == 0 && m == 0) {
      return 1.;
    }

    /// Define zero precision.
    /// CAVEAT: Discretionary choice.
    const double eps = 1.e-15;

    /// Define complex unit.
    std::complex<double> _M_I(0., 1.);

    /// Calculate modulus.
    double xyz_mod_sq = 0.;
    for (int iaxis = 0; iaxis < 3; iaxis++) {
      xyz_mod_sq += pos[iaxis] * pos[iaxis];
    }  // r^2 = x^2 + y^2 + z^2

    double xyz_mod = sqrt(xyz_mod_sq);  // r = √(x^2 + y^2 + z^2)
    double xy_mod = sqrt(pos[0] * pos[0] + pos[1] * pos[1]);
      // r_xy = √(x^2 + y^2)

    /// Calculate the angular variable μ = cos(θ).
    double mu = 0.;
    if (fabs(xyz_mod) < eps) {
      /// Return value in the trivial case.
      return 0.;
    } else {
      mu = pos[2] / xyz_mod;  // μ = z / r
    }

    /// Calculate the angular variable ϕ.
    double phi = 0.;
    if (fabs(xy_mod) < eps) {
      /// Return default value in the special case.
      phi = 0.;
    } else {
      phi = acos(pos[0] / xy_mod);  // ϕ = arccos(x / r_xy)
      if (pos[1] < 0.) {
        phi = - phi + 2.*M_PI;  // ϕ ∈ [π, 2π] if y < 0
      }
    }

    /// Calculate spherical harmonics with m >= 0 via the
    /// normalised associated Legendre polynomial, i.e.
    /// Y_lm = √((2l + 1)/(4π)) √((l - |m|)!/(l + |m|)!) * P_l^|m| * exp(imϕ).
    std::complex<double> ylm = 0.;
    ylm = gsl_sf_legendre_sphPlm(ell, abs(m), mu) * exp(_M_I * double(m) * phi);

    /// Impose parity and conjugation.
    ylm = pow(-1., (m - abs(m))/2.) * std::conj(ylm);

    /// Normalise to the reduced form.
    ylm *= sqrt(4.*M_PI / (2.*ell + 1.));

    return ylm;
  }

  /**
   * Store reduced spherical harmonics computed in Fourier space.
   *
   * @param[in] ell Degree of the spherical harmonic.
   * @param[in] m Order of the spherical harmonic.
   * @param[in] boxsize Box size in each dimension.
   * @param[in] ngrid Grid number in each dimension.
   * @param[out] ylm_out Stored computed values.
   */
  static void store_reduced_spherical_harmonic_in_fourier_space(
    const int ell, const int m,
    const double boxsize[3], const int ngrid[3],
    std::complex<double>* ylm_out
  ) {
    /// Exit in error when no output variable is provided.
    if (ylm_out == NULL) {
      if (currTask == 0) {
        throw InvalidDataOps(
          "[%s ERRO] Cannot store computed spherical harmonics. "
          "Output array is null.\n",
          show_timestamp().c_str()
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

  /**
   * Store reduced spherical harmonics computed in configuration space.
   *
   * @param[in] ell Degree of the spherical harmonic.
   * @param[in] m Order of the spherical harmonic.
   * @param[in] boxsize Box size in each dimension.
   * @param[in] ngrid Grid number in each dimension.
   * @param[out] ylm_out Stored computed values.
   */
  static void store_reduced_spherical_harmonic_in_config_space(
    const int ell, const int m,
    const double boxsize[3], const int ngrid[3],
    std::complex<double>* ylm_out
  ) {
    /// Exit in error when no output variable is provided.
    if (ylm_out == NULL) {
      if (currTask == 0) {
        throw InvalidDataOps(
          "[%s ERRO] Cannot store computed spherical harmonics. "
          "Output array is null.\n",
          show_timestamp().c_str()
        );
      }
    }

    /// Determine the grid size in each dimension.
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
};

#endif  // TRIUMVIRATE_INCLUDE_HARMONIC_HPP_INCLUDED_
