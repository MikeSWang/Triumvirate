#ifndef TRIUM_TOOLS_H_INCLUDED_
#define TRIUM_TOOLS_H_INCLUDED_

#ifndef TRIUM_PARAMETERS_H_INCLUDED_
#include "parameters.hpp"
#endif

/**
 * Collection of tools.
 *
 * This includes the calculation of the 'reduced' (conjugated and
 * suitably normalised) spherical harmonics
 * @f$ \sqrt(\frac{4*\pi}{2\ell + 1}) Y_{\ell m}^\ast @f$.
 *
 * This also includes the storage of calculated values on a regular grid
 * in configuration or Fourier space, where the grid size is determined
 * from an input pameter set.
 *
 * This also includes the set up of isotropic wavenumber or separation
 * bins in configuration or Fourier space.
 *
 */
class ToolCollection {
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
			int ell,
			int m,
			double* pos
		) {
		/// Return the trivial value.
		if (ell == 0 && m == 0) {
			return 1.;
		}

		/// Define zero precision.
		const double eps = 1.e-15;

		/// Calculate magnitudes.
		double xyz_mag_square = 0.;
		for (int axis = 0; axis < 3; axis++) {
			xyz_mag_square += pos[axis] * pos[axis];
		}  // r^2 = x^2 + y^2 + z^2

		double xyz_mag = sqrt(xyz_mag_square);  // r = √(x^2 + y^2 + z^2)

		double xy_mag = sqrt(
			pos[0] * pos[0] + pos[1] * pos[1]
		);  // r_xy = √(x^2 + y^2)

		/// Calculate the angular variable μ = cos(θ).
		double mu = 0.;
		if (fabs(xyz_mag) < eps) {
			return 0.;  // return trivial value
		} else {
			mu = pos[2] / xyz_mag;  // μ = z / r
		}

		/// Calculate the angular variable ϕ.
		double phi = 0.;
		if (fabs(xy_mag) < eps) {
			phi = 0.;  // return default value
		} else {
			phi = acos(pos[0] / xy_mag);  // ϕ = arccos(x / r_xy)
			if (pos[1] < 0.) {
				phi = - phi + 2.*M_PI;  // ϕ ∈ [π, 2π] if y < 0
			}
		}

		/// Calculate spherical harmonics with m >= 0 via the
		/// normalised associated Legendre polynomial, i.e.
		/// Ylm = √((2l + 1)/(4π)) √((l - |m|)!/(l + |m|)!)
		/// * P_l^|m| * exp(imϕ).
		std::complex<double> I_(0., 1.);  // definition of imaginary unit i
		std::complex<double> ylm = 0.;
		ylm = gsl_sf_legendre_sphPlm(ell, abs(m), mu) * exp(I_ * double(m) * phi);

		/// Impose parity and conjugation.
		ylm = pow(-1., (m - abs(m))/2.) * std::conj(ylm);

		/// Normalise to the reduced form.
		ylm *= sqrt(4.*M_PI / (2.*ell + 1.));

		return ylm;
	}

	/**
	 * Store reduced spherical harmonics calculated in Fourier space.
	 *
	 * @param[in] ell_ Degree of the spherical harmonic.
	 * @param[in] m_ Order of the spherical harmonic.
	 * @param[in] params Input parameter set.
	 * @param[out] ylm_out Stored calculated values.
	 * @returns Exit status.
	 */
	static int store_reduced_spherical_harmonic_in_fourier_space(
		int ell_, int m_, ParameterSet& params, std::complex<double>* ylm_out
	) {
		/// Return failure when no storage is provided for calculations.
		if (ylm_out == NULL) {
			return -1;
		}

		/// Determine the fundamental wavenumber in each dimension.
		double dk[3];
		dk[0] = 2.*M_PI / params.boxsize[0];
		dk[1] = 2.*M_PI / params.boxsize[1];
		dk[2] = 2.*M_PI / params.boxsize[2];

		/// Assign a wavevector to each grid cell.
		double kvec[3];
		for (int i = 0; i < params.nmesh[0]; i++) {
			for (int j = 0; j < params.nmesh[1]; j++) {
				for (int k = 0; k < params.nmesh[2]; k++) {
					/// Lay the 'bricks' vertically, then inwards, then to
					/// the right, i.e. along z-axis, y-axis and then x-axis.
					/// The assigned integer coordinate is
					/// (i * nmesh_y * nmesh_z + j * nmesh_z + k)
					/// where nmesh is the mesh number along each axis.
					long long coord_flat = (i * params.nmesh[1] + j) * params.nmesh[2] + k;

					/// Note the origin is at the centre of the mesh grid.
					/// ???: Have the two halves of the grid in each dimension been
					/// rearranged?
					kvec[0] = (i < params.nmesh[0]/2) ? (i * dk[0]) : ((i - params.nmesh[0]) * dk[0]);
					kvec[1] = (j < params.nmesh[1]/2) ? (j * dk[1]) : ((j - params.nmesh[1]) * dk[1]);
					kvec[2] = (k < params.nmesh[2]/2) ? (k * dk[2]) : ((k - params.nmesh[2]) * dk[2]);

					ylm_out[coord_flat] = calc_reduced_spherical_harmonic(ell_, m_, kvec);
				}
			}
		}

		return 0;
	}

	/**
	 * Store reduced spherical harmonics calculated in configuration space.
	 *
	 * @param[in] ell_ Degree of the spherical harmonic.
	 * @param[in] m_ Order of the spherical harmonic.
	 * @param[in] params Input parameter set.
	 * @param[out] ylm_out Stored calculated values.
	 * @returns Exit status.
	 */
	static int store_reduced_spherical_harmonic_in_config_space(
		int ell_, int m_, ParameterSet& params, std::complex<double>* ylm_out
	) {
		/// Return failure when no storage is provided for calculations.
		if (ylm_out == NULL) {
			return -1;
		}

		/// Determine the grid size in each dimension.
		double dr[3];
		dr[0] = params.boxsize[0] / double(params.nmesh[0]);
		dr[1] = params.boxsize[1] / double(params.nmesh[1]);
		dr[2] = params.boxsize[2] / double(params.nmesh[2]);

		/// Assign a position vector to each grid cell.
		double rvec[3];
		for (int i = 0; i < params.nmesh[0]; i++) {
			for (int j = 0; j < params.nmesh[1]; j++) {
				for (int k = 0; k < params.nmesh[2]; k++) {
					/// Lay the 'bricks' vertically, then inwards, then to
					/// the right, i.e. along z-axis, y-axis and then x-axis.
					/// The assigned integer coordinate is
					/// (i * nmesh_y * nmesh_z + j * nmesh_z + k)
					/// where nmesh is the mesh number along each axis.
					long long coord_flat = (i * params.nmesh[1] + j) * params.nmesh[2] + k;

					/// Note the origin is at the centre of the mesh grid.
					/// ???: Have the two halves of the grid in each dimension been
					/// rearranged?
					rvec[0] = (i < params.nmesh[0]/2) ? (i * dr[0]) : ((i - params.nmesh[0]) * dr[0]);
					rvec[1] = (j < params.nmesh[1]/2) ? (j * dr[1]) : ((j - params.nmesh[1]) * dr[1]);
					rvec[2] = (k < params.nmesh[2]/2) ? (k * dr[2]) : ((k - params.nmesh[2]) * dr[2]);

					ylm_out[coord_flat] = calc_reduced_spherical_harmonic(ell_, m_, rvec);
				}
			}
		}

		return 0;
	}

	/**
	 * Set wavenumber bins.
	 *
	 * @param[in] params Input parameter set.
	 * @param[out] kbin_out Set wavenumber bins.
	 * @returns Exit status.
	 */
	static int set_kbin(ParameterSet& params, double* kbin_out) {
		double dk = (params.kmax - params.kmin) / double(params.num_kbin - 1);
		for (int i = 0; i < params.num_kbin; i++) {
			kbin_out[i] = params.kmin + dk * i;
		}
		return 0;
	}

	/**
	 * Set separation bins.
	 *
	 * @param[in] params Input parameter set.
	 * @param[out] rbin_out Set separation bins.
	 * @returns Exit status.
	 */
	static int set_rbin(ParameterSet& params, double* rbin_out) {
		double dr = (params.rmax - params.rmin) / double(params.num_rbin - 1);
		for (int i = 0; i < params.num_rbin; i++) {
			rbin_out[i] = params.rmin + dr * i;
		}
		return 0;
	}
};

#endif
