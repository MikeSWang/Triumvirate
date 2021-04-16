#ifndef TRIUM_TOOLS_H_INCLUDED_
#define TRIUM_TOOLS_H_INCLUDED_

#ifndef TRIUM_PARAMETER_H_INCLUDED_
#include "parameters.hpp"
#endif

/**
 * Collection of tools.
 *
 * This collection includes calculators for the 'reduced' (conjugated
 * and suitably normalised) spherical harmonics
 * @f$ \sqrt(\frac{4*\pi}{2\ell + 1}) Y_{\ell m}^\ast @f$,
 * and storage of calculations performed in configuration/Fourier space.
 *
 */
class ToolClass {
 public:
	/**
	 * Calculate reduced spherical harmonics.
	 *
	 * @param ell Degree of the spherical harmonic.
	 * @param m Order of the spherical harmonic.
	 * @param pos Angular position.
	 * @returns ylm Normalised conjugated value.
	 */
	static std::complex<double> calc_reduced_spherical_harmonic(int ell, int m, double* pos) {
		/// Return the trivial value.
		if (ell == 0 && m == 0) {
			return 1.;
		}

		/// Define constants.
		const int dim = 3;
		const double eps = 1.e-15;

		/// Calculate magnitudes.
		double mag_xyz_square = 0.;
		for (int axis = 0; axis < dim; axis++) {
			mag_xyz_square += pos[axis] * pos[axis];
		}  // r^2 = x^2 + y^2 + z^2

		double mag_xyz = sqrt(mag_xyz_square);  // r = √(x^2 + y^2 + z^2)
		double mag_xy = sqrt(
			pos[0] * pos[0] + pos[1] * pos[1]
		);  // r_xy = √(x^2 + y^2)

		/// Calculate the angular variable μ = cos(θ).
		double mu = 0.;
		if (fabs(mag_xyz) < eps) {
			/// Negligible magnitude, return trivial value.
			return 0.;
		} else {
			/// μ = z / r
			mu = pos[2] / mag_xyz;
		}

		/// Calculate the angular variable ϕ.
		double phi = 0.;
		if (fabs(mag_xy) < eps) {
			/// Negligible partial magnitude, return default value.
			phi = 0.;
		} else {
			/// ϕ = arccos(x / r_xy)
			phi = acos(pos[0] / mag_xy);
		}

		if (pos[1] < 0.) {
			/// ϕ ∈ [π, 2π] if y < 0
			phi = - phi + 2.*M_PI;
		}

		/// Calculate the normal spherical harmonic with m >= 0 via the
		/// normalised associated Legendre polynomial, i.e.
		/// Ylm = √((2l + 1)/(4π)) √((l - |m|)!/(l + |m|)!) * P_l^|m| * exp(imϕ).
		std::complex<double> I_(0., 1.);  // definition of imagary unit i
		std::complex<double> ylm = 0.;
		ylm = gsl_sf_legendre_sphPlm(ell, abs(m), mu) * exp(I_ * double(m) * phi);

		/// Parity and conjugation.
		ylm = pow(-1., (m - abs(m))/2.) * std::conj(ylm);

		/// Renormalise to the reduced form.
		ylm *= sqrt(4.*M_PI / (2.*ell + 1.));

		return ylm;
	}

	/**
	 * Store reduced spherical harmonics calculated in Fourier space.
	 *
	 * @param[in] ell_ Degree of the spherical harmonic.
	 * @param[in] m_ Order of the spherical harmonic.
	 * @param[in] param Input parameter class.
	 * @param[out] ylm_out Stored calculated values.
	 * @returns Exit code.
	 */
	static int store_reduced_spherical_harmonic_in_fourier_space(
		int ell_, int m_, ParameterClass &param, std::complex<double>* ylm_out
	) {
		if (ylm_out == NULL) {
			return -1;
		}

		/// Determine the fundamental wavenumber in each dimension.
		double dk[3];
		dk[0] = 2.*M_PI / param.boxsize[0];
		dk[1] = 2.*M_PI / param.boxsize[1];
		dk[2] = 2.*M_PI / param.boxsize[2];

		/// Assign a wavevector to each grid cell.
		double kvec[3];
		for (int i = 0; i < param.nmesh[0]; i++) {
			for (int j = 0; j < param.nmesh[1]; j++) {
				for (int k = 0; k < param.nmesh[2]; k++) {
					/// Lay the 'bricks' vertically, then inwards, then to the right,
					/// i.e. along z-axis, y-axis and then x-axis.  The assigned
					/// integer coordinate is
					/// 	(i * nmesh_y * nmesh_z + j * nmesh_z + k)
					/// where nmesh is the mesh number along each axis.
					long long coord = (i * param.nmesh[1] + j) * param.nmesh[2] + k;

					/// Note the origin is at the centre of the mesh grid.
					/// ???: Have the two halves of the grid in each dimension been
					/// rearranged?
					kvec[0] = (i < param.nmesh[0]/2) ? (i * dk[0]) : ((i - param.nmesh[0]) * dk[0]);
					kvec[1] = (j < param.nmesh[1]/2) ? (j * dk[1]) : ((j - param.nmesh[1]) * dk[1]);
					kvec[2] = (k < param.nmesh[2]/2) ? (k * dk[2]) : ((k - param.nmesh[2]) * dk[2]);

					ylm_out[coord] = calc_reduced_spherical_harmonic(ell_, m_, kvec);
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
	 * @param[in] param Input parameter class.
	 * @param[out] ylm_out Stored calculated values.
	 * @returns Exit code.
	 */
	static int store_reduced_spherical_harmonic_in_config_space(
		int ell_, int m_, ParameterClass &param, std::complex<double>* ylm_out
	) {
		if (ylm_out == NULL) {
			return -1;
		}

		/// Determine the grid size in each dimension.
		double dr[3];
		dr[0] = param.boxsize[0] / double(param.nmesh[0]);
		dr[1] = param.boxsize[1] / double(param.nmesh[1]);
		dr[2] = param.boxsize[2] / double(param.nmesh[2]);

		/// Assign a position vector to each grid cell.
		double rvec[3];
		for (int i = 0; i < param.nmesh[0]; i++) {
			for (int j = 0; j < param.nmesh[1]; j++) {
				for (int k = 0; k < param.nmesh[2]; k++) {
					/// Lay the 'bricks' vertically, then inwards, then to the right,
					/// i.e. along z-axis, y-axis and then x-axis.  The assigned
					/// integer coordinate is
					/// 	(i * Nmesh_y * Nmesh_z + j * Nmesh_z + k)
					/// where Nmesh is the mesh number along each axis.
					long long coord = (i * param.nmesh[1] + j) * param.nmesh[2] + k;

					/// Note the origin is at the centre of the mesh grid.
					/// ???: Have the two halves of the grid in each dimension been
					/// rearranged?
					rvec[0] = (i < param.nmesh[0]/2) ? (i * dr[0]) : ((i - param.nmesh[0]) * dr[0]);
					rvec[1] = (j < param.nmesh[1]/2) ? (j * dr[1]) : ((j - param.nmesh[1]) * dr[1]);
					rvec[2] = (k < param.nmesh[2]/2) ? (k * dr[2]) : ((k - param.nmesh[2]) * dr[2]);

					ylm_out[coord] = calc_reduced_spherical_harmonic(ell_, m_, rvec);
				}
			}
		}

		return 0;
	}

	/**
	 * Set wavenumber bins.
	 *
	 * @param[in] param Input parameter class.
	 * @param[out] kbin_out Set wavenumber bins.
	 * @returns Exit code.
	 */
	static int set_kbin(ParameterClass &param, double* kbin_out) {
		double dk = (param.kmax - param.kmin) / double(param.num_kbin - 1);
		for (int i = 0; i < param.num_kbin; i++) {
			kbin_out[i] = param.kmin + dk * i;
		}
		return 0;
	}

	/**
	 * Set separation bins.
	 *
	 * @param[in] param Input parameter class.
	 * @param[out] rbin_out Set separation bins.
	 * @returns Exit code.
	 */
	static int set_rbin(ParameterClass &param, double* rbin_out) {
		double dr = (param.rmax - param.rmin) / double(param.num_rbin - 1);
		for (int i = 0; i < param.num_rbin; i++) {
			rbin_out[i] = param.rmin + dr * i;
		}
		return 0;
	}

 private:
};

#endif
