#ifndef TRIUM_BESSEL_H_INCLUDED_
#define TRIUM_BESSEL_H_INCLUDED_

/**
 * Interpolated spherical Bessel function @f$ j_\ell(x) @f$
 * of the first kind.
 *
 */
class SphericalBesselClass {
 public:
	/**
	 * Initialise the interpolated function.
	 *
	 * @param ell Order @f$ \ell @f$ of the spherical Bessel function.
	 */
	SphericalBesselClass(int ell) {
		/// Set up sampling range and number and calculate sample spacing.
		double xmin = 0.;
		double xmax = 10000.;  // NOTE: discretionary
		int nsample = 1000000;  // NOTE: discretionary

		double dx = (xmax - xmin) / double(nsample - 1);

		/// Initialise and evaluate sample points.
		double* x = new double[nsample];
		double* j_ell = new double[nsample];

		for (int i = 0; i < nsample; i++) {
			x[i] = xmin + dx * double(i);
			j_ell[i] = gsl_sf_bessel_jl(ell, x[i]);
		}

		/// Initialise interpolator.
		this->accel = gsl_interp_accel_alloc();
		this->spline = gsl_spline_alloc(gsl_interp_cspline, nsample);

		gsl_spline_init(this->spline, x, j_ell, nsample);

		/// Delete sample points.
		delete[] x; x = NULL;
		delete[] j_ell; j_ell = NULL;
	}

	/**
	 * Deconstruct interpolators.
	 */
	~SphericalBesselClass() {
		if (this->accel != NULL) {
			gsl_interp_accel_free(this->accel); this->accel = NULL;
		}

		if (this->spline != NULL) {
			gsl_spline_free(this->spline); this->spline = NULL;
		}
	}

	/**
	 * Evaluate the interpolated function.
	 *
	 * @param x Argument of the function.
	 * @returns Value of the function.
	 */
	double eval(double x) {
		return gsl_spline_eval(this->spline, x, this->accel);
	}

 private:
	gsl_interp_accel* accel;
	gsl_spline* spline;
};

#endif
