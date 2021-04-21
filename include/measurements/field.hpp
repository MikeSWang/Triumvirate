#ifndef TRIUM_FIELD_H_INCLUDED_
#define TRIUM_FIELD_H_INCLUDED_

/**
 * Density field object.
 *
 * This is based on a particle container
 */
template <class TemplateParticle>
class DensityFieldClass {
public:
	/*****************************************/
	fftw_complex * field;
	/*****************************************/
	const fftw_complex& operator[](int id) { return this->field[id]; }

	DensityFieldClass(ParameterSet & params) {

		this->params = params;

		this->field = NULL;
		this->field = fftw_alloc_complex(this->params.nmesh_tot);
		bytes += double( sizeof(fftw_complex) * this->params.nmesh_tot / 1024.0 / 1024.0 / 1024.0);

		for (int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] = 0.0;
			this->field[i][1] = 0.0;
		}
	}

	~DensityFieldClass() {
		finalise_density_field();
	}

	void finalise_density_field() {
		if(this->field != NULL) {
			fftw_free(this->field); this->field = NULL;
			bytes -= double( sizeof(fftw_complex) * this->params.nmesh_tot / 1024.0 / 1024.0 / 1024.0);
		}
	}

	int calcField(TemplateParticle & particle, fftw_complex * weight) {
		if(0) {
		} else if (this->params.assignment == "NGP") {
			this->calcFieldNGP(particle, weight);
		} else if (this->params.assignment == "CIC") {
			this->calcFieldCIC(particle, weight);
		} else if (this->params.assignment == "TSC") {
			this->calcFieldTSC(particle, weight);
		} else {
			return -1;
		}

		return 0;
	}

	int calcFieldNGP(TemplateParticle & particle, fftw_complex * weight) {

		/* rho = sum_i w_i delta_D(x - x_i) */

		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] = 0.0;
			this->field[i][1] = 0.0;
		}

		/* delta_D = (1/dV) * delta_K, dV = volume/nmesh */
		double dV = this->params.volume / double(this->params.nmesh_tot);
		double factor = 1.0 / dV;

		for(int p = 0; p < particle.n_tot; p++) {
			double w[1][3];
			int iw[1][3];
			for(int axes = 0; axes < 3; axes++) {
				double xp = particle[p].pos[axes] * double(this->params.nmesh[axes]) / this->params.boxsize[axes];
				iw[0][axes] = int(xp + 0.5);
				w[0][axes] = 1.0;
			}

			for(int i = 0; i < 1; i++) {
			for(int j = 0; j < 1; j++) {
			for(int k = 0; k < 1; k++) {
				long long coord = ( iw[i][0] * this->params.nmesh[1] + iw[j][1] ) * this->params.nmesh[2] + iw[k][2];
				if( (coord >= 0) && (coord < this->params.nmesh_tot) ) {
					this->field[coord][0] += weight[p][0] * factor * w[i][0] * w[j][1] * w[k][2];
					this->field[coord][1] += weight[p][1] * factor * w[i][0] * w[j][1] * w[k][2];
				}
			}}}
		}

		return 0;
	}


	int calcFieldCIC(TemplateParticle & particle, fftw_complex * weight) {

		/* rho = sum_i w_i delta_D(x - x_i) */

		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] = 0.0;
			this->field[i][1] = 0.0;
		}

		/* delta_D = (1/dV) * delta_K, dV = volume/nmesh */
		double dV = this->params.volume / double(this->params.nmesh_tot);
		double factor = 1.0 / dV;

		for(int p = 0; p < particle.n_tot; p++) {
			double w[2][3];
			int iw[2][3];
			for(int axes = 0; axes < 3; axes++) {
				double xp = particle[p].pos[axes] * double(this->params.nmesh[axes]) / this->params.boxsize[axes];
				iw[0][axes] = int(xp);
				double dx = xp - double(iw[0][axes]);
				w[0][axes] = 1.0 - dx;
				w[1][axes] = dx;
				iw[1][axes] = iw[0][axes] + 1;
			}

			for(int i = 0; i < 2; i++) {
			for(int j = 0; j < 2; j++) {
			for(int k = 0; k < 2; k++) {
				long long coord = ( iw[i][0] * this->params.nmesh[1] + iw[j][1] ) * this->params.nmesh[2] + iw[k][2];
				if( (coord >= 0) && (coord < this->params.nmesh_tot) ) {
					this->field[coord][0] += weight[p][0] * factor * w[i][0] * w[j][1] * w[k][2];
					this->field[coord][1] += weight[p][1] * factor * w[i][0] * w[j][1] * w[k][2];
				}
			}}}
		}

		return 0;
	}



	int calcFieldTSC(TemplateParticle & particle, fftw_complex * weight) {

		/* rho = sum_i w_i delta_D(x - x_i) */

		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] = 0.0;
			this->field[i][1] = 0.0;
		}

		/* delta_D = (1/dV) * delta_K, dV = volume/nmesh */
		double dV = this->params.volume / double(this->params.nmesh_tot);
		double factor = 1.0 / dV;

		for(int p = 0; p < particle.n_tot; p++) {
			double w[3][3];
			int iw[3][3];
			for(int axes = 0; axes < 3; axes++) {
				double xp = particle[p].pos[axes] * double(this->params.nmesh[axes]) / this->params.boxsize[axes];
				iw[1][axes] = int(xp + 0.5);
				double dx = xp - double(iw[1][axes]);
				w[0][axes] = 0.5 * (0.5 - dx) * (0.5 - dx);
				w[1][axes] = 0.75 - dx * dx;
				w[2][axes] = 0.5 * (0.5 + dx) * (0.5 + dx);
				iw[0][axes] = iw[1][axes] - 1;
				iw[2][axes] = iw[1][axes] + 1;
			}

			for(int i = 0; i < 3; i++) {
			for(int j = 0; j < 3; j++) {
			for(int k = 0; k < 3; k++) {
				long long coord = ( iw[i][0] * this->params.nmesh[1] + iw[j][1] ) * this->params.nmesh[2] + iw[k][2];
				if( (coord >= 0) && (coord < this->params.nmesh_tot) ) {
					this->field[coord][0] += weight[p][0] * factor * w[i][0] * w[j][1] * w[k][2];
					this->field[coord][1] += weight[p][1] * factor * w[i][0] * w[j][1] * w[k][2];
				}
			}}}
		}

		return 0;
	}

	int calcYlmWeightedDensityFluctuation(
            TemplateParticle & P_D, TemplateParticle & P_R,
            LineOfSight* los_D, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		DensityFieldClass<ParticleBOSSClass> n_R(this->params);
		fftw_complex * weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_D.n_tot);
		for (int p = 0; p < P_D.n_tot; p++) {
            double los[3] = {los_D[p].pos[0], los_D[p].pos[1], los_D[p].pos[2]};
            std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
            weight[p][0] = Ylm.real() * P_D[p].w;
            weight[p][1] = Ylm.imag() * P_D[p].w;
		}

		this->calcField(P_D, weight);
		fftw_free(weight); weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_R.n_tot);
		for (int p = 0; p < P_R.n_tot; p++) {
            double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
            std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
            weight[p][0] = Ylm.real() * P_R[p].w;
            weight[p][1] = Ylm.imag() * P_R[p].w;
		}

		n_R.calcField(P_R, weight);
		fftw_free(weight); weight = NULL;

		/****/
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] -= alpha * n_R.field[i][0];
			this->field[i][1] -= alpha * n_R.field[i][1];
		}

		return 0;
	}


	int calcYlmWeightedDensityFluctuationForBispectrumShotnoise(
            TemplateParticle & P_D, TemplateParticle & P_R,
            LineOfSight* los_D, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		DensityFieldClass<ParticleBOSSClass> n_R(this->params);
		fftw_complex * weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_D.n_tot);
		for (int p = 0; p < P_D.n_tot; p++) {
			double los[3] = {los_D[p].pos[0], los_D[p].pos[1], los_D[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			Ylm = std::conj(Ylm); // IMPORTANT! //
			weight[p][0] = Ylm.real() * pow(P_D[p].w, 2);
		       	weight[p][1] = Ylm.imag() * pow(P_D[p].w, 2);
		}
		this->calcField(P_D, weight);
		fftw_free(weight); weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_R.n_tot);
		for (int p = 0; p < P_R.n_tot; p++) {
			double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			Ylm = std::conj(Ylm); // IMPORTANT! //
			weight[p][0] = Ylm.real() * pow(P_R[p].w, 2);
		       	weight[p][1] = Ylm.imag() * pow(P_R[p].w, 2);
		}
		n_R.calcField(P_R, weight);
		fftw_free(weight); weight = NULL;

		/****/
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] += alpha * alpha * n_R.field[i][0];
			this->field[i][1] += alpha * alpha * n_R.field[i][1];
		}

		return 0;
	}

	int calcNormalDensityFluctuation_in_box(TemplateParticle & P_D, ParameterSet & params) {

		fftw_complex * weight = NULL;
		weight = fftw_alloc_complex(P_D.n_tot);
		for (int p = 0; p < P_D.n_tot; p++) {
			weight[p][0] = 1.0;
		       	weight[p][1] = 0.0;
		}
		this->calcField(P_D, weight);
		fftw_free(weight); weight = NULL;

		/****/
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] -= (double(P_D.n_tot)/params.volume);
			this->field[i][1] -= 0.0;
		}

		return 0;
	}

	int calcNormalDensityForBispectrumShotnoise_in_box(TemplateParticle & P_D) {

		fftw_complex * weight = NULL;
		weight = fftw_alloc_complex(P_D.n_tot);
		for (int p = 0; p < P_D.n_tot; p++) {
			weight[p][0] = 1.0;
		       	weight[p][1] = 0.0;
		}
		this->calcField(P_D, weight);
		fftw_free(weight); weight = NULL;

		return 0;
	}

	int calcYlmWeightedMeanDensity(
            TemplateParticle & P_R, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		fftw_complex * weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_R.n_tot);
		for (int p = 0; p < P_R.n_tot; p++) {
			double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			weight[p][0] = Ylm.real() * P_R[p].w;
		       	weight[p][1] = Ylm.imag() * P_R[p].w;
		}
		this->calcField(P_R, weight);
		fftw_free(weight); weight = NULL;

		/****/
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] *= alpha;
			this->field[i][1] *= alpha;
		}

		return 0;
	}

	int calcYlmWeightedMeanDensityForThreePointWindowFunctionShotnoise(
            TemplateParticle & P_R, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		fftw_complex * weight = NULL;
		/****/
		weight = fftw_alloc_complex(P_R.n_tot);
		for (int p = 0; p < P_R.n_tot; p++) {
			double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			Ylm = std::conj(Ylm); // IMPORTANT! //
			weight[p][0] = Ylm.real() * pow(P_R[p].w, 2);
		       	weight[p][1] = Ylm.imag() * pow(P_R[p].w, 2);
		}
		this->calcField(P_R, weight);
		fftw_free(weight); weight = NULL;

		/****/
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] *= alpha * alpha;
			this->field[i][1] *= alpha * alpha;
		}

		return 0;
	}


	int calcFourierTransform() {

		/* d^3x = (dV) \sum_i */
		double dV = this->params.volume / double(this->params.nmesh_tot);

		for(int i = 0;  i < this->params.nmesh_tot; i++) {
			this->field[i][0] *= dV;
			this->field[i][1] *= dV;
		}

		/* FFT */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0],this->params.nmesh[1],this->params.nmesh[2], this->field, this->field, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		return 0;
	}


	int calcAssignmentFunctionCorrection() {

		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];

		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;

			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];

			// assignment function.
			double	wf = this->WindowFunction(kvec);

			this->field[coord][0] /=wf;
			this->field[coord][1] /=wf;
		}}}

		return 0;
	}

	int calcInverseFourierTransform() {

		/* d^3k/(2pi)^3 = (1/V) \sum_i */
		double V = this->params.volume;

		for(int i = 0;  i < this->params.nmesh_tot; i++) {
			this->field[i][0] /= V;
			this->field[i][1] /= V;
		}

		/* Inverse FFT */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0], this->params.nmesh[1], this->params.nmesh[2], this->field, this->field, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		return 0;
	}

	int calcInverseFourierTransformForBispectrum(DensityFieldClass & density, double kmag_in, double dk_in, std::complex<double> * Ylm) {

		/* initialize */
		for (int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] = 0.0;
			this->field[i][1] = 0.0;
		}

		int n_mode = 0;
		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];

		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;

			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];
			double kmag = sqrt( kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2] );

			double lower = ( kmag_in > dk_in*0.5 ) ? (kmag_in - dk_in * 0.5) : 0.0;
			double upper = kmag_in + dk_in*0.5;
			if( (kmag > lower) && (kmag <= upper) ) {

				std::complex<double> dn(density[coord][0], density[coord][1]);

				// the window function.
				double	wf = this->WindowFunction(kvec);
				dn /= wf;

				this->field[coord][0] = ( Ylm[coord] * dn ).real();
				this->field[coord][1] = ( Ylm[coord] * dn ).imag();
				n_mode++;
			} else {
				this->field[coord][0] = 0.0;
				this->field[coord][1] = 0.0;
			}
		}}}

		/* Inverse FFTW */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0], this->params.nmesh[1], this->params.nmesh[2], this->field, this->field, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		for(int coord = 0; coord < this->params.nmesh_tot; coord++) {
			this->field[coord][0] /=  double(n_mode);
			this->field[coord][1] /=  double(n_mode);
		}

		return 0;
	}

	int calcInverseFourierTransformForThreePointFunction(DensityFieldClass & density, double rmag_in, std::complex<double> * Ylm,
			SphericalBesselCalculator & sj ) {

		/* initialize */
		for (int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] = 0.0;
			this->field[i][1] = 0.0;
		}

		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];

		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;

			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];
			double kmag = sqrt( kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2] );

			std::complex<double> dn(density[coord][0], density[coord][1]);

			// the window function.
			double	wf = this->WindowFunction(kvec);
			dn /= wf;

			this->field[coord][0] = sj.eval(kmag * rmag_in) * ( Ylm[coord] * dn ).real() / this->params.volume;
			this->field[coord][1] = sj.eval(kmag * rmag_in) * ( Ylm[coord] * dn ).imag() / this->params.volume;
		}}}

		/* Inverse FFTW */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0], this->params.nmesh[1], this->params.nmesh[2], this->field, this->field, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		return 0;
	}


	double WindowFunction(double * kvec) {
		if(0) {
		} else if (this->params.assignment == "NGP") {
			return this->WindowFunctionNGP(kvec);
		} else if (this->params.assignment == "CIC") {
			return this->WindowFunctionCIC(kvec);
		} else if (this->params.assignment == "TSC") {
			return this->WindowFunctionTSC(kvec);
		} else {
			return 1.0;
		}

		return 1.0;

	}
	double WindowFunctionNGP(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = (int) ( kvec[0] / dk[0] + 1.0e-5);
		int j = (int) ( kvec[1] / dk[1] + 1.0e-5);
		int k = (int) ( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1) / x1 : 1.0;
		double f2 = (j != 0) ? sin(x2) / x2 : 1.0;
		double f3 = (k != 0) ? sin(x3) / x3 : 1.0;

		double ret = f1 * f2 * f3;
		return pow(ret,1);  // W = ret**3
	}


	double WindowFunctionCIC(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = (int) ( kvec[0] / dk[0] + 1.0e-5);
		int j = (int) ( kvec[1] / dk[1] + 1.0e-5);
		int k = (int) ( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1) / x1 : 1.0;
		double f2 = (j != 0) ? sin(x2) / x2 : 1.0;
		double f3 = (k != 0) ? sin(x3) / x3 : 1.0;

		double ret = f1 * f2 * f3;
		return pow(ret,2);  // W = ret**3
	}



	double WindowFunctionTSC(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = (int) ( kvec[0] / dk[0] + 1.0e-5);
		int j = (int) ( kvec[1] / dk[1] + 1.0e-5);
		int k = (int) ( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1) / x1 : 1.0;
		double f2 = (j != 0) ? sin(x2) / x2 : 1.0;
		double f3 = (k != 0) ? sin(x3) / x3 : 1.0;

		double ret = f1 * f2 * f3;
		return pow(ret,3);  // W = ret**3
	}


	double calc_survey_volume(TemplateParticle & P_R) {

		/****/
		fftw_complex * weight = NULL;
		weight = fftw_alloc_complex(P_R.n_tot);
		for (int p = 0; p < P_R.n_tot; p++) {
			weight[p][0] = 1.0;
		       	weight[p][1] = 0.0;
		}
		this->calcField(P_R, weight);
		fftw_free(weight); weight = NULL;

		/****/
		double dV = this->params.volume / double(this->params.nmesh_tot);
		double norm = 0.0;
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			norm += this->field[i][0] * this->field[i][0] * dV;
		}

		double survey_volume = double(P_R.n_tot) * double(P_R.n_tot) / norm;

		return survey_volume;

	}

	int calcNormalDensityFluctuation_in_boxForReconstruction(TemplateParticle & P_D, TemplateParticle & P_R, double alpha) {

		DensityFieldClass<ParticleBOSSClass> n_R(this->params);
		fftw_complex * weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_D.n_tot);
		for (int p = 0; p < P_D.n_tot; p++) {
			weight[p][0] = 1.0;
		       	weight[p][1] = 0.0;
		}
		this->calcField(P_D, weight);
		fftw_free(weight); weight = NULL;

		/****/
		weight = fftw_alloc_complex(P_R.n_tot);
		for (int p = 0; p < P_R.n_tot; p++) {
			weight[p][0] = 1.0;
		       	weight[p][1] = 0.0;
		}
		n_R.calcField(P_R, weight);
		fftw_free(weight); weight = NULL;

		/****/
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			this->field[i][0] -= alpha * n_R.field[i][0];
			this->field[i][1] -= alpha * n_R.field[i][1];
		}

		return 0;
	}

 private:
	ParameterSet params;
};

template <class TemplateParticle>
class TwoPointStatisticsClass {
private:
	ParameterSet params;
public:
	std::complex<double> * pk;
	std::complex<double> * xi;
	int * n_mode_pk;
	int * n_mode_xi;

	TwoPointStatisticsClass(ParameterSet & params){

		this->params = params;
		this->pk = new std::complex<double>[params.num_kbin];
		this->n_mode_pk = new int[params.num_kbin];
	       	for(int i = 0; i < params.num_kbin; i++) {
			this->pk[i] = 0.0;
			this->n_mode_pk[i] = 0.0;
		}

		this->xi = new std::complex<double>[params.num_rbin];
		this->n_mode_xi = new int[params.num_rbin];
	       	for(int i = 0; i < params.num_rbin; i++) {
			this->xi[i] = 0.0;
			this->n_mode_xi[i] = 0.0;
		}


	}

	~TwoPointStatisticsClass(){
		finalizeTwoPointStatistics();
	}

	void finalizeTwoPointStatistics() {
		delete[] this->pk; this->pk = NULL;
		delete[] this->n_mode_pk; this->n_mode_pk = NULL;
		delete[] this->xi; this->xi = NULL;
		delete[] this->n_mode_xi; this->n_mode_xi = NULL;
	}

	std::complex<double> calcShotNoiseForPowerspectrum(
            TemplateParticle & P_D, TemplateParticle & P_R,
            LineOfSight* los_D, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		std::complex<double> sum_D = 0.0;
		std::complex<double> sum_R = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			double los[3] = {los_D[p].pos[0], los_D[p].pos[1], los_D[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			sum_D += pow(P_D[p].w,2) * Ylm;
		}

		for(int p = 0; p < P_R.n_tot; p++) {
			double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			sum_R += pow(P_R[p].w,2) * Ylm;
		}

		return sum_D + pow(alpha,2) * sum_R;

	}

	std::complex<double> calcShotNoiseForTwoPointWindowFunction(
            TemplateParticle & P_R, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		std::complex<double> sum_R = 0.0;
		for(int p = 0; p < P_R.n_tot; p++) {
			double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			sum_R += pow(P_R[p].w,2) * Ylm;
		}

		return pow(alpha,2) * sum_R;

	}

	std::complex<double> calcShotNoiseForBispectrum(
            TemplateParticle & P_D, TemplateParticle & P_R,
            LineOfSight* los_D, LineOfSight* los_R,
            double alpha, int _ELL_, int _M_) {

		std::complex<double> sum_D = 0.0;
		std::complex<double> sum_R = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			double los[3] = {los_D[p].pos[0], los_D[p].pos[1], los_D[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			sum_D += pow(P_D[p].w,3) * Ylm;
		}

		for(int p = 0; p < P_R.n_tot; p++) {
			double los[3] = {los_R[p].pos[0], los_R[p].pos[1], los_R[p].pos[2]};
			std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, los);
			sum_R += pow(P_R[p].w,3) * Ylm;
		}

		return sum_D - pow(alpha,3) * sum_R;

	}


	int calc_power_spec(DensityFieldClass<TemplateParticle> & density1, DensityFieldClass<TemplateParticle> & density2,
		        	   double * kbin, std::complex<double> shotnoise, int _ELL_, int _M_) {

		double dk_bin = kbin[1] - kbin[0];

		int n_temp = 100000;
		double dk_temp = 0.0001;
		std::complex<double> * pk_temp = new std::complex<double>[n_temp];
		int * n_mode_temp = new int[n_temp];
	    for(int i = 0; i < n_temp; i++) {
			pk_temp[i] = 0.0;
			n_mode_temp[i] = 0.0;
		}
	    for(int i = 0; i < this->params.num_kbin; i++) {
			this->pk[i] = 0.0;
			this->n_mode_pk[i] = 0.0;
		}

		/* 場の値の代入*/
		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;
			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];
			double kmag = sqrt( kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2] );
			int kindex = int( kmag / dk_temp + 0.5);
			if(kindex < n_temp) {

				// delta^2 = (2pi)^3 delta_D(k1+k2) P(k), where (2pi)^3 delta_D(k1+k2) = V //
				// P = (V/N^2) |n(k)|^2 //
				std::complex<double> f1(density1[coord][0], density1[coord][1]);
				std::complex<double> f2(density2[coord][0], density2[coord][1]);
				std::complex<double> pk2 = f1 * conj(f2);

				/* subtract shotnoise */
				pk2 -= ( shotnoise * ShotnoiseFunction(kvec) );

				// assignment function //.
				double	wf = WindowFunction(kvec);
				pk2 /= pow(wf,2);

				// spherical harmonics //
				std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, kvec);
				pk2 *= Ylm;

				pk_temp[kindex] += pk2;
				n_mode_temp[kindex]++;
			}

		}}}

		for(int j = 0; j < this->params.num_kbin; j++) {
			double lower = ( kbin[j] > dk_bin*0.5 ) ? (kbin[j] - dk_bin * 0.5) : 0.0;
			double upper = kbin[j] + dk_bin*0.5;
			for(int i = 0; i < n_temp; i++) {
				double kmag_temp = double(i) * dk_temp;
				if( (kmag_temp > lower) && (kmag_temp <= upper) ) {
					this->pk[j] += pk_temp[i];
					this->n_mode_pk[j] += n_mode_temp[i];
				}
			}
		}

		for(int i = 0; i < this->params.num_kbin; i++) {
			if(this->n_mode_pk[i] != 0) {
				this->pk[i] /= double(this->n_mode_pk[i]);
			} else {
				this->pk[i] = 0.0;
			}
		}

		delete[] pk_temp;
		delete[] n_mode_temp;

		return 0;

	}


	int calcCorrelationFunction(DensityFieldClass<TemplateParticle> & density1, DensityFieldClass<TemplateParticle> & density2,
		        	 double * rbin, std::complex<double> shotnoise, int _ELL_, int _M_) {

		fftw_complex * xi3d_temp = fftw_alloc_complex(this->params.nmesh_tot);
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			xi3d_temp[i][0] = 0.0;
			xi3d_temp[i][1] = 0.0;
		}

		/* dk^3 = (1/V) sum_i */
		double factor  = 1.0 / this->params.volume;
		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;
			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];

			// delta^2 = (2pi)^3 delta_D(k1+k2) P(k), where (2pi)^3 delta_D(k1+k2) = V //
			// P = (V/N^2) |n(k)|^2 //
			std::complex<double> f1(density1[coord][0], density1[coord][1]);
			std::complex<double> f2(density2[coord][0], density2[coord][1]);
			std::complex<double> pk2 = f1 * conj(f2);

			/* subtract shotnoise */
			pk2 -= ( shotnoise * ShotnoiseFunction(kvec) );

			/* assiginment function.*/
			double	wf = WindowFunction(kvec);
			pk2 /= pow(wf,2);

			xi3d_temp[coord][0] = factor * pk2.real();
			xi3d_temp[coord][1] = factor * pk2.imag();

		}}}

		/* Inverse FFTW */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0], this->params.nmesh[1], this->params.nmesh[2], xi3d_temp, xi3d_temp, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		int n_temp = 10000;
		double dr_temp = 0.5;
		std::complex<double> * xi_temp = new std::complex<double>[n_temp];
		int * n_mode_temp = new int[n_temp];
	       	for(int i = 0; i < n_temp; i++) {
			xi_temp[i] = 0.0;
			n_mode_temp[i] = 0.0;
		}
	       	for(int i = 0; i < this->params.num_rbin; i++) {
			this->xi[i] = 0.0;
			this->n_mode_xi[i] = 0.0;
		}

		double rvec[3];
		double dr[3];
		dr[0] = this->params.boxsize[0] / double(this->params.nmesh[0]);
		dr[1] = this->params.boxsize[1] / double(this->params.nmesh[1]);
		dr[2] = this->params.boxsize[2] / double(this->params.nmesh[2]);
		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;
			rvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dr[0] : (double) (i - this->params.nmesh[0]) * dr[0];
			rvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dr[1] : (double) (j - this->params.nmesh[1]) * dr[1];
			rvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dr[2] : (double) (k - this->params.nmesh[2]) * dr[2];
			double rmag = sqrt( rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2] );
			int rindex = (int) (rmag / dr_temp + 0.5);
			if(rindex < n_temp) {
				std::complex<double> xi2(xi3d_temp[coord][0], xi3d_temp[coord][1]);

				// spherical harmonics //
				std::complex<double> Ylm = ToolCollection::calc_reduced_spherical_harmonic(_ELL_, _M_, rvec);
				xi2 *= Ylm;

				xi_temp[rindex] += xi2;
				n_mode_temp[rindex]++;

			}
		}}}

		/*****************************************/
		double dr_bin = rbin[1] - rbin[0];

		for(int j = 0; j < this->params.num_rbin; j++) {
			double lower = ( rbin[j] > dr_bin*0.5 ) ? (rbin[j] - dr_bin * 0.5) : 0.0;
			double upper = rbin[j] + dr_bin*0.5;
			for(int i = 0; i < n_temp; i++) {
				double rmag_temp = double(i) * dr_temp;
				if( (rmag_temp > lower) && (rmag_temp <= upper) ) {
					this->xi[j] += xi_temp[i];
					this->n_mode_xi[j] += n_mode_temp[i];
				}
			}
		}

		for(int i = 0; i < this->params.num_rbin; i++) {
			if(this->n_mode_xi[i] != 0) {
				this->xi[i] /= double(this->n_mode_xi[i]);
			} else {
				this->xi[i] = 0.0;
			}
		}


		delete[] xi3d_temp;
		delete[] xi_temp;
		delete[] n_mode_temp;

		return 0;

	}

	int calcShotNoiseForBispectrum_ijk(DensityFieldClass<TemplateParticle> & density1, DensityFieldClass<TemplateParticle> & density2,
		        	           std::complex<double> shotnoise, int _ELL_, int _M_, fftw_complex * xi) {


		/* dk^3 = (1/V) sum_i */
		double factor  = 1.0 / this->params.volume;
		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;
			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];

			// delta^2 = (2pi)^3 delta_D(k1+k2) P(k), where (2pi)^3 delta_D(k1+k2) = V //
			// P = (V/N^2) |n(k)|^2 //
			std::complex<double> f1(density1[coord][0], density1[coord][1]);
			std::complex<double> f2(density2[coord][0], density2[coord][1]);
			std::complex<double> pk2 = f1 * conj(f2);

			/* subtract shotnoise */
			pk2 -= ( shotnoise * ShotnoiseFunction(kvec) );

			// window function.
			double	wf = WindowFunction(kvec);
			pk2 /= pow(wf,2);

			xi[coord][0] = factor * pk2.real();
			xi[coord][1] = factor * pk2.imag();

		}}}

		/* Inverse FFTW */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0], this->params.nmesh[1], this->params.nmesh[2], xi, xi, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		return 0;

	}

	int calcCorrelationFunctionForThreePointFunction(DensityFieldClass<TemplateParticle> & density1, DensityFieldClass<TemplateParticle> & density2,
		        	 double * rbin, std::complex<double> shotnoise, int _ELL_, int _M_,
				 std::complex<double> * Ylm1, std::complex<double> * Ylm2) {

		fftw_complex * xi3d_temp = fftw_alloc_complex(this->params.nmesh_tot);
		for(int i = 0; i < this->params.nmesh_tot; i++) {
			xi3d_temp[i][0] = 0.0;
			xi3d_temp[i][1] = 0.0;
		}

		/* dk^3 = (1/V) sum_i */
		double factor  = 1.0 / this->params.volume;
		double kvec[3];
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;
			kvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - this->params.nmesh[0]) * dk[0];
			kvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - this->params.nmesh[1]) * dk[1];
			kvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - this->params.nmesh[2]) * dk[2];

			// delta^2 = (2pi)^3 delta_D(k1+k2) P(k), where (2pi)^3 delta_D(k1+k2) = V //
			// P = (V/N^2) |n(k)|^2 //
			std::complex<double> f1(density1[coord][0], density1[coord][1]);
			std::complex<double> f2(density2[coord][0], density2[coord][1]);
			std::complex<double> pk2 = f1 * conj(f2);

			/* subtract shotnoise */
			pk2 -= ( shotnoise * ShotnoiseFunction(kvec) );

			/* assiginment function.*/
			double	wf = WindowFunction(kvec);
			pk2 /= pow(wf,2);

			xi3d_temp[coord][0] = factor * pk2.real();
			xi3d_temp[coord][1] = factor * pk2.imag();

		}}}

		/* Inverse FFTW */
		fftw_plan plan = fftw_plan_dft_3d(this->params.nmesh[0], this->params.nmesh[1], this->params.nmesh[2], xi3d_temp, xi3d_temp, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);

		int n_temp = 10000;
		double dr_temp = 0.5;
		std::complex<double> * xi_temp = new std::complex<double>[n_temp];
		int * n_mode_temp = new int[n_temp];
	       	for(int i = 0; i < n_temp; i++) {
			xi_temp[i] = 0.0;
			n_mode_temp[i] = 0.0;
		}
	       	for(int i = 0; i < this->params.num_rbin; i++) {
			this->xi[i] = 0.0;
			this->n_mode_xi[i] = 0.0;
		}

		double rvec[3];
		double dr[3];
		dr[0] = this->params.boxsize[0] / double(this->params.nmesh[0]);
		dr[1] = this->params.boxsize[1] / double(this->params.nmesh[1]);
		dr[2] = this->params.boxsize[2] / double(this->params.nmesh[2]);
		for(int i = 0; i < this->params.nmesh[0]; i++) {
		for(int j = 0; j < this->params.nmesh[1]; j++) {
		for(int k = 0; k < this->params.nmesh[2]; k++) {
			long long coord = ( i * this->params.nmesh[1] + j ) * this->params.nmesh[2] + k;
			rvec[0] = (i < this->params.nmesh[0]/2) ? (double) i * dr[0] : (double) (i - this->params.nmesh[0]) * dr[0];
			rvec[1] = (j < this->params.nmesh[1]/2) ? (double) j * dr[1] : (double) (j - this->params.nmesh[1]) * dr[1];
			rvec[2] = (k < this->params.nmesh[2]/2) ? (double) k * dr[2] : (double) (k - this->params.nmesh[2]) * dr[2];
			double rmag = sqrt( rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2] );
			int rindex = (int) (rmag / dr_temp + 0.5);
			if(rindex < n_temp) {
				std::complex<double> xi2(xi3d_temp[coord][0], xi3d_temp[coord][1]);

				// spherical harmonics //
				xi2 *= Ylm1[coord] * Ylm2[coord];

				xi_temp[rindex] += xi2;
				n_mode_temp[rindex]++;

			}
		}}}

		/*****************************************/

		double dr_bin[this->params.num_rbin-1];
		for(int j = 0; j < this->params.num_rbin-1; j++) {
			dr_bin[j] = rbin[j+1] - rbin[j];
		}
		for(int j = 0; j < this->params.num_rbin; j++) {

			double lower = 0.0;
			if(j ==0) {
				lower = ( rbin[j] > dr_bin[j] * 0.5) ? (rbin[j] - dr_bin[j] * 0.5) : 0.0;
			} else {
				lower = rbin[j] - dr_bin[j-1] * 0.5;
			}
			double upper = 0.0;
			if(j == this->params.num_rbin-1) {
				upper = (rbin[j] + dr_bin[j-1] * 0.5);
			} else {
				upper = (rbin[j] + dr_bin[j] * 0.5);
			}

			for(int i = 0; i < n_temp; i++) {
				double rmag_temp = double(i) * dr_temp;
				if( (rmag_temp > lower) && (rmag_temp <= upper) ) {
					this->xi[j] += xi_temp[i];
					this->n_mode_xi[j] += n_mode_temp[i];
				}
			}
		}

		double dV = this->params.volume / double(this->params.nmesh_tot);
		for(int i = 0; i < this->params.num_rbin; i++) {
			if(this->n_mode_xi[i] != 0) {
				this->xi[i] *= ( pow(-1.0, this->params.ell1+this->params.ell2)
						/ double(this->n_mode_xi[i]) / double(this->n_mode_xi[i]) / dV );
			} else {
				this->xi[i] = 0.0;
			}
		}

		delete[] xi3d_temp;
		delete[] xi_temp;
		delete[] n_mode_temp;

		return 0;

	}


	double ShotnoiseFunction(double * kvec) {
		if(0) {
		} else if (this->params.assignment == "NGP") {
			return this->ShotnoiseFunctionNGP(kvec);
		} else if (this->params.assignment == "CIC") {
			return this->ShotnoiseFunctionCIC(kvec);
		} else if (this->params.assignment == "TSC") {
			return this->ShotnoiseFunctionTSC(kvec);
		} else {
			return 0.0;
		}
		return 0.0;
	}

	double ShotnoiseFunctionNGP(const double * kvec) {
		return 1.0;
	}

	double ShotnoiseFunctionCIC(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = int( kvec[0] / dk[0] + 1.0e-5);
		int j = int( kvec[1] / dk[1] + 1.0e-5);
		int k = int( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1): 0.0;
		double f2 = (j != 0) ? sin(x2): 0.0;
		double f3 = (k != 0) ? sin(x3): 0.0;
		double ret = ( 1.0 - (2.0/3.0) * f1 * f1 ) * ( 1.0 - (2.0/3.0) * f2 * f2 ) * ( 1.0 - (2.0/3.0) * f3 * f3 );
		return ret;
	}

	double ShotnoiseFunctionTSC(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = int( kvec[0] / dk[0] + 1.0e-5);
		int j = int( kvec[1] / dk[1] + 1.0e-5);
		int k = int( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1): 0.0;
		double f2 = (j != 0) ? sin(x2): 0.0;
		double f3 = (k != 0) ? sin(x3): 0.0;
		double ret = ( 1.0 - f1 * f1 + 2.0 * pow(f1,4) / 15.0 )
			     * ( 1.0 - f2 * f2 + 2.0 * pow(f2,4) / 15.0 ) * ( 1.0 - f3 * f3 + 2.0 * pow(f3,4) / 15.0 );
		return ret;
	}

	double WindowFunction(double * kvec) {
		if(0) {
		} else if (this->params.assignment == "NGP") {
			return this->WindowFunctionNGP(kvec);
		} else if (this->params.assignment == "CIC") {
			return this->WindowFunctionCIC(kvec);
		} else if (this->params.assignment == "TSC") {
			return this->WindowFunctionTSC(kvec);
		} else {
			return 1.0;
		}

		return 1.0;
	}

	double WindowFunctionNGP(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = (int) ( kvec[0] / dk[0] + 1.0e-5);
		int j = (int) ( kvec[1] / dk[1] + 1.0e-5);
		int k = (int) ( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1) / x1 : 1.0;
		double f2 = (j != 0) ? sin(x2) / x2 : 1.0;
		double f3 = (k != 0) ? sin(x3) / x3 : 1.0;

		double ret = f1 * f2 * f3;
		return pow(ret,1);  // W = ret**3
	}


	double WindowFunctionCIC(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = (int) ( kvec[0] / dk[0] + 1.0e-5);
		int j = (int) ( kvec[1] / dk[1] + 1.0e-5);
		int k = (int) ( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1) / x1 : 1.0;
		double f2 = (j != 0) ? sin(x2) / x2 : 1.0;
		double f3 = (k != 0) ? sin(x3) / x3 : 1.0;

		double ret = f1 * f2 * f3;
		return pow(ret,2);  // W = ret**3
	}



	double WindowFunctionTSC(const double * kvec) {
		double dk[3];
		dk[0] = 2.0 * M_PI / this->params.boxsize[0];
		dk[1] = 2.0 * M_PI / this->params.boxsize[1];
		dk[2] = 2.0 * M_PI / this->params.boxsize[2];
		int i = (int) ( kvec[0] / dk[0] + 1.0e-5);
		int j = (int) ( kvec[1] / dk[1] + 1.0e-5);
		int k = (int) ( kvec[2] / dk[2] + 1.0e-5);
		double x1 = M_PI * double(i) / double(this->params.nmesh[0]);
		double x2 = M_PI * double(j) / double(this->params.nmesh[1]);
		double x3 = M_PI * double(k) / double(this->params.nmesh[2]);
		double f1 = (i != 0) ? sin(x1) / x1 : 1.0;
		double f2 = (j != 0) ? sin(x2) / x2 : 1.0;
		double f3 = (k != 0) ? sin(x3) / x3 : 1.0;

		double ret = f1 * f2 * f3;
		return pow(ret,3);  // W = ret**3
	}



	std::complex<double> calcShotNoiseForPowerspectrum_in_boxForReconstruction(TemplateParticle & P_D, TemplateParticle & P_R, double alpha) {

		std::complex<double> sum_D = double(P_D.n_tot);
		std::complex<double> sum_R = double(P_R.n_tot);

		return sum_D + pow(alpha,2) * sum_R;

	}
};

#endif
