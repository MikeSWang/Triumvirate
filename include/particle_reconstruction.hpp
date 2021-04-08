#ifndef __particle_reconstruction__
#define __particle_reconstruction__

#ifndef __parameter__
#include "parameter.hpp"
#endif

#ifndef __field__
#include "field.hpp"
#endif

int calcReconstructionParticles(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        losStruct * los_D, losStruct * los_R,
        ParameterClass & param, double alpha,
        double Vsurvey, double b1_fid, double RG, int nmesh[3]) {

	double nmesh_save[3] = {param.nmesh[0], param.nmesh[1], param.nmesh[2]};
	double nmesh_tot_save = param.nmesh_tot;

	param.nmesh[0] = nmesh[0];
	param.nmesh[1] = nmesh[1];
	param.nmesh[2] = nmesh[2];
	param.nmesh_tot = nmesh[0] * nmesh[1] * nmesh[2];

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> D(param);
	D.calcYlmWeightedDensityFluctuation(P_D, P_R, los_D, los_R, alpha, 0, 0);
	/* Fourier transform*/
	D.calcFourierTransform();
	D.calcAssignmentFunctionCorrection();

	double kvec[3];
	double dk[3];
	dk[0] = 2.0 * M_PI / param.boxsize[0];
	dk[1] = 2.0 * M_PI / param.boxsize[1];
	dk[2] = 2.0 * M_PI / param.boxsize[2];

	fftw_complex * Psi_0 = NULL;
	fftw_complex * Psi_1 = NULL;
	fftw_complex * Psi_2 = NULL;
	Psi_0 = fftw_alloc_complex(param.nmesh_tot);
	Psi_1 = fftw_alloc_complex(param.nmesh_tot);
	Psi_2 = fftw_alloc_complex(param.nmesh_tot);
	std::complex<double> _I_(0.0, 1.0);

	/****************************************/
	/* set Psi(kvec) */
	for(int i = 0; i < param.nmesh[0]; i++) {
	for(int j = 0; j < param.nmesh[1]; j++) {
	for(int k = 0; k < param.nmesh[2]; k++) {
		long long coord = ( i * param.nmesh[1] + j ) * param.nmesh[2] + k;
		kvec[0] = (i < param.nmesh[0]/2) ? (double) i * dk[0] : (double) (i - param.nmesh[0]) * dk[0];
		kvec[1] = (j < param.nmesh[1]/2) ? (double) j * dk[1] : (double) (j - param.nmesh[1]) * dk[1];
		kvec[2] = (k < param.nmesh[2]/2) ? (double) k * dk[2] : (double) (k - param.nmesh[2]) * dk[2];
		double kmag = sqrt( kvec[0] * kvec[0] + kvec[1] * kvec[1] + kvec[2] * kvec[2] );
		/*** IMPORTANT ***/
		if(kmag < 1.0e-10) {
		    continue;
		}
		/*****************/
		double kmag2 = kmag * kmag;
		double W = exp( - pow(kmag * RG, 2) / 2.0);
		std::complex<double> fD(D[coord][0], D[coord][1]);

		std::complex<double> Psi_temp[3];
		for(int axis = 0; axis < 3; axis++) {
		    Psi_temp[axis] = - _I_ * (W / b1_fid) * ( kvec[axis] / kmag2 ) * fD;
		}

		Psi_0[coord][0] = Psi_temp[0].real();
		Psi_0[coord][1] = Psi_temp[0].imag();

		Psi_1[coord][0] = Psi_temp[1].real();
		Psi_1[coord][1] = Psi_temp[1].imag();

		Psi_2[coord][0] = Psi_temp[2].real();
		Psi_2[coord][1] = Psi_temp[2].imag();

	}}}
	/****************************************/

	/****************************************/
	/* compute Psi(x) */
	for(int i = 0;  i < param.nmesh_tot; i++) {
		Psi_0[i][0] /= param.volume ;
		Psi_0[i][1] /= param.volume;
		Psi_1[i][0] /= param.volume ;
		Psi_1[i][1] /= param.volume;
		Psi_2[i][0] /= param.volume ;
		Psi_2[i][1] /= param.volume;
	}

	/* Inverse FFT */
	fftw_plan plan_0 = fftw_plan_dft_3d(param.nmesh[0], param.nmesh[1], param.nmesh[2], Psi_0, Psi_0, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_1 = fftw_plan_dft_3d(param.nmesh[0], param.nmesh[1], param.nmesh[2], Psi_1, Psi_1, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_plan plan_2 = fftw_plan_dft_3d(param.nmesh[0], param.nmesh[1], param.nmesh[2], Psi_2, Psi_2, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(plan_0);
	fftw_execute(plan_1);
	fftw_execute(plan_2);
	fftw_destroy_plan(plan_0);
	fftw_destroy_plan(plan_1);
	fftw_destroy_plan(plan_2);
	/****************************************/

	double num_D_weight = 0.0;
	for(int p = 0; p < P_D.n_tot; p++) {
		num_D_weight += P_D[p].w;
	}

	/** normalization **/
	for(int i = 0; i < param.nmesh_tot; i++) {
	    Psi_0[i][0] *= (Vsurvey / num_D_weight);
	    Psi_0[i][1] *= (Vsurvey / num_D_weight);

	    Psi_1[i][0] *= (Vsurvey / num_D_weight);
	    Psi_1[i][1] *= (Vsurvey / num_D_weight);

	    Psi_2[i][0] *= (Vsurvey / num_D_weight);
	    Psi_2[i][1] *= (Vsurvey / num_D_weight);

	}

	/****************************************/
	/* add Psi to particle positions */
	for(int p = 0; p < P_D.n_tot; p++) {
		double w[3][3];
		int iw[3][3];
		for(int axes = 0; axes < 3; axes++) {
			double xp = P_D[p].pos[axes] * double(param.nmesh[axes]) / param.boxsize[axes];
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
			long long coord = ( iw[i][0] * param.nmesh[1] + iw[j][1] ) * param.nmesh[2] + iw[k][2];
			if( (coord >= 0) && (coord < param.nmesh_tot) ) {
//				this->field[coord][0] += weight[p][0] * factor * w[i][0] * w[j][1] * w[k][2];
//				this->field[coord][1] += weight[p][1] * factor * w[i][0] * w[j][1] * w[k][2];
				P_D[p].pos[0] += w[i][0] * w[j][1] * w[k][2] * Psi_0[coord][0];
				P_D[p].pos[1] += w[i][0] * w[j][1] * w[k][2] * Psi_1[coord][0];
				P_D[p].pos[2] += w[i][0] * w[j][1] * w[k][2] * Psi_2[coord][0];
			}
		}}}

	}


	/****************************************/
	/* add Psi to particle positions */
	for(int p = 0; p < P_R.n_tot; p++) {
		double w[3][3];
		int iw[3][3];
		for(int axes = 0; axes < 3; axes++) {
			double xp = P_R[p].pos[axes] * double(param.nmesh[axes]) / param.boxsize[axes];
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
			long long coord = ( iw[i][0] * param.nmesh[1] + iw[j][1] ) * param.nmesh[2] + iw[k][2];
			if( (coord >= 0) && (coord < param.nmesh_tot) ) {
//				this->field[coord][0] += weight[p][0] * factor * w[i][0] * w[j][1] * w[k][2];
//				this->field[coord][1] += weight[p][1] * factor * w[i][0] * w[j][1] * w[k][2];
				P_R[p].pos[0] += w[i][0] * w[j][1] * w[k][2] * Psi_0[coord][0];
				P_R[p].pos[1] += w[i][0] * w[j][1] * w[k][2] * Psi_1[coord][0];
				P_R[p].pos[2] += w[i][0] * w[j][1] * w[k][2] * Psi_2[coord][0];
			}
		}}}
	}

	fftw_free(Psi_0); Psi_0 = NULL;
	fftw_free(Psi_1); Psi_1 = NULL;
	fftw_free(Psi_2); Psi_2 = NULL;
	D.finalizeDensityField();

	/* reset param. */
	param.nmesh[0] = nmesh_save[0];
	param.nmesh[1] = nmesh_save[1];
	param.nmesh[2] = nmesh_save[2];
	param.nmesh_tot = nmesh_tot_save;

	return 0;
}

#endif

