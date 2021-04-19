/* standard libraries */
#include <iostream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <sys/stat.h>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <complex>
#include <cmath>
#include <ctime>
#include <vector>

/* tooles */
#include <fftw3.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_coupling.h>

#define wigner_3j(j1,j2,j3,m1,m2,m3) (gsl_sf_coupling_3j(2*j1,2*j2,2*j3,2*m1,2*m2,2*m3))

#include "common.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "bessel.hpp"
#include "particles.hpp"
#include "field.hpp"
#include "powerspec.hpp"
#include "bispec.hpp"

#include "particle_reconstruction.hpp"

int main(int argc, char *argv[]) {

	/* time */
	timeStart = clock();

	/***********************************************/
	/* if there is no parameter file, exit program */
	/***********************************************/
	if(argc != 2){
		printf("Parameters are missing.\n");
		printf("Call with <ParameterFile>\n");
		exit(1);
	}

	/******************************/
	/* read "param.ini"           */
	/******************************/
	ParameterSet param;
	if( param.read_parameters(argv) ) {
		printf("ERROR\n");
		exit(1);
	}
	if(thisTask == 0) { printf("done\n");}

	/* set files */
	param.set_io_files();

	/* set kbin */
	double kbin[param.num_kbin];
	ToolCollection::set_kbin(param, kbin);
	/* set rbin */
	double rbin[param.num_rbin];
	ToolCollection::set_rbin(param, rbin);

	/***************************/
	/* read particles */
	/***************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if(thisTask == 0) { printf("reading particles | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}
	/********/

	ParticleBOSSClass P_D, P_R;
	if ( P_D.read_particles_BOSS(param.data_catalogue_file) ) {
		printf("ERROR2\n");
		exit(1);
	}
	if ( P_R.read_particles_BOSS(param.rand_catalogue_file) ) {
		printf("ERROR3\n");
		exit(1);
	}

    /*************/
    /* calc. LOS */
    /*************/
    LineOfSight* los_D = new LineOfSight[P_D.n_tot];
    LineOfSight* los_R = new LineOfSight[P_R.n_tot];
    for (int p = 0; p < P_D.n_tot; p++) {
        double pos_mag = sqrt( P_D[p].pos[0] * P_D[p].pos[0] + P_D[p].pos[1] * P_D[p].pos[1] + P_D[p].pos[2] * P_D[p].pos[2] );
        los_D[p].pos[0] = P_D[p].pos[0] / pos_mag;
        los_D[p].pos[1] = P_D[p].pos[1] / pos_mag;
        los_D[p].pos[2] = P_D[p].pos[2] / pos_mag;
    }
    for (int p = 0; p < P_R.n_tot; p++) {
        double pos_mag = sqrt( P_R[p].pos[0] * P_R[p].pos[0] + P_R[p].pos[1] * P_R[p].pos[1] + P_R[p].pos[2] * P_R[p].pos[2] );
        los_R[p].pos[0] = P_R[p].pos[0] / pos_mag;
        los_R[p].pos[1] = P_R[p].pos[1] / pos_mag;
        los_R[p].pos[2] = P_R[p].pos[2] / pos_mag;
    }

    /********************/
	double alpha = 0.0;
	if(param.catalogue_type== "mock" || param.catalogue_type== "survey") {

		/****************************************************/
		/* calc the ratio between the numbers of galaxies and random particles, namely "alpha" */
		alpha = ParticleBOSSClass::calcAlpha(P_D, P_R);

		/***************************************************************/
		/* place particles in a box for FFT */
		ParticleBOSSClass::resetParticleForFFT(P_D, P_R, param);

	}
	if(param.catalogue_type== "sim") {
		P_D.calcPeriodicBoundary(param);
	}

	/******** calc. survey volume *********/
//	int save_nmesh[3] = {param.nmesh[0], param.nmesh[1], param.nmesh[2]};
//	int save_nmesh_tot = param.nmesh_tot;
//	param.nmesh[0] = 512;
//	param.nmesh[1] = 512;
//	param.nmesh[2] = 512;
//	param.nmesh_tot = param.nmesh[0] * param.nmesh[1] * param.nmesh[2];
	DensityFieldClass<ParticleBOSSClass> Vol(param);
	double Vsurvey = Vol.calcSurveyVolume(P_R);
	Vol.finalizeDensityField();
//	param.nmesh[0] = save_nmesh[0];
//	param.nmesh[1] = save_nmesh[1];
//	param.nmesh[2] = save_nmesh[2];
//	param.nmesh_tot = save_nmesh_tot;
//
//	std::cout << param.nmesh[0] << std::endl;
//	std::cout << param.nmesh[1] << std::endl;
//	std::cout << param.nmesh[2] << std::endl;

	/* time */
	durationInSec = double(clock() - timeStart);
	if(thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}
	/********/

	if(thisTask == 0) {
		std::cout << "AAAAAAAAAAAAAA" << std::endl;
		std::cout << "N = " << P_D.n_tot << std::endl;
		std::cout << "V = " << param.volume << std::endl;
		std::cout << "L_x = " << param.boxsize[0] << std::endl;
		std::cout << "L_y = " << param.boxsize[1] << std::endl;
		std::cout << "L_z = " << param.boxsize[2] << std::endl;
		std::cout << "AAAAAAAAAAAAAA" << std::endl;
		std::cout << "           " << std::endl;
		std::cout << "BBBBBBBBBBBBBB" << std::endl;
		std::cout <<  P_D.pos_min[0] << "   " << P_R.pos_min[0] << std::endl;
		std::cout <<  P_D.pos_min[1] << "   " << P_R.pos_min[1] << std::endl;
		std::cout <<  P_D.pos_min[2] << "   " << P_R.pos_min[2] << std::endl;
		std::cout << "              " << std::endl;
		std::cout <<  P_D.pos_max[0] << "   " << P_R.pos_max[0] << std::endl;
		std::cout <<  P_D.pos_max[1] << "   " << P_R.pos_max[1] << std::endl;
		std::cout <<  P_D.pos_max[2] << "   " << P_R.pos_max[2] << std::endl;
		std::cout << "BBBBBBBBBBBBBB" << std::endl;
		std::cout << "              " << std::endl;
		std::cout << "CCCCCCCCCCCCCC" << std::endl;
		std::cout << "Survey_Volume = " << Vsurvey << std::endl;
		std::cout << "Survey_L = " << pow(Vsurvey,1.0/3.0) << std::endl;
	}

	if(param.reconstruction == "true") {
		/****************************/
		/*****  reconstruction  *****/
		/****************************/
		if(thisTask == 0) {
		    std::cout << "computing recostruction" << std::endl;
		}
		double b1_fid = param.b1_fid;
		double RG = param.RG;
		int nmesh[3] = {512, 512, 512};
		calcReconstructionParticles(P_D, P_R, los_D, los_R, param, alpha, Vsurvey, b1_fid, RG, nmesh);
		ParticleBOSSClass::resetParticleForFFT(P_D, P_R, param);

		/* time */
		durationInSec = double(clock() - timeStart);
		if(thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}
		/********/
	}

	/******************************/
	/* ith bin   */
	/******************************/
//   if(param.catalogue_type=="survey") {

//        if(thisTask >= param.num_rbin) {
//            param.ith_rbin = 0;
//        } else {
//            param.ith_rbin = thisTask;
//        }

//        if(thisTask >= param.num_kbin) {
//            param.ith_kbin = 0;
//        } else {
//            param.ith_kbin = thisTask;
//        }

//    }

	/******************************/
	/* calculate power spectrum   */
	/******************************/
	if (param.catalogue_type == "mock" || param.catalogue_type =="survey") {

//        calcPowerSpectrum(P_D, P_R, los_D, los_R, param, alpha, kbin, Vsurvey);
//        calcTwoPointFunction(P_D, P_R, los_D, los_R, param, alpha, rbin, Vsurvey);
//        calcTwoPointWindowFunction(P_R, los_R, param, alpha, rbin, Vsurvey);

        calcBiSpectrum(P_D, P_R, los_D, los_R, param, alpha, kbin, Vsurvey);
//        calcThreePointFunction(P_D, P_R, los_D, los_R, param, alpha, rbin, Vsurvey);
//        calcThreePointWindowFunction(P_R, los_R, param, alpha, rbin, Vsurvey);
//        calcThreePointWindowFunctionFor3PCF(P_R, los_R, param, alpha, rbin, Vsurvey);

	}

	if(param.catalogue_type== "sim") {
//		calcPowerSpectrumForBOX(P_D, param, kbin);
//		calcTwoPointFunctionForBOX(P_D, param, rbin);
//	    	calcBiSpectrumForBOX(P_D, param, kbin);
//		calcThreePointFunctionForBOX(P_D, param, rbin);
	}

    delete[] los_D; los_D = NULL;
    delete[] los_R; los_R = NULL;
	P_D.finalise_particles();
	P_R.finalise_particles();
	if(thisTask == 0 ) {
		std::cout << "bytes = " << bytes << std::endl;
	}



	/* time */
	durationInSec = double(clock() - timeStart);
	if(thisTask == 0) { printf("end program | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}
	/********/

	return 0;
}

