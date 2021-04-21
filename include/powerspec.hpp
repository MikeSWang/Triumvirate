#ifndef __powerspec__
#define __powerspec__

#ifndef TRIUM_FIELD_H_INCLUDED_
#include "field.hpp"
#endif

int calc_power_spec(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        LineOfSight* los_D, LineOfSight* los_R,
        ParameterSet & param, double alpha, double * kbin, double vol_survey) {

	if(thisTask == 0) { printf("start to compute power spectrum...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensityFluctuation(P_D, P_R, los_D, los_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define pk_save */
	std::complex<double> * pk_save = new std::complex<double>[param.num_kbin];
	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc power spectrum */
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensityFluctuation(P_D, P_R, los_D, los_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();

		/* calc shotnoise term */
		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrum(P_D, P_R, los_D, los_R, alpha, param.ELL, _M_);

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {

		/*****/
		/* (-1)**m1 delta_{m1, -M} */
		double w = wigner_3j(param.ell1, 0, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, 0, param.ELL, _m1_, 0, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		stat.calc_power_spec(dn_LM, dn_00, kbin, shotnoise, param.ell1, _m1_);

		for(int i =0; i < param.num_kbin; i++) {
			pk_save[i] += w * stat.pk[i];
	       	}

	}
		durationInSec = double(clock() - timeStart);
		if(thisTask == 0) {
			printf("M = %d | %.3f sec\n", _M_, durationInSec / CLOCKS_PER_SEC);
		}
	}

	double norm = ParticleBOSSClass::calc_norm_for_power_spectrum(P_D, vol_survey);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(fp);

	delete[] pk_save;
	return 0;
}

int calc_2pt_func(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        LineOfSight* los_D, LineOfSight* los_R,
        ParameterSet & param, double alpha, double * rbin, double vol_survey) {

	if(thisTask == 0) { printf("start to compute two-point correlation function...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedDensityFluctuation(P_D, P_R, los_D, los_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.num_rbin];
	for(int i =0; i < param.num_rbin; i++) {
		xi_save[i] = 0.0;
	}

	/* calc power spectrum */
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedDensityFluctuation(P_D, P_R, los_D, los_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();

		/* calc shotnoise term */
		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrum(P_D, P_R, los_D, los_R, alpha, param.ELL, _M_);

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {

		/*****/
		/* (-1)**m1 delta_{m1, -M} */
		double w = wigner_3j(param.ell1, 0, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, 0, param.ELL, _m1_, 0, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		stat.calcCorrelationFunction(dn_LM, dn_00, rbin, shotnoise, param.ell1, _m1_);

		for(int i =0; i < param.num_rbin; i++) {
			xi_save[i] += w * stat.xi[i];
	       	}

	}
		durationInSec = double(clock() - timeStart);
		if(thisTask == 0) {
			printf("M = %d | %.3f sec\n", _M_, durationInSec / CLOCKS_PER_SEC);
		}
	}

	double norm = ParticleBOSSClass::calc_norm_for_power_spectrum(P_D, vol_survey);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete[] xi_save;
	return 0;
}


int calc_power_specWindowFunction(
        ParticleBOSSClass & P_R,
        LineOfSight* los_R,
        ParameterSet & param, double alpha, double * kbin, double vol_survey) {

	if(thisTask == 0) { printf("start to compute two-point window function...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedMeanDensity(P_R, los_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define pk_save */
	std::complex<double> * pk_save = new std::complex<double>[param.num_kbin];

	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc power spectrum */
	std::cout << "BITE = " << bytes << std::endl;

	/* calc shotnoise term */
	TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
	std::complex<double> shotnoise = stat.calcShotNoiseForTwoPointWindowFunction(P_R, los_R, alpha, param.ELL, 0);

	std::cout << "BITE = " << bytes << std::endl;

	stat.calc_power_spec(dn_00, dn_00, kbin, shotnoise, param.ell1, 0);

	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] += stat.pk[i];
	}

	std::cout << "BITE = " << bytes << std::endl;

	double norm = ParticleBOSSClass::calc_norm_for_power_spectrum(P_R, vol_survey);
	norm /= (alpha * alpha);

	/***********************************************/
	/***********************************************/
	/********** IMPORTANT *************************/
	norm /= param.volume;
	/***********************************************/
	/***********************************************/
	/***********************************************/

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d_window", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(fp);

	std::cout << "XIXIXIXIXIXIXXI" << std::endl;
	std::cout << norm * pk_save[0].real() << std::endl;
	std::cout << "XIXIXIXIXIXIXXI" << std::endl;

	delete[] pk_save;
	return 0;
}


int calc_2pt_func_window(
        ParticleBOSSClass & P_R,
        LineOfSight* los_R,
        ParameterSet & param, double alpha, double * rbin, double vol_survey) {
	if(thisTask == 0) { printf("start to compute two-point window function...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn_00(param);
	dn_00.calcYlmWeightedMeanDensity(P_R, los_R, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.num_rbin];


	for(int i =0; i < param.num_rbin; i++) {
		xi_save[i] = 0.0;
	}

	/* calc power spectrum */
	for(int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityFieldClass<ParticleBOSSClass> dn_LM(param);
		dn_LM.calcYlmWeightedMeanDensity(P_R, los_R, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calcFourierTransform();

		/* calc shotnoise term */
		TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calcShotNoiseForTwoPointWindowFunction(P_R, los_R, alpha, param.ELL, _M_);

	for(int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {

		/*****/
		/* (-1)**m1 delta_{m1, -M} */
		double w = wigner_3j(param.ell1, 0, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, 0, param.ELL, _m1_, 0, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1);
		if(fabs(w) < 1.0e-10) {
			continue;
		}

		stat.calcCorrelationFunction(dn_LM, dn_00, rbin, shotnoise, param.ell1, _m1_);

		for(int i =0; i < param.num_rbin; i++) {
			xi_save[i] += w * stat.xi[i];
	       	}

	}
		durationInSec = double(clock() - timeStart);
		if(thisTask == 0) {
			printf("M = %d | %.3f sec\n", _M_, durationInSec / CLOCKS_PER_SEC);
		}
	}

	double norm = ParticleBOSSClass::calc_norm_for_power_spectrum(P_R, vol_survey);
	norm /= (alpha * alpha);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d_window", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete[] xi_save;
	return 0;
}

int calc_power_spec_in_box(ParticleBOSSClass & P_D, ParameterSet & param, double * kbin) {
	if(thisTask == 0) { printf("start to compute power spectrum...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn(param);
	dn.calcNormalDensityFluctuation_in_box(P_D, param);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define pk_save */
//	std::complex<double> pk_save[param.num_kbin];
	std::complex<double> * pk_save = new std::complex<double>[param.num_kbin];
	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc shotnoise term */
	TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
	std::complex<double> shotnoise = double(P_D.n_tot);

	stat.calc_power_spec(dn, dn, kbin, shotnoise, param.ELL, 0);

	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] += double(2*param.ELL+1) * stat.pk[i];
       	}

	durationInSec = double(clock() - timeStart);
	if(thisTask == 0) {
		printf("%.3f sec\n", durationInSec / CLOCKS_PER_SEC);
	}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(fp);

	delete[] pk_save;

	return 0;
}

int calc_2pt_func_in_box(ParticleBOSSClass & P_D, ParameterSet & param, double * rbin) {
	if(thisTask == 0) { printf("start to compute two-point correlation function...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn(param);
	dn.calcNormalDensityFluctuation_in_box(P_D, param);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[param.num_rbin];
	for(int i =0; i < param.num_rbin; i++) {
		xi_save[i] = 0.0;
	}

	/* calc power spectrum */
	TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
	std::complex<double> shotnoise = double(P_D.n_tot);

	stat.calcCorrelationFunction(dn, dn, rbin, shotnoise, param.ELL, 0);

	for(int i =0; i < param.num_rbin; i++) {
		xi_save[i] += double(2*param.ELL+1) * stat.xi[i];
	}

	durationInSec = double(clock() - timeStart);
	if(thisTask == 0) {
		printf("%.3f sec\n", durationInSec / CLOCKS_PER_SEC);
	}

	double norm = param.volume/double(P_D.n_tot)/double(P_D.n_tot);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete[] xi_save;

	return 0;
}


int calc_power_spec_in_boxForReconstruction(ParticleBOSSClass & P_D, ParticleBOSSClass & P_R, ParameterSet & param, double alpha, double * kbin) {
	if(thisTask == 0) { printf("start to compute power spectrum...\n");}

	if( !( (param.ELL == param.ell1) && (param.ell2 == 0) ) ) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityFieldClass<ParticleBOSSClass> dn(param);
	dn.calcNormalDensityFluctuation_in_boxForReconstruction(P_D, P_R, alpha);
	/* Fourier transform*/
	dn.calcFourierTransform();

	/* define pk_save */
//	std::complex<double> pk_save[param.num_kbin];
	std::complex<double> * pk_save = new std::complex<double>[param.num_kbin];
	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc shotnoise term */
	TwoPointStatisticsClass<ParticleBOSSClass> stat(param);
	std::complex<double> shotnoise = stat.calcShotNoiseForPowerspectrum_in_boxForReconstruction(P_D, P_R, alpha);

	stat.calc_power_spec(dn, dn, kbin, shotnoise, param.ELL, 0);

	for(int i =0; i < param.num_kbin; i++) {
		pk_save[i] += double(2*param.ELL+1) * stat.pk[i];
       	}

	durationInSec = double(clock() - timeStart);
	if(thisTask == 0) {
		printf("%.3f sec\n", durationInSec / CLOCKS_PER_SEC);
	}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d", param.output_dir.c_str(), param.ELL);
	fp = fopen(buf, "w");
	for(int i =0; i < param.num_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(fp);

	delete[] pk_save;

	return 0;
}


#endif

