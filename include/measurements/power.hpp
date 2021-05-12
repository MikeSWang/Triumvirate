#ifndef TRIUM_POWER_H_INCLUDED_
#define TRIUM_POWER_H_INCLUDED_

#ifndef TRIUM_FIELD_H_INCLUDED_
#include "field.hpp"
#endif

/**
 * Calculate power spectrum from catalogues.
 *
 * @param particles_data Data-source particle container.
 * @param particles_rand Random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params Input parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param vol_survey Survey volume.
 * @returns Exit status.
 */
int calc_power_spec(
	ParticlesCatalogue& particles_data, ParticlesCatalogue & particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
  double* kbin,
  double vol_survey
) {
	if (thisTask == 0) {
		printf("[Info] :: Measuring power spectrum.\n");
	}

	if (!(params.ell1 == params.ELL && params.ell2 == 0)) {
		printf("[Error] :: Disallowed multipole degree combination for power spectrum measurements. ");
		printf("Please set `ell1 = ELL` and `ell2 = 0`.\n");
		exit(1);
	}

	/// Compute the density fluctuation.
	DensityField<ParticlesCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_overdensity(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);

	/// Fourier transform.
	dn_00.calc_fourier_transform();

	/// Compute power spectrum.
	std::complex<double>* pk_save = new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		pk_save[i] = 0.;
	}

	for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
		/// Compute spherical harmonic transform of the density fluctuation.
		DensityField<ParticlesCatalogue> dn_lm(params);
		dn_lm.calc_ylm_weighted_overdensity(
			particles_data, particles_rand,
			los_data, los_rand,
			alpha,
			params.ELL, M_
		);

		/// Fourier transform.
		dn_lm.calc_fourier_transform();

		/// Compute shot noise.
		TwoPointStatistics<ParticlesCatalogue> stats(params);
		std::complex<double> shotnoise = stats.calc_shotnoise_for_power_spec(
			particles_data, particles_rand,
			los_data, los_rand,
			alpha,
			params.ELL, M_
		);

		for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
			/// ???: (-1)^m_1 \delta_{m_1, -M}.
			double coupling = (2*params.ELL + 1) * (2*params.ell1 + 1)
				* wigner_3j(params.ell1, 0, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, 0, params.ELL, m1_, 0, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			stats.calc_power_spec(dn_lm, dn_00, kbin, shotnoise, params.ell1, m1_);

			for (int i = 0; i < params.num_kbin; i++) {
				pk_save[i] += coupling * stats.pk[i];
			}
		}

		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) {
			printf(
				"[Info] Computed order M = %d (... %.3f seconds elapsed in total).\n",
				M_, durationInSec / CLOCKS_PER_SEC
			);
		}
	}

	double norm = ParticlesCatalogue::calc_norm_for_power_spec(particles_data, vol_survey);

	char buf[1024];
	sprintf(buf, "%s/pk%d", params.output_dir.c_str(), params.ELL);

	FILE* ptr_file;
	ptr_file = fopen(buf, "w");
	for (int i = 0; i < params.num_kbin; i++) {
		fprintf(ptr_file, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(ptr_file);

	delete[] pk_save;

	return 0;
}

/**
 * Calculate two-point correlation function.
 *
 * @param particles_data Data-source particle container.
 * @param particles_rand Random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params Input parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param vol_survey Survey volume.
 * @returns Exit status.
 */
int calc_2pt_func(
	ParticlesCatalogue& particles_data, ParticlesCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double vol_survey
) {

	if (thisTask == 0) { printf("start to compute two-point correlation function...\n");}

	if (!((params.ELL == params.ell1) && (params.ell2 == 0))) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticlesCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_overdensity(particles_data, particles_rand, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[params.num_rbin];
	for (int i =0; i < params.num_rbin; i++) {
		xi_save[i] = 0.0;
	}

	/* calc power spectrum */
	for (int m_ = - params.ELL; m_ <= params.ELL; m_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityField<ParticlesCatalogue> dn_lm(params);
		dn_lm.calc_ylm_weighted_overdensity(particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, m_);
		/* Fourier transform*/
		dn_lm.calc_fourier_transform();

		/* calc shotnoise term */
		TwoPointStatistics<ParticlesCatalogue> stat(params);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_power_spec(particles_data, particles_rand, los_data, los_rand, alpha, params.ELL, m_);

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {

		/*****/
		/* (-1)**m1 delta_{m1, -M} */
		double w = wigner_3j(params.ell1, 0, params.ELL, 0, 0, 0) * wigner_3j(params.ell1, 0, params.ELL, m1_, 0, m_);
		w *= double(2*params.ELL+1) * double(2*params.ell1+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		stat.calc_corr_func(dn_lm, dn_00, rbin, shotnoise, params.ell1, m1_);

		for (int i =0; i < params.num_rbin; i++) {
			xi_save[i] += w * stat.xi[i];
	       	}

	}
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) {
			printf("M = %d | %.3f sec\n", m_, durationInSec / CLOCKS_PER_SEC);
		}
	}

	double norm = ParticlesCatalogue::calc_norm_for_power_spec(particles_data, vol_survey);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d", params.output_dir.c_str(), params.ELL);
	fp = fopen(buf, "w");
	for (int i =0; i < params.num_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete[] xi_save;
	return 0;
}


int calc_power_specWindowFunction(
        ParticlesCatalogue & particles_rand,
        LineOfSight* los_rand,
        ParameterSet & params, double alpha, double * kbin, double vol_survey) {

	if (thisTask == 0) { printf("start to compute two-point window function...\n");}

	if (!((params.ELL == params.ell1) && (params.ell2 == 0))) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticlesCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* define pk_save */
	std::complex<double> * pk_save = new std::complex<double>[params.num_kbin];

	for (int i =0; i < params.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc power spectrum */
	std::cout << "BITE = " << bytes << std::endl;

	/* calc shotnoise term */
	TwoPointStatistics<ParticlesCatalogue> stat(params);
	std::complex<double> shotnoise = stat.calc_shotnoise_for_2pt_func_window(particles_rand, los_rand, alpha, params.ELL, 0);

	std::cout << "BITE = " << bytes << std::endl;

	stat.calc_power_spec(dn_00, dn_00, kbin, shotnoise, params.ell1, 0);

	for (int i =0; i < params.num_kbin; i++) {
		pk_save[i] += stat.pk[i];
	}

	std::cout << "BITE = " << bytes << std::endl;

	double norm = ParticlesCatalogue::calc_norm_for_power_spec(particles_rand, vol_survey);
	norm /= (alpha * alpha);

	/***********************************************/
	/***********************************************/
	/********** IMPORTANT *************************/
	norm /= params.volume;
	/***********************************************/
	/***********************************************/
	/***********************************************/

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d_window", params.output_dir.c_str(), params.ELL);
	fp = fopen(buf, "w");
	for (int i =0; i < params.num_kbin; i++) {
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
        ParticlesCatalogue & particles_rand,
        LineOfSight* los_rand,
        ParameterSet & params, double alpha, double * rbin, double vol_survey) {
	if (thisTask == 0) { printf("start to compute two-point window function...\n");}

	if (!((params.ELL == params.ell1) && (params.ell2 == 0))) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticlesCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[params.num_rbin];


	for (int i =0; i < params.num_rbin; i++) {
		xi_save[i] = 0.0;
	}

	/* calc power spectrum */
	for (int m_ = - params.ELL; m_ <= params.ELL; m_++) {

		/****************************************/
		/* calc the ylm-weighted density fluctuation */
		DensityField<ParticlesCatalogue> dn_lm(params);
		dn_lm.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, params.ELL, m_);
		/* Fourier transform*/
		dn_lm.calc_fourier_transform();

		/* calc shotnoise term */
		TwoPointStatistics<ParticlesCatalogue> stat(params);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_2pt_func_window(particles_rand, los_rand, alpha, params.ELL, m_);

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {

		/*****/
		/* (-1)**m1 delta_{m1, -M} */
		double w = wigner_3j(params.ell1, 0, params.ELL, 0, 0, 0) * wigner_3j(params.ell1, 0, params.ELL, m1_, 0, m_);
		w *= double(2*params.ELL+1) * double(2*params.ell1+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		stat.calc_corr_func(dn_lm, dn_00, rbin, shotnoise, params.ell1, m1_);

		for (int i =0; i < params.num_rbin; i++) {
			xi_save[i] += w * stat.xi[i];
	       	}

	}
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) {
			printf("M = %d | %.3f sec\n", m_, durationInSec / CLOCKS_PER_SEC);
		}
	}

	double norm = ParticlesCatalogue::calc_norm_for_power_spec(particles_rand, vol_survey);
	norm /= (alpha * alpha);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d_window", params.output_dir.c_str(), params.ELL);
	fp = fopen(buf, "w");
	for (int i =0; i < params.num_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete[] xi_save;
	return 0;
}

int calc_power_spec_in_box(ParticlesCatalogue & particles_data, ParameterSet & params, double * kbin) {
	if (thisTask == 0) { printf("start to compute power spectrum...\n");}

	if (!((params.ELL == params.ell1) && (params.ell2 == 0))) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticlesCatalogue> dn(params);
	dn.calc_density_field_in_box(particles_data, params);
	/* Fourier transform*/
	dn.calc_fourier_transform();

	/* define pk_save */
//	std::complex<double> pk_save[params.num_kbin];
	std::complex<double> * pk_save = new std::complex<double>[params.num_kbin];
	for (int i =0; i < params.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc shotnoise term */
	TwoPointStatistics<ParticlesCatalogue> stat(params);
	std::complex<double> shotnoise = double(particles_data.n_tot);

	stat.calc_power_spec(dn, dn, kbin, shotnoise, params.ELL, 0);

	for (int i =0; i < params.num_kbin; i++) {
		pk_save[i] += double(2*params.ELL+1) * stat.pk[i];
       	}

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf("%.3f sec\n", durationInSec / CLOCKS_PER_SEC);
	}

	double norm = params.volume / double(particles_data.n_tot) / double(particles_data.n_tot);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d", params.output_dir.c_str(), params.ELL);
	fp = fopen(buf, "w");
	for (int i =0; i < params.num_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(fp);

	delete[] pk_save;

	return 0;
}

int calc_2pt_func_in_box(ParticlesCatalogue & particles_data, ParameterSet & params, double * rbin) {
	if (thisTask == 0) { printf("start to compute two-point correlation function...\n");}

	if (!((params.ELL == params.ell1) && (params.ell2 == 0))) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticlesCatalogue> dn(params);
	dn.calc_density_field_in_box(particles_data, params);
	/* Fourier transform*/
	dn.calc_fourier_transform();

	/* define xi_save */
	std::complex<double> * xi_save = new std::complex<double>[params.num_rbin];
	for (int i =0; i < params.num_rbin; i++) {
		xi_save[i] = 0.0;
	}

	/* calc power spectrum */
	TwoPointStatistics<ParticlesCatalogue> stat(params);
	std::complex<double> shotnoise = double(particles_data.n_tot);

	stat.calc_corr_func(dn, dn, rbin, shotnoise, params.ELL, 0);

	for (int i =0; i < params.num_rbin; i++) {
		xi_save[i] += double(2*params.ELL+1) * stat.xi[i];
	}

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf("%.3f sec\n", durationInSec / CLOCKS_PER_SEC);
	}

	double norm = params.volume/double(particles_data.n_tot)/double(particles_data.n_tot);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/xi%d", params.output_dir.c_str(), params.ELL);
	fp = fopen(buf, "w");
	for (int i =0; i < params.num_rbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", rbin[i], norm * xi_save[i].real());
	}
	fclose(fp);

	delete[] xi_save;

	return 0;
}


int calc_power_spec_in_boxForReconstruction(ParticlesCatalogue & particles_data, ParticlesCatalogue & particles_rand, ParameterSet & params, double alpha, double * kbin) {
	if (thisTask == 0) { printf("start to compute power spectrum...\n");}

	if (!((params.ELL == params.ell1) && (params.ell2 == 0))) {
		printf("This multipole combination is not allowed.\n");
		printf("It should be ELL == ell1 and ell2 == 0\n");
		exit(1);
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticlesCatalogue> dn(params);
	dn.calc_density_field_in_box_for_recon(particles_data, particles_rand, alpha);
	/* Fourier transform*/
	dn.calc_fourier_transform();

	/* define pk_save */
//	std::complex<double> pk_save[params.num_kbin];
	std::complex<double> * pk_save = new std::complex<double>[params.num_kbin];
	for (int i =0; i < params.num_kbin; i++) {
		pk_save[i] = 0.0;
	}

	/* calc shotnoise term */
	TwoPointStatistics<ParticlesCatalogue> stat(params);
	std::complex<double> shotnoise = stat.calc_shotnoise_for_power_spec_in_box_for_recon(particles_data, particles_rand, alpha);

	stat.calc_power_spec(dn, dn, kbin, shotnoise, params.ELL, 0);

	for (int i =0; i < params.num_kbin; i++) {
		pk_save[i] += double(2*params.ELL+1) * stat.pk[i];
       	}

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf("%.3f sec\n", durationInSec / CLOCKS_PER_SEC);
	}

	double norm = params.volume / double(particles_data.n_tot) / double(particles_data.n_tot);

	FILE * fp;
	char buf[1024];
	sprintf(buf, "%s/pk%d", params.output_dir.c_str(), params.ELL);
	fp = fopen(buf, "w");
	for (int i =0; i < params.num_kbin; i++) {
		fprintf(fp, "%.5f \t %.7e\n", kbin[i], norm * pk_save[i].real());
	}
	fclose(fp);

	delete[] pk_save;

	return 0;
}


#endif

