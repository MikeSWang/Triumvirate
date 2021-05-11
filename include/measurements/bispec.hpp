#ifndef __bispec__
#define __bispec__

int calc_bispec(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        LineOfSight* los_data, LineOfSight* los_rand,
        ParameterSet & param, double alpha, double * kbin, double vol_survey) {

	if (thisTask == 0) { printf("start to compute bispectrum...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.num_kbin];
	for (int i = 0; i < param.num_kbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00_sn(param);
	dn_00_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00_sn.calc_fourier_transform();

	/* calc shot noise terms */
	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> N_LM_sn(param);
		N_LM_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		N_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		if ((param.ell1 == 0) && (param.ell2 == 0)) {
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[i] += w * shotnoise;
	       		}
		}

		if ((param.ell2 == 0)) {
			stat.calc_power_spec(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell1, _m1_);
			if (param.form == "diag") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[i] += w * stat.pk[i];
			       	}
			} else if (param.form == "full") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[i] += w * stat.pk[param.ith_kbin];
			       	}
			}
		}

		if ((param.ell1 == 0)) {
			stat.calc_power_spec(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell2, _m2_);
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[i] += w * stat.pk[i];
		       	}
		}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}}}

	/****/
	dn_00_sn.finalise_density_field();
	/****************************************/

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);
		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		fftw_complex * xi = fftw_alloc_complex(param.nmesh_tot);
		bytes += double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);
		for (int i = 0; i < (param.nmesh_tot); i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calc_shotnoise_for_bispec_ijk(dn_LM_sn, N_00_sn, shotnoise, param.ELL, _M_, xi);

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];
			double kmag1;
			if (param.form == "diag") {
				kmag1 = kmag2;
			} else if (param.form == "full") {
				kmag1 = kbin[param.ith_kbin];
			}

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
			double rvec[3];
			double dr[3];
			dr[0] = param.boxsize[0] / double(param.nmesh[0]);
			dr[1] = param.boxsize[1] / double(param.nmesh[1]);
			dr[2] = param.boxsize[2] / double(param.nmesh[2]);
			for (int i = 0; i < param.nmesh[0]; i++) {
			for (int j = 0; j < param.nmesh[1]; j++) {
			for (int k = 0; k < param.nmesh[2]; k++) {
				long long coord_flat = (i * param.nmesh[1] + j) * param.nmesh[2] + k;
				rvec[0] = (i < param.nmesh[0]/2) ? (double) i * dr[0] : (double) (i - param.nmesh[0]) * dr[0];
				rvec[1] = (j < param.nmesh[1]/2) ? (double) j * dr[1] : (double) (j - param.nmesh[1]) * dr[1];
				rvec[2] = (k < param.nmesh[2]/2) ? (double) k * dr[2] : (double) (k - param.nmesh[2]) * dr[2];
				double rmag = sqrt(rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[i][0], xi[i][1]);
				double j1 = sj1.eval(kmag1 * rmag);
				double j2 = sj2.eval(kmag2 * rmag);

				sn_sum += (j1 * j2 * ff * Ylm1[i] * Ylm2[i]);

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.nmesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[ik] += (w * sn_sum);

			/* time */
			durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		bytes -= double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}

	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing bispectrum...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();


	/* define bk_save:: bispectrum */
	std::complex<double> * bk_save = new std::complex<double>[param.num_kbin];
	for (int i = 0; i < param.num_kbin; i++) {
		bk_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);
		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		if (param.form == "full") {
			kmag1 = kbin[param.ith_kbin];
			dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag1, dk, Ylm1);
		}

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				kmag1 = kmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag1, dk, Ylm1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_bispec(dn_00, kmag2, dk, Ylm2);

			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[ik] += (w * bk_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = ParticleBOSSClass::calc_norm_for_power_spec(P_D, vol_survey);

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/bk%d%d%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_kbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n", kbin[i], kbin[i],
					norm * (bk_save[i].real() - sn_save[i].real()), norm * (bk_save[i].imag() - sn_save[i].imag()),
					norm * sn_save[i].real(), norm * sn_save[i].imag());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/bk%d%d%d_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_kbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i],
					norm * (bk_save[i].real() - sn_save[i].real()), norm * (bk_save[i].imag() - sn_save[i].imag()),
					norm * sn_save[i].real(), norm * sn_save[i].imag());
		}
	}
	fclose(fp);

	delete[] sn_save;
	delete[] bk_save;

	return 0;
}


int calc_3pt_func(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        LineOfSight* los_data, LineOfSight* los_rand,
        ParameterSet & param, double alpha, double * rbin, double vol_survey) {


	if (thisTask == 0) { printf("start to compute three-point function...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.num_rbin];
	for (int i = 0; i < param.num_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		stat.calc_corr_func_for_3pt_func(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for (int i =0; i < param.num_rbin; i++) {
			if (param.form == "diag") {
				sn_save[i] += w * stat.xi[i];
			} else if (param.form == "full") {
				if (i == param.ith_rbin) {
					sn_save[i] += w * stat.xi[i];
				} else {
					sn_save[i] += 0.0;
				}
			}
	       	}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}


	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing three-point function...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);


	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.num_rbin];
	for (int i = 0; i < param.num_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {


		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double rmag1;
		if (param.form == "full") {
			rmag1 = rbin[param.ith_rbin];
			dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
		}

		for (int ir = 0; ir < param.num_rbin; ir++) {

			double rmag2 = rbin[ir];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				rmag1 = rmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag2, Ylm2, sj2);

			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("r2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}


	double norm = ParticleBOSSClass::calc_norm_for_bispec(P_D, vol_survey);

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/zeta%d%d%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_rbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", rbin[i], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/zeta%d%d%d_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_rbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	}
	fclose(fp);

	delete[] sn_save;
	delete[] zeta_save;

	return 0;
}



int calc_3pt_func_window(
        ParticleBOSSClass & P_R,
        LineOfSight* los_rand,
        ParameterSet & param, double alpha, double * rbin, double vol_survey) {

	if (thisTask == 0) { printf("start to compute three-point window...\n");}

//	int n_temp = param.num_rbin;
	int n_temp = 10;
	int NR = 3;
	/*************************************************************************/
	/* reset rbin */
	rbin[0] = 0.0;
	rbin[1] = 1.0;
	rbin[2] = 10.0;
	rbin[3] = 20.0;
	rbin[4] = 30.0;
	rbin[5] = 40.0;
	rbin[6] = 50.0;
	rbin[7] = 60.0;
	double rmin = 70.0;
	double dlnr = (log(param.rmax) - log(rmin)) / double((param.num_rbin-8) - 1);

	param.ith_rbin = thisTask;

	for (int i = 8; i < param.num_rbin; i++) {
		rbin[i] = rmin * exp(dlnr * double(i-8));
	}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[n_temp];
	for (int i = 0; i < n_temp; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_ylm_weighted_mean_density_for_3pt_window_shotnoise(P_R, los_rand, alpha, 0, 0);

	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_ylm_weighted_mean_density(P_R, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_2pt_func_window(P_R, los_rand, alpha, param.ELL, _M_);

		stat.calc_corr_func_for_3pt_func(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for (int i = 0; i < n_temp; i++) {
			if (param.form == "diag") {
				sn_save[i] += w * stat.xi[i + NR * n_temp];
			} else if (param.form == "full") {
				if (i + NR*n_temp == param.ith_rbin) {
					sn_save[i] += w * stat.xi[i + NR * n_temp];
				} else {
					sn_save[i] += 0.0;
				}
			}
	       	}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}

	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing three-point function...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_ylm_weighted_mean_density(P_R, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);


	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[n_temp];
	for (int i = 0; i < n_temp; i++) {
		zeta_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_ylm_weighted_mean_density(P_R, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double rmag1;
		if (param.form == "full") {
			rmag1 = rbin[param.ith_rbin];
			dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
		}

		for (int ir = 0; ir < n_temp; ir++) {

			double rmag2 = rbin[ir + NR * n_temp];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				rmag1 = rmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag2, Ylm2, sj2);

			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("r2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = ParticleBOSSClass::calc_norm_for_bispec(P_R, vol_survey);
	norm /= (alpha * alpha * alpha);

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/zeta_window_%d%d%d_%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, NR);
		fp = fopen(buf, "w");
		for (int i = 0; i < n_temp; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", rbin[i+NR*n_temp], rbin[i+NR*n_temp], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/zeta_window_%d%d%d_%02d_%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin, NR);
		fp = fopen(buf, "w");
		for (int i = 0; i < n_temp; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", rbin[param.ith_rbin], rbin[i+NR*n_temp], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	}

	fclose(fp);

	delete[] sn_save;
	delete[] zeta_save;

	return 0;
}

int calc_3pt_func_window_for_3pcf(
        ParticleBOSSClass & P_R,
        LineOfSight* los_rand,
        ParameterSet & param, double alpha, double * rbin, double vol_survey) {

	if (thisTask == 0) { printf("start to compute three-point window function...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.num_rbin];
	for (int i = 0; i < param.num_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_ylm_weighted_mean_density_for_3pt_window_shotnoise(P_R, los_rand, alpha, 0, 0);

	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_ylm_weighted_mean_density(P_R, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_2pt_func_window(P_R, los_rand, alpha, param.ELL, _M_);

		stat.calc_corr_func_for_3pt_func(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for (int i =0; i < param.num_rbin; i++) {
			if (param.form == "diag") {
				sn_save[i] += w * stat.xi[i];
			} else if (param.form == "full") {
				if (i == param.ith_rbin) {
					sn_save[i] += w * stat.xi[i];
				} else {
					sn_save[i] += 0.0;
				}
			}
	       	}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}

	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing three-point function...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_ylm_weighted_mean_density(P_R, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);


	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.num_rbin];
	for (int i = 0; i < param.num_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_ylm_weighted_mean_density(P_R, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double rmag1;
		if (param.form == "full") {
			rmag1 = rbin[param.ith_rbin];
			dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
		}

		for (int ir = 0; ir < param.num_rbin; ir++) {

			double rmag2 = rbin[ir];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				rmag1 = rmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag2, Ylm2, sj2);

			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("r2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = ParticleBOSSClass::calc_norm_for_bispec(P_R, vol_survey);
	norm /= (alpha * alpha * alpha);

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/zeta%d%d%d_window", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_rbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", rbin[i], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/zeta%d%d%d_window_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_rbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	}

	fclose(fp);

	delete[] sn_save;
	delete[] zeta_save;

	return 0;
}




int calc_bispec_in_box(ParticleBOSSClass & P_D, ParameterSet & param, double * kbin) {
	if (thisTask == 0) { printf("start to compute bispectrum...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.num_kbin];
	for (int i = 0; i < param.num_kbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00_sn(param);
	dn_00_sn.calc_density_field_in_box(P_D, param);
	/* Fourier transform*/
	dn_00_sn.calc_fourier_transform();

	/* calc shot noise terms */
	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> N_LM_sn(param);
		N_LM_sn.calc_density_field_in_box_for_bispec(P_D);
		/* Fourier transform*/
		N_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = double(P_D.n_tot);

		if ((param.ell1 == 0) && (param.ell2 == 0)) {
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[i] += w * shotnoise;
	       		}
		}

		if ((param.ell2 == 0)) {
			stat.calc_power_spec(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell1, _m1_);
			if (param.form == "diag") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[i] += w * stat.pk[i];
			       	}
			} else if (param.form == "full") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[i] += w * stat.pk[param.ith_kbin];
			       	}
			}
		}

		if ((param.ell1 == 0)) {
			stat.calc_power_spec(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell2, _m2_);
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[i] += w * stat.pk[i];
		       	}
		}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}}

	/****/
	dn_00_sn.finalise_density_field();
	/****************************************/

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_density_field_in_box_for_bispec(P_D);
	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
		ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_density_field_in_box(P_D, param);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = double(P_D.n_tot);

		fftw_complex * xi = fftw_alloc_complex(param.nmesh_tot);
		bytes += double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);
		for (int i = 0; i < (param.nmesh_tot); i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calc_shotnoise_for_bispec_ijk(dn_LM_sn, N_00_sn, shotnoise, param.ELL, _M_, xi);

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];
			double kmag1;
			if (param.form == "diag") {
				kmag1 = kmag2;
			} else if (param.form == "full") {
				kmag1 = kbin[param.ith_kbin];
			}

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
			double rvec[3];
			double dr[3];
			dr[0] = param.boxsize[0] / double(param.nmesh[0]);
			dr[1] = param.boxsize[1] / double(param.nmesh[1]);
			dr[2] = param.boxsize[2] / double(param.nmesh[2]);
			for (int i = 0; i < param.nmesh[0]; i++) {
			for (int j = 0; j < param.nmesh[1]; j++) {
			for (int k = 0; k < param.nmesh[2]; k++) {
				long long coord_flat = (i * param.nmesh[1] + j) * param.nmesh[2] + k;
				rvec[0] = (i < param.nmesh[0]/2) ? (double) i * dr[0] : (double) (i - param.nmesh[0]) * dr[0];
				rvec[1] = (j < param.nmesh[1]/2) ? (double) j * dr[1] : (double) (j - param.nmesh[1]) * dr[1];
				rvec[2] = (k < param.nmesh[2]/2) ? (double) k * dr[2] : (double) (k - param.nmesh[2]) * dr[2];
				double rmag = sqrt(rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[i][0], xi[i][1]);
				double j1 = sj1.eval(kmag1 * rmag);
				double j2 = sj2.eval(kmag2 * rmag);

				sn_sum += (j1 * j2 * ff * Ylm1[i] * Ylm2[i]);

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.nmesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[ik] += (w * sn_sum);

			/* time */
			durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		bytes -= double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}

	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing bispectrum...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_density_field_in_box(P_D, param);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();


	/* define bk_save:: bispectrum */
	std::complex<double> * bk_save = new std::complex<double>[param.num_kbin];
	for (int i = 0; i < param.num_kbin; i++) {
		bk_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;


		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
		ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_density_field_in_box(P_D, param);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		if (param.form == "full") {
			kmag1 = kbin[param.ith_kbin];
			dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag1, dk, Ylm1);
		}

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				kmag1 = kmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag1, dk, Ylm1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_bispec(dn_00, kmag2, dk, Ylm2);

			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[ik] += (w * bk_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}


		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);
	norm *= (param.volume/double(P_D.n_tot));

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/bk%d%d%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_kbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[i], kbin[i], norm * (bk_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/bk%d%d%d_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_kbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i], norm * (bk_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	}
	fclose(fp);

	delete[] sn_save;
	delete[] bk_save;

	return 0;
}

int calc_3pt_func_in_box(ParticleBOSSClass & P_D, ParameterSet & param, double * rbin) {
	if (thisTask == 0) { printf("start to compute three-point function...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.num_rbin];
	for (int i = 0; i < param.num_rbin; i++) {
		sn_save[i] = 0.0;
	}

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_density_field_in_box_for_bispec(P_D);
	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
		ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);


		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_density_field_in_box(P_D, param);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = double(P_D.n_tot);

		stat.calc_corr_func_for_3pt_func(dn_LM_sn, N_00_sn, rbin, shotnoise, param.ell1, _m1_, Ylm1, Ylm2);

		for (int i =0; i < param.num_rbin; i++) {
			if (param.form == "diag") {
				sn_save[i] += w * stat.xi[i];
			} else if (param.form == "full") {
				if (i == param.ith_rbin) {
					sn_save[i] += w * stat.xi[i];
				} else {
					sn_save[i] += 0.0;
				}
			}
	       	}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}


		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}

	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* three-point function */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing three-point function...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_density_field_in_box(P_D, param);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);

	/* define zeta_save:: bispectrum */
	std::complex<double> * zeta_save = new std::complex<double>[param.num_rbin];
	for (int i = 0; i < param.num_rbin; i++) {
		zeta_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		int _M_ = 0;

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
		ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);


		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_density_field_in_box(P_D, param);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double rmag1;
		if (param.form == "full") {
			rmag1 = rbin[param.ith_rbin];
			dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
		}

		for (int ir = 0; ir < param.num_rbin; ir++) {

			double rmag2 = rbin[ir];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				rmag1 = rmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag1, Ylm1, sj1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag2, Ylm2, sj2);

			/* calc bispectrum */
			std::complex<double> zeta_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			std::complex<double> _I_(0.0, 1.0);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				zeta_sum += (pow(_I_, param.ell1+param.ell2) * fac * f1 * f2 * f3);
			}

			zeta_save[ir] += (w * zeta_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("r2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", rmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}


		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = param.volume / double(P_D.n_tot) / double(P_D.n_tot);
	norm *= (param.volume / double(P_D.n_tot));

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/zeta%d%d%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_rbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", rbin[i], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/zeta%d%d%d_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_rbin);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_rbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n",  rbin[param.ith_rbin], rbin[i], norm * (zeta_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	}
	fclose(fp);

	delete[] sn_save;
	delete[] zeta_save;

	return 0;
}


int calc_bispecChoiceOfLOS(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        LineOfSight* los_data, LineOfSight* los_rand,
        ParameterSet & param, double alpha, double * kbin, int los, double vol_survey
       ) {

	if (thisTask == 0) { printf("start to compute bispectrum...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::complex<double> * sn_save = new std::complex<double>[param.num_kbin];
	for (int i = 0; i < param.num_kbin; i++) {
		sn_save[i] = 0.0;
	}

	/* calc shot noise terms */
	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************/
		/* calc the normal density fluctuation */
		/* dn = n - bar{n} */
		DensityField<ParticleBOSSClass> dn_sn(param);
		if (los == 0) {
			dn_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		} else {
			dn_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		}
		/* Fourier transform*/
		dn_sn.calc_fourier_transform();

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> N_sn(param);
		if (los == 0) {
			N_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		} else {
			N_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		}
		/* Fourier transform*/
		N_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		if ((param.ell1 == 0) && (param.ell2 == 0)) {
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[i] += w * shotnoise;
	       		}
		}

		if ((param.ell2 == 0)) {
			stat.calc_power_spec(dn_sn, N_sn, kbin, shotnoise, param.ell1, _m1_);
			if (param.form == "diag") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[i] += w * stat.pk[i];
			       	}
			} else if (param.form == "full") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[i] += w * stat.pk[param.ith_kbin];
			       	}
			}
		}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}}}


	/* calc shot noise terms */
	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************/
		/* calc the normal density fluctuation */
		/* dn = n - bar{n} */
		DensityField<ParticleBOSSClass> dn_sn(param);
		if (los == 1) {
			dn_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		} else {
			dn_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		}
		/* Fourier transform*/
		dn_sn.calc_fourier_transform();

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> N_sn(param);
		if (los == 1) {
			N_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		} else {
			N_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		}
		/* Fourier transform*/
		N_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		if ((param.ell1 == 0)) {
			stat.calc_power_spec(dn_sn, N_sn, kbin, shotnoise, param.ell2, _m2_);
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[i] += w * stat.pk[i];
		       	}
		}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}}}
















	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);

		}
	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************/
		/* calc N00 */
		DensityField<ParticleBOSSClass> N_sn(param);
		if (los == 2) {
			N_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		} else {
			N_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		}
		/* Fourier transform*/
		N_sn.calc_fourier_transform();

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_sn(param);
		if (los == 2) {
			dn_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		} else {
			dn_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		}
		/* Fourier transform*/
		dn_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		fftw_complex * xi = fftw_alloc_complex(param.nmesh_tot);
		bytes += double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);
		for (int i = 0; i < (param.nmesh_tot); i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calc_shotnoise_for_bispec_ijk(dn_sn, N_sn, shotnoise, param.ELL, _M_, xi);

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];
			double kmag1;
			if (param.form == "diag") {
				kmag1 = kmag2;
			} else if (param.form == "full") {
				kmag1 = kbin[param.ith_kbin];
			}

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
			double rvec[3];
			double dr[3];
			dr[0] = param.boxsize[0] / double(param.nmesh[0]);
			dr[1] = param.boxsize[1] / double(param.nmesh[1]);
			dr[2] = param.boxsize[2] / double(param.nmesh[2]);
			for (int i = 0; i < param.nmesh[0]; i++) {
			for (int j = 0; j < param.nmesh[1]; j++) {
			for (int k = 0; k < param.nmesh[2]; k++) {
				long long coord_flat = (i * param.nmesh[1] + j) * param.nmesh[2] + k;
				rvec[0] = (i < param.nmesh[0]/2) ? (double) i * dr[0] : (double) (i - param.nmesh[0]) * dr[0];
				rvec[1] = (j < param.nmesh[1]/2) ? (double) j * dr[1] : (double) (j - param.nmesh[1]) * dr[1];
				rvec[2] = (k < param.nmesh[2]/2) ? (double) k * dr[2] : (double) (k - param.nmesh[2]) * dr[2];
				double rmag = sqrt(rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[i][0], xi[i][1]);
				double j1 = sj1.eval(kmag1 * rmag);
				double j2 = sj2.eval(kmag2 * rmag);

				sn_sum += (j1 * j2 * ff * Ylm1[i] * Ylm2[i]);

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.nmesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[ik] += (w * sn_sum);

			/* time */
			durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		bytes -= double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}


	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing bispectrum...\n");}

	/* define bk_save:: bispectrum */
	std::complex<double> * bk_save = new std::complex<double>[param.num_kbin];
	for (int i = 0; i < param.num_kbin; i++) {
		bk_save[i] = 0.0;
	}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);

		}
	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/****************************************/
		/* calc the normal density fluctuation */
		/* dn = n - bar{n} */
		DensityField<ParticleBOSSClass> dn1(param);
		if (los == 0) {
			dn1.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		} else {
			dn1.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		}
		/* Fourier transform*/
		dn1.calc_fourier_transform();

		DensityField<ParticleBOSSClass> dn2(param);
		if (los == 1) {
			dn2.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		} else {
			dn2.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		}
		/* Fourier transform*/
		dn2.calc_fourier_transform();

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn3(param);
		if (los == 2) {
			dn3.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		} else {
			dn3.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
		}
		/* Fourier transform*/
		dn3.calc_fourier_transform();
		/* divided by a assignment function */
		dn3.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn3.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		if (param.form == "full") {
			kmag1 = kbin[param.ith_kbin];
				dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn1, kmag1, dk, Ylm1);
		}

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				kmag1 = kmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn1, kmag1, dk, Ylm1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_bispec(dn2, kmag2, dk, Ylm2);

			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn3[i][0], dn3[i][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[ik] += (w * bk_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = ParticleBOSSClass::calc_norm_for_bispec(P_D, vol_survey);

	FILE * fp;
	char buf[1024];
	if (param.form == "diag") {
		sprintf(buf, "%s/bk%d%d%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_kbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[i], kbin[i], norm * (bk_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	} else if (param.form == "full") {
		sprintf(buf, "%s/bk%d%d%d_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, param.ith_kbin);
		fp = fopen(buf, "w");
		for (int i = 0; i < param.num_kbin; i++) {
			fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i], norm * (bk_save[i].real() - sn_save[i].real()), norm * sn_save[i].real());
		}
	}
	fclose(fp);

	delete[] sn_save;
	delete[] bk_save;

	return 0;
}




int calc_bispecMmode(
        ParticleBOSSClass & P_D, ParticleBOSSClass & P_R,
        LineOfSight* los_data, LineOfSight* los_rand,
        ParameterSet & param, double alpha, double * kbin, double vol_survey) {

    if (thisTask == 0) { printf("start to compute bispectrum...\n");}

	if (fabs(wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0)) < 1.0e-10) {
		if (thisTask == 0) {
			printf("This multipole combination is not allowed.\n");
			printf("It should be wigner_3j(ell1, ell2, ELL, 0,0,0) != 0\n");
		}
		exit(1);
	}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* shotnoise */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("computing shotnoise terms...| %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/* define sn_save:: shotnoise term */
	std::vector< std::vector< std::complex<double> > > sn_save;

	sn_save.resize(2*param.ELL+1);
	for (int i = 0; i < param.num_kbin; i++) {
		sn_save[i].resize(param.num_kbin);
	}

	for (int i = 0; i < param.num_kbin; i++) {
	for (int m = 0; m < 2*param.ELL+1; m++) {
		sn_save[m][i] = 0.0;
	}}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00_sn(param);
	dn_00_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00_sn.calc_fourier_transform();

	/* calc shot noise terms */
	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {
	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> N_LM_sn(param);
		N_LM_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		N_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		if ((param.ell1 == 0) && (param.ell2 == 0)) {
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[_M_+param.ELL][i] += shotnoise;
	       		}
		}

		if ((param.ell2 == 0)) {
			stat.calc_power_spec(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell1, _m1_);
			if (param.form == "diag") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[_M_+param.ELL][i] += stat.pk[i];
			       	}
			} else if (param.form == "full") {
				for (int i =0; i < param.num_kbin; i++) {
					sn_save[_M_+param.ELL][i] += stat.pk[param.ith_kbin];
			       	}
			}
		}

		if ((param.ell1 == 0)) {
			stat.calc_power_spec(dn_00_sn, N_LM_sn, kbin, shotnoise, param.ell2, _m2_);
			for (int i =0; i < param.num_kbin; i++) {
				sn_save[_M_+param.ELL][i] += stat.pk[i];
		       	}
		}

		/* time */
		durationInSec = double(clock() - timeStart);
		if (thisTask == 0) { printf("m1 = %d, m2 = %d, M = %d| %.3f sec\n", _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

	}}}

	/****/
	dn_00_sn.finalise_density_field();
	/****************************************/

	/****************************************/
	/* calc F00 */
	DensityField<ParticleBOSSClass> N_00_sn(param);
	N_00_sn.calc_ylm_weighted_overdensity_for_bispec_shotnoise(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	N_00_sn.calc_fourier_transform();

	/* store spherical bessel functions */
	SphericalBesselCalculator sj1(param.ell1);
	SphericalBesselCalculator sj2(param.ell2);

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/**********************************************/
		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM_sn(param);
		dn_LM_sn.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM_sn.calc_fourier_transform();

		TwoPointStatistics<ParticleBOSSClass> stat(param);
		std::complex<double> shotnoise = stat.calc_shotnoise_for_bispec(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);

		fftw_complex * xi = fftw_alloc_complex(param.nmesh_tot);
		bytes += double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);
		for (int i = 0; i < (param.nmesh_tot); i++) {
			xi[i][0] = 0.0;
			xi[i][1] = 0.0;
		}

		stat.calc_shotnoise_for_bispec_ijk(dn_LM_sn, N_00_sn, shotnoise, param.ELL, _M_, xi);

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];
			double kmag1;
			if (param.form == "diag") {
				kmag1 = kmag2;
			} else if (param.form == "full") {
				kmag1 = kbin[param.ith_kbin];
			}

			/* calc shotnoise */
			std::complex<double> sn_sum = 0.0;
			double rvec[3];
			double dr[3];
			dr[0] = param.boxsize[0] / double(param.nmesh[0]);
			dr[1] = param.boxsize[1] / double(param.nmesh[1]);
			dr[2] = param.boxsize[2] / double(param.nmesh[2]);
			for (int i = 0; i < param.nmesh[0]; i++) {
			for (int j = 0; j < param.nmesh[1]; j++) {
			for (int k = 0; k < param.nmesh[2]; k++) {
				long long coord_flat = (i * param.nmesh[1] + j) * param.nmesh[2] + k;
				rvec[0] = (i < param.nmesh[0]/2) ? (double) i * dr[0] : (double) (i - param.nmesh[0]) * dr[0];
				rvec[1] = (j < param.nmesh[1]/2) ? (double) j * dr[1] : (double) (j - param.nmesh[1]) * dr[1];
				rvec[2] = (k < param.nmesh[2]/2) ? (double) k * dr[2] : (double) (k - param.nmesh[2]) * dr[2];
				double rmag = sqrt(rvec[0] * rvec[0] +  rvec[1] * rvec[1] +  rvec[2] * rvec[2]);
				std::complex<double> ff(xi[i][0], xi[i][1]);
				double j1 = sj1.eval(kmag1 * rmag);
				double j2 = sj2.eval(kmag2 * rmag);

				sn_sum += (j1 * j2 * ff * Ylm1[i] * Ylm2[i]);

			}}}

			std::complex<double> _I_(0.0,1.0);
			double fac = param.volume / double(param.nmesh_tot);
			sn_sum *= pow(_I_, param.ell1 + param.ell2) * fac;

			sn_save[_M_+param.ELL][ik] += (sn_sum);

			/* time */
			durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}
		}

		fftw_free(xi); xi = NULL;
		bytes -= double(sizeof(fftw_complex) * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

	}}

	/****/
	N_00_sn.finalise_density_field();
	/****************************************/

	/* time */
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) { printf("done | %.3f sec\n", durationInSec / CLOCKS_PER_SEC);}

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/

	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* bispectrum */
	/*****************************************************************************/
	/*****************************************************************************/
	/*****************************************************************************/
	/* time */
	if (thisTask == 0) { printf("computing bispectrum...\n");}

	/****************************************/
	/* calc the normal density fluctuation */
	/* dn = n - bar{n} */
	DensityField<ParticleBOSSClass> dn_00(param);
	dn_00.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, 0, 0);
	/* Fourier transform*/
	dn_00.calc_fourier_transform();


	/* define bk_save:: bispectrum */
//	std::complex<double> bk_save[2*param.ELL+1][param.num_kbin];
	std::vector< std::vector< std::complex<double> > > bk_save;

	bk_save.resize(2*param.ELL+1);
	for (int i = 0; i < param.num_kbin; i++) {
		bk_save[i].resize(param.num_kbin);
	}

	for (int i = 0; i < param.num_kbin; i++) {
	for (int m = 0; m < 2*param.ELL+1; m++) {
		bk_save[m][i] = 0.0;
	}}

	for (int _m1_ = - param.ell1; _m1_ <= param.ell1; _m1_++) {
	for (int _m2_ = - param.ell2; _m2_ <= param.ell2; _m2_++) {

		/****************************************************/
		/* store spherical harmonics */
		std::complex<double> * Ylm1 = new std::complex<double> [param.nmesh_tot];
		std::complex<double> * Ylm2 = new std::complex<double> [param.nmesh_tot];
		bytes += double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);

		std::string flag = "FALSE";
		for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {
			/*****/
			double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
			w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);

			if (fabs(w) > 1.0e-10) {
				flag = "TRUE";
			}
		}

		if (flag == "TRUE") {

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell1, _m1_, param, Ylm1);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(param.ell2, _m2_, param, Ylm2);

		}

	for (int _M_ = - param.ELL; _M_ <= param.ELL; _M_++) {

		/*****/
		double w = wigner_3j(param.ell1, param.ell2, param.ELL, 0, 0, 0) * wigner_3j(param.ell1, param.ell2, param.ELL, _m1_, _m2_, _M_);
		w *= double(2*param.ELL+1) * double(2*param.ell1+1) * double(2*param.ell2+1);
		if (fabs(w) < 1.0e-10) {
			continue;
		}

		/* calc the yLM-weighted density perturbation */
		DensityField<ParticleBOSSClass> dn_LM(param);
		dn_LM.calc_ylm_weighted_overdensity(P_D, P_R, los_data, los_rand, alpha, param.ELL, _M_);
		/* Fourier transform*/
		dn_LM.calc_fourier_transform();
		/* divided by a assignment function */
		dn_LM.calc_assignment_compensation();
		/* Inverse Fourier transform*/
		dn_LM.calc_inverse_fourier_transform();

		/* calc dn_tilde1 */
		DensityField<ParticleBOSSClass> dn_tilde1(param);
		double kmag1;
		double dk = kbin[1] - kbin[0];
		if (param.form == "full") {
			kmag1 = kbin[param.ith_kbin];
			dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag1, dk, Ylm1);
		}

		for (int ik = 0; ik < param.num_kbin; ik++) {

			double kmag2 = kbin[ik];

			/* calc dn_tilde1 */
			if (param.form == "diag") {
				kmag1 = kmag2;
				dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag1, dk, Ylm1);
			}

			/* calc dn_tilde2 */
			DensityField<ParticleBOSSClass> dn_tilde2(param);
			dn_tilde2.calc_inverse_fourier_transform_for_bispec(dn_00, kmag2, dk, Ylm2);

			/* calc bispectrum */
			std::complex<double> bk_sum = 0.0;
			double fac = param.volume / double(param.nmesh_tot);
			for (int i = 0; i < param.nmesh_tot; i++) {
				std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);
				std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);
				std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);
				bk_sum += (fac * f1 * f2 * f3);
			}

			bk_save[_M_+param.ELL][ik] += (bk_sum);

			double durationInSec = double(clock() - timeStart);
			if (thisTask == 0) { printf("k2 = %.3f, m1 = %d, m2 = %d, M = %d| %.3f sec\n", kmag2, _m1_, _m2_, _M_, durationInSec / CLOCKS_PER_SEC);}

		}

	}

		delete[] Ylm1; Ylm1 = NULL;
		delete[] Ylm2; Ylm2 = NULL;
		bytes -= double(sizeof(std::complex<double>) * 2 * (param.nmesh_tot) / 1024.0 / 1024.0 / 1024.0);


	}}

	double norm = ParticleBOSSClass::calc_norm_for_bispec(P_D, vol_survey);

	for (int _M_ = 0; _M_<2*param.ELL+1;_M_++) {
		FILE * fp;
		char buf[1024];
		if (param.form == "diag") {
			sprintf(buf, "%s/bk%d%d%d_M%d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, _M_);
			fp = fopen(buf, "w");
			for (int i = 0; i < param.num_kbin; i++) {
				fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[i], kbin[i], norm * (bk_save[_M_][i].real() - sn_save[_M_][i].real()), norm * sn_save[_M_][i].real());
			}
		} else if (param.form == "full") {
			sprintf(buf, "%s/bk%d%d%d_M%d_%02d", param.output_dir.c_str(), param.ell1, param.ell2, param.ELL, _M_, param.ith_kbin);
			fp = fopen(buf, "w");
			for (int i = 0; i < param.num_kbin; i++) {
				fprintf(fp, "%.5f \t %.5f \t %.7e \t %.7e\n", kbin[param.ith_kbin], kbin[i], norm * (bk_save[_M_][i].real() - sn_save[_M_][i].real()), norm * sn_save[_M_][i].real());
			}
		}
		fclose(fp);
	}

	return 0;
}

#endif
