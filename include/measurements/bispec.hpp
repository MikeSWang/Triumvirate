#ifndef TRIUM_BISPEC_H_INCLUDED_
#define TRIUM_BISPEC_H_INCLUDED_

/**
 * Calculate bispectrum from catalogues.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params (Reference to) the input parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  ParameterSet& params,
  double alpha,
  double* kbin,
  double survey_vol_norm
) {
  if (thisTask == 0) {
    printf("[Status] :: Measuring bispectrum.\n");
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
  ) {
    if (thisTask == 0) {
      printf(
        "[Error] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  durationInSec = double(clock() - timeStart);
  if (thisTask == 0) {
    printf(
      "[Status] :: Computing shot noise terms "
      "(... %.3f seconds elapsed in total).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    shotnoise_save[i] = 0.;
  }

  /// Compute shot noise terms.
  DensityField<ParticleCatalogue> dn_00_shotnoise(params);
  dn_00_shotnoise.calc_ylm_weighted_fluctuation(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  dn_00_shotnoise.calc_fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < 1.e-10) {
          continue;
        }

        DensityField<ParticleCatalogue> shotnoise_LM(params);
        shotnoise_LM.calc_ylm_weighted_fields_for_bispec_shotnoise(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        shotnoise_LM.calc_fourier_transform();

        TwoPointStatistics<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );

        if (params.ell1 == 0 && params.ell2 == 0) {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * shotnoise;
          }
        }

        if (params.ell2 == 0) {
          stats.calc_power_spec(
            dn_00_shotnoise, shotnoise_LM, kbin, shotnoise, params.ell1, m1_
          );
          if (params.form == "diag") {
            for (int i = 0; i < params.num_kbin; i++) {
              shotnoise_save[i] += coupling * stats.pk[i];
            }
          } else if (params.form == "full") {
            for (int i = 0; i < params.num_kbin; i++) {
              shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
            }
          }
        }

        if (params.ell1 == 0) {
          stats.calc_power_spec(
            dn_00_shotnoise, shotnoise_LM, kbin, shotnoise, params.ell2, m2_
          );
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[i];
          }
        }

        durationInSec = double(clock() - timeStart);
        if (thisTask == 0) {
          printf(
            "[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
            "(... %.3f seconds elapsed in total).\n",
            m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
          );
        }
      }
    }
  }

  dn_00_shotnoise.finalise_density_field();

  DensityField<ParticleCatalogue> shotnoise_00(params);
  shotnoise_00.calc_ylm_weighted_fields_for_bispec_shotnoise(
    particles_data, particles_rand,
    los_data, los_rand,
    alpha,
    0, 0
  );
  shotnoise_00.calc_fourier_transform();

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
      bytes += 2 * sizeof(std::complex<double>)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > 1.e-10) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        ToolCollection::store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params, ylm_a
        );
        ToolCollection::store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < 1.e-10) {
          continue;
        }

        DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
        dn_LM_shotnoise.calc_ylm_weighted_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM_shotnoise.calc_fourier_transform();

        TwoPointStatistics<ParticleCatalogue> stats(params);
        std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );

        fftw_complex* xi = fftw_alloc_complex(params.nmesh_tot);
        bytes += sizeof(fftw_complex) *
          double(params.nmesh_tot) / 1024. / 1024. / 1024.;
        for (int i = 0; i < params.nmesh_tot; i++) {
          xi[i][0] = 0.;
          xi[i][1] = 0.;
        }

        stats.calc_shotnoise_for_bispec_ijk(
          dn_LM_shotnoise, shotnoise_00, shotnoise, params.ELL, M_, xi
        );

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_a;
          double kmag_b = kbin[i_kbin];
          if (params.form == "diag") {
            kmag_a = kmag_b;
          } else if (params.form == "full") {
            kmag_a = kbin[params.ith_kbin];
          }

          std::complex<double> shotnoise_sum = 0.;
          double rvec[3];
          double dr[3];
          dr[0] = params.boxsize[0] / double(params.nmesh[0]);
          dr[1] = params.boxsize[1] / double(params.nmesh[1]);
          dr[2] = params.boxsize[2] / double(params.nmesh[2]);
          for (int i = 0; i < params.nmesh[0]; i++) {
            for (int j = 0; j < params.nmesh[1]; j++) {
              for (int k = 0; k < params.nmesh[2]; k++) {
                long long coord_flat =
                  (i * params.nmesh[1] + j) * params.nmesh[2] + k;
                rvec[0] = (i < params.nmesh[0]/2) ?
                  i * dr[0] : (i - params.nmesh[0]) * dr[0];
                rvec[1] = (j < params.nmesh[1]/2) ?
                  j * dr[1] : (j - params.nmesh[1]) * dr[1];
                rvec[2] = (k < params.nmesh[2]/2) ?
                  k * dr[2] : (k - params.nmesh[2]) * dr[2];
                double rmag = sqrt(
                  rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
                );

                std::complex<double> ff(xi[coord_flat][0], xi[coord_flat][1]);  // ???
                double j1 = sj1.eval(kmag_a * rmag);
                double j2 = sj2.eval(kmag_b * rmag);

                shotnoise_sum += j1 * j2 * ylm_a[i] * ylm_b[i] * ff;
              }
            }
          }

          std::complex<double> I_(0., 1.);
          double factor = params.volume / double(params.nmesh_tot);
          shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

          shotnoise_save[i_kbin] += coupling * shotnoise_sum;

          durationInSec = double(clock() - timeStart);
          if (thisTask == 0) {
            printf(
              "[Status] :: Computed for wavenumber and orders "
              "k2 = %.3f, m1 = %d, m2 = %d, M = %d "
              "(... %.3f seconds elapsed in total).\n",
              kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
            );
          }
        }

        fftw_free(xi); xi = NULL;
        bytes -= sizeof(fftw_complex)
          * double(params.nmesh_tot) / 1024. / 1024. / 1024.;
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      bytes -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;
    }
  }

  shotnoise_00.finalise_density_field();

  durationInSec = double(clock() - timeStart);
  if (thisTask == 0) {
    printf(
      "[Status] :: Computed shot noise terms "
      "(... %.3f seconds elapsed in total).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// Initialise output bispectrum.
  if (thisTask == 0) {
    printf("[Status] :: Measuring bispectrum.\n");
  }

  std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    bk_save[i] = 0.;
  }

  /// Compute bispectrum.
  DensityField<ParticleCatalogue> dn_00(params);
  dn_00.calc_ylm_weighted_fluctuation(
    particles_data, particles_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.calc_fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      std::complex<double>* ylm_a = new std::complex<double> [params.nmesh_tot];
      std::complex<double>* ylm_b = new std::complex<double> [params.nmesh_tot];
      bytes += 2 * sizeof(std::complex<double>)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;

      std::string flag_nontrivial = "FALSE";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) > 1.e-10) {
          flag_nontrivial = "TRUE";
        }
      }

      if (flag_nontrivial == "TRUE") {
        ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params, ylm_a
        );
        ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params, ylm_b
        );
      }

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = double(2*params.ELL + 1)
          * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
          * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
          * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
        if (fabs(coupling) < 1.e-10) {
          continue;
        }

        DensityField<ParticleCatalogue> dn_LM(params);
        dn_LM.calc_ylm_weighted_fluctuation(
          particles_data, particles_rand,
          los_data, los_rand,
          alpha,
          params.ELL, M_
        );
        dn_LM.calc_fourier_transform();
        dn_LM.calc_assignment_compensation();
        dn_LM.calc_inverse_fourier_transform();

        DensityField<ParticleCatalogue> dn_tilde1(params);
        double kmag_a;
        double dk = kbin[1] - kbin[0];
        if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
          /// Compute ``\delta\tilde{n}_1``.
          dn_tilde1.calc_inverse_fourier_transform_for_bispec(
            dn_00, kmag_a, dk, ylm_a
          );
        }

        for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
          double kmag_b = kbin[i_kbin];

          /// Compute ``\delta\tilde{n}_1``.
          if (params.form == "diag") {
            kmag_a = kmag_b;
            dn_tilde1.calc_inverse_fourier_transform_for_bispec(
              dn_00, kmag_a, dk, ylm_a
            );
          }

          /// Compute ``\delta\tilde{n}_2``.
          DensityField<ParticleCatalogue> dn_tilde2(params);
          dn_tilde2.calc_inverse_fourier_transform_for_bispec(
            dn_00, kmag_b, dk, ylm_b
          );

          double factor = params.volume / double(params.nmesh_tot);
          std::complex<double> bk_sum = 0.;
          for (int i = 0; i < params.nmesh_tot; i++) {
            std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
            std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
            std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
            bk_sum += factor * f1 * f2 * f3;
          }

          bk_save[i_kbin] += coupling * bk_sum;

          double durationInSec = double(clock() - timeStart);
          if (thisTask == 0) {
            printf(
              "[Status] :: Computed for wavenumber and order "
              "k2 = %.3f, m1 = %d, m2 = %d, M = %d "
              "(... %.3f seconds elapsed in total).\n",
              kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
            );
          }
        }
      }

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
      bytes -= 2 * sizeof(std::complex<double>)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;
    }
  }

  /// Normalise and then save the output.
  double norm = ParticleCatalogue::calc_norm_for_bispec(
    particles_data, survey_vol_norm
  );

  FILE* saved_file_ptr;
  char buf[1024];
  if (params.form == "diag") {
    sprintf(
      buf, "%s/bk%d%d%d",
      params.output_dir.c_str(), params.ell1, params.ell2, params.ELL
    );
    saved_file_ptr = fopen(buf, "w");
    for (int i = 0; i < params.num_kbin; i++) {
      fprintf(
        saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[i], kbin[i],
        norm * (bk_save[i].real() - shotnoise_save[i].real()),
        norm * (bk_save[i].imag() - shotnoise_save[i].imag()),
        norm * shotnoise_save[i].real(),
        norm * shotnoise_save[i].imag()
      );
    }
  } else if (params.form == "full") {
    sprintf(
      buf, "%s/bk%d%d%d_%02d",
      params.output_dir.c_str(),
      params.ell1, params.ell2, params.ELL,
      params.ith_kbin
    );
    saved_file_ptr = fopen(buf, "w");
    for (int i = 0; i < params.num_kbin; i++) {
      fprintf(
        saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e \t %.7e \t %.7e\n",
        kbin[params.ith_kbin], kbin[i],
        norm * (bk_save[i].real() - shotnoise_save[i].real()),
        norm * (bk_save[i].imag() - shotnoise_save[i].imag()),
        norm * shotnoise_save[i].real(),
        norm * shotnoise_save[i].imag()
      );
    }
  }
  fclose(saved_file_ptr);

  delete[] shotnoise_save;
  delete[] bk_save;

  return 0;
}

/**
 * Calculate bispectrum in a periodic box.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param params (Reference to) the input parameter set.
 * @param kbin Wavenumber bins.
 * @returns Exit status.
 */
int calc_bispec_in_box(
  ParticleCatalogue& particles_data,
  ParameterSet& params,
  double* kbin
) {
  if (thisTask == 0) {
    printf("[Status] :: Measuring bispectrum in box.\n");
  }

  if (
    fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
  ) {
    if (thisTask == 0) {
      printf(
        "[Error] :: Disallowed multipole degree combination "
        "for bispectrum measurements. "
        "Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
      );
    }
    exit(1);
  }

  /// Initialise output shot noise terms.
  durationInSec = double(clock() - timeStart);
  if (thisTask == 0) {
    printf(
      "[Status] :: Computing shot noise terms "
      "(... %.3f seconds elapsed in total).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  std::complex<double>* shotnoise_save =
    new std::complex<double>[params.num_kbin];
  for (int i = 0; i < params.num_kbin; i++) {
    shotnoise_save[i] = 0.;
  }

  /// Compute shot noise terms.
  DensityField<ParticleCatalogue> dn_00_shotnoise(params);
  dn_00_shotnoise.calc_fluctuation_in_box(particles_data, params);
  dn_00_shotnoise.calc_fourier_transform();

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;

      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (fabs(coupling) < 1.e-10) {
        continue;
      }

      DensityField<ParticleCatalogue> shotnoise_LM(params);
      shotnoise_LM.calc_density_field_in_box_for_bispec(particles_data);
      shotnoise_LM.calc_fourier_transform();

      TwoPointStatistics<ParticleCatalogue> stats(params);
      std::complex<double> shotnoise = double(particles_data.nparticles);
        // NOTE: ``double`` conversion essential here

      if (params.ell1 == 0 && params.ell2 == 0) {
        for (int i = 0; i < params.num_kbin; i++) {
          shotnoise_save[i] += coupling * shotnoise;
        }
      }

      if (params.ell2 == 0) {
        stats.calc_power_spec(
          dn_00_shotnoise, shotnoise_LM,
          kbin,
          shotnoise,
          params.ell1, m1_
        );
        if (params.form == "diag") {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[i];
          }
        } else if (params.form == "full") {
          for (int i = 0; i < params.num_kbin; i++) {
            shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
          }
        }
      }

      if (params.ell1 == 0) {
        stats.calc_power_spec(
          dn_00_shotnoise, shotnoise_LM,
          kbin,
          shotnoise,
          params.ell2, m2_
        );
        for (int i = 0; i < params.num_kbin; i++) {
          shotnoise_save[i] += coupling * stats.pk[i];
        }
      }

      durationInSec = double(clock() - timeStart);
      if (thisTask == 0) {
        printf(
          "[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
          "(... %.3f seconds elapsed in total).\n",
          m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
        );
      }
    }
  }

  dn_00_shotnoise.finalise_density_field();

  DensityField<ParticleCatalogue> shotnoise_00(params);
  shotnoise_00.calc_density_field_in_box_for_bispec(particles_data);
  shotnoise_00.calc_fourier_transform();

  SphericalBesselCalculator sj1(params.ell1);
  SphericalBesselCalculator sj2(params.ell2);

  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      int M_ = 0;

      double coupling = double(2*params.ELL + 1)
        * double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
        * wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
        * wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
      if (fabs(coupling) < 1.e-10) {
        continue;
      }

      std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
      std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
      bytes += 2 * sizeof(std::complex<double>)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;

      ToolCollection::store_reduced_spherical_harmonic_in_config_space(
        params.ell1, m1_, params, ylm_a
      );
      ToolCollection::store_reduced_spherical_harmonic_in_config_space(
        params.ell2, m2_, params, ylm_b
      );

      DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
      dn_LM_shotnoise.calc_fluctuation_in_box(particles_data, params);
      dn_LM_shotnoise.calc_fourier_transform();

      TwoPointStatistics<ParticleCatalogue> stats(params);
      std::complex<double> shotnoise = double(particles_data.nparticles);

      fftw_complex* xi = fftw_alloc_complex(params.nmesh_tot);
      bytes += sizeof(fftw_complex)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;
      for (int i = 0; i < params.nmesh_tot; i++) {
        xi[i][0] = 0.;
        xi[i][1] = 0.;
      }

      stats.calc_shotnoise_for_bispec_ijk(
        dn_LM_shotnoise, shotnoise_00, shotnoise,
        params.ELL, M_,
        xi
      );

      for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
        double kmag_b = kbin[i_kbin];
        double kmag_a;
        if (params.form == "diag") {
          kmag_a = kmag_b;
        } else if (params.form == "full") {
          kmag_a = kbin[params.ith_kbin];
        }

        double dr[3];
        dr[0] = params.boxsize[0] / double(params.nmesh[0]);
        dr[1] = params.boxsize[1] / double(params.nmesh[1]);
        dr[2] = params.boxsize[2] / double(params.nmesh[2]);

        double rvec[3];
        std::complex<double> shotnoise_sum = 0.;
        for (int i = 0; i < params.nmesh[0]; i++) {
          for (int j = 0; j < params.nmesh[1]; j++) {
            for (int k = 0; k < params.nmesh[2]; k++) {
              long long coord_flat =
                (i * params.nmesh[1] + j) * params.nmesh[2] + k;
              rvec[0] = (i < params.nmesh[0]/2) ?
                i * dr[0] : (i - params.nmesh[0]) * dr[0];
              rvec[1] = (j < params.nmesh[1]/2) ?
                j * dr[1] : (j - params.nmesh[1]) * dr[1];
              rvec[2] = (k < params.nmesh[2]/2) ?
                k * dr[2] : (k - params.nmesh[2]) * dr[2];

              double rmag = sqrt(
                rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
              );

              std::complex<double> ff(xi[coord_flat][0], xi[coord_flat][1]);  // ???
              double j1 = sj1.eval(kmag_a * rmag);
              double j2 = sj2.eval(kmag_b * rmag);

              shotnoise_sum += j1 * j2 * ff * ylm_a[i] * ylm_b[i];
            }
          }
        }

        std::complex<double> I_(0., 1.);
        double factor = params.volume / double(params.nmesh_tot);
        shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

        shotnoise_save[i_kbin] += coupling * shotnoise_sum;

        durationInSec = double(clock() - timeStart);
        if (thisTask == 0) {
          printf(
            "[Status] :: Computed for wavenumber and orders "
            "k2 = %.3f, m1 = %d, m2 = %d, M = %d "
            "(... %.3f seconds elapsed in total).\n",
            kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
          );
        }
      }

      fftw_free(xi); xi = NULL;
      bytes -= sizeof(fftw_complex)
        * double(params.nmesh_tot) / 1024. / 1024. / 1024.;

      delete[] ylm_a; ylm_a = NULL;
      delete[] ylm_b; ylm_b = NULL;
  		bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	if (thisTask == 0) {
		printf("[Status] :: Measuring bispectrum.\n");
	}

	std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		bk_save[i] = 0.;
	}

	/// Compute bispectrum.
	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_fluctuation_in_box(particles_data, params);
	dn_00.calc_fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			int M_ = 0;

			double coupling = double(2*params.ELL + 1)
				* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell1, m1_, params, ylm_a
			);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell2, m2_, params, ylm_b
			);

			DensityField<ParticleCatalogue> dn_LM(params);
			dn_LM.calc_fluctuation_in_box(particles_data, params);
			dn_LM.calc_fourier_transform();
			dn_LM.calc_assignment_compensation();
			dn_LM.calc_inverse_fourier_transform();

			DensityField<ParticleCatalogue> dn_tilde1(params);
			double kmag_a;
			double dk = kbin[1] - kbin[0];
			if (params.form == "full") {
				/// Compute ``\delta\tilde{n}_1``.
				kmag_a = kbin[params.ith_kbin];
				dn_tilde1.calc_inverse_fourier_transform_for_bispec(
					dn_00, kmag_a, dk, ylm_a
				);
			}

			for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
				double kmag_b = kbin[i_kbin];

				/// Compute ``\delta\tilde{n}_1``.
				if (params.form == "diag") {
					kmag_a = kmag_b;
					dn_tilde1.calc_inverse_fourier_transform_for_bispec(
						dn_00, kmag_a, dk, ylm_a
					);
				}

				/// Compute ``\delta\tilde{n}_2``.
				DensityField<ParticleCatalogue> dn_tilde2(params);
				dn_tilde2.calc_inverse_fourier_transform_for_bispec(
					dn_00, kmag_b, dk, ylm_b
				);

				double factor = params.volume / double(params.nmesh_tot);
				std::complex<double> bk_sum = 0.;
				for (int i = 0; i < params.nmesh_tot; i++) {
					std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
					std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
					std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
					bk_sum += factor * f1 * f2 * f3;
				}

				bk_save[i_kbin] += coupling * bk_sum;

				double durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for wavenumber and order "
						"k2 = %.3f, m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = params.volume
		/ double(particles_data.nparticles) / double(particles_data.nparticles);
	norm *= params.volume / double(particles_data.nparticles);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/bk%d%d%d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[i], kbin[i],
				norm * (bk_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/bk%d%d%d_%02d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_kbin
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[params.ith_kbin], kbin[i],
				norm * (bk_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] bk_save;

	return 0;
}

/**
 * Calculate three-point function from catalogues.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params (Reference to) the input parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_func(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (thisTask == 0) {
		printf("[Status] :: Measuring three-point function.\n");
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (thisTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point function measurements. "
				"Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_00(params);
	shotnoise_00.calc_ylm_weighted_fields_for_bispec_shotnoise(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	shotnoise_00.calc_fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
				dn_LM_shotnoise.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);

				stats.calc_corr_func_for_3pt_func(
					dn_LM_shotnoise, shotnoise_00,
					rbin, shotnoise,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < params.num_rbin; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else if (params.form == "full") {
						if (i == params.ith_rbin) {
							shotnoise_save[i] += coupling * stats.xi[i];
						} else {
							shotnoise_save[i] += 0.;
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point function.
	if (thisTask == 0) {
		printf("[Status] :: Measuring three-point function.\n");
	}

	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00.calc_fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double> [params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double> [params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.calc_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				DensityField<ParticleCatalogue> dn_tilde1(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];
					/// Compute ``\delta\tilde{n}_1``.
					if (params.form == "diag") {
						rmag_a = rmag_b;
						dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					/// Compute ``\delta\tilde{n}_2``.
					DensityField<ParticleCatalogue> dn_tilde2(params);
					dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(
						dn_00, rmag_b, ylm_b, sj2
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh_tot);
					for (int i = 0; i < params.nmesh_tot; i++) {
						std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
						std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
						std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
						zeta_sum += pow(I_, params.ell1+params.ell2) * factor * f1 * f2 * f3;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double durationInSec = double(clock() - timeStart);
					if (thisTask == 0) {
						printf(
							"[Status] :: Computed for separation and order "
							"r2 = %.3f, m1 = %d, m2 = %d, M = %d "
							"(... %.3f seconds elapsed in total).\n",
							rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
						);
					}
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = ParticleCatalogue::calc_norm_for_bispec(
		particles_data, survey_vol_norm
	);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/zeta%d%d%d_%02d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate three-point function window from catalogues.
 *
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_rand Random-source particle lines of sight.
 * @param params (Reference to) the input parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_func_window(
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (thisTask == 0) {
		printf("[Status] :: Measuring three-point function window.\n");
	}

	int n_temp = 10;  // ???: NOTE: discretionary
	int NR = 3;  // ???

	rbin[0] = 0.;  // ???
	rbin[1] = 1.;
	rbin[2] = 10.;
	rbin[3] = 20.;
	rbin[4] = 30.;
	rbin[5] = 40.;
	rbin[6] = 50.;
	rbin[7] = 60.;
	double rmin = 70.;
	double dlnr = (log(params.rmax) - log(rmin))
		/ double((params.num_rbin - 8) - 1);

	params.ith_rbin = thisTask;  // ???

	for (int i = 8; i < params.num_rbin; i++) {
		rbin[i] = rmin * exp(dlnr * (i - 8));
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (thisTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point function window measurements. "
				"Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save = new std::complex<double>[n_temp];
	for (int i = 0; i < n_temp; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_00(params);
	shotnoise_00.calc_ylm_weighted_mean_density_for_3pt_window_shotnoise(
		particles_rand, los_rand, alpha, 0, 0
	);
	shotnoise_00.calc_fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double> [params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double> [params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
				dn_LM_shotnoise.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_)
				;
				dn_LM_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise =
					stats.calc_shotnoise_for_2pt_func_window(
						particles_rand, los_rand, alpha, params.ELL, M_
					);

				stats.calc_corr_func_for_3pt_func(
					dn_LM_shotnoise, shotnoise_00,
					rbin,
					shotnoise,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < n_temp; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i + NR * n_temp];
					} else if (params.form == "full") {
						if (i + NR * n_temp == params.ith_rbin) {  // ???
							shotnoise_save[i] += coupling * stats.xi[i + NR * n_temp];
						} else {
							shotnoise_save[i] += 0.;  // ??? redundant
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point function window.
	if (thisTask == 0) {
		printf("[Status] :: Measuring three-point function window.\n");
	}

	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	dn_00.calc_fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);

	std::complex<double>* zeta_save = new std::complex<double>[n_temp];
	for (int i = 0; i < n_temp; i++) {
		zeta_save[i] = 0.;
	}

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.calc_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				DensityField<ParticleCatalogue> dn_tilde1(params);
				double rmag_a;
				if (params.form == "full") {
					rmag_a = rbin[params.ith_rbin];
					dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				for (int i_rbin = 0; i_rbin < n_temp; i_rbin++) {
					double rmag_b = rbin[i_rbin + NR * n_temp];
					/// Compute ``\delta\tilde{n}_1``.
					if (params.form == "diag") {
						rmag_a = rmag_b;
						dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					/// Compute ``\delta\tilde{n}_2``.
					DensityField<ParticleCatalogue> dn_tilde2(params);
					dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(
						dn_00, rmag_b, ylm_b, sj2
					);

					std::complex<double> I_(0., 1.);
					std::complex<double> zeta_sum = 0.;
					double factor = params.volume / double(params.nmesh_tot);
					for (int i = 0; i < params.nmesh_tot; i++) {
						std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
						std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
						std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
						zeta_sum += pow(I_, params.ell1+params.ell2) * factor * f1 * f2 * f3;
					}

					zeta_save[i_rbin] += coupling * zeta_sum;

					double durationInSec = double(clock() - timeStart);
					if (thisTask == 0) {
						printf(
							"[Status] :: Computed for separation and order "
							"r2 = %.3f, m1 = %d, m2 = %d, M = %d "
							"(... %.3f seconds elapsed in total).\n",
							rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
						);
					}
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = ParticleCatalogue::calc_norm_for_bispec(
			particles_rand, survey_vol_norm
	);
	norm /= alpha * alpha * alpha;

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta_window_%d%d%d_%d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			NR
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < n_temp; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i+NR*n_temp], rbin[i+NR*n_temp],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/zeta_window_%d%d%d_%02d_%d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin, NR
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < n_temp; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i+NR*n_temp],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate three-point function window for three-point window function
 * from catalogues.
 *
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_rand Random-source particle lines of sight.
 * @param params (Reference to) the input parameter set.
 * @param alpha Alpha ratio.
 * @param rbin Separation bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_3pt_func_window_for_3pcf(  // ???
	ParticleCatalogue& particles_rand,
	LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* rbin,
	double survey_vol_norm
) {
	if (thisTask == 0) {
		printf(
			"[Status] :: Measuring three-point function window "
			"for three-point function.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (thisTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point function window measurements. "
				"Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_00(params);
	shotnoise_00.calc_ylm_weighted_mean_density_for_3pt_window_shotnoise(
		particles_rand, los_rand, alpha, 0, 0
	);
	shotnoise_00.calc_fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
				dn_LM_shotnoise.calc_ylm_weighted_mean_density(
					particles_rand, los_rand, alpha, params.ELL, M_
				);
				dn_LM_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise =
					stats.calc_shotnoise_for_2pt_func_window(
						particles_rand, los_rand, alpha, params.ELL, M_
					);

				stats.calc_corr_func_for_3pt_func(
					dn_LM_shotnoise, shotnoise_00,
					rbin, shotnoise,
					params.ell1, m1_,
					ylm_a, ylm_b
				);

				for (int i = 0; i < params.num_rbin; i++) {
					if (params.form == "diag") {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else if (params.form == "full") {
						if (i == params.ith_rbin) {
							shotnoise_save[i] += coupling * stats.xi[i];
						} else {
							shotnoise_save[i] += 0.;
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
					* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output three-point function window for three-point function.
	if (thisTask == 0) {
		printf(
			"[Status] :: Measuring three-point function window "
			"for three-point function.\n"
		);
	}

	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, 0, 0);
	dn_00.calc_fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double> [params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double> [params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_mean_density(particles_rand, los_rand, alpha, params.ELL, M_);
				dn_LM.calc_fourier_transform();
				dn_LM.calc_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				DensityField<ParticleCatalogue> dn_tilde1(params);
				double rmag_a;
				if (params.form == "full") {
					/// Compute ``\delta\tilde{n}_1``.
					rmag_a = rbin[params.ith_rbin];
					dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(dn_00, rmag_a, ylm_a, sj1);
				}

				for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
					double rmag_b = rbin[i_rbin];
					/// Compute ``\delta\tilde{n}_1``.
					if (params.form == "diag") {
						rmag_a = rmag_b;
						dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
							dn_00, rmag_a, ylm_a, sj1
						);
					}

					/// Compute ``\delta\tilde{n}_2``.
					DensityField<ParticleCatalogue> dn_tilde2(params);
					dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(
						dn_00, rmag_b, ylm_b, sj2
					);

				std::complex<double> I_(0., 1.);
				std::complex<double> zeta_sum = 0.;
				double factor = params.volume / double(params.nmesh_tot);
				for (int i = 0; i < params.nmesh_tot; i++) {
					std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
					std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
					std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
					zeta_sum += pow(I_, params.ell1+params.ell2) * factor * f1 * f2 * f3;
				}

				zeta_save[i_rbin] += (coupling * zeta_sum);

				double durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for separation and order "
						"r2 = %.3f, m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}
		}

		delete[] ylm_a; ylm_a = NULL;
		delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = ParticleCatalogue::calc_norm_for_bispec(
		particles_rand, survey_vol_norm
	);
	norm /= alpha * alpha * alpha;

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d_window",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/zeta%d%d%d_window_%02d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate three-point function in a periodic box.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param params (Reference to) the input parameter set.
 * @param rbin Separation bins.
 * @returns Exit status.
 */
int calc_3pt_func_in_box(
	ParticleCatalogue& particles_data,
	ParameterSet& params,
	double* rbin
) {
	if (thisTask == 0) {
		printf("[Status] :: Measuring three-point function in box.\n");
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (thisTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for three-point function measurements. "
				"Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> shotnoise_00(params);
	shotnoise_00.calc_density_field_in_box_for_bispec(particles_data);
	shotnoise_00.calc_fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			int M_ = 0;

			double coupling = double(2*params.ELL + 1)
				* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			std::complex<double>* ylm_a = new std::complex<double> [params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double> [params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			ToolCollection::store_reduced_spherical_harmonic_in_config_space(
				params.ell1, m1_, params, ylm_a
			);
			ToolCollection::store_reduced_spherical_harmonic_in_config_space(
				params.ell2, m2_, params, ylm_b
			);

			DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
			dn_LM_shotnoise.calc_fluctuation_in_box(particles_data, params);
			dn_LM_shotnoise.calc_fourier_transform();

			TwoPointStatistics<ParticleCatalogue> stats(params);
			std::complex<double> shotnoise = double(particles_data.nparticles);

			stats.calc_corr_func_for_3pt_func(
				dn_LM_shotnoise, shotnoise_00,
				rbin,
				shotnoise,
				params.ell1, m1_,
				ylm_a, ylm_b
			);

			for (int i = 0; i < params.num_rbin; i++) {
				if (params.form == "diag") {
					shotnoise_save[i] += coupling * stats.xi[i];
				} else if (params.form == "full") {
					if (i == params.ith_rbin) {
						shotnoise_save[i] += coupling * stats.xi[i];
					} else {
						shotnoise_save[i] += 0.;
					}
				}
			}

			durationInSec = double(clock() - timeStart);
			if (thisTask == 0) {
				printf(
					"[Status] :: Computed for wavenumber and orders "
					"m1 = %d, m2 = %d, M = %d "
					"(... %.3f seconds elapsed in total).\n",
					m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
				);
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
					* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	if (thisTask == 0) {
		printf("[Status] :: Measuring three-point function.\n");
	}

	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_fluctuation_in_box(particles_data, params);
	dn_00.calc_fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);

	std::complex<double>* zeta_save = new std::complex<double>[params.num_rbin];
	for (int i = 0; i < params.num_rbin; i++) {
		zeta_save[i] = 0.;
	}

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			int M_ = 0;

			double coupling = double(2*params.ELL + 1)
				* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
				* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
			if (fabs(coupling) < 1.e-10) {
				continue;
			}

			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell1, m1_, params, ylm_a
			);
			ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
				params.ell2, m2_, params, ylm_b
			);

			DensityField<ParticleCatalogue> dn_LM(params);
			dn_LM.calc_fluctuation_in_box(particles_data, params);
			dn_LM.calc_fourier_transform();
			dn_LM.calc_assignment_compensation();
			dn_LM.calc_inverse_fourier_transform();

			DensityField<ParticleCatalogue> dn_tilde1(params);
			double rmag_a;
			if (params.form == "full") {
				/// Compute ``\delta\tilde{n}_1``.
				rmag_a = rbin[params.ith_rbin];
				dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
					dn_00, rmag_a, ylm_a, sj1
				);
			}

			for (int i_rbin = 0; i_rbin < params.num_rbin; i_rbin++) {
				double rmag_b = rbin[i_rbin];

				/// Compute ``\delta\tilde{n}_1``.
				if (params.form == "diag") {
					rmag_a = rmag_b;
					dn_tilde1.calc_inverse_fourier_transform_for_3pt_func(
						dn_00, rmag_a, ylm_a, sj1
					);
				}

				/// Compute ``\delta\tilde{n}_2``.
				DensityField<ParticleCatalogue> dn_tilde2(params);
				dn_tilde2.calc_inverse_fourier_transform_for_3pt_func(
					dn_00, rmag_b, ylm_b, sj2
				);

				double factor = params.volume / double(params.nmesh_tot);
				std::complex<double> I_(0., 1.);
				std::complex<double> zeta_sum = 0.;
				for (int i = 0; i < params.nmesh_tot; i++) {
					std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
					std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
					std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
					zeta_sum += factor * pow(I_, params.ell1+params.ell2) * f1 * f2 * f3;
				}

				zeta_save[i_rbin] += coupling * zeta_sum;

				double durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for separation and order "
						"r2 = %.3f, m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						rmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = params.volume
		/ double(particles_data.nparticles) / double(particles_data.nparticles);
	norm *= params.volume / double(particles_data.nparticles);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/zeta%d%d%d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[i], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/zeta%d%d%d_%02d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_rbin
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_rbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				rbin[params.ith_rbin], rbin[i],
				norm * (zeta_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] zeta_save;

	return 0;
}

/**
 * Calculate bispectrum from catalogues with respect to a choice of
 * line of sight.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params (Reference to) the input parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param los Choice of line of sight.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispecChoiceOfLOS(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* kbin,
	int los,
	double survey_vol_norm
) {
	if (thisTask == 0) {
		printf(
			"[Status] :: Measuring bispectrum for the choice of line of sight.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (thisTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for bispectrum measurements. "
				"Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::complex<double>* shotnoise_save =
		new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		shotnoise_save[i] = 0.;
	}

	/// Compute shot noise terms.
	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_shotnoise(params);
				if (los == 0) {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_shotnoise.calc_fourier_transform();

				DensityField<ParticleCatalogue> N_shotnoise(params);
					// NOTE: standard naming convention overriden
				if (los == 0) {
					N_shotnoise.calc_ylm_weighted_fields_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				} else {
					N_shotnoise.calc_ylm_weighted_fields_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				}
				N_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);

				if (params.ell1 == 0 && params.ell2 == 0) {
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[i] += coupling * shotnoise;
					}
				}

				if (params.ell2 == 0) {
					stats.calc_power_spec(
						dn_shotnoise, N_shotnoise, kbin, shotnoise, params.ell1, m1_
					);
					if (params.form == "diag") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[i] += coupling * stats.pk[i];
						}
					} else if (params.form == "full") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[i] += coupling * stats.pk[params.ith_kbin];
						}
					}
				}

				durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC);
				}
			}
		}
	}

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_shotnoise(params);
				if (los == 1) {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_shotnoise.calc_fourier_transform();

				DensityField<ParticleCatalogue> N_shotnoise(params);
				if (los == 1) {
					N_shotnoise.calc_ylm_weighted_fields_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				} else {
					N_shotnoise.calc_ylm_weighted_fields_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				}
				N_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);

				if (params.ell1 == 0) {
					stats.calc_power_spec(
						dn_shotnoise, N_shotnoise, kbin, shotnoise, params.ell2, m2_
					);
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[i] += coupling * stats.pk[i];
					}
				}

				durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC);
				}
			}
		}
	}

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double> [params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double> [params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					*	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> N_shotnoise(params);
				if (los == 2) {
					N_shotnoise.calc_ylm_weighted_fields_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				} else {
					N_shotnoise.calc_ylm_weighted_fields_for_bispec_shotnoise(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				}
				N_shotnoise.calc_fourier_transform();

				DensityField<ParticleCatalogue> dn_shotnoise(params);
				if (los == 2) {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn_shotnoise.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);

				fftw_complex* xi = fftw_alloc_complex(params.nmesh_tot);
				bytes += sizeof(fftw_complex)
					* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
				for (int i = 0; i < params.nmesh_tot; i++) {
					xi[i][0] = 0.;
					xi[i][1] = 0.;
				}

				stats.calc_shotnoise_for_bispec_ijk(
					dn_shotnoise, N_shotnoise, shotnoise,
					params.ELL, M_,
					xi
				);

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];
					double kmag_a;
					if (params.form == "diag") {
						kmag_a = kmag_b;
					} else if (params.form == "full") {
						kmag_a = kbin[params.ith_kbin];
					}

					double dr[3];
					dr[0] = params.boxsize[0] / double(params.nmesh[0]);
					dr[1] = params.boxsize[1] / double(params.nmesh[1]);
					dr[2] = params.boxsize[2] / double(params.nmesh[2]);

					double rvec[3];
					std::complex<double> shotnoise_sum = 0.;
					for (int i = 0; i < params.nmesh[0]; i++) {
						for (int j = 0; j < params.nmesh[1]; j++) {
							for (int k = 0; k < params.nmesh[2]; k++) {
								long long coord_flat =
									(i * params.nmesh[1] + j) * params.nmesh[2] + k;
								rvec[0] = (i < params.nmesh[0]/2) ?
									i * dr[0] : (i - params.nmesh[0]) * dr[0];
								rvec[1] = (j < params.nmesh[1]/2) ?
									j * dr[1] : (j - params.nmesh[1]) * dr[1];
								rvec[2] = (k < params.nmesh[2]/2) ?
									k * dr[2] : (k - params.nmesh[2]) * dr[2];
								double rmag = sqrt(
									rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
								);

								std::complex<double> ff(xi[coord_flat][0], xi[coord_flat][1]);  // ???
								double j1 = sj1.eval(kmag_a * rmag);
								double j2 = sj2.eval(kmag_b * rmag);

								shotnoise_sum += (j1 * j2 * ff * ylm_a[i] * ylm_b[i]);
							}
						}
					}

					std::complex<double> I_(0., 1.);
					double factor = params.volume / double(params.nmesh_tot);
					shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

					shotnoise_save[i_kbin] += coupling * shotnoise_sum;

					durationInSec = double(clock() - timeStart);
					if (thisTask == 0) {
						printf(
							"[Status] :: Computed for wavenumber and orders "
							"k2 = %.3f, m1 = %d, m2 = %d, M = %d "
							"(... %.3f seconds elapsed in total).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
						);
					}
				}

				fftw_free(xi); xi = NULL;
				bytes -= sizeof(fftw_complex)
					* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	if (thisTask == 0) {
		printf("[Status] :: Measuring bispectrum.\n");
	}

	std::complex<double>* bk_save = new std::complex<double>[params.num_kbin];
	for (int i = 0; i < params.num_kbin; i++) {
		bk_save[i] = 0.;
	}

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				  *	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn1(params);
				if (los == 0) {
					dn1.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
					  los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn1.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn1.calc_fourier_transform();

				DensityField<ParticleCatalogue> dn2(params);
				if (los == 1) {
					dn2.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn2.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn2.calc_fourier_transform();

				DensityField<ParticleCatalogue> dn3(params);
				if (los == 2) {
					dn3.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						params.ELL, M_
					);
				} else {
					dn3.calc_ylm_weighted_fluctuation(
						particles_data, particles_rand,
						los_data, los_rand,
						alpha,
						0, 0
					);
				}
				dn3.calc_fourier_transform();
				dn3.calc_assignment_compensation();
				dn3.calc_inverse_fourier_transform();

				DensityField<ParticleCatalogue> dn_tilde1(params);
				double kmag_a;
				double dk = kbin[1] - kbin[0];
				if (params.form == "full") {
					kmag_a = kbin[params.ith_kbin];
					dn_tilde1.calc_inverse_fourier_transform_for_bispec(
						dn1, kmag_a, dk, ylm_a
					);
				}

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];

					if (params.form == "diag") {
						kmag_a = kmag_b;
						dn_tilde1.calc_inverse_fourier_transform_for_bispec(
							dn1, kmag_a, dk, ylm_a
						);
					}

					/// Compute ``\delta\tilde{n}_2``.
					DensityField<ParticleCatalogue> dn_tilde2(params);
					dn_tilde2.calc_inverse_fourier_transform_for_bispec(
						dn2, kmag_b, dk, ylm_b
					);

					double factor = params.volume / double(params.nmesh_tot);
					std::complex<double> bk_sum = 0.;
					for (int i = 0; i < params.nmesh_tot; i++) {
						std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
						std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
						std::complex<double> f3(dn3[i][0], dn3[i][1]);  // ???
						bk_sum += factor * f1 * f2 * f3;
					}

					bk_save[i_kbin] += coupling * bk_sum;

					double durationInSec = double(clock() - timeStart);
					if (thisTask == 0) {
						printf(
							"[Status] :: Computed for wavenumber and order "
							"k2 = %.3f, m1 = %d, m2 = %d, M = %d "
							"(... %.3f seconds elapsed in total).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
						);
					}
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = ParticleCatalogue::calc_norm_for_bispec(
		particles_data, survey_vol_norm
	);

	FILE* saved_file_ptr;
	char buf[1024];
	if (params.form == "diag") {
		sprintf(
			buf, "%s/bk%d%d%d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[i], kbin[i],
				norm * (bk_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	} else if (params.form == "full") {
		sprintf(
			buf, "%s/bk%d%d%d_%02d",
			params.output_dir.c_str(),
			params.ell1, params.ell2, params.ELL,
			params.ith_kbin
		);
		saved_file_ptr = fopen(buf, "w");
		for (int i = 0; i < params.num_kbin; i++) {
			fprintf(
				saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
				kbin[params.ith_kbin], kbin[i],
				norm * (bk_save[i].real() - shotnoise_save[i].real()),
				norm * shotnoise_save[i].real()
			);
		}
	}
	fclose(saved_file_ptr);

	delete[] shotnoise_save;
	delete[] bk_save;

	return 0;
}

/**
 * Calculate bispectrum from catalogues for modes with individual
 * orders @f$M@f$.
 *
 * @param particles_data (Reference to) the data-source particle container.
 * @param particles_rand (Reference to) the random-source particle container.
 * @param los_data Data-source particle lines of sight.
 * @param los_rand Random-source particle lines of sight.
 * @param params (Reference to) the input parameter set.
 * @param alpha Alpha ratio.
 * @param kbin Wavenumber bins.
 * @param survey_vol_norm Survey volume normalisation constant.
 * @returns Exit status.
 */
int calc_bispec_for_M_mode(
	ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
	LineOfSight* los_data, LineOfSight* los_rand,
	ParameterSet& params,
	double alpha,
	double* kbin,
	double survey_vol_norm
) {
	if (thisTask == 0) {
		printf(
			"[Status] :: Measuring bispectrum for individual modes with order `m`.\n"
		);
	}

	if (
		fabs(wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)) < 1.e-10
	) {
		if (thisTask == 0) {
			printf(
				"[Error] :: Disallowed multipole degree combination "
				"for bispectrum measurements. "
				"Please ensure wigner_3j(ell1, ell2, ELL, 0, 0, 0) != 0.\n"
			);
		}
		exit(1);
	}

	/// Initialise output shot noise terms.
	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computing shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	std::vector<
		std::vector<
			std::complex<double>
		>
	> shotnoise_save;

	shotnoise_save.resize(2*params.ELL + 1);
	for (int i = 0; i < params.num_kbin; i++) {
		shotnoise_save[i].resize(params.num_kbin);
	}

	for (int i = 0; i < params.num_kbin; i++) {
		for (int m = 0; m < 2*params.ELL + 1; m++) {
			shotnoise_save[m][i] = 0.;
		}
	}

	/// Compute shot noise terms.
	DensityField<ParticleCatalogue> dn_00_shotnoise(params);
	dn_00_shotnoise.calc_ylm_weighted_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00_shotnoise.calc_fourier_transform();

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> shotnoise_LM(params);
				shotnoise_LM.calc_ylm_weighted_fields_for_bispec_shotnoise(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				shotnoise_LM.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);

				if (params.ell1 == 0 && params.ell2 == 0) {
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[M_ + params.ELL][i] += shotnoise;
					}
				}

				if (params.ell2 == 0) {
					stats.calc_power_spec(
						dn_00_shotnoise, shotnoise_LM, kbin, shotnoise, params.ell1, m1_
					);
					if (params.form == "diag") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[M_+params.ELL][i] += stats.pk[i];
						}
					} else if (params.form == "full") {
						for (int i = 0; i < params.num_kbin; i++) {
							shotnoise_save[M_+params.ELL][i] += stats.pk[params.ith_kbin];
						}
					}
				}

				if ((params.ell1 == 0)) {
					stats.calc_power_spec(
						dn_00_shotnoise, shotnoise_LM, kbin, shotnoise, params.ell2, m2_
					);
					for (int i = 0; i < params.num_kbin; i++) {
						shotnoise_save[M_+params.ELL][i] += stats.pk[i];
					}
				}

				durationInSec = double(clock() - timeStart);
				if (thisTask == 0) {
					printf(
						"[Status] :: Computed for orders m1 = %d, m2 = %d, M = %d "
						"(... %.3f seconds elapsed in total).\n",
						m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
					);
				}
			}
		}
	}

	dn_00_shotnoise.finalise_density_field();

	DensityField<ParticleCatalogue> shotnoise_00(params);
	shotnoise_00.calc_ylm_weighted_fields_for_bispec_shotnoise(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	shotnoise_00.calc_fourier_transform();

	SphericalBesselCalculator sj1(params.ell1);
	SphericalBesselCalculator sj2(params.ell2);

	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
				  *	wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_config_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				};

				DensityField<ParticleCatalogue> dn_LM_shotnoise(params);
				dn_LM_shotnoise.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM_shotnoise.calc_fourier_transform();

				TwoPointStatistics<ParticleCatalogue> stats(params);
				std::complex<double> shotnoise = stats.calc_shotnoise_for_bispec(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);

				fftw_complex* xi = fftw_alloc_complex(params.nmesh_tot);
				bytes += sizeof(fftw_complex) *
					double(params.nmesh_tot) / 1024. / 1024. / 1024.;
				for (int i = 0; i < params.nmesh_tot; i++) {
					xi[i][0] = 0.;
					xi[i][1] = 0.;
				}

				stats.calc_shotnoise_for_bispec_ijk(
					dn_LM_shotnoise, shotnoise_00, shotnoise, params.ELL, M_, xi
				);

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];
					double kmag_a;
					if (params.form == "diag") {
						kmag_a = kmag_b;
					} else if (params.form == "full") {
						kmag_a = kbin[params.ith_kbin];
					}

					double dr[3];
					dr[0] = params.boxsize[0] / double(params.nmesh[0]);
					dr[1] = params.boxsize[1] / double(params.nmesh[1]);
					dr[2] = params.boxsize[2] / double(params.nmesh[2]);

					double rvec[3];
					std::complex<double> shotnoise_sum = 0.;
					for (int i = 0; i < params.nmesh[0]; i++) {
						for (int j = 0; j < params.nmesh[1]; j++) {
							for (int k = 0; k < params.nmesh[2]; k++) {
								long long coord_flat =
									(i * params.nmesh[1] + j) * params.nmesh[2] + k;
								rvec[0] = (i < params.nmesh[0]/2) ?
									i * dr[0] : (i - params.nmesh[0]) * dr[0];
								rvec[1] = (j < params.nmesh[1]/2) ?
									j * dr[1] : (j - params.nmesh[1]) * dr[1];
								rvec[2] = (k < params.nmesh[2]/2) ?
									k * dr[2] : (k - params.nmesh[2]) * dr[2];
								double rmag = sqrt(
									rvec[0] * rvec[0] + rvec[1] * rvec[1] + rvec[2] * rvec[2]
								);

								std::complex<double> ff(xi[coord_flat][0], xi[coord_flat][1]);  // ???
								double j1 = sj1.eval(kmag_a * rmag);
								double j2 = sj2.eval(kmag_b * rmag);

								shotnoise_sum += (j1 * j2 * ff * ylm_a[i] * ylm_b[i]);
							}
						}
					}

					std::complex<double> I_(0., 1.);
					double factor = params.volume / double(params.nmesh_tot);
					shotnoise_sum *= factor * pow(I_, params.ell1 + params.ell2);

					shotnoise_save[M_ + params.ELL][i_kbin] += shotnoise_sum;

					durationInSec = double(clock() - timeStart);
					if (thisTask == 0) {
						printf(
							"[Status] :: Computed for wavenumber and orders "
							"k2 = %.3f, m1 = %d, m2 = %d, M = %d "
							"(... %.3f seconds elapsed in total).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
						);
					}
				}

				fftw_free(xi); xi = NULL;
				bytes -= sizeof(fftw_complex)
					* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	shotnoise_00.finalise_density_field();

	durationInSec = double(clock() - timeStart);
	if (thisTask == 0) {
		printf(
			"[Status] :: Computed shot noise terms "
			"(... %.3f seconds elapsed in total).\n",
			durationInSec / CLOCKS_PER_SEC
		);
	}

	/// Initialise output bispectrum.
	if (thisTask == 0) {
		printf("[Status] :: Measuring bispectrum.\n");
	}

	DensityField<ParticleCatalogue> dn_00(params);
	dn_00.calc_ylm_weighted_fluctuation(
		particles_data, particles_rand,
		los_data, los_rand,
		alpha,
		0, 0
	);
	dn_00.calc_fourier_transform();

	std::vector<
		std::vector<
			std::complex<double>
		>
	> bk_save;

	bk_save.resize(2*params.ELL + 1);
	for (int i = 0; i < params.num_kbin; i++) {
		bk_save[i].resize(params.num_kbin);
	}

	for (int i = 0; i < params.num_kbin; i++) {
		for (int m = 0; m < 2*params.ELL + 1; m++) {
			bk_save[m][i] = 0.;
		}
	}

	/// Compute bispectrum.
	for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
		for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
			std::complex<double>* ylm_a = new std::complex<double>[params.nmesh_tot];
			std::complex<double>* ylm_b = new std::complex<double>[params.nmesh_tot];
			bytes += 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;

			std::string flag_nontrivial = "FALSE";
			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) > 1.e-10) {
					flag_nontrivial = "TRUE";
				}
			}

			if (flag_nontrivial == "TRUE") {
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell1, m1_, params, ylm_a
				);
				ToolCollection::store_reduced_spherical_harmonic_in_fourier_space(
					params.ell2, m2_, params, ylm_b
				);
			}

			for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
				double coupling = double(2*params.ELL + 1)
					* double(2*params.ell1 + 1) * double(2*params.ell2 + 1)
					* wigner_3j(params.ell1, params.ell2, params.ELL, 0, 0, 0)
					* wigner_3j(params.ell1, params.ell2, params.ELL, m1_, m2_, M_);
				if (fabs(coupling) < 1.e-10) {
					continue;
				}

				DensityField<ParticleCatalogue> dn_LM(params);
				dn_LM.calc_ylm_weighted_fluctuation(
					particles_data, particles_rand,
					los_data, los_rand,
					alpha,
					params.ELL, M_
				);
				dn_LM.calc_fourier_transform();
				dn_LM.calc_assignment_compensation();
				dn_LM.calc_inverse_fourier_transform();

				DensityField<ParticleCatalogue> dn_tilde1(params);
				double kmag_a;
				double dk = kbin[1] - kbin[0];
				if (params.form == "full") {
					/// Compute ``\delta\tilde{n}_1``.
					kmag_a = kbin[params.ith_kbin];
					dn_tilde1.calc_inverse_fourier_transform_for_bispec(dn_00, kmag_a, dk, ylm_a);
				}

				for (int i_kbin = 0; i_kbin < params.num_kbin; i_kbin++) {
					double kmag_b = kbin[i_kbin];

					/// Compute ``\delta\tilde{n}_1``.
					if (params.form == "diag") {
						kmag_a = kmag_b;
						dn_tilde1.calc_inverse_fourier_transform_for_bispec(
							dn_00, kmag_a, dk, ylm_a
						);
					}

					/// Compute ``\delta\tilde{n}_2``.
					DensityField<ParticleCatalogue> dn_tilde2(params);
					dn_tilde2.calc_inverse_fourier_transform_for_bispec(
						dn_00, kmag_b, dk, ylm_b
					);

					double factor = params.volume / double(params.nmesh_tot);
					std::complex<double> bk_sum = 0.;
					for (int i = 0; i < params.nmesh_tot; i++) {
						std::complex<double> f1(dn_tilde1[i][0], dn_tilde1[i][1]);  // ???
						std::complex<double> f2(dn_tilde2[i][0], dn_tilde2[i][1]);  // ???
						std::complex<double> f3(dn_LM[i][0], dn_LM[i][1]);  // ???
						bk_sum += factor * f1 * f2 * f3;
					}

					bk_save[M_ + params.ELL][i_kbin] += bk_sum;

					double durationInSec = double(clock() - timeStart);
					if (thisTask == 0) {
						printf(
							"[Status] :: Computed for wavenumber and order "
							"k2 = %.3f, m1 = %d, m2 = %d, M = %d "
							"(... %.3f seconds elapsed in total).\n",
							kmag_b, m1_, m2_, M_, durationInSec / CLOCKS_PER_SEC
						);
					}
				}
			}

			delete[] ylm_a; ylm_a = NULL;
			delete[] ylm_b; ylm_b = NULL;
			bytes -= 2 * sizeof(std::complex<double>)
				* double(params.nmesh_tot) / 1024. / 1024. / 1024.;
		}
	}

	/// Normalise and then save the output.
	double norm = ParticleCatalogue::calc_norm_for_bispec(
		particles_data, survey_vol_norm
	);

	for (int M_ = 0; M_<2*params.ELL + 1; M_++) {
		FILE* saved_file_ptr;
		char buf[1024];
		if (params.form == "diag") {
			sprintf(
				buf, "%s/bk%d%d%d_M%d",
				params.output_dir.c_str(),
				params.ell1, params.ell2,
				params.ELL, M_
			);
			saved_file_ptr = fopen(buf, "w");
			for (int i = 0; i < params.num_kbin; i++) {
				fprintf(
					saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
					kbin[i], kbin[i],
					norm * (bk_save[M_][i].real() - shotnoise_save[M_][i].real()),
					norm * shotnoise_save[M_][i].real()
				);
			}
		} else if (params.form == "full") {
			sprintf(
				buf, "%s/bk%d%d%d_M%d_%02d",
				params.output_dir.c_str(),
				params.ell1, params.ell2,
				params.ELL, M_,
				params.ith_kbin
			);
			saved_file_ptr = fopen(buf, "w");
			for (int i = 0; i < params.num_kbin; i++) {
				fprintf(
					saved_file_ptr, "%.5f \t %.5f \t %.7e \t %.7e\n",
					kbin[params.ith_kbin], kbin[i],
					norm * (bk_save[M_][i].real() - shotnoise_save[M_][i].real()),
					norm * shotnoise_save[M_][i].real()
				);
			}
		}
		fclose(saved_file_ptr);
	}

	return 0;
}

#endif
