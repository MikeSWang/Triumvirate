// Copyright (C) [GPLv3 Licence]
//
// This file is part of the Triumvirate program. See the COPYRIGHT
// and LICENCE files at the top-level directory of this distribution
// for details of copyright and licensing.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <https://www.gnu.org/licenses/>.

/**
 * @file threept.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "threept.hpp"

namespace trvs = trv::sys;
namespace trvm = trv::maths;

namespace trv {

// ***********************************************************************
// Coupling coefficients
// ***********************************************************************

double calc_coupling_coeff_3pt(
  int ell1, int ell2, int ELL, int m1, int m2, int M
) {
  return double(2*ell1 + 1) * double(2*ell2 + 1) * double(2*ELL + 1)
    * trvm::wigner_3j(ell1, ell2, ELL, 0, 0, 0)
    * trvm::wigner_3j(ell1, ell2, ELL, m1, m2, M);
}

void validate_multipole_coupling(trv::ParameterSet& params) {
  double coupling_ = trvm::wigner_3j(
    params.ell1, params.ell2, params.ELL, 0, 0, 0
  );
  if (std::fabs(coupling_) < trvm::eps_coupling) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Specified three-point correlator multipole "
        "vanishes identically owing to zero-valued Wigner 3-j symbol."
      );
      throw trvs::InvalidParameterError(
        "Specified three-point correlator multipole "
        "vanishes identically owing to zero-valued Wigner 3-j symbol.\n"
      );
    }
  }
}


// ***********************************************************************
// Normalisation
// ***********************************************************************

double calc_bispec_normalisation_from_particles(
  ParticleCatalogue& particles, double alpha
) {
  if (particles.pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
      throw trvs::InvalidDataError("Particle data are uninitialised.\n");
    }
  }

  double norm = 0.;  // I₃

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:norm)
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    norm += particles[pid].ws
      * std::pow(particles[pid].nz, 2) * std::pow(particles[pid].wc, 3);
  }

  if (norm == 0.) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Particle 'nz' values appear to be all zeros. "
        "Check the input catalogue contains valid 'nz' field."
      );
      throw trvs::InvalidDataError(
        "Particle 'nz' values appear to be all zeros. "
        "Check the input catalogue contains valid 'nz' field.\n"
      );
    }
  }

  double norm_factor = 1. / (alpha * norm);  // 1/I₃

  return norm_factor;
}

double calc_bispec_normalisation_from_mesh(
  ParticleCatalogue& particles, trv::ParameterSet& params, double alpha
) {
  MeshField catalogue_mesh(params);

  double norm_factor =
    catalogue_mesh.calc_grid_based_powlaw_norm(particles, 3);

  catalogue_mesh.finalise_density_field();  // likely redundant but safe

  norm_factor /= std::pow(alpha, 3);

  return norm_factor;
}


// ***********************************************************************
// Shot noise
// ***********************************************************************

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_bispec(
  ParticleCatalogue& particles_data, ParticleCatalogue& particles_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  double alpha, int ell, int m
) {
  double sn_data_real = 0., sn_data_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:sn_data_real, sn_data_imag)
#endif
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    double los_[3] = {
      los_data[pid].pos[0], los_data[pid].pos[1], los_data[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    std::complex<double> sn_part = ylm * std::pow(particles_data[pid].w, 3);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_data_real += sn_part_real;
    sn_data_imag += sn_part_imag;
  }

  std::complex<double> sn_data(sn_data_real, sn_data_imag);

  double sn_rand_real = 0., sn_rand_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:sn_rand_real, sn_rand_imag)
#endif
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    double los_[3] = {
      los_rand[pid].pos[0], los_rand[pid].pos[1], los_rand[pid].pos[2]
    };

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    std::complex<double> sn_part = ylm * std::pow(particles_rand[pid].w, 3);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_rand_real += sn_part_real;
    sn_rand_imag += sn_part_imag;
  }

  std::complex<double> sn_rand(sn_rand_real, sn_rand_imag);

  return sn_data + std::pow(alpha, 3) * sn_rand;
}

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_bispec(
  ParticleCatalogue& particles, LineOfSight* los,
  double alpha, int ell, int m
) {
  double sn_real = 0., sn_imag = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for reduction(+:sn_real, sn_imag)
#endif
  for (int pid = 0; pid < particles.ntotal; pid++) {
    double los_[3] = {los[pid].pos[0], los[pid].pos[1], los[pid].pos[2]};

    std::complex<double> ylm = trvm::SphericalHarmonicCalculator::
      calc_reduced_spherical_harmonic(ell, m, los_);

    std::complex<double> sn_part = ylm * std::pow(particles[pid].w, 3);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_real += sn_part_real;
    sn_imag += sn_part_imag;
  }

  std::complex<double> sn(sn_real, sn_imag);

  return std::pow(alpha, 3) * sn;
}


// ***********************************************************************
// Full statistics
// ***********************************************************************

// STYLE: Standard naming convention is not always followed for
// intermediary quantities in the functions below.

// Hereafter 'the Paper' refers to Sugiyama et al. (2019) [1803.02132].

trv::BispecMeasurements compute_bispec(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing bispectrum from paired survey-type catalogues..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int* nmodes_save = new int[kbinning.num_bins];
  double* k1_save = new double[kbinning.num_bins];
  double* k2_save = new double[kbinning.num_bins];
  std::complex<double>* bk_save = new std::complex<double>[kbinning.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[kbinning.num_bins];
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k1_save[ibin] = 0.;
    k2_save[ibin] = 0.;
    bk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params);  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  MeshField& dn_00_for_sn = dn_00;  // δn_00(k) (for shot noise)

  double vol_cell = dn_00.vol_cell;

  MeshField N_00(params);  // N_00(k)
  N_00.compute_ylm_wgtd_quad_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  // Compute bispectrum terms including shot noise.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      // Initialise reduced-spherical-harmonic weights on mesh grids.
      std::vector< std::complex<double> > ylm_k_a(params.nmesh);
      std::vector< std::complex<double> > ylm_k_b(params.nmesh);
      std::vector< std::complex<double> > ylm_r_a(params.nmesh);
      std::vector< std::complex<double> > ylm_r_b(params.nmesh);
      trvs::gbytesMem +=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
      trvs::update_maxmem();

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Raw bispectrum
        // ·······························································

        // Compute bispectrum components in eqs. (41) & (42) in the Paper.
        MeshField G_LM(params);  // G_LM
        G_LM.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        MeshField F_lm_a(params);  // F_lm_a
        MeshField F_lm_b(params);  // F_lm_b
        if (params.form == "full") {
          double k_lower = kbinning.bin_edges[params.idx_bin];
          double k_upper = kbinning.bin_edges[params.idx_bin + 1];

          double k_eff_a_;
          int nmodes_a_;

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
          );

          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            // nmodes1_save[ibin] = nmodes_a_ inferred from *_b_
            k1_save[ibin] = k_eff_a_;
          }
        }

        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          double k_lower = kbinning.bin_edges[ibin];
          double k_upper = kbinning.bin_edges[ibin + 1];

          double k_eff_b_;
          int nmodes_b_;

          F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_b, k_lower, k_upper, k_eff_b_, nmodes_b_
          );

          k2_save[ibin] = k_eff_b_;
          nmodes_save[ibin] = nmodes_b_;

          if (params.form == "diag") {
            double k_eff_a_;  // redundant as this is `k_eff_b_`
            int nmodes_a_;  // redundant as this is `nmodes_b_`

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_00, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
            );

            k1_save[ibin] = k_eff_a_;
          }

          std::complex<double> bk_component = 0.;  // B_{l₁ l₂ L}^{m₁ m₂ M}
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
            bk_component += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          bk_save[ibin] += coupling * vol_cell * bk_component;
        }

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eqs. (45) & (46) in the Paper.
        MeshField dn_LM_for_sn(params);  // δn_LM(k) (for shot noise)
        dn_LM_for_sn.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM_for_sn.fourier_transform();

        MeshField N_LM(params);  // N_LM(k)
        N_LM.compute_ylm_wgtd_quad_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        N_LM.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );  // \bar{S}_LM

        if (params.ell1 == 0 && params.ell2 == 0) {
          // When l₁ = l₂ = 0, the Wigner 3-j symbol enforces L = 0
          // and the pre-factors involving degrees and orders become 1.
          std::complex<double> S_ijk = coupling * Sbar_LM;  // S|{i = j = k}
          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            sn_save[ibin] += S_ijk;
          }
        }

        if (params.ell2 == 0) {
          // When l₂ = 0, the Wigner 3-j symbol enforces L = l₁.
          FieldStats stats_sn(params);  // S|{i ≠ j = k}
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_sn, N_LM, Sbar_LM, params.ell1, m1_, kbinning
          );
          if (params.form == "diag") {
            for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          } else
          if (params.form == "full") {
            for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[params.idx_bin] - stats_sn.sn[params.idx_bin]
              );
            }
          }
        }

        if (params.ell1 == 0) {
          // When l₁ = 0, the Wigner 3-j symbol enforces L = l₂.
          FieldStats stats_sn(params);  // S|{j ≠ i = k}
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_00_for_sn, N_LM, Sbar_LM, params.ell2, m2_, kbinning
          );
          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            sn_save[ibin] +=
              coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
          }
        }

        FieldStats stats_sn(params);
        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          double k_a = k1_save[ibin];
          double k_b = k2_save[ibin];

          std::complex<double> S_ij_k = parity
            * stats_sn.compute_uncoupled_shotnoise_for_bispec_per_bin(
              dn_LM_for_sn, N_00, ylm_r_a, ylm_r_b, sj_a, sj_b,
              Sbar_LM, k_a, k_b
            );  // S|{i = j ≠ k}

          sn_save[ibin] += coupling * S_ij_k;
        }

        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Bispectrum term at orders (m1, m2, M) = (%d, %d, %d) computed.",
            m1_, m2_, M_
          );
        }
      }

      trvs::gbytesMem -=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
    }
  }

  dn_00.finalise_density_field();  // ~dn_00 (likely redundant but safe)
  N_00.finalise_density_field();  // ~N_00 (likely redundant but safe)

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    if (params.form == "diag") {
      bispec_out.k1_bin.push_back(kbinning.bin_centres[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.k1_bin.push_back(kbinning.bin_centres[params.idx_bin]);
    }
    bispec_out.k1_eff.push_back(k1_save[ibin]);
    bispec_out.k2_bin.push_back(kbinning.bin_centres[ibin]);
    bispec_out.k2_eff.push_back(k2_save[ibin]);
    bispec_out.nmodes.push_back(nmodes_save[ibin]);
    bispec_out.bk_raw.push_back(norm_factor * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] nmodes_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed bispectrum from paired survey-type catalogues."
    );
  }

  return bispec_out;
}

trv::ThreePCFMeasurements compute_3pcf(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing three-point correlation function "
      "from paired survey-type catalogues..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;

  double parity = std::pow(-1, params.ell1 + params.ell2);

  // Set up output.
  int* npairs_save = new int[rbinning.num_bins];
  double* r1_save = new double[rbinning.num_bins];
  double* r2_save = new double[rbinning.num_bins];
  std::complex<double>* zeta_save =
    new std::complex<double>[rbinning.num_bins];
  std::complex<double>* sn_save =
    new std::complex<double>[rbinning.num_bins];
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r1_save[ibin] = 0.;
    r2_save[ibin] = 0.;
    zeta_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params);  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  double vol_cell = dn_00.vol_cell;

  MeshField N_00(params);  // N_00(k)
  N_00.compute_ylm_wgtd_quad_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  // Compute 3PCF terms including shot noise.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for vanishing cases where all Wigner-3j symbols are zero.
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      // Initialise reduced-spherical-harmonic weights on mesh grids.
      std::vector< std::complex<double> > ylm_r_a(params.nmesh);
      std::vector< std::complex<double> > ylm_r_b(params.nmesh);
      std::vector< std::complex<double> > ylm_k_a(params.nmesh);
      std::vector< std::complex<double> > ylm_k_b(params.nmesh);
      trvs::gbytesMem +=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
      trvs::update_maxmem();

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eq. (51) in the Paper.
        MeshField dn_LM_for_sn(params);  // δn_LM(k) (for shot noise)
        dn_LM_for_sn.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        dn_LM_for_sn.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );  // \bar{S}_LM

        FieldStats stats_sn(params);  // S|{i = j ≠ k}
        stats_sn.compute_uncoupled_shotnoise_for_3pcf(
          dn_LM_for_sn, N_00, ylm_r_a, ylm_r_b, Sbar_LM, rbinning
        );

        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          if (params.form == "diag") {
            sn_save[ibin] += coupling * stats_sn.xi[ibin];
          } else
          if (params.form == "full") {
            // Enforce the Kronecker delta in eq. (51) in the Paper.
            if (ibin == params.idx_bin) {
              sn_save[ibin] += coupling * stats_sn.xi[ibin];
            }
            // else {
            //   sn_save[ibin] += 0.;
            // }
          }
        }

        // Only record the binned coordinates and counts once.
        if (M_ == 0 && m1_ == 0 && m2_ == 0) {
          for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
            npairs_save[ibin] = stats_sn.npairs[ibin];
            r2_save[ibin] = stats_sn.r[ibin];
          }
          if (params.form == "diag") {
            for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
              r1_save[ibin] = stats_sn.r[ibin];
            }
          } else
          if (params.form == "full") {
            for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
              r1_save[ibin] = stats_sn.r[params.idx_bin];
            }
          }
        }

        // ·······························································
        // Raw 3PCF
        // ·······························································

        // Compute 3PCF components in eqs. (42), (48) & (49) in the Paper.
        MeshField G_LM(params);  // G_LM
        G_LM.compute_ylm_wgtd_field(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        MeshField F_lm_a(params);  // F_lm_a
        MeshField F_lm_b(params);  // F_lm_b

        if (params.form == "full") {
          double r_a = r1_save[params.idx_bin];
          F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
            dn_00, ylm_k_a, sj_a, r_a
          );
        }

        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          double r_b = r2_save[ibin];

          F_lm_b.inv_fourier_transform_sjl_ylm_wgtd_field(
            dn_00, ylm_k_b, sj_b, r_b
          );

          if (params.form == "diag") {
            double r_a = r_b;
            F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
              dn_00, ylm_k_a, sj_a, r_a
            );
          }

          std::complex<double> zeta_component = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
            zeta_component += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          zeta_save[ibin] += parity * coupling * vol_cell * zeta_component;
        }

        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Three-point correlation function term at orders "
            "(m1, m2, M) = (%d, %d, %d) computed.",
            m1_, m2_, M_
          );
        }
      }

      trvs::gbytesMem -=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
    }
  }

  dn_00.finalise_density_field();  // ~dn_00 (likely redundant but safe)
  N_00.finalise_density_field();  // ~N_00 (likely redundant but safe)

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::ThreePCFMeasurements threepcf_out;
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    if (params.form == "diag") {
      threepcf_out.r1_bin.push_back(rbinning.bin_centres[ibin]);
    } else
    if (params.form == "full") {
      threepcf_out.r1_bin.push_back(rbinning.bin_centres[params.idx_bin]);
    }
    threepcf_out.r1_eff.push_back(r1_save[ibin]);
    threepcf_out.r2_bin.push_back(rbinning.bin_centres[ibin]);
    threepcf_out.r2_eff.push_back(r2_save[ibin]);
    threepcf_out.npairs.push_back(npairs_save[ibin]);
    threepcf_out.zeta_raw.push_back(norm_factor * zeta_save[ibin]);
    threepcf_out.zeta_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] npairs_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed 3-point correlation function "
      "from paired survey-type catalogues."
    );
  }

  return threepcf_out;
}

trv::BispecMeasurements compute_bispec_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing bispectrum from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int* nmodes_save = new int[kbinning.num_bins];
  double* k1_save = new double[kbinning.num_bins];
  double* k2_save = new double[kbinning.num_bins];
  std::complex<double>* bk_save = new std::complex<double>[kbinning.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[kbinning.num_bins];
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k1_save[ibin] = 0.;
    k2_save[ibin] = 0.;
    bk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params);  // δn_00(k)
  dn_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn_00.fourier_transform();

  MeshField& dn_00_for_sn = dn_00;  // δn_00(k) (for shot noise)
  MeshField& dn_L0_for_sn = dn_00;  // δn_L0(k) (for shot noise)

  double vol_cell = dn_00.vol_cell;

  // Under the global plane-parallel approximation, y_{LM} = δᴰ_{M0}
  // (L-invariant) for the line-of-sight spherical harmonic.
  MeshField N_L0(params);  // N_L0(k)
  N_L0.compute_unweighted_field(catalogue_data);
  N_L0.fourier_transform();

  MeshField& N_00 = N_L0;  // N_00(k)

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  // Compute bispectrum terms including shot noise.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Under the global plane-parallel approximation, δᴰ_{M0} enforces
      // M = 0 for any spherical-harmonic-weighted field fluctuations.
      int M_ = 0;

      // Calculate the coupling coefficient.
      double coupling = trv::calc_coupling_coeff_3pt(
        params.ell1, params.ell2, params.ELL, m1_, m2_, M_
      );  // Wigner 3-j's
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      // Initialise/reset spherical harmonic mesh grids.
      std::vector< std::complex<double> > ylm_k_a(params.nmesh);
      std::vector< std::complex<double> > ylm_k_b(params.nmesh);
      std::vector< std::complex<double> > ylm_r_a(params.nmesh);
      std::vector< std::complex<double> > ylm_r_b(params.nmesh);
      trvs::gbytesMem +=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
      trvs::update_maxmem();

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator
        ::store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );

      // ·································································
      // Raw bispectrum
      // ·································································

      MeshField G_00(params);  // G_00
      G_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
      G_00.fourier_transform();
      G_00.apply_assignment_compensation();
      G_00.inv_fourier_transform();

      MeshField F_lm_a(params);  // F_lm_a
      MeshField F_lm_b(params);  // F_lm_b
      if (params.form == "full") {
        double k_lower = kbinning.bin_edges[params.idx_bin];
        double k_upper = kbinning.bin_edges[params.idx_bin + 1];

        double k_eff_a_;
        int nmodes_a_;

        F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
          dn_00, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
        );

        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          // nmodes1_save[ibin] = nmodes_a_ inferred from *_b_
          k1_save[ibin] = k_eff_a_;
        }
      }

      for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
        double k_lower = kbinning.bin_edges[ibin];
        double k_upper = kbinning.bin_edges[ibin + 1];

        double k_eff_b_;
        int nmodes_b_;

        F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
          dn_00, ylm_k_b, k_lower, k_upper, k_eff_b_, nmodes_b_
        );

        k2_save[ibin] = k_eff_b_;
        nmodes_save[ibin] = nmodes_b_;

        if (params.form == "diag") {
          double k_eff_a_;  // redundant as this is `k_eff_b_`
          int nmodes_a_;  // redundant as this is `nmodes_b_`

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_00, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
          );

          k1_save[ibin] = k_eff_a_;
        }

        std::complex<double> bk_component = 0.;  // B_{l₁ l₂ L}^{m₁ m₂ M}
        for (int gid = 0; gid < params.nmesh; gid++) {
          std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
          std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
          std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
          bk_component += F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;
        }

        bk_save[ibin] += coupling * vol_cell * bk_component;
      }

      // ·································································
      // Shot noise
      // ·································································

      // Under the global plane-parallel approximation, y_{LM} = δᴰ_{M0}
      // (L-invariant) for the line-of-sight spherical harmonic.
      // Also note the field is unweighted from simulation sources.
      std::complex<double> Sbar_LM =
        double(catalogue_data.ntotal);  // \bar{S}_LM
      std::complex<double> Sbar_L0 = Sbar_LM;  // \bar{S}_L0

      if (params.ell1 == 0 && params.ell2 == 0) {
        std::complex<double> S_ijk = coupling * Sbar_LM;  // S|{i = j = k}
        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          sn_save[ibin] += coupling * S_ijk;
        }
      }

      if (params.ell2 == 0) {
        // When l₂ = 0, the Wigner 3-j symbol enforces L = l₁.
        FieldStats stats_sn(params);  // S|{i ≠ j = k}
        stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_sn, N_L0, Sbar_LM, params.ell1, m1_, kbinning
        );
        if (params.form == "diag") {
          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            sn_save[ibin] +=
              coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
          }
        } else
        if (params.form == "full") {
          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            sn_save[ibin] += coupling * (
              stats_sn.pk[params.idx_bin] - stats_sn.sn[params.idx_bin]
            );
          }
        }
      }

      if (params.ell1 == 0) {
        // When l₁ = 0, the Wigner 3-j symbol enforces L = l₂.
        FieldStats stats_sn(params);  // S|{j ≠ i = k}
        stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
          dn_00_for_sn, N_L0, Sbar_LM, params.ell2, m2_, kbinning
        );
        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          sn_save[ibin] += coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
        }
      }

      FieldStats stats_sn(params);
      for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
        double k_a = k1_save[ibin];
        double k_b = k2_save[ibin];

        std::complex<double> S_ij_k = parity *
          stats_sn.compute_uncoupled_shotnoise_for_bispec_per_bin(
            dn_L0_for_sn, N_00, ylm_r_a, ylm_r_b, sj_a, sj_b,
            Sbar_L0, k_a, k_b
          );  // S|{i = j ≠ k}

        sn_save[ibin] += coupling * S_ij_k;
      }

      if (trvs::currTask == 0) {
        trvs::logger.stat(
          "Bispectrum term at orders (m1, m2, M) = (%d, %d, 0) computed.",
          m1_, m2_
        );
      }

      trvs::gbytesMem -=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
    }
  }

  dn_00.finalise_density_field();  // ~dn_00 (likely redundant but safe)
  N_L0.finalise_density_field();  // ~N_L0 (likely redundant but safe)

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    if (params.form == "diag") {
      bispec_out.k1_bin.push_back(kbinning.bin_centres[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.k1_bin.push_back(kbinning.bin_centres[params.idx_bin]);
    }
    bispec_out.k1_eff.push_back(k1_save[ibin]);
    bispec_out.k2_bin.push_back(kbinning.bin_centres[ibin]);
    bispec_out.k2_eff.push_back(k2_save[ibin]);
    bispec_out.nmodes.push_back(nmodes_save[ibin]);
    bispec_out.bk_raw.push_back(norm_factor * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] nmodes_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed bispectrum from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  return bispec_out;
}

trv::ThreePCFMeasurements compute_3pcf_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing three-point correlation function "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double parity = std::pow(-1, params.ell1 + params.ell2);

  // Set up output.
  int* npairs_save = new int[rbinning.num_bins];
  double* r1_save = new double[rbinning.num_bins];
  double* r2_save = new double[rbinning.num_bins];
  std::complex<double>* zeta_save =
    new std::complex<double>[rbinning.num_bins];
  std::complex<double>* sn_save =
    new std::complex<double>[rbinning.num_bins];
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r1_save[ibin] = 0.;
    r2_save[ibin] = 0.;
    zeta_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params);  // δn_00(k)
  dn_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn_00.fourier_transform();

  MeshField& dn_L0_for_sn = dn_00;  // δn_L0(k)

  double vol_cell = dn_00.vol_cell;

  MeshField N_00(params);  // N_00(k)
  N_00.compute_unweighted_field(catalogue_data);
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  // Compute 3PCF terms including shot noise.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Under the global plane-parallel approximation, δᴰ_{M0} enforces
      // M = 0 for any spherical-harmonic-weighted field fluctuations.
      int M_ = 0;

      // Calculate the coupling coefficient.
      double coupling = trv::calc_coupling_coeff_3pt(
        params.ell1, params.ell2, params.ELL, m1_, m2_, M_
      );  // Wigner 3-j's
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      // Initialise/reset spherical harmonic mesh grids.
      std::vector< std::complex<double> > ylm_r_a(params.nmesh);
      std::vector< std::complex<double> > ylm_r_b(params.nmesh);
      std::vector< std::complex<double> > ylm_k_a(params.nmesh);
      std::vector< std::complex<double> > ylm_k_b(params.nmesh);
      trvs::gbytesMem +=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
      trvs::update_maxmem();

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );

      // ·································································
      // Shot noise
      // ·································································

      // Under the global plane-parallel approximation, y_{LM} = δᴰ_{M0}
      // (L-invariant) for the line-of-sight spherical harmonic.
      // Also note the field is unweighted from simulation sources.
      std::complex<double> Sbar_L0 =
        double(catalogue_data.ntotal);  // \bar{S}_L0

      FieldStats stats_sn(params);  // S|{i = j ≠ k}
      stats_sn.compute_uncoupled_shotnoise_for_3pcf(
        dn_L0_for_sn, N_00, ylm_r_a, ylm_r_b, Sbar_L0, rbinning
      );

      for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
        if (params.form == "diag") {
          sn_save[ibin] += coupling * stats_sn.xi[ibin];
        } else
        if (params.form == "full") {
            // Enforce the Kronecker delta in eq. (51) in the Paper.
          if (ibin == params.idx_bin) {
            sn_save[ibin] += coupling * stats_sn.xi[ibin];
          }
          // else {
          //   sn_save[ibin] += 0.;
          // }
        }
      }

      // Only record the binned coordinates and counts once.
      if (m1_ == 0 && m2_ == 0) {
        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          npairs_save[ibin] = stats_sn.npairs[ibin];
          r2_save[ibin] = stats_sn.r[ibin];
        }
        if (params.form == "diag") {
          for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
            r1_save[ibin] = stats_sn.r[ibin];
          }
        } else
        if (params.form == "full") {
          for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
            r1_save[ibin] = stats_sn.r[params.idx_bin];
          }
        }
      }

      // ·································································
      // Raw 3PCF
      // ·································································

      // Compute 3PCF components in eqs. (42), (48) & (49) in the Paper.
      MeshField G_00(params);  // G_00
      G_00.compute_unweighted_field_fluctuations_insitu(catalogue_data);
      G_00.fourier_transform();
      G_00.apply_assignment_compensation();
      G_00.inv_fourier_transform();

      MeshField F_lm_a(params);  // F_lm_a
      MeshField F_lm_b(params);  // F_lm_b

      if (params.form == "full") {
        double r_a = r1_save[params.idx_bin];
        F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
          dn_00, ylm_k_a, sj_a, r_a
        );
      }

      for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
        double r_b = r2_save[ibin];

        F_lm_b.inv_fourier_transform_sjl_ylm_wgtd_field(
          dn_00, ylm_k_b, sj_b, r_b
        );

        if (params.form == "diag") {
          double r_a = r_b;
          F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
            dn_00, ylm_k_a, sj_a, r_a
          );
        }

        std::complex<double> zeta_component = 0.;
        for (int gid = 0; gid < params.nmesh; gid++) {
          std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
          std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
          std::complex<double> G_00_gridpt(G_00[gid][0], G_00[gid][1]);
          zeta_component += F_lm_a_gridpt * F_lm_b_gridpt * G_00_gridpt;
        }

        zeta_save[ibin] += parity * coupling * vol_cell * zeta_component;
      }

      if (trvs::currTask == 0) {
        trvs::logger.stat(
          "Three-point correlation function term at orders "
          "(m1, m2, M) = (%d, %d, 0) computed.",
          m1_, m2_
        );
      }

      trvs::gbytesMem -=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
    }
  }

  dn_00.finalise_density_field();  // ~dn_00 (likely redundant but safe)
  N_00.finalise_density_field();  // ~N_00 (likely redundant but safe)

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::ThreePCFMeasurements threepcf_out;
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    if (params.form == "diag") {
      threepcf_out.r1_bin.push_back(rbinning.bin_centres[ibin]);
    } else
    if (params.form == "full") {
      threepcf_out.r1_bin.push_back(rbinning.bin_centres[params.idx_bin]);
    }
    threepcf_out.r1_eff.push_back(r1_save[ibin]);
    threepcf_out.r2_bin.push_back(rbinning.bin_centres[ibin]);
    threepcf_out.r2_eff.push_back(r2_save[ibin]);
    threepcf_out.npairs.push_back(npairs_save[ibin]);
    threepcf_out.zeta_raw.push_back(norm_factor * zeta_save[ibin]);
    threepcf_out.zeta_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] npairs_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed 3-point correlation function "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  return threepcf_out;
}

trv::ThreePCFWindowMeasurements compute_3pcf_window(
  ParticleCatalogue& catalogue_rand, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double alpha, double norm_factor, bool wide_angle
) {
  std::string msg_tag = wide_angle ? "wide-angle corrections " : "";

  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing three-point correlation function window %s"
      "from random catalogue...",
      msg_tag.c_str()
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double parity = std::pow(-1, params.ell1 + params.ell2);

  // Set up output.
  int* npairs_save = new int[rbinning.num_bins];
  double* r1_save = new double[rbinning.num_bins];
  double* r2_save = new double[rbinning.num_bins];
  std::complex<double>* zeta_save =
    new std::complex<double>[rbinning.num_bins];
  std::complex<double>* sn_save =
    new std::complex<double>[rbinning.num_bins];
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r1_save[ibin] = 0.;
    r2_save[ibin] = 0.;
    zeta_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField n_00(params);  // n_00(k)
  n_00.compute_ylm_wgtd_field(catalogue_rand, los_rand, alpha, 0, 0);
  n_00.fourier_transform();

  double vol_cell = n_00.vol_cell;

  MeshField N_00(params);  // N_00(k)
  N_00.compute_ylm_wgtd_quad_field(catalogue_rand, los_rand, alpha, 0, 0);
  N_00.fourier_transform();

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  // Compute 3PCF window terms including shot noise.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for vanishing cases where all Wigner-3j symbols are zero.
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      // Initialise reduced-spherical-harmonic weights on mesh grids.
      std::vector< std::complex<double> > ylm_r_a(params.nmesh);
      std::vector< std::complex<double> > ylm_r_b(params.nmesh);
      std::vector< std::complex<double> > ylm_k_a(params.nmesh);
      std::vector< std::complex<double> > ylm_k_b(params.nmesh);
      trvs::gbytesMem +=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
      trvs::update_maxmem();

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eq. (51) in the Paper.
        MeshField n_LM_for_sn(params);  // n_LM(k) (for shot noise)
        n_LM_for_sn.compute_ylm_wgtd_field(
          catalogue_rand, los_rand, alpha, params.ELL, M_
        );
        n_LM_for_sn.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_rand, los_rand, alpha, params.ELL, M_
        );  // \bar{S}_LM

        FieldStats stats_sn(params);  // S|{i = j ≠ k}
        stats_sn.compute_uncoupled_shotnoise_for_3pcf(
          n_LM_for_sn, N_00, ylm_r_a, ylm_r_b, Sbar_LM, rbinning
        );

        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          if (params.form == "diag") {
            sn_save[ibin] += coupling * stats_sn.xi[ibin];
          } else if (params.form == "full") {
            // Enforce the Kronecker delta in eq. (51) in the Paper.
            if (ibin == params.idx_bin) {
              sn_save[ibin] += coupling * stats_sn.xi[ibin];
            }
            // else {
            //   sn_save[ibin] += 0.;
            // }
          }
        }

        // Only record the binned coordinates and counts once.
        if (M_ == 0 && m1_ == 0 && m2_ == 0) {
          for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
            npairs_save[ibin] = stats_sn.npairs[ibin];
            r2_save[ibin] = stats_sn.r[ibin];
          }
          if (params.form == "diag") {
            for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
              r1_save[ibin] = stats_sn.r[ibin];
            }
          } else
          if (params.form == "full") {
            for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
              r1_save[ibin] = stats_sn.r[params.idx_bin];
            }
          }
        }

        // ·······························································
        // Raw 3PCF
        // ·······························································

        // Compute 3PCF components in eqs. (42), (48) & (49) in the Paper.
        MeshField G_LM(params);  // G_LM
        G_LM.compute_ylm_wgtd_field(
          catalogue_rand, los_rand, alpha, params.ELL, M_
        );
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        // Perform wide-angle corrections if required.
        if (wide_angle) {
          G_LM.apply_wide_angle_pow_law_kernel();
        }

        MeshField F_lm_a(params);  // F_lm_a
        MeshField F_lm_b(params);  // F_lm_b

        if (params.form == "full") {
          double r_a = r1_save[params.idx_bin];
          F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
            n_00, ylm_k_a, sj_a, r_a
          );
        }

        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          double r_b = r2_save[ibin];

          F_lm_b.inv_fourier_transform_sjl_ylm_wgtd_field(
            n_00, ylm_k_b, sj_b, r_b
          );

          if (params.form == "diag") {
            double r_a = r_b;
            F_lm_a.inv_fourier_transform_sjl_ylm_wgtd_field(
              n_00, ylm_k_a, sj_a, r_a
            );
          }

          std::complex<double> zeta_component = 0.;
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
            zeta_component += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          zeta_save[ibin] += parity * coupling * vol_cell * zeta_component;
        }

        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Three-point correlation function window term at orders "
            "(m1, m2, M) = (%d, %d, %d) computed.",
            m1_, m2_, M_
          );
        }
      }

      trvs::gbytesMem -=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
    }
  }

  n_00.finalise_density_field();  // ~n_00 (likely redundant but safe)
  N_00.finalise_density_field();  // ~N_00 (likely redundant but safe)

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::ThreePCFWindowMeasurements threepcfwin_out;
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    if (params.form == "diag") {
      threepcfwin_out.r1_bin.push_back(rbinning.bin_centres[ibin]);
    } else
    if (params.form == "full") {
      threepcfwin_out.r1_bin.push_back(rbinning.bin_centres[params.idx_bin]);
    }
    threepcfwin_out.r1_eff.push_back(r1_save[ibin]);
    threepcfwin_out.r2_bin.push_back(rbinning.bin_centres[ibin]);
    threepcfwin_out.r2_eff.push_back(r2_save[ibin]);
    threepcfwin_out.npairs.push_back(npairs_save[ibin]);
    threepcfwin_out.zeta_raw.push_back(norm_factor * zeta_save[ibin]);
    threepcfwin_out.zeta_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] npairs_save; delete[] r1_save; delete[] r2_save;
  delete[] zeta_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed three-point correlation function window %s"
      "from random catalogue.",
      msg_tag.c_str()
    );
  }

  return threepcfwin_out;
}

#ifdef TRV_USE_LEGACY_CODE
trv::BispecMeasurements compute_bispec_for_los_choice(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  int los_choice,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing bispectrum from paired survey-type catalogues "
      "for line-of-sight choice %d...",
      los_choice
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up/check input.
  validate_multipole_coupling(params);

  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;

  std::complex<double> parity = std::pow(trvm::M_I, params.ell1 + params.ell2);

  // Set up output.
  int* nmodes_save = new int[kbinning.num_bins];
  double* k1_save = new double[kbinning.num_bins];
  double* k2_save = new double[kbinning.num_bins];
  std::complex<double>* bk_save = new std::complex<double>[kbinning.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[kbinning.num_bins];
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k1_save[ibin] = 0.;
    k2_save[ibin] = 0.;
    bk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute common field quantities.
  MeshField dn_00(params);  // δn_00(r)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );

  double vol_cell = dn_00.vol_cell;
  double dr[3] = {dn_00.dr[0], dn_00.dr[1], dn_00.dr[2]};

  MeshField N_00(params);  // N_00(r)
  N_00.compute_ylm_wgtd_quad_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );

  trvm::SphericalBesselCalculator sj_a(params.ell1);  // j_l_a
  trvm::SphericalBesselCalculator sj_b(params.ell2);  // j_l_b

  // Compute bispectrum terms.
  for (int m1_ = - params.ell1; m1_ <= params.ell1; m1_++) {
    for (int m2_ = - params.ell2; m2_ <= params.ell2; m2_++) {
      // Check for if all Wigner-3j symbols are zero.
      std::string flag_vanishing = "true";
      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );
        if (std::fabs(coupling) > trvm::eps_coupling) {
          flag_vanishing = "false";
          break;
        }
      }
      if (flag_vanishing == "true") {continue;}

      // Initialise reduced-spherical-harmonic weights on mesh grids.
      std::vector< std::complex<double> > ylm_k_a(params.nmesh);
      std::vector< std::complex<double> > ylm_k_b(params.nmesh);
      std::vector< std::complex<double> > ylm_r_a(params.nmesh);
      std::vector< std::complex<double> > ylm_r_b(params.nmesh);
      trvs::gbytesMem +=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
      trvs::update_maxmem();

      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_k_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_fourier_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_k_b
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell1, m1_, params.boxsize, params.ngrid, ylm_r_a
        );
      trvm::SphericalHarmonicCalculator::
        store_reduced_spherical_harmonic_in_config_space(
          params.ell2, m2_, params.boxsize, params.ngrid, ylm_r_b
        );

      for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
        // Calculate the coupling coefficient.
        double coupling = trv::calc_coupling_coeff_3pt(
          params.ell1, params.ell2, params.ELL, m1_, m2_, M_
        );  // Wigner 3-j's
        if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

        // ·······························································
        // Raw bispectrum
        // ·······························································

        // Compute bispectrum components in eqs. (41) & (42) in the Paper.
        MeshField dn_LM_a(params);  // δn_LM_a
        if (los_choice == 0) {
          dn_LM_a.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_a = dn_00;
        }
        dn_LM_a.fourier_transform();

        MeshField dn_LM_b(params);  // δn_LM_b
        if (los_choice == 1) {
          dn_LM_b.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_b = dn_00;
        }
        dn_LM_b.fourier_transform();

        MeshField G_LM(params);  // G_LM
        if (los_choice == 2) {
          G_LM.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          G_LM = dn_00;
        }
        G_LM.fourier_transform();
        G_LM.apply_assignment_compensation();
        G_LM.inv_fourier_transform();

        MeshField F_lm_a(params);  // F_lm_a
        MeshField F_lm_b(params);  // F_lm_b
        if (params.form == "full") {
          double k_lower = kbinning.bin_edges[params.idx_bin];
          double k_upper = kbinning.bin_edges[params.idx_bin + 1];

          double k_eff_a_;
          int nmodes_a_;

          F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_LM_a, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
          );

          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            // nmodes1_save[ibin] = nmodes_a_ inferred from *_b_
            k1_save[ibin] = k_eff_a_;
          }
        }

        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          double k_cen = kbinning.bin_centres[ibin];
          double k_lower = kbinning.bin_edges[ibin];
          double k_upper = kbinning.bin_edges[ibin + 1];

          double k_eff_b_;
          int nmodes_b_;

          F_lm_b.inv_fourier_transform_ylm_wgtd_field_band_limited(
            dn_LM_b, ylm_k_b, k_lower, k_upper, k_eff_b_, nmodes_b_
          );

          k2_save[ibin] = k_eff_b_;
          nmodes_save[ibin] = nmodes_b_;

          if (params.form == "diag") {
            double k_eff_a_;  // redundant as this is `k_eff_b_`
            int nmodes_a_;  // redundant as this is `nmodes_b_`

            F_lm_a.inv_fourier_transform_ylm_wgtd_field_band_limited(
              dn_LM_a, ylm_k_a, k_lower, k_upper, k_eff_a_, nmodes_a_
            );

            k1_save[ibin] = k_eff_a_;
          }

          std::complex<double> bk_component = 0.;  // B_{l₁ l₂ L}^{m₁ m₂ M}
          for (int gid = 0; gid < params.nmesh; gid++) {
            std::complex<double> F_lm_a_gridpt(F_lm_a[gid][0], F_lm_a[gid][1]);
            std::complex<double> F_lm_b_gridpt(F_lm_b[gid][0], F_lm_b[gid][1]);
            std::complex<double> G_LM_gridpt(G_LM[gid][0], G_LM[gid][1]);
            bk_component += F_lm_a_gridpt * F_lm_b_gridpt * G_LM_gridpt;
          }

          bk_save[ibin] += coupling * vol_cell * bk_component;
        }

        // ·······························································
        // Shot noise
        // ·······························································

        // Compute shot noise components in eqs. (45) & (46) in the Paper.
        MeshField dn_LM_a_for_sn(params);  // δn_LM_a(k) (for shot noise)
        if (los_choice == 0) {
          dn_LM_a_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_a_for_sn = dn_00;
        }
        dn_LM_a_for_sn.fourier_transform();

        MeshField dn_LM_b_for_sn(params);  // δn_LM_b(k) (for shot noise)
        if (los_choice == 1) {
          dn_LM_b_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_b_for_sn = dn_00;
        }
        dn_LM_b_for_sn.fourier_transform();

        MeshField dn_LM_c_for_sn(params);  // δn_LM_c(k) (for shot noise)
        if (los_choice == 2) {
          dn_LM_c_for_sn.compute_ylm_wgtd_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        } else {
          dn_LM_c_for_sn = dn_00;
        }
        dn_LM_c_for_sn.fourier_transform();

        MeshField N_LM_a(params);  // N_LM_a(k)
        if (los_choice == 0) {
          N_LM_a = N_00;
        } else {
          N_LM_a.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand,
            alpha, params.ELL, M_
          );
        }
        N_LM_a.fourier_transform();

        MeshField N_LM_b(params);  // N_LM_b(k)
        if (los_choice == 1) {
          N_LM_b = N_00;
        } else {
          N_LM_b.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM_b.fourier_transform();

        MeshField N_LM_c(params);  // N_LM_c(k)
        if (los_choice == 2) {
          N_LM_c = N_00;
        } else {
          N_LM_c.compute_ylm_wgtd_quad_field(
            catalogue_data, catalogue_rand, los_data, los_rand, alpha,
            params.ELL, M_
          );
        }
        N_LM_c.fourier_transform();

        std::complex<double> Sbar_LM = calc_ylm_wgtd_shotnoise_amp_for_bispec(
          catalogue_data, catalogue_rand, los_data, los_rand, alpha,
          params.ELL, M_
        );  // \bar{S}_LM

        if (params.ell1 == 0 && params.ell2 == 0) {
          // When l₁ = l₂ = 0, the Wigner 3-j symbol enforces L = 0
          // and the pre-factors involving degrees and orders become 1.
          std::complex<double> S_ijk = coupling * Sbar_LM;  // S|{i = j = k}
          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            sn_save[ibin] += S_ijk;
          }
        }

        if (params.ell2 == 0) {
          // When l₂ = 0, the Wigner 3-j symbol enforces L = l₁.
          FieldStats stats_sn(params);  // S|{i ≠ j = k}
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_LM_a_for_sn, N_LM_a, Sbar_LM, params.ell1, m1_, kbinning
          );
          if (params.form == "diag") {
            for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[ibin] - stats_sn.sn[ibin]
              );
            }
          } else if (params.form == "full") {
            for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
              sn_save[ibin] += coupling * (
                stats_sn.pk[params.idx_bin] - stats_sn.sn[params.idx_bin]
              );
            }
          }
        }

        if (params.ell1 == 0) {
          // When l₁ = 0, the Wigner 3-j symbol enforces L = l₂.
            FieldStats stats_sn(params);  // S|{j ≠ i = k}
          stats_sn.compute_ylm_wgtd_2pt_stats_in_fourier(
            dn_LM_b_for_sn, N_LM_b, Sbar_LM,
            params.ell2, m2_, kbinning
          );
          for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
            sn_save[ibin] +=
              coupling * (stats_sn.pk[ibin] - stats_sn.sn[ibin]);
          }
        }

        FieldStats stats_sn(params);
        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          double k_a = k1_save[ibin];
          double k_b = k2_save[ibin];

          std::complex<double> S_ij_k = parity
            * stats_sn.compute_uncoupled_shotnoise_for_bispec_per_bin(
              dn_LM_c_for_sn, N_LM_c, ylm_r_a, ylm_r_b, sj_a, sj_b,
              Sbar_LM, k_a, k_b
            );  // S|{i = j ≠ k}

          sn_save[ibin] += coupling * S_ij_k;
        }

        if (trvs::currTask == 0) {
          trvs::logger.stat(
            "Bispectrum term at orders (m1, m2, M) = (%d, %d, %d) computed.",
            m1_, m2_, M_
          );
        }
      }

      delete[] ylm_k_a; ylm_k_a = nullptr;
      delete[] ylm_k_b; ylm_k_b = nullptr;
      delete[] ylm_r_a; ylm_r_a = nullptr;
      delete[] ylm_r_b; ylm_r_b = nullptr;
      trvs::gbytesMem -=
        trvs::size_in_gb< std::complex<double> >(4 * params.nmesh);
    }
  }

  dn_00.finalise_density_field();  // ~dn_00 (likely redundant but safe)
  N_00.finalise_density_field();  // ~N_00 (likely redundant but safe)

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::BispecMeasurements bispec_out;
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    if (params.form == "diag") {
      bispec_out.k1_bin.push_back(kbinning.bin_centres[ibin]);
    } else
    if (params.form == "full") {
      bispec_out.k1_bin.push_back(kbinning.bin_centres[params.idx_bin]);
    }
    bispec_out.k1_eff.push_back(k1_save[ibin]);
    bispec_out.k2_bin.push_back(kbinning.bin_centres[ibin]);
    bispec_out.k2_eff.push_back(k2_save[ibin]);
    bispec_out.nmodes.push_back(nmodes_save[ibin]);
    bispec_out.bk_raw.push_back(norm_factor * bk_save[ibin]);
    bispec_out.bk_shot.push_back(norm_factor * sn_save[ibin]);
  }

  delete[] nmodes_save; delete[] k1_save; delete[] k2_save;
  delete[] bk_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed bispectrum from paired survey-type catalogues "
      "for line-of-sight choice %d.",
      los_choice
    );
  }

  return bispec_out;
}
#endif  // TRV_USE_LEGACY_CODE

}  // namespace trv
