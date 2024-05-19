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
 * @file twopt.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori Sugiyama (https://github.com/naonori)
 *
 */

#include "twopt.hpp"

/// @cond DOXYGEN_DOC_CONST
/// CAVEAT: Discretionary choice such that the eps_norm/norm is of order 0.01%.
const double eps_norm = 1.e-5;
/// @endcond

namespace trvs = trv::sys;
namespace trvm = trv::maths;

namespace trv {

// ***********************************************************************
// Coupling coefficients
// ***********************************************************************

double calc_coupling_coeff_2pt(int ell, int ELL, int m, int M) {
  return (2*ell + 1) * (2*ELL + 1)
    * trvm::wigner_3j(ell, 0, ELL, 0, 0, 0)
    * trvm::wigner_3j(ell, 0, ELL, m, 0, M);
}


// ***********************************************************************
// Normalisation
// ***********************************************************************

double calc_powspec_normalisation_from_particles(
  ParticleCatalogue& particles, double alpha
) {
  if (particles.pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
    }
    throw trvs::InvalidDataError("Particle data are uninitialised.\n");
  }

  double norm = 0.;  // I₂

#ifdef TRV_USE_OMP
#pragma omp parallel for simd reduction(+:norm)
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    norm += particles[pid].ws
      * particles[pid].nz * std::pow(particles[pid].wc, 2);
  }

  if (norm == 0.) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Particle 'nz' values appear to be all zeros. "
        "Check the input catalogue contains valid 'nz' field."
      );
    }
    throw trvs::InvalidDataError(
      "Particle 'nz' values appear to be all zeros. "
      "Check the input catalogue contains valid 'nz' field.\n"
    );
  }

  double norm_factor = 1. / (alpha * norm);  // 1/I₂

  return norm_factor;
}

double calc_powspec_normalisation_from_mesh(
  trv::ParticleCatalogue& particles, trv::ParameterSet& params, double alpha
) {
  trv::MeshField catalogue_mesh(params, false, "`catalogue_mesh`");

  double norm_factor =
    catalogue_mesh.calc_grid_based_powlaw_norm(particles, 2);

  norm_factor /= std::pow(alpha, 2);

  return norm_factor;
}

double calc_powspec_normalisation_from_meshes(
  trv::ParticleCatalogue& particles_data,
  trv::ParticleCatalogue& particles_rand,
  trv::ParameterSet& params, double alpha
) {
  // Assign particles to mesh.
  trv::MeshField mesh_data(params, false, "`mesh_data`");
  trv::MeshField mesh_rand(params, false, "`mesh_rand`");

  fftw_complex* weight_data = nullptr;
  weight_data = fftw_alloc_complex(particles_data.ntotal);

  fftw_complex* weight_rand = nullptr;
  weight_rand = fftw_alloc_complex(particles_rand.ntotal);

  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles_data.ntotal);
  trvs::gbytesMem += trvs::size_in_gb<fftw_complex>(particles_rand.ntotal);
  trvs::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for simd
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    weight_data[pid][0] = particles_data[pid].w;
    weight_data[pid][1] = 0.;
  }

#ifdef TRV_USE_OMP
#pragma omp parallel for simd
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    weight_rand[pid][0] = particles_rand[pid].w;
    weight_rand[pid][1] = 0.;
  }

  mesh_data.assign_weighted_field_to_mesh(particles_data, weight_data);
  mesh_rand.assign_weighted_field_to_mesh(particles_rand, weight_rand);

  fftw_free(weight_data); weight_data = nullptr;
  fftw_free(weight_rand); weight_rand = nullptr;

  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles_data.ntotal);
  trvs::gbytesMem -= trvs::size_in_gb<fftw_complex>(particles_rand.ntotal);

  // Calculate normalisation.
  double norm = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for simd reduction(+:norm)
#endif  // TRV_USE_OMP
  for (long long gid = 0; gid < params.nmesh; gid++) {
    norm += mesh_data.field[gid][0] * mesh_rand.field[gid][0];
  }

  double vol_cell = params.volume / double(params.nmesh);

  double norm_factor = 1. / (alpha * vol_cell * norm);  // 1/I₂

  return norm_factor;
}

double calc_powspec_normalisation_from_meshes(
  trv::ParticleCatalogue& particles_data,
  trv::ParticleCatalogue& particles_rand,
  trv::ParameterSet& params, double alpha,
  double padding, double cellsize, std::string assignment
) {
  // Modify parameter set for normalisation.
  trv::ParameterSet params_norm(params);

  double boxsize_norm = (1. + padding) * std::max(
    *std::max_element(particles_data.pos_span, particles_data.pos_span + 3),
    *std::max_element(particles_rand.pos_span, particles_rand.pos_span + 3)
  );

  int ngrid_norm = std::ceil(boxsize_norm / cellsize);

  ngrid_norm += ngrid_norm % 2;  // ensure even
  boxsize_norm = ngrid_norm * cellsize;  // enforce cell size

  for (int iaxis = 0; iaxis < 3; iaxis++) {
    params_norm.boxsize[iaxis] = boxsize_norm;
    params_norm.ngrid[iaxis] = ngrid_norm;
  }

  params_norm.assignment = assignment;

  params_norm.validate();

  // Reuse existing overloaded method.
  double norm_factor = calc_powspec_normalisation_from_meshes(
    particles_data, particles_rand, params_norm, alpha
  );

  return norm_factor;
}


// ***********************************************************************
// Shot noise
// ***********************************************************************

double calc_powspec_shotnoise_from_particles(
  ParticleCatalogue& particles, double alpha
) {
  if (particles.pdata == nullptr) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Particle data are uninitialised.");
    }
    throw trvs::InvalidDataError("Particle data are uninitialised.\n");
  }

  double shotnoise = 0.;

#ifdef TRV_USE_OMP
#pragma omp parallel for simd reduction(+:shotnoise)
#endif  // TRV_USE_OMP
  for (int pid = 0; pid < particles.ntotal; pid++) {
    shotnoise +=
      std::pow(particles[pid].ws, 2) * std::pow(particles[pid].wc, 2);
  }

  return shotnoise;
}

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_powspec(
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

    std::complex<double> sn_part = ylm * std::pow(particles_data[pid].w, 2);
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

    std::complex<double> sn_part = ylm * std::pow(particles_rand[pid].w, 2);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_rand_real += sn_part_real;
    sn_rand_imag += sn_part_imag;
  }

  std::complex<double> sn_rand(sn_rand_real, sn_rand_imag);

  return sn_data + std::pow(alpha, 2) * sn_rand;
}

std::complex<double> calc_ylm_wgtd_shotnoise_amp_for_powspec(
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

    std::complex<double> sn_part = ylm * std::pow(particles[pid].w, 2);
    double sn_part_real = sn_part.real();
    double sn_part_imag = sn_part.imag();

    sn_real += sn_part_real;
    sn_imag += sn_part_imag;
  }

  std::complex<double> sn(sn_real, sn_imag);

  return std::pow(alpha, 2) * sn;
}


// ***********************************************************************
// Full statistics
// ***********************************************************************

// STYLE: Standard naming convention is not always followed for
// intermediary quantities in the functions below.

trv::PowspecMeasurements compute_powspec(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing power spectrum from paired survey-type catalogues..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up input.
  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;
  int ell1 = params.ELL;

  // Set up output.
  int* nmodes_save = new int[kbinning.num_bins];
  double* k_save = new double[kbinning.num_bins];
  std::complex<double>* pk_save = new std::complex<double>[kbinning.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[kbinning.num_bins];
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  MeshField dn_00(params, true, "`dn_00`");  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  FieldStats stats_2pt(params);

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params, true, "`dn_LM`");  // δn_LM(k)
    dn_LM.compute_ylm_wgtd_field(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    std::complex<double> sn_amp = trv::calc_ylm_wgtd_shotnoise_amp_for_powspec(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );  // \bar{N}_LM(k)

    // Compute quantity equivalent to (-1)^m₁ δᴰ_{m₁, -M} which, after
    // being summed over m₁, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = calc_coupling_coeff_2pt(ell1, params.ELL, m1, M_);
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      stats_2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
        dn_LM, dn_00, sn_amp, ell1, m1, kbinning
      );

      for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
        pk_save[ibin] += coupling * stats_2pt.pk[ibin];
        sn_save[ibin] += coupling * stats_2pt.sn[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
          nmodes_save[ibin] = stats_2pt.nmodes[ibin];
          k_save[ibin] = stats_2pt.k[ibin];
        }
      }
    }

    if (trvs::currTask == 0) {
      trvs::logger.stat("Power spectrum term computed at order M = %d.", M_);
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    powspec_out.kbin.push_back(kbinning.bin_centres[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmodes.push_back(nmodes_save[ibin]);
    powspec_out.pk_raw.push_back(norm_factor * pk_save[ibin]);
    powspec_out.pk_shot.push_back(norm_factor * sn_save[ibin]);
  }
  powspec_out.dim = kbinning.num_bins;

  delete[] nmodes_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed power spectrum from paired survey-type catalogues."
    );
  }

  return powspec_out;
}

trv::TwoPCFMeasurements compute_corrfunc(
  ParticleCatalogue& catalogue_data, ParticleCatalogue& catalogue_rand,
  LineOfSight* los_data, LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing two-point correlation function from "
      "paired survey-type catalogues..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up input.
  double alpha = catalogue_data.wstotal / catalogue_rand.wstotal;
  int ell1 = params.ELL;

  // Set up output.
  int* npairs_save = new int[rbinning.num_bins];
  double* r_save = new double[rbinning.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[rbinning.num_bins];
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  MeshField dn_00(params, true, "`dn_00`");  // δn_00(k)
  dn_00.compute_ylm_wgtd_field(
    catalogue_data, catalogue_rand, los_data, los_rand, alpha, 0, 0
  );
  dn_00.fourier_transform();

  FieldStats stats_2pt(params);

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params, true, "`dn_LM`");  // δn_LM(k)
    dn_LM.compute_ylm_wgtd_field(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();

    std::complex<double> sn_amp = trv::calc_ylm_wgtd_shotnoise_amp_for_powspec(
      catalogue_data, catalogue_rand, los_data, los_rand, alpha, params.ELL, M_
    );  // \bar{N}_LM(k)

    // Compute quantity equivalent to (-1)^m₁ δᴰ_{m₁, -M} which, after
    // being summed over m₁, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = calc_coupling_coeff_2pt(ell1, params.ELL, m1, M_);
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      stats_2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, ell1, m1, rbinning
      );

      for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
        xi_save[ibin] += coupling * stats_2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          npairs_save[ibin] = stats_2pt.npairs[ibin];
          r_save[ibin] = stats_2pt.r[ibin];
        }
      }
    }

    if (trvs::currTask == 0) {
      trvs::logger.stat(
        "Two-point correlation function term computed at order M = %d.", M_
      );
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::TwoPCFMeasurements corrfunc_out;
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    corrfunc_out.rbin.push_back(rbinning.bin_centres[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npairs.push_back(npairs_save[ibin]);
    corrfunc_out.xi.push_back(norm_factor * xi_save[ibin]);
  }
  corrfunc_out.dim = rbinning.num_bins;

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed two-point correlation function "
      "from paired survey-type catalogues."
    );
  }

  return corrfunc_out;
}

trv::PowspecMeasurements compute_powspec_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning kbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing power spectrum from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  int* nmodes_save = new int[kbinning.num_bins];
  double* k_save = new double[kbinning.num_bins];
  std::complex<double>* pk_save = new std::complex<double>[kbinning.num_bins];
  std::complex<double>* sn_save = new std::complex<double>[kbinning.num_bins];
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    nmodes_save[ibin] = 0;
    k_save[ibin] = 0.;
    pk_save[ibin] = 0.;
    sn_save[ibin] = 0.;
  }  // likely redundant but safe

  // Check input normalisation matches expectation.
  double norm = double(catalogue_data.ntotal) * double(catalogue_data.ntotal)
    / params.volume;
  if (std::fabs(1 - norm * norm_factor) > eps_norm) {
    if (trvs::currTask == 0) {
      trvs::logger.warn(
        "Power spectrum normalisation input differs from "
        "expected value for an unweight field in a periodic box."
      );
    }
  }

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute power spectrum.
  MeshField dn(params, true, "`dn`");  // δn(k)
  dn.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(catalogue_data.ntotal);  // \bar{N}

  // Under the global plane-parallel approximation, δᴰ_{M0} enforces
  // M = 0 for any spherical-harmonic-weighted field fluctuations.
  FieldStats stats_2pt(params);
  stats_2pt.compute_ylm_wgtd_2pt_stats_in_fourier(
    dn, dn, sn_amp, params.ELL, 0, kbinning
  );

  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    nmodes_save[ibin] = stats_2pt.nmodes[ibin];
    k_save[ibin] = stats_2pt.k[ibin];
    pk_save[ibin] += double(2*params.ELL + 1) * stats_2pt.pk[ibin];
    sn_save[ibin] += double(2*params.ELL + 1) * stats_2pt.sn[ibin];
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  // Fill in output struct.
  trv::PowspecMeasurements powspec_out;
  for (int ibin = 0; ibin < kbinning.num_bins; ibin++) {
    powspec_out.kbin.push_back(kbinning.bin_centres[ibin]);
    powspec_out.keff.push_back(k_save[ibin]);
    powspec_out.nmodes.push_back(nmodes_save[ibin]);
    powspec_out.pk_raw.push_back(norm_factor * pk_save[ibin]);
    powspec_out.pk_shot.push_back(norm_factor * sn_save[ibin]);
  }
  powspec_out.dim = kbinning.num_bins;

  delete[] nmodes_save; delete[] k_save; delete[] pk_save; delete[] sn_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed power spectrum "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  return powspec_out;
}

trv::TwoPCFMeasurements compute_corrfunc_in_gpp_box(
  ParticleCatalogue& catalogue_data,
  trv::ParameterSet& params, trv::Binning& rbinning,
  double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing two-point correlation function "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  int* npairs_save = new int[rbinning.num_bins];
  double* r_save = new double[rbinning.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[rbinning.num_bins];
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }  // likely redundant but safe

  // Check input normalisation matches expectation.
  double norm = double(catalogue_data.ntotal) * double(catalogue_data.ntotal)
    / params.volume;
  if (std::fabs(1 - norm * norm_factor) > eps_norm) {
    if (trvs::currTask == 0) {
      trvs::logger.warn(
        "Two-point correlation function normalisation input differs from "
        "expected value for an unweight field in a periodic box."
      );
    }
  }

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  // Compute 2PCF.
  MeshField dn(params, true, "`dn`");  // δn(k)
  dn.compute_unweighted_field_fluctuations_insitu(catalogue_data);
  dn.fourier_transform();

  std::complex<double> sn_amp = double(catalogue_data.ntotal);  // \bar{N}

  // Under the global plane-parallel approximation, δᴰ_{M0} enforces
  // M = 0 for any spherical-harmonic-weighted field fluctuations.
  FieldStats stats_2pt(params);
  stats_2pt.compute_ylm_wgtd_2pt_stats_in_config(
    dn, dn, sn_amp, params.ELL, 0, rbinning
  );

  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = stats_2pt.npairs[ibin];
    r_save[ibin] = stats_2pt.r[ibin];
    xi_save[ibin] += double(2*params.ELL + 1) * stats_2pt.xi[ibin];
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  // Fill in output struct.
  trv::TwoPCFMeasurements corrfunc_out;
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    corrfunc_out.rbin.push_back(rbinning.bin_centres[ibin]);
    corrfunc_out.reff.push_back(r_save[ibin]);
    corrfunc_out.npairs.push_back(npairs_save[ibin]);
    corrfunc_out.xi.push_back(norm_factor * xi_save[ibin]);
  }
  corrfunc_out.dim = rbinning.num_bins;

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed two-point correlation function "
      "from a periodic-box simulation-type catalogue "
      "in the global plane-parallel approximation."
    );
  }

  return corrfunc_out;
}

trv::TwoPCFWindowMeasurements compute_corrfunc_window(
  trv::ParticleCatalogue& catalogue_rand, trv::LineOfSight* los_rand,
  trv::ParameterSet& params, trv::Binning rbinning,
  double alpha, double norm_factor
) {
  trvs::logger.reset_level(params.verbose);

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "Computing two-point correlation function window "
      "from a random catalogue..."
    );
  }

  // ---------------------------------------------------------------------
  // Set-up
  // ---------------------------------------------------------------------

  // Set up input.
  int ell1 = params.ELL;

  // Set up output.
  int* npairs_save = new int[rbinning.num_bins];
  double* r_save = new double[rbinning.num_bins];
  std::complex<double>* xi_save = new std::complex<double>[rbinning.num_bins];
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    npairs_save[ibin] = 0;
    r_save[ibin] = 0.;
    xi_save[ibin] = 0.;
  }  // likely redundant but safe

  // ---------------------------------------------------------------------
  // Measurement
  // ---------------------------------------------------------------------

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_init_threads();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  MeshField dn_00(params, true, "`dn_00`");
  dn_00.compute_ylm_wgtd_field(catalogue_rand, los_rand, alpha, 0, 0);
  dn_00.fourier_transform();  // δn_00(k)

  FieldStats stats_2pt(params);

  for (int M_ = - params.ELL; M_ <= params.ELL; M_++) {
    MeshField dn_LM(params, true, "`dn_LM`");
    dn_LM.compute_ylm_wgtd_field(
      catalogue_rand, los_rand, alpha, params.ELL, M_
    );
    dn_LM.fourier_transform();  // δn_LM(k)

    std::complex<double> sn_amp = trv::calc_ylm_wgtd_shotnoise_amp_for_powspec(
      catalogue_rand, los_rand, alpha, params.ELL, M_
    );  // \bar{N}_LM(k)

    // Calculate equivalent to (-1)^m₁ δᴰ_{m₁, -M} which, after being
    // summed over m₁, agrees with Hand et al. (2017) [1704.02357].
    for (int m1 = - ell1; m1 <= ell1; m1++) {
      double coupling = trv::calc_coupling_coeff_2pt(ell1, params.ELL, m1, M_);
      if (std::fabs(coupling) < trvm::eps_coupling) {continue;}

      stats_2pt.compute_ylm_wgtd_2pt_stats_in_config(
        dn_LM, dn_00, sn_amp, ell1, m1, rbinning
      );

      for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
        xi_save[ibin] += coupling * stats_2pt.xi[ibin];
      }

      if (M_ == 0 && m1 == 0) {
        for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
          npairs_save[ibin] = stats_2pt.npairs[ibin];
          r_save[ibin] = stats_2pt.r[ibin];
        }
      }
    }

    if (trvs::currTask == 0) {
      trvs::logger.stat(
        "Two-point correlation function window term computed at order M = %d.",
        M_
      );
    }
  }

  // ---------------------------------------------------------------------
  // Results
  // ---------------------------------------------------------------------

  trv::TwoPCFWindowMeasurements corrfunc_win_out;
  for (int ibin = 0; ibin < rbinning.num_bins; ibin++) {
    corrfunc_win_out.rbin.push_back(rbinning.bin_centres[ibin]);
    corrfunc_win_out.reff.push_back(r_save[ibin]);
    corrfunc_win_out.npairs.push_back(npairs_save[ibin]);
    corrfunc_win_out.xi.push_back(norm_factor * xi_save[ibin]);
  }
  corrfunc_win_out.dim = rbinning.num_bins;

  delete[] npairs_save; delete[] r_save; delete[] xi_save;

  if (trvs::currTask == 0) {
    trvs::logger.stat(
      "... computed two-point correlation function window "
      "from a random catalogue."
    );
  }

  return corrfunc_win_out;
}

}  // namespace trv
