/**
 * @file measurements.cpp
 * @brief Perform N-point correlator clustering measurements.
 *
 */

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <string>

#include "common.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "particles.hpp"
#include "twopt.hpp"
#include "threept.hpp"

/**
 * Main program performing measurements.
 */
int main(int argc, char* argv[]) {
  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "%s\n[STAT] (+%s) Program has started.\n",
      std::string(80, '>').c_str(),
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /* * Initialisation ****************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Initialising program...\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /// Configure parameters.
  if (argc != 2) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[ERRO] (+%s) Missing parameter file in call.\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
      );
    }
    exit(1);
  }

  ParameterSet params;  // program parameters
  if (params.read_from_file(argv)) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[ERRO] (+%s) Failed to initialise parameters.\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
      );
    }
    exit(1);
  }

  if (!(params.printout())) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[INFO] (+%s) Check 'parameters_used' file in your "
        "measurement output directory for reference.\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
      );
    }
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) ... initialised program.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /** Data processing ****************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Reading catalogues...\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /// Read catalogue files.
  ParticleCatalogue particles_data, particles_rand;  // catalogues
  if (particles_data.read_particle_data_from_file(
    params.data_catalogue_file, params.catalogue_header
  )) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[ERRO] (+%s) Failed to load data-source catalogue file.\n",
        calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
      );
    }
    exit(1);
  }

  std::string flag_rand = "false";  // random catalogue status
  if (is_set(params.rand_catalogue_file)) {
    if (particles_rand.read_particle_data_from_file(
      params.rand_catalogue_file, params.catalogue_header
    )) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Failed to load random-source catalogue file.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      exit(1);
    }
    flag_rand = "true";
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) ... read catalogues.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /* * Measurements ******************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Making measurements...\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /// Set up measurements.
  std::string flag_npoint;  // n-point order of the measured statistics
  if (
    params.measurement_type == "powspec"
    || params.measurement_type == "2pcf"
    || params.measurement_type == "2pcf-win"
  ) {
    flag_npoint = "2pt";
  }
  if (
    params.measurement_type == "bispec"
    || params.measurement_type == "3pcf"
    || params.measurement_type == "3pcf-win"
    || params.measurement_type == "3pcf-win-wa"
  ) {
    flag_npoint = "3pt";
  }

  double kbin[params.num_kbin];  // wavenumber bins
  BinScheme::set_kbin(params, kbin);

  double rbin[params.num_rbin];  // separation bins
  BinScheme::set_rbin(params, rbin);

  bool save = true;  // whether to save the results or not

  /// Compute line of sight.
  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Computing lines of sight...\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  LineOfSight* los_data = new LineOfSight[particles_data.ntotal];  // data LoS
  for (int pid = 0; pid < particles_data.ntotal; pid++) {
    double los_mag = sqrt(
      particles_data[pid].pos[0] * particles_data[pid].pos[0]
      + particles_data[pid].pos[1] * particles_data[pid].pos[1]
      + particles_data[pid].pos[2] * particles_data[pid].pos[2]
    );

    los_data[pid].pos[0] = particles_data[pid].pos[0] / los_mag;
    los_data[pid].pos[1] = particles_data[pid].pos[1] / los_mag;
    los_data[pid].pos[2] = particles_data[pid].pos[2] / los_mag;
  }

  LineOfSight* los_rand = new LineOfSight[particles_rand.ntotal];  // random LoS
  for (int pid = 0; pid < particles_rand.ntotal; pid++) {
    double los_mag = sqrt(
      particles_rand[pid].pos[0] * particles_rand[pid].pos[0]
      + particles_rand[pid].pos[1] * particles_rand[pid].pos[1]
      + particles_rand[pid].pos[2] * particles_rand[pid].pos[2]
    );

    los_rand[pid].pos[0] = particles_rand[pid].pos[0] / los_mag;
    los_rand[pid].pos[1] = particles_rand[pid].pos[1] / los_mag;
    los_rand[pid].pos[2] = particles_rand[pid].pos[2] / los_mag;
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) ... computed lines of sight.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /// Offset particle positions for measurements.
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    double ngrid_pad[3] = {3., 3., 3.};
    ParticleCatalogue::pad_pair_in_box(
      particles_data, particles_rand,
      params.boxsize, params.ngrid,
      ngrid_pad
    );
  }
  if (params.catalogue_type == "sim") {
    particles_data.offset_coords_for_periodicity(params.boxsize);
  }

  /// Compute number density alpha ratio.
  double alpha;  // alpha ratio
  if (flag_rand == "true") {
    alpha = particles_data.wtotal / particles_rand.wtotal;
  } else {
    alpha = 1.;
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[INFO] (+%s) Alpha contrast: %.6e.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str(), alpha
    );
  }

  /// Compute normalisation factor for clustering statistics.
  double norm;  // normalisation factor
  if (
    params.norm_convention == "mesh"
  ) {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      }
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      }
    } else {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      }
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      }
    }
  } else if (
    params.norm_convention == "particle"
  ) {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_particles(particles_rand, alpha);
      }
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_particles(particles_rand, alpha);
      }
    } else {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_particles(particles_data, alpha);
      }
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_particles(particles_data, alpha);
      }
    }
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[INFO] (+%s) Normalisation constant: %.6e.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str(), norm
    );
  }

  /// Perform clustering-statistics algorithms.
  if (params.measurement_type == "powspec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      compute_powspec(
        particles_data, particles_rand, los_data, los_rand,
        params, kbin, alpha, norm, save
      );
    }
    if (params.catalogue_type == "sim") {
      compute_powspec_in_box(particles_data, params, kbin, save);
    }
  }
  if (params.measurement_type == "2pcf") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      compute_corrfunc(
        particles_data, particles_rand, los_data, los_rand,
        params, rbin, alpha, norm, save
      );
    }
    if (params.catalogue_type == "sim") {
      compute_corrfunc_in_box(particles_data, params, rbin, save);
    }
  }
  if (params.measurement_type == "2pcf-win") {
    compute_corrfunc_window(
      particles_rand, los_rand, params, rbin, alpha, norm, save
    );
  }

  if (params.measurement_type == "bispec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      calc_bispec(
        particles_data, particles_rand, los_data, los_rand,
        params, kbin, alpha, norm
      );
    }
    if (params.catalogue_type == "sim") {
      calc_bispec_in_box(particles_data, params, kbin);
    }
  }
  if (params.measurement_type == "3pcf") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      calc_3pcf(
        particles_data, particles_rand, los_data, los_rand,
        params, alpha, rbin, norm
      );
    }
    if (params.catalogue_type == "sim") {
      calc_3pcf_in_box(particles_data, params, rbin);
    }
  }
  if (params.measurement_type == "3pcf-win") {
    calc_3pcf_window(
      particles_rand, los_rand, params, alpha, rbin, norm
    );
  }
  if (params.measurement_type == "3pcf-win-wa") {
    calc_3pcf_window_for_wide_angle(
      particles_rand, los_rand, params, alpha, rbin, norm
    );
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) ... made measurements.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /* * Finalisation ******************************************************** */

  particles_data.finalise_particles();
  particles_rand.finalise_particles();

  delete[] los_data; los_data = NULL;
  delete[] los_rand; los_rand = NULL;

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Persistent memory usage: %.0f bytes.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
      gbytesMem
    );
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Program has completed.\n%s\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
      std::string(80, '<').c_str()
    );
  }

  return 0;
}
