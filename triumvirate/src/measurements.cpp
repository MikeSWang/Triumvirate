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

  ParameterSet params;
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

  /// Set up measurements.
  double kbin[params.num_kbin];
  BinScheme::set_kbin(params, kbin);

  double rbin[params.num_rbin];
  BinScheme::set_rbin(params, rbin);

  bool save = true;

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) ... program initialised.\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /** Data processing ****************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) Reading and processing catalogues...\n",
      calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
    );
  }

  /// Read catalogue files.
  ParticleCatalogue particles_data, particles_rand;
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

  /// Compute line of sight.
  LineOfSight* los_data = new LineOfSight[particles_data.ntotal];
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

  LineOfSight* los_rand = new LineOfSight[particles_rand.ntotal];
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

  /// Offset particle positions for measurements.
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    ParticleCatalogue::boxify_catalogues_for_fft(
      particles_data, particles_rand, params.boxsize, params.ngrid
    );
  }
  if (params.catalogue_type == "sim") {
    particles_data.offset_coords_for_periodicity(params.boxsize);
  }

  /// Compute number density alpha ratio.
  double alpha = 0.;
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    alpha = particles_data.wtotal / particles_rand.wtotal;
  }

  /// Compute inverse-effective-volume normalisation for clustering statistics.
  double norm;
  if (params.norm_convention == "mesh") {
    norm = calc_powspec_normalisation_from_mesh(particles_rand, params);
  } else if (params.norm_convention == "particle") {
    norm = calc_powspec_normalisation_from_particles(particles_rand, alpha);
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[STAT] (+%s) ... catalogues read and processed.\n",
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

  if (params.measurement_type == "powspec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      compute_powspec(
        particles_data, particles_rand, los_data, los_rand,
        params, kbin,
        alpha, norm / particles_data.wtotal / particles_data.wtotal,
        save
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
        params, alpha, kbin, norm
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
      "[STAT] (+%s) ... measurements made.\n",
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
      bytesMem
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
