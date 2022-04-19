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
      "%s\n[%s STAT] Program has started.\n",
      std::string(80, '>').c_str(),
      show_timestamp()
    );
  }

  /* * Initialisation ****************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] Initialising program...\n",
      show_timestamp()
    );
  }

  /// Configure parameters.
  if (argc != 2) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[%s ERRO] Missing parameter file in call.\n",
        show_timestamp()
      );
    }
    exit(1);
  }

  ParameterSet params;  // program parameters
  if (params.read_from_file(argv)) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[%s ERRO] Failed to initialise parameters.\n",
        show_timestamp()
      );
    }
    exit(1);
  }

  if (!(params.printout())) {
    if (currTask == 0) {
      clockElapsed = double(clock() - clockStart);
      printf(
        "[%s INFO] Check 'parameters_used' file in your "
        "measurement output directory for reference.\n",
        show_timestamp()
      );
    }
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] ... initialised program.\n",
      show_timestamp()
    );
  }

  /** Data processing ****************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] Reading catalogues...\n",
      show_timestamp()
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
        "[%s ERRO] Failed to load data-source catalogue file.\n",
        show_timestamp()
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
          "[%s ERRO] Failed to load random-source catalogue file.\n",
          show_timestamp()
        );
      }
      exit(1);
    }
    flag_rand = "true";
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] ... read catalogues.\n",
      show_timestamp()
    );
  }

  /* * Measurements ******************************************************** */

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] Making measurements...\n",
      show_timestamp()
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
      "[%s STAT] Computing lines of sight...\n",
      show_timestamp()
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
      "[%s STAT] ... computed lines of sight.\n",
      show_timestamp()
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
  } else
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
      "[%s INFO] Alpha contrast: %.6e.\n",
      show_timestamp(), alpha
    );
  }

  /// Compute normalisation factor for clustering statistics.
  double norm;  // normalisation factor
  if (params.norm_convention == "mesh") {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      } else
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
      } else
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      }
    }
  } else
  if (params.norm_convention == "particle") {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_particles(particles_rand, alpha);
      } else
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_particles(particles_rand, alpha);
      }
    } else {
      if (flag_npoint == "2pt") {
        norm = calc_powspec_normalisation_from_particles(particles_data, alpha);
      } else
      if (flag_npoint == "3pt") {
        norm = calc_bispec_normalisation_from_particles(particles_data, alpha);
      }
    }
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s INFO] Normalisation constant: %.6e.\n",
      show_timestamp(), norm
    );
  }

  /// Perform clustering-statistics algorithms.
  if (params.measurement_type == "powspec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      compute_powspec(
        particles_data, particles_rand, los_data, los_rand,
        params, kbin, alpha, norm, save
      );
    } else
    if (params.catalogue_type == "sim") {
      compute_powspec_in_box(particles_data, params, kbin, save);
    }
  } else
  if (params.measurement_type == "2pcf") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      compute_corrfunc(
        particles_data, particles_rand, los_data, los_rand,
        params, rbin, alpha, norm, save
      );
    } else
    if (params.catalogue_type == "sim") {
      compute_corrfunc_in_box(particles_data, params, rbin, save);
    }
  } else
  if (params.measurement_type == "2pcf-win") {
    compute_corrfunc_window(
      particles_rand, los_rand, params, rbin, alpha, norm, save
    );
  } else

  if (params.measurement_type == "bispec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      calc_bispec(
        particles_data, particles_rand, los_data, los_rand,
        params, kbin, alpha, norm
      );
    } else
    if (params.catalogue_type == "sim") {
      calc_bispec_in_box(particles_data, params, kbin, norm);
    }
  } else
  if (params.measurement_type == "3pcf") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      calc_3pcf(
        particles_data, particles_rand, los_data, los_rand,
        params, rbin, alpha, norm
      );
    } else
    if (params.catalogue_type == "sim") {
      calc_3pcf_in_box(particles_data, params, rbin, norm);
    }
  } else
  if (params.measurement_type == "3pcf-win") {
    calc_3pcf_window(
      particles_rand, los_rand, params, rbin, alpha, norm
    );
  } else
  if (params.measurement_type == "3pcf-win-wa") {
    calc_3pcf_window_for_wide_angle(
      particles_rand, los_rand, params, rbin, alpha, norm
    );
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] ... made measurements.\n",
      show_timestamp()
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
      "[%s STAT] Persistent memory usage: %.0f bytes.\n",
      show_timestamp(),
      gbytesMem
    );
  }

  if (currTask == 0) {
    clockElapsed = double(clock() - clockStart);
    printf(
      "[%s STAT] Program has completed.\n%s\n",
      show_timestamp(),
      std::string(80, '<').c_str()
    );
  }

  return 0;
}
