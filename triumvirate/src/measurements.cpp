#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <cmath>
#include <string>

#include "common.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "harmonic.hpp"
#include "bessel.hpp"
#include "particles.hpp"
#include "field.hpp"
#include "twopt.hpp"
#include "threept.hpp"

/**
 * Main program performing measurements.
 */
int main(int argc, char* argv[]) {
  if (currTask == 0) {
    printf("%s\n", std::string(80, '>').c_str());
    printf("[Status] :: Program has started.\n");
  }

  /// ******************
  ///   Initialisation
  /// ******************

  timeStart = clock();

  if (currTask == 0) {
    printf("[Status] :: Initialising program...\n");
  }

  /// Initialise parameters.
  if (argc != 2 && currTask == 0) {
    printf(
      "[Error] :: Missing parameter file. Add parameter file to the call.\n"
    );
    exit(1);
  }

  ParameterSet params;
  if (params.read_from_file(argv) && currTask == 0) {
    printf("[Error] :: Failed to initialise parameters.\n");
    exit(1);
  }

  params.printout();

  if (currTask == 0) {
    printf("[Status] :: Parameters have been initialised.\n");
    printf(
      "[Info] :: Check 'parameters_used' file in your "
      "'measurement_dir' directory for reference.\n"
    );
  }

  /// TODO: Internalise binning.
  /// Set up measurements.
  double kbin[params.num_kbin];
  BinScheme::set_kbin(params, kbin);

  double rbin[params.num_rbin];
  BinScheme::set_rbin(params, rbin);

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: ... program initialised (%.3f seconds elapsed in total).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// *******************
  ///   Data Processing
  /// *******************

  if (currTask == 0) {
    printf("[Status] :: Reading and processing catalogues...\n");
  }

  /// Read catalogue files.
  ParticleCatalogue particles_data, particles_rand;
  if (particles_data.read_particle_data_from_file(params.data_catalogue_file)) {
    printf("[Error] :: Failed to load data-source catalogue file.\n");
    exit(1);
  }
  if (particles_rand.read_particle_data_from_file(params.rand_catalogue_file)) {
    printf("[Error] :: Failed to load random-source catalogue file.\n");
    exit(1);
  }

  /// Compute line of sight.
  LineOfSight* los_data = new LineOfSight[particles_data.ntotal];
  LineOfSight* los_rand = new LineOfSight[particles_rand.ntotal];
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

  /// TODO: Internalise `alpha`.
  /// Compute number density alpha ratio.
  double alpha = 0.;
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    alpha = particles_data.wtotal / particles_rand.wtotal;
  }

  /// Offset particle positions for measurements.
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    ParticleCatalogue::align_catalogues_for_fft(
      particles_data, particles_rand, params.boxsize, params.ngrid
    );
  }
  if (params.catalogue_type == "sim") {
    particles_data.offset_coords_for_periodicity(params.boxsize);
  }

  /// TODO: Reimplement survey volume normalisation.
  /// Compute catalogue field volume.
  DensityField<ParticleCatalogue> field_rand(params);

  double vol_norm = field_rand.calc_volume_normalisation(particles_rand);

  field_rand.finalise_density_field();

  /// TODO: Reimplement printout.
  /// Print out catalogue field properties.
  if (currTask == 0) {
		std::cout << "[Status] :: See printout below." << std::endl;
    std::cout << "-- Data Catalogue --" << std::endl;
    std::cout << "Number count : N = "
      << particles_data.ntotal << std::endl;
    std::cout << "Volume : V = "
      << params.volume << std::endl;
    std::cout << "Size along x-axis : L_x = "
      << params.boxsize[0] << std::endl;
    std::cout << "Size along y-axis : L_y = "
      << params.boxsize[1] << std::endl;
    std::cout << "Size along z-axis : L_z = "
      << params.boxsize[2] << std::endl;
    std::cout << "" << std::endl;

    std::cout << "-- Data Density Field --" << std::endl;
    std::cout << particles_data.pos_min[0] << " <= x_data <= "
      << particles_data.pos_max[0] << std::endl;
    std::cout << particles_data.pos_min[1] << " <= y_data <= "
      << particles_data.pos_max[1] << std::endl;
    std::cout << particles_data.pos_min[2] << " <= z_data <= "
      << particles_data.pos_max[2] << std::endl;
    std::cout << "" << std::endl;

    std::cout << "-- Random Density Field --" << std::endl;
    std::cout << particles_rand.pos_min[0] << " <= x_rand <= "
      << particles_rand.pos_max[0] << std::endl;
    std::cout << particles_rand.pos_min[1] << " <= y_rand <= "
      << particles_rand.pos_max[1] << std::endl;
    std::cout << particles_rand.pos_min[2] << " <= z_rand <= "
      << particles_rand.pos_max[2] << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Volume normalisation : V_norm = "
      << survey_vol_norm << std::endl;
    std::cout << "Size of effective volume (geometric mean) : L_norm = "
      << pow(survey_vol_norm, 1./3.) << std::endl;
  }

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: ... catalogues read and processed "
      "(%.3f seconds elapsed in total).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// ****************
  ///   Measurements
  /// ****************

  if (currTask == 0) {
    printf("[Status] :: Making measurements...\n");
  }

  /// TODO: Internalsie arguments etc.

  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    if (params.measurement_type == "powspec") {
      calc_powspec(
        particles_data, particles_rand,
        los_data, los_rand,
        params,
        alpha,
        kbin,
        vol_norm
      );
    }
    if (params.measurement_type == "2pcf") {
      calc_corrfunc(
        particles_data, particles_rand,
        los_data, los_rand,
        params,
        alpha,
        rbin,
        vol_norm
      );
    }
    if (params.measurement_type == "2pcf-win") {
      calc_corrfunc_window(
        particles_rand, los_rand, params, alpha, rbin, vol_norm
      );
    }

    if (params.measurement_type == "bispec") {
      calc_bispec(
        particles_data, particles_rand,
        los_data, los_rand,
        params,
        alpha,
        kbin,
        vol_norm
      );
    }
    if (params.measurement_type == "3pcf") {
      calc_3pt_corrfunc(
        particles_data, particles_rand,
        los_data, los_rand,
        params,
        alpha,
        rbin,
        vol_norm
      );
    }
    if (params.measurement_type == "3pcf-win") {
      calc_3pt_corrfunc_window(
        particles_rand, los_rand, params, alpha, rbin, vol_norm
      );
    }
    if (params.measurement_type == "3pcf-win-wa") {
      calc_3pt_corrfunc_window_for_wide_angle(
        particles_rand, los_rand, params, alpha, rbin, vol_norm
      );
    }
  }

  if (params.catalogue_type == "sim") {
    calc_powspec_in_box(particles_data, params, kbin);
    calc_corrfunc_in_box(particles_data, params, rbin);

    calc_bispec_in_box(particles_data, params, kbin);
    calc_3pt_corrfunc_in_box(particles_data, params, rbin);
  }

  durationInSec = double(clock() - timeStart);
  if (currTask == 0) {
    printf(
      "[Status] :: ... measurements made (%.3f seconds elapsed in total).\n",
      durationInSec / CLOCKS_PER_SEC
    );
  }

  /// ****************
  ///   Finalisation
  /// ****************

  delete[] los_data; los_data = NULL;
  delete[] los_rand; los_rand = NULL;
  particles_data.finalise_particles();
  particles_rand.finalise_particles();

  durationInSec = double(clock() - timeStart);

  if (currTask == 0) {
    printf("[Info] :: Persistent memory usage: %.0f bytes.\n", bytesMem);
    printf(
      "[Info] :: Total time elapsed: %.3f seconds.\n",
      durationInSec / CLOCKS_PER_SEC
    );
    printf("[Status] :: Program has completed.\n");
    printf("%s\n", std::string(80, '<').c_str());
  }

  return 0;
}
