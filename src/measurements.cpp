#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#include <fftw3.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_coupling.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_spline.h>

/// Provide standalone functions.
#define wigner_3j(j1,j2,j3,m1,m2,m3) ( \
  gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3) \
)

#include "measurements/common.hpp"
#include "measurements/parameters.hpp"
#include "measurements/tools.hpp"
#include "measurements/bessel.hpp"
#include "measurements/particles.hpp"
#include "measurements/field.hpp"
#include "measurements/power.hpp"
#include "measurements/bispec.hpp"

/**
 * Main measurement program.
 */
int main(int argc, char *argv[]) {
  if (thisTask == 0) {
    printf("Program started.\n\n");
  }

  /// ******************
  ///   Initialisation
  /// ******************

  timeStart = clock();

  if (thisTask == 0) {
    printf("%s\n", std::string(80, '>').c_str());
    printf("Initialising...\n\n");
  }

  /// Initialise parameters.
  if (argc != 2 && thisTask == 0) {
    printf("[Error] :: Missing parameter file. Add parameter file to the call.\n");
    exit(1);
  }

  ParameterSet params;
  if (params.read_parameters(argv) && thisTask == 0) {
    printf("[Error] :: Failed to initialise parameters.\n");
    exit(1);
  }

  params.set_io_files();

  if (thisTask == 0) {
    printf("[Status] :: Parameters have been initialised.\n");
    printf("[Info] :: Check 'used_parameters' file in your 'output_dir' directory for reference.\n");
  }

  /// Set up measurements.
  double kbin[params.num_kbin];
  ToolCollection::set_kbin(params, kbin);

  double rbin[params.num_rbin];
  ToolCollection::set_rbin(params, rbin);

  durationInSec = double(clock() - timeStart);
  if (thisTask == 0) {
    printf("\n... total time elapsed: %.3f seconds.\n", durationInSec / CLOCKS_PER_SEC);
    printf("%s\n\n", std::string(80, '<').c_str());
  }

  /// *******************
  ///   Data Processing
  /// *******************
  if (thisTask == 0) {
    printf("%s\n", std::string(80, '>').c_str());
    printf("Reading and processing catalogues...\n\n");
  }

  /// Read catalogue files.
  ParticlesCatalogue particles_data, particles_rand;
  if (particles_data.read_particles_catalogue(params.data_catalogue_file)) {
    printf("[Error] :: Failed to load data-source catalogue file.\n");
    exit(1);
  }
  if (particles_rand.read_particles_catalogue(params.rand_catalogue_file)) {
    printf("[Error] :: Failed to load random-source catalogue file.\n");
    exit(1);
  }

  /// Compute line of sight.
  LineOfSight* los_data = new LineOfSight[particles_data.nparticles];
  LineOfSight* los_rand = new LineOfSight[particles_rand.nparticles];
  for (int id = 0; id < particles_data.nparticles; id++) {
      double los_mag = sqrt(
        particles_data[id].pos[0] * particles_data[id].pos[0]
        + particles_data[id].pos[1] * particles_data[id].pos[1]
        + particles_data[id].pos[2] * particles_data[id].pos[2]
      );
      los_data[id].pos[0] = particles_data[id].pos[0] / los_mag;
      los_data[id].pos[1] = particles_data[id].pos[1] / los_mag;
      los_data[id].pos[2] = particles_data[id].pos[2] / los_mag;
  }
  for (int id = 0; id < particles_rand.nparticles; id++) {
      double los_mag = sqrt(
        particles_rand[id].pos[0] * particles_rand[id].pos[0]
        + particles_rand[id].pos[1] * particles_rand[id].pos[1]
        + particles_rand[id].pos[2] * particles_rand[id].pos[2]
      );
      los_rand[id].pos[0] = particles_rand[id].pos[0] / los_mag;
      los_rand[id].pos[1] = particles_rand[id].pos[1] / los_mag;
      los_rand[id].pos[2] = particles_rand[id].pos[2] / los_mag;
  }

  /// Compute number density alpha ratio.
  double alpha = 0.;
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    alpha = ParticlesCatalogue::calc_alpha_ratio(particles_data, particles_rand);
  }

  /// Offset particle positions for measurement standardisation.
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    ParticlesCatalogue::offset_particles_for_fft(particles_data, particles_rand, params);
  }
  if (params.catalogue_type == "sim") {
    particles_data.offset_particles_for_periodicity(params);
  }

  /// Compute catalogue field volume.
  DensityField<ParticlesCatalogue> field_rand(params);
  double survey_vol_norm = field_rand.calc_survey_volume_norm(particles_rand);
  field_rand.finalise_density_field();

  /// Print out catalogue field properties.
  if (thisTask == 0) {
		std::cout << "[Status] :: See printout below." << std::endl;
    std::cout << "-- Data Catalogue --" << std::endl;
    std::cout << "Number count : N = " << particles_data.nparticles << std::endl;
    std::cout << "Volume : V = " << params.volume << std::endl;
    std::cout << "Size along x-axis : L_x = " << params.boxsize[0] << std::endl;
    std::cout << "Size along y-axis : L_y = " << params.boxsize[1] << std::endl;
    std::cout << "Size along z-axis : L_z = " << params.boxsize[2] << std::endl;
    std::cout << "" << std::endl;

    std::cout << "-- Data Density Field --" << std::endl;
    std::cout << particles_data.pos_min[0] << " <= x_data <= " << particles_data.pos_max[0] << std::endl;
    std::cout << particles_data.pos_min[1] << " <= y_data <= " << particles_data.pos_max[1] << std::endl;
    std::cout << particles_data.pos_min[2] << " <= z_data <= " << particles_data.pos_max[2] << std::endl;
    std::cout << "" << std::endl;

    std::cout << "-- Random Density Field --" << std::endl;
    std::cout << particles_rand.pos_min[0] << " <= x_rand <= " << particles_rand.pos_max[0] << std::endl;
    std::cout << particles_rand.pos_min[1] << " <= y_rand <= " << particles_rand.pos_max[1] << std::endl;
    std::cout << particles_rand.pos_min[2] << " <= z_rand <= " << particles_rand.pos_max[2] << std::endl;
    std::cout << "" << std::endl;

    std::cout << "Survey volume normalisation : V_survey = " << survey_vol_norm << std::endl;
    std::cout << "Size of field (geometric mean) : L_survey = " << pow(survey_vol_norm, 1./3.) << std::endl;
  }

  durationInSec = double(clock() - timeStart);
  if (thisTask == 0) {
    printf("\n... total time elapsed: %.3f seconds.\n", durationInSec / CLOCKS_PER_SEC);
    printf("%s\n\n", std::string(80, '<').c_str());
  }

  /// ****************
  ///   Measurements
  /// ****************
  if (thisTask == 0) {
    printf("%s\n", std::string(80, '>').c_str());
    printf("Making measurements...\n\n");
  }

  /// ???: What's happening here.
  if (params.catalogue_type == "survey") {
    if (thisTask >= params.num_rbin) {
      params.ith_rbin = 0;
    } else {
      params.ith_rbin = thisTask;
    }

    if (thisTask >= params.num_kbin) {
      params.ith_kbin = 0;
    } else {
      params.ith_kbin = thisTask;
    }
  }

  /// IDEA: Formalise the choosing of measurement programs by parsing
  /// an additional command-line argument.
  if (params.catalogue_type == "mock" || params.catalogue_type == "survey") {
    /*
    calc_power_spec(particles_data, particles_rand, los_data, los_rand, params, alpha, kbin, survey_vol_norm);
    calc_2pt_func(particles_data, particles_rand, los_data, los_rand, params, alpha, rbin, survey_vol_norm);
    calc_2pt_func_window(particles_rand, los_rand, params, alpha, rbin, survey_vol_norm);
    */

    calc_bispec(particles_data, particles_rand, los_data, los_rand, params, alpha, kbin, survey_vol_norm);
    /*
    calc_3pt_func(particles_data, particles_rand, los_data, los_rand, params, alpha, rbin, survey_vol_norm);
    calc_3pt_func_window(particles_rand, los_rand, params, alpha, rbin, survey_vol_norm);
    calc_3pt_func_window_for_3pcf(particles_rand, los_rand, params, alpha, rbin, survey_vol_norm);
    */
  }

  if (params.catalogue_type == "sim") {
    /*
    calc_power_spec_in_box(particles_data, params, kbin);
    calc_2pt_func_in_box(particles_data, params, rbin);
    */

    calc_bispec_in_box(particles_data, params, kbin);
    /*
    calc_3pt_func_in_box(particles_data, params, rbin);
    */
  }

  durationInSec = double(clock() - timeStart);
  if (thisTask == 0) {
    printf("\n... total time elapsed: %.3f seconds.\n", durationInSec / CLOCKS_PER_SEC);
    printf("%s\n\n", std::string(80, '<').c_str());
  }

  /// ****************
  ///   Finalisation
  /// ****************
  delete[] los_data; los_data = NULL;
  delete[] los_rand; los_rand = NULL;
  particles_data.finalise_particles();
  particles_rand.finalise_particles();

  durationInSec = double(clock() - timeStart);

  if (thisTask == 0) {
    printf("Persistent memory usage: %.0f bytes.\n", bytes);
    printf("Total time elapsed: %.3f seconds.\n\n", durationInSec / CLOCKS_PER_SEC);

    printf("Program ended.\n");
  }

  return 0;
}
