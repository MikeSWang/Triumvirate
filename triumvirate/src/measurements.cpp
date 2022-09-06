/**
 * @file measurements.cpp
 * @brief Perform N-point correlator clustering measurements.
 *
 */

#include "monitor.hpp"
#include "parameters.hpp"
#include "tools.hpp"
#include "particles.hpp"
#include "twopt.hpp"
#include "threept.hpp"

/**
 * Main program performing measurements.
 */
int main(int argc, char* argv[]) {
  if (trv::mon::currTask == 0) {
    std::printf(
      "%s\n[%s STAT] Program started.\n",
      std::string(80, '>').c_str(), trv::mon::show_timestamp().c_str()
    );
  }

  /* * Initialisation ****************************************************** */

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Initialising program...\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  /// Configure parameters.
  if (argc != 2) {
    if (trv::mon::currTask == 0) {
      throw trv::mon::IOError(
        "[%s ERRO] Missing parameter file in call.\n",
        trv::mon::show_timestamp().c_str()
      );
    }

  }

  trv::ParameterSet params;  // program parameters
  if (params.read_from_file(argv[1])) {
    if (trv::mon::currTask == 0) {
      throw trv::mon::IOError(
        "[%s ERRO] Failed to initialise parameters.\n",
        trv::mon::show_timestamp().c_str()
      );
    }
  }

  if (!(params.printout())) {
    if (trv::mon::currTask == 0) {
      std::printf(
        "[%s INFO] Check 'parameters_used*' file in your "
        "measurement output directory for reference.\n",
        trv::mon::show_timestamp().c_str()
      );
    }
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] ... initialised program.\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  /** Data processing ****************************************************** */

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Reading catalogues...\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  /// Read catalogue files.
  trv::obj::ParticleCatalogue particles_data, particles_rand;  // catalogues

  std::string flag_data = "false";  // data catalogue status
  if (trv::mon::if_filepath_is_set(params.data_catalogue_file)) {
    if (particles_data.read_particle_data_from_file(
      params.data_catalogue_file, params.catalogue_columns, params.volume
    )) {
      if (trv::mon::currTask == 0) {
        throw trv::mon::IOError(
          "[%s ERRO] Failed to load data-source catalogue file.\n",
          trv::mon::show_timestamp().c_str()
        );
      }
    }
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      flag_data = "true";
    }
  }

  std::string flag_rand = "false";  // random catalogue status
  if (trv::mon::if_filepath_is_set(params.rand_catalogue_file)) {
    if (particles_rand.read_particle_data_from_file(
      params.rand_catalogue_file, params.catalogue_columns, params.volume
    )) {
      if (trv::mon::currTask == 0) {
        throw trv::mon::IOError(
          "[%s ERRO] Failed to load random-source catalogue file.\n",
          trv::mon::show_timestamp().c_str()
        );
      }
    }
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      flag_rand = "true";
    }
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] ... read catalogues.\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  /* * Measurements ******************************************************** */

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Making measurements...\n",
      trv::mon::show_timestamp().c_str()
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
  } else
  if (
    params.measurement_type == "bispec"
    || params.measurement_type == "3pcf"
    || params.measurement_type == "3pcf-win"
    || params.measurement_type == "3pcf-win-wa"
  ) {
    flag_npoint = "3pt";
  }

  trv::Binning binning(params);
  binning.set_bins();

  bool save = true;  // whether to save the results or not

  /// Compute line of sight.
  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Computing lines of sight...\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  trv::obj::LineOfSight* los_data =
    new trv::obj::LineOfSight[particles_data.ntotal];  // data LoS
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

  trv::obj::LineOfSight* los_rand =
    new trv::obj::LineOfSight[particles_rand.ntotal];  // random LoS
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

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] ... computed lines of sight.\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  /// Offset particle positions for measurements.
  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Offset particle coordinates for box alignment.\n",
      trv::mon::show_timestamp().c_str()
    );
  }
  if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
    if (params.alignment == "pad") {
      if (params.padscale == "grid") {
        double ngrid_pad[3] = {
          params.padfactor, params.padfactor, params.padfactor
        };
        if (flag_data == "true") {
          trv::obj::ParticleCatalogue::pad_pair_in_box(
            particles_data, particles_rand,
            params.boxsize, params.ngrid,
            ngrid_pad
          );
        } else {
          trv::obj::ParticleCatalogue::pad_in_box(
            particles_rand, params.boxsize, params.ngrid, ngrid_pad
          );
        }
      } else
      if (params.padscale == "box") {
        double boxsize_pad[3] = {
          params.padfactor, params.padfactor, params.padfactor
        };
        if (flag_data == "true") {
          trv::obj::ParticleCatalogue::pad_pair_in_box(
            particles_data, particles_rand,
            params.boxsize, boxsize_pad
          );
        } else {
          trv::obj::ParticleCatalogue::pad_in_box(
            particles_rand, params.boxsize, boxsize_pad
          );
        }
      }
    } else
    if (params.alignment == "centre") {
      if (flag_data == "true") {
        trv::obj::ParticleCatalogue::centre_pair_in_box(
          particles_data, particles_rand, params.boxsize
        );
      } else {
        trv::obj::ParticleCatalogue::centre_in_box(
          particles_rand, params.boxsize
        );
      }
    }
  } else
  if (params.catalogue_type == "sim") {
    particles_data.offset_coords_for_periodicity(params.boxsize);
  }

  /// Compute number density alpha ratio.
  double alpha;  // alpha ratio
  if (flag_data == "true" && flag_rand == "true") {
    alpha = particles_data.wtotal / particles_rand.wtotal;
  } else {
    alpha = 1.;
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s INFO] Alpha contrast: %.6e.\n",
      trv::mon::show_timestamp().c_str(), alpha
    );
  }

  /// Compute normalisation factor for clustering statistics.
  double norm;  // normalisation factor
  if (params.norm_convention == "mesh") {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm = trv::algo::calc_powspec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      } else
      if (flag_npoint == "3pt") {
        norm = trv::algo::calc_bispec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      }
    } else {
      if (flag_npoint == "2pt") {
        norm = trv::algo::calc_powspec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      } else
      if (flag_npoint == "3pt") {
        norm = trv::algo::calc_bispec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      }
    }
  } else
  if (params.norm_convention == "particle") {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm = trv::algo::
          calc_powspec_normalisation_from_particles(particles_rand, alpha);
      } else
      if (flag_npoint == "3pt") {
        norm = trv::algo::
          calc_bispec_normalisation_from_particles(particles_rand, alpha);
      }
    } else {
      if (flag_npoint == "2pt") {
        norm = trv::algo::
          calc_powspec_normalisation_from_particles(particles_data, alpha);
      } else
      if (flag_npoint == "3pt") {
        norm = trv::algo::
          calc_bispec_normalisation_from_particles(particles_data, alpha);
      }
    }
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s INFO] Normalisation constant: %.6e.\n",
      trv::mon::show_timestamp().c_str(),
      norm
    );
  }

  /// Compute, optionally, the alternative normalisation factor
  /// for clustering statistics.
  double norm_alt;  // alternative normalisation factor
  #ifdef DBG_NORMALT
  if (params.norm_convention == "mesh") {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm_alt = trv::algo::
          calc_powspec_normalisation_from_particles(particles_rand, alpha);
      } else
      if (flag_npoint == "3pt") {
        norm_alt = trv::algo::
          calc_bispec_normalisation_from_particles(particles_rand, alpha);
      }
    } else {
      if (flag_npoint == "2pt") {
        norm_alt = trv::algo::
          calc_powspec_normalisation_from_particles(particles_data, alpha);
      } else
      if (flag_npoint == "3pt") {
        norm_alt = trv::algo::
          calc_bispec_normalisation_from_particles(particles_data, alpha);
      }
    }
  } else
  if (params.norm_convention == "particle") {
    if (flag_rand == "true") {
      if (flag_npoint == "2pt") {
        norm_alt = trv::algo::calc_powspec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      } else
      if (flag_npoint == "3pt") {
        norm_alt = trv::algo::calc_bispec_normalisation_from_mesh(
          particles_rand, params, alpha
        );
      }
    } else {
      if (flag_npoint == "2pt") {
        norm_alt = trv::algo::calc_powspec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      } else
      if (flag_npoint == "3pt") {
        norm_alt = trv::algo::calc_bispec_normalisation_from_mesh(
          particles_data, params, alpha
        );
      }
    }
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s INFO] Alternative normalisation constant: %.6e.\n",
      trv::mon::show_timestamp().c_str(),
      norm_alt
    );
  }
  #endif  // DBG_NORMALT

  /// Perform clustering-statistics algorithms.
  if (params.measurement_type == "powspec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      trv::algo::compute_powspec(
        particles_data, particles_rand, los_data, los_rand,
        params, binning.bin_centres, alpha, norm, norm_alt, save
      );
    } else
    if (params.catalogue_type == "sim") {
      trv::algo::compute_powspec_in_box(
        particles_data, params, binning.bin_centres, norm, norm_alt, save
      );
    }
  } else
  if (params.measurement_type == "2pcf") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      trv::algo::compute_corrfunc(
        particles_data, particles_rand, los_data, los_rand,
        params, binning.bin_centres, alpha, norm, norm_alt, save
      );
    } else
    if (params.catalogue_type == "sim") {
      trv::algo::compute_corrfunc_in_box(
        particles_data, params, binning.bin_centres, norm, norm_alt, save
      );
    }
  } else
  if (params.measurement_type == "2pcf-win") {
    trv::algo::compute_corrfunc_window(
      particles_rand, los_rand, params, binning.bin_centres, alpha, norm, norm_alt, save
    );
  } else

  if (params.measurement_type == "bispec") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      trv::algo::compute_bispec(
        particles_data, particles_rand, los_data, los_rand,
        params, binning.bin_centres, alpha, norm, norm_alt, save
      );
    } else
    if (params.catalogue_type == "sim") {
      trv::algo::compute_bispec_in_box(
        particles_data, params, binning.bin_centres, norm, norm_alt, save
      );
    }
  } else
  if (params.measurement_type == "3pcf") {
    if (params.catalogue_type == "survey" || params.catalogue_type == "mock") {
      trv::algo::compute_3pcf(
        particles_data, particles_rand, los_data, los_rand,
        params, binning.bin_centres, alpha, norm, norm_alt, save
      );
    } else
    if (params.catalogue_type == "sim") {
      trv::algo::compute_3pcf_in_box(
        particles_data, params, binning.bin_centres, norm, norm_alt, save
      );
    }
  } else
  if (params.measurement_type == "3pcf-win") {
    bool wa = false;
    trv::algo::compute_3pcf_window(
      particles_rand, los_rand, params, binning.bin_centres, alpha, norm, norm_alt, wa, save
    );
  } else
  if (params.measurement_type == "3pcf-win-wa") {
    bool wa = true;
    trv::algo::compute_3pcf_window(
      particles_rand, los_rand, params, binning.bin_centres, alpha, norm, norm_alt, wa, save
    );
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] ... made measurements.\n",
      trv::mon::show_timestamp().c_str()
    );
  }

  /* * Finalisation ******************************************************** */

  particles_data.finalise_particles();
  particles_rand.finalise_particles();

  delete[] los_data; los_data = NULL;
  delete[] los_rand; los_rand = NULL;

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Persistent memory usage: %.0f gigabytes.\n",
      trv::mon::show_timestamp().c_str(), trv::mon::gbytesMem
    );
  }

  if (trv::mon::currTask == 0) {
    std::printf(
      "[%s STAT] Program has completed.\n%s\n",
      trv::mon::show_timestamp().c_str(), std::string(80, '<').c_str()
    );
  }

  return 0;
}
