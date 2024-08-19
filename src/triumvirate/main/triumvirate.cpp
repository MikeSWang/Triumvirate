// Triumvirate: Three-Point Clustering Measurements in LSS
//
// Copyright (C) 2023 Mike S Wang & Naonori S Sugiyama [GPL-3.0]
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
 * @file triumvirate.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori S Sugiyama (https://github.com/naonori)
 * @brief Perform two- and three-point clustering statistic measurements.
 *
 */

#include <cstdio>
#include <stdexcept>
#include <string>

#include "monitor.hpp"
#include "parameters.hpp"
#include "particles.hpp"
#include "dataobjs.hpp"
#include "io.hpp"
#include "twopt.hpp"
#include "threept.hpp"

/**
 * @brief Set 'custom' binning.
 *
 * @todo This method is intended as a temporary alternative to
 *       trv::Binning::set_bins(std::vector) used for custom binning
 *       schemes specified by the external parameter file.
 *
 * @param _binning Binning.
 */
void _set_custom_bins(trv::Binning& _binning) {
  // Clear any existing binning.
  _binning.bin_edges.clear();
  _binning.bin_centres.clear();
  _binning.bin_widths.clear();

  // Implement 'custom' binning scheme for 'fourier' space.
  if (_binning.space == "fourier") {
    throw std::domain_error("No custom binning specified in Fourier space.");
  }

  // Implement 'custom' binning scheme for 'config' space.
  double dbin_pad = 10.;
  int nbin_pad = 5;

  _binning.bin_edges.push_back(0.);
  _binning.bin_centres.push_back(0.5);
  _binning.bin_widths.push_back(1.);

  _binning.bin_edges.push_back(1.);
  _binning.bin_centres.push_back(5.);
  _binning.bin_widths.push_back(9.);

  for (int ibin = 1; ibin < nbin_pad; ibin++) {
    double edge_left = dbin_pad * ibin;
    double centre = edge_left + dbin_pad / 2.;

    _binning.bin_edges.push_back(edge_left);
    _binning.bin_centres.push_back(centre);
    _binning.bin_widths.push_back(dbin_pad);
  }

  double bin_min = dbin_pad * nbin_pad;

  double dlnbin = (std::log(_binning.bin_max) - std::log(bin_min))
    / double(_binning.num_bins - nbin_pad);

  for (int ibin = nbin_pad; ibin < _binning.num_bins; ibin++) {
    double edge_left = bin_min * std::exp(dlnbin * (ibin - nbin_pad));
    double edge_right = bin_min * std::exp(dlnbin * (ibin - nbin_pad + 1));
    double centre = (edge_left + edge_right) / 2.;

    _binning.bin_edges.push_back(edge_left);
    _binning.bin_centres.push_back(centre);
    _binning.bin_widths.push_back(edge_right - edge_left);
  }

  _binning.bin_edges.push_back(_binning.bin_max);
}

/**
 * @brief A 'black-box' program for measuring two- and three-point
 *        clustering statistics.
 *
 * @param argc Number of command-line arguments.
 * @param argv Command-line arguments.
 * @returns Exit status.
 */
int main(int argc, char* argv[]) {

  for (int idx_arg = 1; idx_arg < argc; idx_arg++) {
    std::string arg = argv[idx_arg];
    if (arg == "-h" || arg == "--help") {
      trv::sys::display_help();
      return 0;
    }
    if (arg == "-V" || arg == "--version") {
      trv::sys::display_prog_logo();
      trv::sys::display_prog_licence();
      trv::sys::display_prog_info();
      return 0;
    }
    if (argc > 2 && arg.rfind("-", 0) == 0) {
      std::fprintf(stderr, "Unknown option: %s\n\n", arg.c_str());
      trv::sys::display_help();
      return 1;
    }
  }

#ifdef TRV_DISP
  if (trv::sys::currTask == 0) {
    trv::sys::display_prog_logo();
    trv::sys::display_prog_licence(true);
    trv::sys::display_prog_info(true);
  }
#endif  // TRV_DISP

  if (trv::sys::currTask == 0) {
    trv::sys::display_prog_logbars(0);
  }

  // =====================================================================
  // A Initialisation
  // =====================================================================

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat(
      "[MAIN:TRV:A] Parameters and source data are being initialised."
    );
  }

  // ---------------------------------------------------------------------
  // A.1 Parameter I/O
  // ---------------------------------------------------------------------

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat("[MAIN:TRV:A] Reading parameters...");
  }

  if (argc < 2) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.error(
        "Failed to initialise program: missing parameter file."
      );
      throw trv::sys::IOError(
        "Failed to initialise program: missing parameter file.\n"
      );
    }
  }

  trv::ParameterSet params;  // program parameters
  if (params.read_from_file(argv[1])) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.error(
        "Failed to initialise program: invalidated parameters."
      );
      throw trv::sys::IOError(
        "Failed to initialise program: invalidated parameters.\n"
      );
    }
  }

  trv::override_paramset_by_envvars(params);

  trv::sys::make_write_dir(params.measurement_dir);
  if (params.print_to_file()) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.warn(
        "Failed to print used parameters to file "
        "in the measurement output directory."
      );
    }
  }

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat("[MAIN:TRV:A] ... read parameters.");
  }

  trv::sys::logger.reset_level(params.verbose);

  // ---------------------------------------------------------------------
  // A.2 Data I/O
  // ---------------------------------------------------------------------

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.stat("[MAIN:TRV:A] Reading catalogues...");
    }
  }

  trv::ParticleCatalogue catalogue_data; // data-source catalogue
  std::string flag_data = "false";       // data-source catalogue status
  if (params.catalogue_type == "survey" || params.catalogue_type == "sim") {
    if (!(trv::sys::if_filepath_is_set(params.data_catalogue_file))) {
      if (trv::sys::currTask == 0) {
        trv::sys::logger.error(
          "Failed to initialise program: "
          "unspecified data-source catalogue file."
        );
        throw trv::sys::IOError(
          "Failed to initialise program: "
          "unspecified data-source catalogue file.\n"
        );
      }
    }
    if (catalogue_data.load_catalogue_file(
      params.data_catalogue_file, params.catalogue_columns, params.volume
    )) {
      if (trv::sys::currTask == 0) {
        trv::sys::logger.error(
          "Failed to initialise program: "
          "unloadable data-source catalogue file."
        );
        throw trv::sys::IOError(
          "Failed to initialise program: "
          "unloadable data-source catalogue file.\n"
        );
      }
    }
    flag_data = "true";
  }

  trv::ParticleCatalogue catalogue_rand; // random-source catalogue
  std::string flag_rand = "false";       // random-source catalogue status
  if (params.catalogue_type == "survey" || params.catalogue_type == "random") {
    if (!(trv::sys::if_filepath_is_set(params.rand_catalogue_file))) {
      if (trv::sys::currTask == 0) {
        trv::sys::logger.error(
          "Failed to initialise program: "
          "unspecified random-source catalogue file."
        );
        throw trv::sys::IOError(
          "Failed to initialise program: "
          "unspecified random-source catalogue file.\n"
        );
      }
    }
    if (catalogue_rand.load_catalogue_file(
      params.rand_catalogue_file, params.catalogue_columns, params.volume
    )) {
      if (trv::sys::currTask == 0) {
        trv::sys::logger.error(
          "Failed to initialise program: "
          "unloadable random-source catalogue file."
        );
        throw trv::sys::IOError(
          "Failed to initialise program: "
          "unloadable random-source catalogue file.\n"
        );
      }
    }
    flag_rand = "true";
  }

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.stat("[MAIN:TRV:A] ... read catalogues.");
    }
  }

  if (params.use_fftw_wisdom != "") {
    trv::sys::make_write_dir(params.use_fftw_wisdom);
  }

  // =====================================================================
  // B Measurements
  // =====================================================================

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat(
      "[MAIN:TRV:B] Clustering statistics are being measured."
    );
  }

  // ---------------------------------------------------------------------
  // B.1 Binning
  // ---------------------------------------------------------------------

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat("[MAIN:TRV:B] Setting up binning...");
  }

  trv::Binning binning(params);  // binning
  if (params.binning == "custom") {
    _set_custom_bins(binning);
  } else {
    binning.set_bins();
  }

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat("[MAIN:TRV:B] ... set up binning.");
  }

  // ---------------------------------------------------------------------
  // B.2 Line of sight
  // ---------------------------------------------------------------------

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.stat("[MAIN:TRV:B] Computing lines of sight...");
    }
  }

  trv::LineOfSight* los_data = nullptr;
  if (flag_data == "true") {
    // data-source LoS
    los_data = new trv::LineOfSight[catalogue_data.ntotal];
    trv::sys::gbytesMem +=
      trv::sys::size_in_gb<struct trv::LineOfSight>(catalogue_data.ntotal);
    trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int pid = 0; pid < catalogue_data.ntotal; pid++) {
      double los_mag =
        trv::maths::get_vec3d_magnitude(catalogue_data[pid].pos);

      if (los_mag == 0.) {
        trv::sys::logger.warn(
          "A data-catalogue particle coincides with the origin."
        );
        los_mag = 1.;
      }

      los_data[pid].pos[0] = catalogue_data[pid].pos[0] / los_mag;
      los_data[pid].pos[1] = catalogue_data[pid].pos[1] / los_mag;
      los_data[pid].pos[2] = catalogue_data[pid].pos[2] / los_mag;
    }
  }

  trv::LineOfSight* los_rand = nullptr;
  if (flag_rand == "true") {
    // random-source LoS
    los_rand = new trv::LineOfSight[catalogue_rand.ntotal];
    trv::sys::gbytesMem +=
      trv::sys::size_in_gb<struct trv::LineOfSight>(catalogue_rand.ntotal);
    trv::sys::update_maxmem();

#ifdef TRV_USE_OMP
#pragma omp parallel for
#endif  // TRV_USE_OMP
    for (int pid = 0; pid < catalogue_rand.ntotal; pid++) {
      double los_mag =
        trv::maths::get_vec3d_magnitude(catalogue_rand[pid].pos);

      if (los_mag == 0.) {
        trv::sys::logger.warn(
          "A random-catalogue particle coincides with the origin."
        );
        los_mag = 1.;
      }

      los_rand[pid].pos[0] = catalogue_rand[pid].pos[0] / los_mag;
      los_rand[pid].pos[1] = catalogue_rand[pid].pos[1] / los_mag;
      los_rand[pid].pos[2] = catalogue_rand[pid].pos[2] / los_mag;
    }
  }

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.stat("[MAIN:TRV:B] ... computed lines of sight.");
    }
  }

  // ---------------------------------------------------------------------
  // B.3 Box alignment
  // ---------------------------------------------------------------------

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.stat(
        "[MAIN:TRV:B] Aligning catalogues inside measurement box..."
      );
    }
  }

  if (params.catalogue_type == "survey") {
    if (params.alignment == "pad") {
      if (params.padscale == "grid") {
        double ngrid_pad[3] = {
          params.padfactor, params.padfactor, params.padfactor
        };
        trv::ParticleCatalogue::pad_grids(
          catalogue_data, catalogue_rand,
          params.boxsize, params.ngrid,
          ngrid_pad
        );
      } else
      if (params.padscale == "box") {
        double boxsize_pad[3] = {
          params.padfactor, params.padfactor, params.padfactor
        };
        trv::ParticleCatalogue::pad_in_box(
          catalogue_data, catalogue_rand,
          params.boxsize, boxsize_pad
        );
      }
    } else
    if (params.alignment == "centre") {
      trv::ParticleCatalogue::centre_in_box(
        catalogue_data, catalogue_rand, params.boxsize
      );
    }
  } else
  if (params.catalogue_type == "sim") {
    catalogue_data.offset_coords_for_periodicity(params.boxsize);
  } else
  if (params.catalogue_type == "random") {
    if (params.alignment == "pad") {
      if (params.padscale == "grid") {
        double ngrid_pad[3] = {
          params.padfactor, params.padfactor, params.padfactor
        };
        trv::ParticleCatalogue::pad_grids(
          catalogue_rand, params.boxsize, params.ngrid, ngrid_pad
        );
      } else
      if (params.padscale == "box") {
        double boxsize_pad[3] = {
          params.padfactor, params.padfactor, params.padfactor
        };
        trv::ParticleCatalogue::pad_in_box(
          catalogue_rand, params.boxsize, boxsize_pad
        );
      }
    } else
    if (params.alignment == "centre") {
      trv::ParticleCatalogue::centre_in_box(
        catalogue_rand, params.boxsize
      );
      catalogue_rand.offset_coords_for_periodicity(params.boxsize);
    }
  }

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.stat(
        "[MAIN:TRV:B] ... aligned catalogues inside measurement box."
      );
    }
  }

  // ---------------------------------------------------------------------
  // B.4 Constants
  // ---------------------------------------------------------------------

  double alpha;  // alpha contrast
  if (flag_data == "true" && flag_rand == "true") {
    alpha = catalogue_data.wstotal / catalogue_rand.wstotal;
  } else {
    alpha = 1.;
  }

  if (params.catalogue_type != "none") {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.info("Alpha contrast: %.6e.", alpha);
    }
  }

  trv::ParticleCatalogue& catalogue_for_norm =
    (flag_rand == "true") ? catalogue_rand : catalogue_data;
  double alpha_for_norm = (flag_rand == "true") ? alpha : 1.;
  double norm_factor_part = 0., norm_factor_mesh = 0., norm_factor_meshes = 0.;

  if (params.npoint == "2pt") {
    norm_factor_part = trv::calc_powspec_normalisation_from_particles(
      catalogue_for_norm, alpha_for_norm
    );
    norm_factor_mesh = trv::calc_powspec_normalisation_from_mesh(
      catalogue_for_norm, params, alpha_for_norm
    );
    // Mixed-mesh normalisation is only implemented for
    // paired survey-like catalogues.
    if (params.catalogue_type == "survey") {
      // Use default parameters for mixed-mesh normalisation in `pypower`.
      const double PADDING = 0.1;
      const double CELLSIZE = 10.;
      const std::string ASSIGNMENT = "cic";
      // Box size for normalisation is internally set and as such,
      // the current alignment of the catalogues is not applicable, but
      // this should have no effect on the normalisation.
      norm_factor_meshes = trv::calc_powspec_normalisation_from_meshes(
        catalogue_data, catalogue_rand, params, alpha,
        PADDING, CELLSIZE, ASSIGNMENT
      );
    }
  } else
  if (params.npoint == "3pt") {
    norm_factor_part = trv::calc_bispec_normalisation_from_particles(
      catalogue_for_norm, alpha_for_norm
    );
    norm_factor_mesh = trv::calc_bispec_normalisation_from_mesh(
      catalogue_for_norm, params, alpha_for_norm
    );
  }

  double norm_factor = 0.;
  if (params.npoint != "none") {
    if (params.norm_convention == "none") {
      norm_factor = 1.;
      if (trv::sys::currTask == 0) {
        trv::sys::logger.info(
          "Normalisation factors: "
          "%.6e (particle), %.6e (mesh), %.6e (mesh-mixed) (none used).",
          norm_factor_part, norm_factor_mesh, norm_factor_meshes
        );
      }
    } else
    if (params.norm_convention == "particle") {
      norm_factor = norm_factor_part;
      if (trv::sys::currTask == 0) {
        trv::sys::logger.info(
          "Normalisation factors: "
          "%.6e (particle; used), %.6e (mesh), %.6e (mesh-mixed).",
          norm_factor, norm_factor_mesh, norm_factor_meshes
        );
      }
    } else
    if (params.norm_convention == "mesh") {
      norm_factor = norm_factor_mesh;
      if (trv::sys::currTask == 0) {
        trv::sys::logger.info(
          "Normalisation factors: "
          "%.6e (particle), %.6e (mesh; used), %.6e (mesh-mixed).",
          norm_factor_part, norm_factor, norm_factor_meshes
        );
      }
    } else
    if (params.norm_convention == "mesh-mixed") {
      norm_factor = norm_factor_meshes;
      if (trv::sys::currTask == 0) {
        trv::sys::logger.info(
          "Normalisation factors: "
          "%.6e (particle), %.6e (mesh), %.6e (mesh-mixed; used).",
          norm_factor_part, norm_factor_mesh, norm_factor
        );
      }
    }
  }

  // ---------------------------------------------------------------------
  // B.5 Clustering algorithms
  // ---------------------------------------------------------------------

  char save_filepath[1024];
  if (params.statistic_type == "powspec") {
    std::snprintf(
      save_filepath, sizeof(save_filepath), "%s/pk%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );
    std::FILE* save_fileptr = nullptr;
    trv::PowspecMeasurements meas_powspec;  // power spectrum
    if (params.catalogue_type == "survey") {
      meas_powspec = trv::compute_powspec(
        catalogue_data, catalogue_rand, los_data, los_rand,
        params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data, catalogue_rand,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    } else
    if (params.catalogue_type == "sim") {
      meas_powspec = trv::compute_powspec_in_gpp_box(
        catalogue_data, params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    }
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_powspec
    );
    std::fclose(save_fileptr);
  } else
  if (params.statistic_type == "2pcf") {
    std::snprintf(
      save_filepath, sizeof(save_filepath), "%s/xi%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );
    std::FILE* save_fileptr = nullptr;
    trv::TwoPCFMeasurements meas_2pcf;  // two-point correlation function
    if (params.catalogue_type == "survey") {
      meas_2pcf = trv::compute_corrfunc(
        catalogue_data, catalogue_rand, los_data, los_rand,
        params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data, catalogue_rand,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    } else
    if (params.catalogue_type == "sim") {
      meas_2pcf = trv::compute_corrfunc_in_gpp_box(
        catalogue_data, params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    }
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_2pcf
    );
    std::fclose(save_fileptr);
  } else
  if (params.statistic_type == "2pcf-win") {
    std::snprintf(
      save_filepath, sizeof(save_filepath), "%s/xiw%d%s",
      params.measurement_dir.c_str(), params.ELL, params.output_tag.c_str()
    );

    trv::TwoPCFWindowMeasurements meas_2pcf_win = trv::compute_corrfunc_window(
      catalogue_rand, los_rand, params, binning, alpha, norm_factor
    );  // two-point correlation function window
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    trv::io::print_measurement_header_to_file(
      save_fileptr, params, catalogue_rand,
      norm_factor_part, norm_factor_mesh, norm_factor_meshes
    );
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_2pcf_win
    );
    std::fclose(save_fileptr);
  } else
  if (params.statistic_type == "bispec") {
    if (params.form == "full" || params.form == "diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/bk%d%d%d_%s%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.form.c_str(),
        params.output_tag.c_str()
      );
    } else
    if (params.form == "off-diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/bk%d%d%d_offdiag%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.idx_bin,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "row") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/bk%d%d%d_row%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.idx_bin,
        params.output_tag.c_str()
      );
    }
    std::FILE* save_fileptr = nullptr;
    trv::BispecMeasurements meas_bispec;  // bispectrum
    if (params.catalogue_type == "survey") {
      meas_bispec = trv::compute_bispec(
        catalogue_data, catalogue_rand, los_data, los_rand,
        params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data, catalogue_rand,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    } else
    if (params.catalogue_type == "sim") {
      meas_bispec = trv::compute_bispec_in_gpp_box(
        catalogue_data, params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    }
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_bispec
    );
    std::fclose(save_fileptr);
  } else
  if (params.statistic_type == "3pcf") {
    if (params.form == "full" || params.form == "diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zeta%d%d%d_%s%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.form.c_str(),
        params.output_tag.c_str()
      );
    } else
    if (params.form == "off-diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zeta%d%d%d_offdiag%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.idx_bin,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "row") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zeta%d%d%d_row%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.idx_bin,
        params.output_tag.c_str()
      );
    }
    std::FILE* save_fileptr = nullptr;
    trv::ThreePCFMeasurements meas_3pcf;  // three-point correlation function
    if (params.catalogue_type == "survey") {
      meas_3pcf = trv::compute_3pcf(
        catalogue_data, catalogue_rand, los_data, los_rand,
        params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data, catalogue_rand,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    } else
    if (params.catalogue_type == "sim") {
      meas_3pcf = trv::compute_3pcf_in_gpp_box(
        catalogue_data, params, binning, norm_factor
      );
      save_fileptr = std::fopen(save_filepath, "w");
      trv::io::print_measurement_header_to_file(
        save_fileptr, params, catalogue_data,
        norm_factor_part, norm_factor_mesh, norm_factor_meshes
      );
    }
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_3pcf
    );
    std::fclose(save_fileptr);
  } else
  if (params.statistic_type == "3pcf-win") {
    if (params.form == "full" || params.form == "diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zetaw%d%d%d_%s%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL,
        params.form.c_str(),
        params.output_tag.c_str()
      );
    } else
    if (params.form == "off-diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zetaw%d%d%d_offdiag%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.idx_bin,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "row") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zetaw%d%d%d_row%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.idx_bin,
        params.output_tag.c_str()
      );
    }
    bool wa = false;

    trv::ThreePCFWindowMeasurements meas_3pcf_win = trv::compute_3pcf_window(
      catalogue_rand, los_rand, params, binning, alpha, norm_factor, wa
    );  // three-point correlation function window
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    trv::io::print_measurement_header_to_file(
      save_fileptr, params, catalogue_rand,
      norm_factor_part, norm_factor_mesh, norm_factor_meshes
    );
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_3pcf_win
    );
    std::fclose(save_fileptr);
  } else
  if (params.statistic_type == "3pcf-win-wa") {
    if (params.form == "full" || params.form == "diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zetaw%d%d%d_wa%d%d_%s%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.i_wa, params.j_wa,
        params.form.c_str(),
        params.output_tag.c_str()
      );
    } else
    if (params.form == "off-diag") {
      std::snprintf(
        save_filepath, sizeof(save_filepath),
        "%s/zetaw%d%d%d_wa%d%d_offdiag%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.i_wa, params.j_wa,
        params.idx_bin,
        params.output_tag.c_str()
      );
    } else
    if (params.form == "row") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s/zetaw%d%d%d_wa%d%d_row%d%s",
        params.measurement_dir.c_str(),
        params.ell1, params.ell2, params.ELL, params.i_wa, params.j_wa,
        params.idx_bin,
        params.output_tag.c_str()
      );
    }

    bool wa = true;

    trv::ThreePCFWindowMeasurements meas_3pcf_win_wa =
      trv::compute_3pcf_window(
        catalogue_rand, los_rand, params, binning, alpha, norm_factor, wa
      );  // three-point correlation function window wide-angle corrections
    std::FILE* save_fileptr = std::fopen(save_filepath, "w");
    trv::io::print_measurement_header_to_file(
      save_fileptr, params, catalogue_rand,
      norm_factor_part, norm_factor_mesh, norm_factor_meshes
    );
    trv::io::print_measurement_datatab_to_file(
      save_fileptr, params, meas_3pcf_win_wa
    );
    std::fclose(save_fileptr);
  }

  if (params.save_binned_vectors != "") {
    trv::FieldStats binning_meshgrid(params, false);
    trv::BinnedVectors binned_vectors = binning_meshgrid.record_binned_vectors(
      binning, params.save_binned_vectors
    );
    if (params.statistic_type == "modes" || params.statistic_type == "pairs") {
      std::snprintf(
        save_filepath, sizeof(save_filepath), "%s",
        params.save_binned_vectors.c_str()
      );
    }
  }

  if (trv::sys::currTask == 0) {
    trv::sys::logger.info("Measurements saved to: %s", save_filepath);
  }

  // =====================================================================
  // C Finalisation
  // =====================================================================

  if (trv::sys::currTask == 0) {
    trv::sys::logger.stat("[MAIN:TRV:C] Data objects are being cleared.");
  }

  // Clear persistent and dynamic memory.
#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
  fftw_cleanup_threads();
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
  fftw_cleanup();
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

  catalogue_data.finalise_particles();
  catalogue_rand.finalise_particles();

  if (los_data != nullptr) {
    delete[] los_data; los_data = nullptr;
  }
  if (los_rand != nullptr) {
    delete[] los_rand; los_rand = nullptr;
  }
  trv::sys::gbytesMem -= trv::sys::size_in_gb<struct trv::LineOfSight>(
    (catalogue_data.ntotal + catalogue_rand.ntotal)
  );

  if (trv::sys::count_fft > 0 || trv::sys::count_ifft > 0) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.info(
        "Number of FFTs: %d forward, %d backward.",
        trv::sys::count_fft, trv::sys::count_ifft
      );
    }
  }

  if (
    trv::sys::max_count_grid > 0 ||
    trv::sys::max_count_cgrid > 0 ||
    trv::sys::max_count_rgrid > 0
  ) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.info(
        "Maximum number of concurrent 3-d grids: "
        "%.1f complex-equivalent, %d complex, %d real.",
        trv::sys::max_count_grid,
        trv::sys::max_count_cgrid, trv::sys::max_count_rgrid
      );
    }
  }

  if (trv::sys::currTask == 0) {
    trv::sys::logger.info(
      "Minimal estimate of peak memory usage: %.1f gibibytes.",
      trv::sys::gbytesMaxMem
    );
  }
  if (trv::sys::gbytesMem > 0.) {
    if (trv::sys::currTask == 0) {
      trv::sys::logger.warn(
        "Uncleared dynamically allocated memory: %.1f gibibytes.",
        trv::sys::gbytesMem
      );
    }
  }

  if (trv::sys::currTask == 0) {
    trv::sys::display_prog_logbars(1);
  }

  return 0;
}
