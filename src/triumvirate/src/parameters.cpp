// Triumvirate: Three-Point Clustering Measurements in LSS
//
// Copyright (C) 2023 Mike S Wang & Naonori S Sugiyama [GPL-3.0-or-later]
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
 * @file parameters.cpp
 * @authors Mike S Wang (https://github.com/MikeSWang),
 *          Naonori S Sugiyama (https://github.com/naonori)
 *
 */

#include "parameters.hpp"

namespace trvs = trv::sys;

namespace trv {

ParameterSet::ParameterSet(const ParameterSet& other) {
  // Copy I/O parameters.
  this->catalogue_dir = other.catalogue_dir;
  this->measurement_dir = other.measurement_dir;
  this->data_catalogue_file = other.data_catalogue_file;
  this->rand_catalogue_file = other.rand_catalogue_file;
  this->catalogue_columns = other.catalogue_columns;
  this->catalogue_dataset = other.catalogue_dataset;
  this->output_tag = other.output_tag;

  // Copy mesh sampling parameters.
  for (int i = 0; i < 3; i++) {
    this->boxsize[i] = other.boxsize[i];
    this->ngrid[i] = other.ngrid[i];
  }
  this->expand = other.expand;
  this->alignment = other.alignment;
  this->padscale = other.padscale;
  this->padfactor = other.padfactor;
  this->assignment = other.assignment;
  this->interlace = other.interlace;
  this->volume = other.volume;
  this->nmesh = other.nmesh;
  this->assignment_order = other.assignment_order;

  // Copy measurement parameters.
  this->catalogue_type = other.catalogue_type;
  this->statistic_type = other.statistic_type;
  this->npoint = other.npoint;
  this->space = other.space;
  this->ell1 = other.ell1;
  this->ell2 = other.ell2;
  this->ELL = other.ELL;
  this->i_wa = other.i_wa;
  this->j_wa = other.j_wa;
  this->form = other.form;
  this->norm_convention = other.norm_convention;
  this->binning = other.binning;
  this->bin_min = other.bin_min;
  this->bin_max = other.bin_max;
  this->num_bins = other.num_bins;
  this->idx_bin = other.idx_bin;
  this->cutoff_nyq = other.cutoff_nyq;

  // Copy misc parameters.
  this->fftw_scheme = other.fftw_scheme;
  this->fftw_planner_flag = other.fftw_planner_flag;
  this->use_fftw_wisdom = other.use_fftw_wisdom;
  this->fftw_wisdom_file_f = other.fftw_wisdom_file_f;
  this->fftw_wisdom_file_b = other.fftw_wisdom_file_b;
  this->save_binned_vectors = other.save_binned_vectors;
  this->verbose = other.verbose;
  this->progbar = other.progbar;

  // Copy private parameters.
  this->data_catalogue_files = other.data_catalogue_files;
  this->rand_catalogue_files = other.rand_catalogue_files;
}

int ParameterSet::read_from_file(char* parameter_filepath) {
  // ---------------------------------------------------------------------
  // Initialisation
  // ---------------------------------------------------------------------

  // Load parameter file.
  std::string param_filepath = parameter_filepath;

  std::ifstream fin(param_filepath.c_str());

  // Initialise temporary variables to hold the extracted parameters.
  // NOTE: Assume file systems allow path character lengths up to 4096,
  // and up to 64 files of this maximum lengths are stored in the
  // char arrays.
  char catalogue_dir_[1024] = "";
  char measurement_dir_[1024] = "";
  char data_catalogue_file_[262144] = "";
  char rand_catalogue_file_[262144] = "";
  char catalogue_columns_[1024] = "";
  char catalogue_dataset_[1024] = "";
  char output_tag_[1024] = "";

  double boxsize_x, boxsize_y, boxsize_z;
  int ngrid_x, ngrid_y, ngrid_z;

  char alignment_[16] = "";
  char padscale_[16] = "";
  char assignment_[16] = "";
  char interlace_[16] = "";

  char catalogue_type_[16] = "";
  char statistic_type_[16] = "";
  char form_[16] = "";
  char norm_convention_[16] = "";
  char binning_[16] = "";

  char fftw_scheme_[16] = "";
  char use_fftw_wisdom_[1024] = "";
  char save_binned_vectors_[1024] = "";
  char progbar_[16] = "";

  // ---------------------------------------------------------------------
  // Extraction
  // ---------------------------------------------------------------------

  std::string line_str;
  char dummy_str[1024], dummy_equal[1024];
  while (std::getline(fin, line_str)) {
    // Check if the line is a parameter assignment.
    if (line_str.find("#") == 0) {
      continue;
    }  // skip comment lines
    if (
      std::sscanf(
        line_str.data(), "%1023s %1023s %1023s",
        dummy_str, dummy_equal, dummy_str
      ) != 3
    ) {
      continue;
    }  // skip non-relation lines
    if (std::strcmp(dummy_equal, "=") != 0) {
      continue;
    }  // skip non-assignment lines

    // Define convenience function for scanning string parameters.
    auto scan_par_str = [line_str, dummy_str, dummy_equal](
      const char* par_name, const char* fmt, const char* par_value
    ) {
      if (line_str.find(par_name) != std::string::npos) {
        std::sscanf(
          line_str.data(), fmt, dummy_str, dummy_equal, par_value
        );
      }
    };

    // -- I/O ------------------------------------------------------------

    scan_par_str("catalogue_dir", "%s %s %s", catalogue_dir_);
    scan_par_str("measurement_dir", "%s %s %s", measurement_dir_);
    scan_par_str("data_catalogue_file", "%s %s %s", data_catalogue_file_);
    scan_par_str("rand_catalogue_file", "%s %s %s", rand_catalogue_file_);
    scan_par_str("catalogue_columns", "%s %s %s", catalogue_columns_);
    scan_par_str("catalogue_dataset", "%s %s %s", catalogue_dataset_);
    scan_par_str("output_tag", "%s %s %s", output_tag_);

    // -- Mesh sampling --------------------------------------------------

    if (line_str.find("boxsize_x") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &boxsize_x
      );
    }
    if (line_str.find("boxsize_y") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &boxsize_y
      );
    }
    if (line_str.find("boxsize_z") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &boxsize_z
      );
    }

    if (line_str.find("ngrid_x") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &ngrid_x
      );
    }
    if (line_str.find("ngrid_y") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &ngrid_y
      );
    }
    if (line_str.find("ngrid_z") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d", dummy_str, dummy_equal, &ngrid_z
      );
    }

    if (line_str.find("expand") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &this->expand
      );
    }

    scan_par_str("alignment", "%1023s %1023s %1023s", alignment_);
    scan_par_str("padscale", "%1023s %1023s %1023s", padscale_);

    if (line_str.find("padfactor") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &this->padfactor
      );
    }

    scan_par_str("assignment", "%1023s %1023s %1023s", assignment_);
    scan_par_str("interlace", "%1023s %1023s %1023s", interlace_);

    // -- Measurement ----------------------------------------------------

    scan_par_str("catalogue_type", "%1023s %1023s %1023s", catalogue_type_);
    scan_par_str("statistic_type", "%1023s %1023s %1023s", statistic_type_);
    scan_par_str("form", "%1023s %1023s %1023s", form_);
    scan_par_str("norm_convention", "%1023s %1023s %1023s", norm_convention_);
    scan_par_str("binning", "%1023s %1023s %1023s", binning_);

    if (line_str.find("ell1") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->ell1
      );
    }
    if (line_str.find("ell2") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->ell2
      );
    }
    if (line_str.find("ELL") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->ELL
      );
    }

    if (line_str.find("i_wa") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->i_wa
      );
    }
    if (line_str.find("j_wa") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->j_wa
      );
    }

    if (line_str.find("bin_min") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &this->bin_min
      );
    }
    if (line_str.find("bin_max") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &this->bin_max
      );
    }

    if (line_str.find("num_bins") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->num_bins
      );
    }
    if (line_str.find("idx_bin") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->idx_bin
      );
    }

    if (line_str.find("cutoff_nyq") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %lg",
        dummy_str, dummy_equal, &this->cutoff_nyq
      );
    }

    // -- Misc -------------------------------------------------------------

    scan_par_str("fftw_scheme", "%1023s %1023s %1023s", fftw_scheme_);
    scan_par_str("use_fftw_wisdom", "%1023s %1023s %1023s", use_fftw_wisdom_);
    scan_par_str(
      "save_binned_vectors", "%1023s %1023s %1023s", save_binned_vectors_
    );
    scan_par_str("progbar", "%1023s %1023s %1023s", progbar_);

    if (line_str.find("verbose") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%1023s %1023s %d",
        dummy_str, dummy_equal, &this->verbose
      );
    }
  }

  // ---------------------------------------------------------------------
  // Attribution
  // ---------------------------------------------------------------------

  // Attribute numerical parameters (directly extracted above).

  // Attribute string parameters.
  this->catalogue_dir = catalogue_dir_;
  this->measurement_dir = measurement_dir_;
  this->data_catalogue_file = data_catalogue_file_;
  this->rand_catalogue_file = rand_catalogue_file_;
  this->catalogue_columns = catalogue_columns_;
  this->catalogue_dataset = catalogue_dataset_;
  this->output_tag = output_tag_;

  if (strlen(alignment_) > 0) {
    this->alignment = alignment_;
  }
  if (strlen(padscale_) > 0) {
    this->padscale = padscale_;
  }
  if (strlen(assignment_) > 0) {
    this->assignment = assignment_;
  }
  if (strlen(interlace_) > 0) {
    this->interlace = interlace_;
  }

  this->catalogue_type = catalogue_type_;
  this->statistic_type = statistic_type_;
  if (strlen(form_) > 0) {
    this->form = form_;
  }
  if (strlen(norm_convention_) > 0) {
    this->norm_convention = norm_convention_;
  }
  if (strlen(binning_) > 0) {
    this->binning = binning_;
  }

  if (strlen(fftw_scheme_) > 0) {
    this->fftw_scheme = fftw_scheme_;
  }
  if (strlen(use_fftw_wisdom_) > 0) {
    this->use_fftw_wisdom = use_fftw_wisdom_;
  }
  if (strlen(save_binned_vectors_) > 0) {
    this->save_binned_vectors = save_binned_vectors_;
  }
  if (strlen(progbar_) > 0) {
    this->progbar = progbar_;
  }

  // Attribute derived parameters.
  this->boxsize[0] = boxsize_x;
  this->boxsize[1] = boxsize_y;
  this->boxsize[2] = boxsize_z;

  this->ngrid[0] = ngrid_x;
  this->ngrid[1] = ngrid_y;
  this->ngrid[2] = ngrid_z;

  this->volume = boxsize_x * boxsize_y * boxsize_z;
  this->nmesh = ngrid_x * ngrid_y * ngrid_z;

  // ---------------------------------------------------------------------
  // Debugging mode
  // ---------------------------------------------------------------------

#ifdef DBG_PARS
  // Define convenience function for displaying debugged parameters.
  auto debug_par_str = [](const std::string& name, const std::string& value) {
    std::cout << name << ": " << value << std::endl;
  };
  auto debug_par_int = [](const std::string& name, int value) {
    std::cout << name << ": " << value << std::endl;
  };
  auto debug_par_longlong = [](const std::string& name, long long value) {
    std::cout << name << ": " << value << std::endl;
  };
  auto debug_par_double = [](const std::string& name, double value) {
    std::cout << name << ": " << value << std::endl;
  };

  // Display debugged parameters.
  debug_par_str("catalogue_dir", this->catalogue_dir);
  debug_par_str("measurement_dir", this->measurement_dir);
  debug_par_str("data_catalogue_file", this->data_catalogue_file);
  debug_par_str("rand_catalogue_file", this->rand_catalogue_file);
  debug_par_str("catalogue_columns", this->catalogue_columns);
  debug_par_str("catalogue_dataset", this->catalogue_dataset);
  debug_par_str("output_tag", this->output_tag);

  debug_par_str("alignment", this->alignment);
  debug_par_str("padscale", this->padscale);
  debug_par_str("assignment", this->assignment);
  debug_par_str("interlace", this->interlace);

  debug_par_str("catalogue_type", this->catalogue_type);
  debug_par_str("statistic_type", this->statistic_type);
  debug_par_str("form", this->form);
  debug_par_str("norm_convention", this->norm_convention);
  debug_par_str("binning", this->binning);

  debug_par_str("fftw_scheme", this->fftw_scheme);
  debug_par_str("use_fftw_wisdom", this->use_fftw_wisdom);
  debug_par_str("save_binned_vectors", this->save_binned_vectors);
  debug_par_str("progbar", this->progbar);

  debug_par_int("ngrid[0]", this->ngrid[0]);
  debug_par_int("ngrid[1]", this->ngrid[1]);
  debug_par_int("ngrid[2]", this->ngrid[2]);

  debug_par_longlong("nmesh", this->nmesh);

  debug_par_int("ell1", this->ell1);
  debug_par_int("ell2", this->ell2);
  debug_par_int("ELL", this->ELL);
  debug_par_int("i_wa", this->i_wa);
  debug_par_int("j_wa", this->j_wa);

  debug_par_int("num_bins", this->num_bins);
  debug_par_int("idx_bin", this->idx_bin);

  debug_par_double("boxsize[0]", this->boxsize[0]);
  debug_par_double("boxsize[1]", this->boxsize[1]);
  debug_par_double("boxsize[2]", this->boxsize[2]);
  debug_par_double("expand", this->expand);
  debug_par_double("volume", this->volume);
  debug_par_double("padfactor", this->padfactor);
  debug_par_double("bin_min", this->bin_min);
  debug_par_double("bin_max", this->bin_max);
  debug_par_double("cutoff_nyq", this->cutoff_nyq);
#endif  // DBG_PARS

  return this->validate(true);
}

int ParameterSet::validate(bool init) {
  trvs::logger.reset_level(this->verbose);

  // Validate and derive string parameters.
  // Any duplicate '/' has no effect.
  if (
    this->catalogue_dir.find_first_not_of(" \t\n\r\v\f") != std::string::npos
    && init
  ) {
    this->catalogue_dir += "/";  // transmutation
  }
  if (
    this->measurement_dir.find_first_not_of(" \t\n\r\v\f") == std::string::npos
  ) {
    this->measurement_dir = "./";  // transmutation
  } else
  if (init) {
    this->measurement_dir += "/";  // transmutation
  }
  trvs::expand_envar_in_path(this->catalogue_dir);
  trvs::expand_envar_in_path(this->measurement_dir);

  // Parse catalogue filenames by delimiter.
  if (
    this->data_catalogue_file.find_first_not_of(" \t\n\r\v\f")
    != std::string::npos
  ) {
    if (!trvs::has_extension(this->data_catalogue_file, trvs::fn_delimiter)) {
      this->data_catalogue_file += trvs::fn_delimiter;  // transmutation
    }
  }
  if (
    this->rand_catalogue_file.find_first_not_of(" \t\n\r\v\f")
    != std::string::npos
  ) {
    if (!trvs::has_extension(this->rand_catalogue_file, trvs::fn_delimiter)) {
      this->rand_catalogue_file += trvs::fn_delimiter;  // transmutation
    }
  }

  this->data_catalogue_files.clear();
  this->rand_catalogue_files.clear();

  std::string data_ctlg_files, rand_ctlg_files;
  std::string data_ctlg_file_, rand_ctlg_file_;
  size_t dlpos = 0;

  if (this->catalogue_type == "survey") {
    if (this->data_catalogue_file != "") {
      data_ctlg_files = this->data_catalogue_file;
      dlpos = 0;
      while (
        (dlpos = data_ctlg_files.find(trvs::fn_delimiter)) != std::string::npos
      ) {
        data_ctlg_file_ = data_ctlg_files.substr(0, dlpos);
        data_ctlg_files.erase(0, dlpos + trvs::fn_delimiter.length());
        if (data_ctlg_file_.rfind("/", 0) != 0 && init) {
          data_ctlg_file_ = this->catalogue_dir + data_ctlg_file_;
        }  // transmutation
        trvs::expand_envar_in_path(data_ctlg_file_);
        this->data_catalogue_files.push_back(data_ctlg_file_);
      }
      // if (this->data_catalogue_files.empty()) {
      //   if (this->data_catalogue_file.rfind("/", 0) != 0 && init) {
      //     this->data_catalogue_file =
      //       this->catalogue_dir + this->data_catalogue_file;
      //   }  // transmutation
      //   this->data_catalogue_files.push_back(this->data_catalogue_file);
      // }
    }
    if (this->rand_catalogue_file != "") {
      rand_ctlg_files = this->rand_catalogue_file;
      dlpos = 0;
      while (
        (dlpos = rand_ctlg_files.find(trvs::fn_delimiter)) != std::string::npos
      ) {
        rand_ctlg_file_ = rand_ctlg_files.substr(0, dlpos);
        rand_ctlg_files.erase(0, dlpos + trvs::fn_delimiter.length());
        if (rand_ctlg_file_.rfind("/", 0) != 0 && init) {
          rand_ctlg_file_ = this->catalogue_dir + rand_ctlg_file_;
        }  // transmutation
        trvs::expand_envar_in_path(rand_ctlg_file_);
        this->rand_catalogue_files.push_back(rand_ctlg_file_);
      }
      // if (this->rand_catalogue_files.empty()) {
      //   if (this->rand_catalogue_file.rfind("/", 0) != 0 && init) {
      //     this->rand_catalogue_file =
      //       this->catalogue_dir + this->rand_catalogue_file;
      //   }  // transmutation
      //   this->rand_catalogue_files.push_back(this->rand_catalogue_file);
      // }
    }
  } else
  if (this->catalogue_type == "random") {
    this->data_catalogue_file = "";  // transmutation
    if (this->rand_catalogue_file != "") {
      rand_ctlg_files = this->rand_catalogue_file;
      dlpos = 0;
      while (
        (dlpos = rand_ctlg_files.find(trvs::fn_delimiter)) != std::string::npos
      ) {
        rand_ctlg_file_ = rand_ctlg_files.substr(0, dlpos);
        rand_ctlg_files.erase(0, dlpos + trvs::fn_delimiter.length());
        if (rand_ctlg_file_.rfind("/", 0) != 0 && init) {
          rand_ctlg_file_ = this->catalogue_dir + rand_ctlg_file_;
        }  // transmutation
        trvs::expand_envar_in_path(rand_ctlg_file_);
        this->rand_catalogue_files.push_back(rand_ctlg_file_);
      }
      // if (this->rand_catalogue_files.empty()) {
      //   if (this->rand_catalogue_file.rfind("/", 0) != 0 && init) {
      //     this->rand_catalogue_file =
      //       this->catalogue_dir + this->rand_catalogue_file;
      //   }  // transmutation
      //   this->rand_catalogue_files.push_back(this->rand_catalogue_file);
      // }
    }
  } else
  if (this->catalogue_type == "sim") {
    if (this->data_catalogue_file != "") {
      data_ctlg_files = this->data_catalogue_file;
      dlpos = 0;
      while (
        (dlpos = data_ctlg_files.find(trvs::fn_delimiter)) != std::string::npos
      ) {
        data_ctlg_file_ = data_ctlg_files.substr(0, dlpos);
        data_ctlg_files.erase(0, dlpos + trvs::fn_delimiter.length());
        if (data_ctlg_file_.rfind("/", 0) != 0 && init) {
          data_ctlg_file_ = this->catalogue_dir + data_ctlg_file_;
        }  // transmutation
        trvs::expand_envar_in_path(data_ctlg_file_);
        this->data_catalogue_files.push_back(data_ctlg_file_);
      }
      // if (this->data_catalogue_files.empty()) {
      //   if (this->data_catalogue_file.rfind("/", 0) != 0 && init) {
      //     this->data_catalogue_file =
      //       this->catalogue_dir + this->data_catalogue_file;
      //   }  // transmutation
      //   this->data_catalogue_files.push_back(this->data_catalogue_file);
      // }
    }
    this->rand_catalogue_file = "";  // transmutation
  } else
  if (this->catalogue_type == "none") {
    // Nothing should happen.
  } else {
#ifndef TRV_EXTCALL
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Catalogue type must be 'survey', 'random', 'sim' or 'none: "
        "`catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Catalogue type must be 'survey', 'random', 'sim' or 'none': "
      "`catalogue_type` = '%s'.",
      this->catalogue_type.c_str()
    );
#endif  // !TRV_EXTCALL
  }
  this->data_catalogue_file = trvs::join_strings(
    this->data_catalogue_files, trvs::fn_delimiter
  );
  this->rand_catalogue_file = trvs::join_strings(
    this->rand_catalogue_files, trvs::fn_delimiter
  );

  if (!(this->alignment == "centre" || this->alignment == "pad")) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Box alignment must be 'centre' or 'pad': `alignment` = '%s'.",
        this->alignment.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Box alignment must be 'centre' or 'pad': `alignment` = '%s'.",
      this->alignment.c_str()
    );
  }
  if (!(this->padscale == "box" || this->padscale == "grid")) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Pad scale must be 'box' or 'grid': `padscale` = '%s'.",
        this->padscale.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Pad scale must be 'box' or 'grid': `padscale` = '%s'.",
      this->padscale.c_str()
    );
  }

  if (this->assignment == "ngp") {
    this->assignment_order = 1;
  } else
  if (this->assignment == "cic") {
    this->assignment_order = 2;
  } else
  if (this->assignment == "tsc") {
    this->assignment_order = 3;
  } else
  if (this->assignment == "pcs") {
    this->assignment_order = 4;
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Mesh assignment scheme must be "
        "'ngp', 'cic', 'tsc' or 'pcs': `assignment` = '%s'.",
        this->assignment.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Mesh assignment scheme must be "
      "'ngp', 'cic', 'tsc' or 'pcs': `assignment` = '%s'.",
      this->assignment.c_str()
    );
  }
  if (this->interlace == "true" || this->interlace == "on") {
    this->interlace = "true";  // transmutation
  } else
  if (this->interlace == "false" || this->interlace == "off") {
    this->interlace = "false";  // transmutation
  } else {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Interlacing must be 'true'/'on' or 'false'/'off': "
        "`interlace` = '%s'.",
        this->interlace.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Interlacing must be 'true'/'on' or 'false'/'off': "
      "`interlace` = '%s'.",
      this->interlace.c_str()
    );
  }

  if (this->statistic_type == "powspec") {
    this->npoint = "2pt"; this->space = "fourier";  // derivation
    if (this->catalogue_type == "random" || this->catalogue_type == "none") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Power spectrum requires 'sim' or 'survey' catalogue(s): "
          "`catalogue_type` = '%s'.",
          this->catalogue_type.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "Power spectrum requires 'sim' or 'survey' catalogue(s): "
        "`catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
  } else
  if (this->statistic_type == "2pcf") {
    this->npoint = "2pt"; this->space = "config";  // derivation
    if (this->catalogue_type == "random" || this->catalogue_type == "none") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Two-point correlation function requires 'sim' or 'survey' "
          "catalogue(s): `catalogue_type` = '%s'.",
          this->catalogue_type.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "Two-point correlation function requires 'sim' or 'survey' "
        "catalogue(s): `catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
  } else
  if (this->statistic_type == "2pcf-win") {
    this->npoint = "2pt"; this->space = "config";  // derivation
    if (this->catalogue_type != "random") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Two-point correlation function window requires 'random' catalogue: "
          "`catalogue_type` = '%s'.",
          this->catalogue_type.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "Two-point correlation function window requires 'random' catalogue: "
        "`catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
  } else
  if (this->statistic_type == "bispec") {
    this->npoint = "3pt"; this->space = "fourier";  // derivation
    if (this->catalogue_type == "random" || this->catalogue_type == "none") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Bispectrum requires 'sim' or 'survey' catalogue(s): "
          "`catalogue_type` = '%s'.",
          this->catalogue_type.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "Bispectrum requires 'sim' or 'survey' catalogue(s): "
        "`catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
  } else
  if (this->statistic_type == "3pcf") {
    this->npoint = "3pt"; this->space = "config";  // derivation
    if (this->catalogue_type == "random" || this->catalogue_type == "none") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Three-point correlation function requires 'sim' or 'survey' "
          "catalogue(s): `catalogue_type` = '%s'.",
          this->catalogue_type.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "Three-point correlation function requires 'sim' or 'survey' "
        "catalogue(s): `catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
  } else
  if (
    this->statistic_type == "3pcf-win"
    || this->statistic_type == "3pcf-win-wa"
  ) {
    this->npoint = "3pt"; this->space = "config";  // derivation
    if (this->catalogue_type != "random") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Three-point correlation function window requires 'random' "
          "catalogue: `catalogue_type` = '%s'.",
          this->catalogue_type.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "Three-point correlation function window requires 'random' "
        "catalogue: `catalogue_type` = '%s'.",
        this->catalogue_type.c_str()
      );
    }
  } else
  if (this->statistic_type == "modes") {
    this->npoint = "none"; this->space = "fourier";  // derivation
  } else
  if (this->statistic_type == "pairs") {
    this->npoint = "none"; this->space = "config";  // derivation
  } else {
#ifndef TRV_EXTCALL
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Statistic type is not recognised: `statistic_type` = '%s'.",
        this->statistic_type.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Statistic type is not recognised: `statistic_type` = '%s'.",
      this->statistic_type.c_str()
    );
#endif  // !TRV_EXTCALL
  }
  if (!(
    this->form == "full"
    || this->form == "diag"
    || this->form == "off-diag"
    || this->form == "row"
  )) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Three-point statistic form is not recognised: `form` = '%s'.",
        this->form.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Three-point statistic form is not recognised: `form` = '%s'.",
      this->form.c_str()
    );
  } else {
    if (this->form == "full" && this->ell1 == this->ell2) {
      this->shape = "triu";  // derivation
    } else {
      this->shape = this->form;  // derivation
    }
  }
  if (!(
    this->norm_convention == "none"
    || this->norm_convention == "particle"
    || this->norm_convention == "mesh"
    || this->norm_convention == "mesh-mixed"
  )) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Normalisation convention must be 'mesh', 'particle', "
        "'mesh' or 'mesh-mixed': `norm_convention` = '%s'.",
        this->norm_convention.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Normalisation convention must be 'mesh', 'particle', "
      "'mesh' or 'mesh-mixed': `norm_convention` = '%s'.",
      this->norm_convention.c_str()
    );
  }
  if (this->norm_convention == "mesh-mixed" && this->npoint != "2pt") {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Normalisation convention 'mesh-mixed' only applies to "
        "two-point statistics: `npoint` = '%s'.",
        this->npoint.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Normalisation convention 'mesh-mixed' only applies to "
      "two-point statistics: `npoint` = '%s'.",
      this->npoint.c_str()
    );
  }
  if (!(
    this->binning == "lin"
    || this->binning == "log"
    || this->binning == "linpad"
    || this->binning == "logpad"
    || this->binning == "custom"
  )) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Binning scheme is unrecognised: `binning` = '%s'.",
        this->binning.c_str()
      );
    }
    throw trvs::InvalidParameterError(
      "Binning scheme is unrecognised: `binning` = '%s'.",
      this->binning.c_str()
    );
  }

  if (!trvs::is_gpu_enabled()) {
    if (this->fftw_scheme == "estimate") {
      this->fftw_planner_flag = FFTW_ESTIMATE;  // derivation
    } else
    if (this->fftw_scheme == "measure") {
        this->fftw_planner_flag = FFTW_MEASURE;  // derivation
    } else
    if (this->fftw_scheme == "patient") {
      this->fftw_planner_flag = FFTW_PATIENT;  // derivation
    } else {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "FFTW planner scheme is not supported: `fftw_scheme` = '%s'.",
          this->fftw_scheme.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "FFTW planner scheme is not supported: `fftw_scheme` = '%s'.",
        this->fftw_scheme.c_str()
      );
    }
  } else {
    if (!(this->fftw_scheme == "measure" || this->fftw_scheme == "")) {
      if (trvs::currTask == 0) {
        trvs::logger.info(
          "FFTW planner scheme is ignored in GPU mode: `fftw_scheme` = '%s'.",
          this->fftw_scheme.c_str()
        );
      }
    }
    this->fftw_scheme = "";  // transmutation
  }

  if (trvs::is_gpu_enabled()) {
    if (!(this->use_fftw_wisdom == "false")) {
      if (trvs::currTask == 0) {
        trvs::logger.info("FFTW wisdom is disabled in GPU mode.");
      }
    }
    this->use_fftw_wisdom = "";  // transmutation
  } else {
    if (this->use_fftw_wisdom == "false") {
      this->use_fftw_wisdom = "";  // transmutation
    } else
    if (init) {
      this->use_fftw_wisdom += "/";  // transmutation
    }
  }
  trvs::expand_envar_in_path(this->use_fftw_wisdom);

  if (this->use_fftw_wisdom != "") {
    if (this->fftw_scheme != "measure" && this->fftw_scheme != "patient") {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "FFTW wisdom is enabled but the planner scheme "
          "is not 'measure' or 'patient': "
          "`fftw_scheme` = '%s'.",
          this->fftw_scheme.c_str()
        );
      }
      throw trvs::InvalidParameterError(
        "FFTW wisdom is enabled but the planner scheme "
        "is not 'measure' or 'patient': "
        "`fftw_scheme` = '%s'.",
        this->fftw_scheme.c_str()
      );
    }

    char fftw_wisdom_file_f_[1024];
    char fftw_wisdom_file_b_[1024];

#if defined(TRV_USE_OMP) && defined(TRV_USE_FFTWOMP)
    std::snprintf(
      fftw_wisdom_file_f_, sizeof(fftw_wisdom_file_f_),
      "%sfftw_omp_cif_%dx%dx%d.wisdom",
      this->use_fftw_wisdom.c_str(),
      this->ngrid[0], this->ngrid[1], this->ngrid[2]
    );
    std::snprintf(
      fftw_wisdom_file_b_, sizeof(fftw_wisdom_file_b_),
      "%sfftw_omp_cib_%dx%dx%d.wisdom",
      this->use_fftw_wisdom.c_str(),
      this->ngrid[0], this->ngrid[1], this->ngrid[2]
    );
#else  // !TRV_USE_OMP || !TRV_USE_FFTWOMP
    std::snprintf(
      fftw_wisdom_file_f_, sizeof(fftw_wisdom_file_f_),
      "%sfftw_cif_%dx%dx%d.wisdom",
      this->use_fftw_wisdom.c_str(),
      this->ngrid[0], this->ngrid[1], this->ngrid[2]
    );
    std::snprintf(
      fftw_wisdom_file_b_, sizeof(fftw_wisdom_file_b_),
      "%sfftw_cib_%dx%dx%d.wisdom",
      this->use_fftw_wisdom.c_str(),
      this->ngrid[0], this->ngrid[1], this->ngrid[2]
    );
#endif  // TRV_USE_OMP && TRV_USE_FFTWOMP

    this->fftw_wisdom_file_f = fftw_wisdom_file_f_;
    this->fftw_wisdom_file_b = fftw_wisdom_file_b_;
  }

  char default_bvec_sfilepath[1024];
  std::snprintf(
    default_bvec_sfilepath, sizeof(default_bvec_sfilepath),
    "%sbinned_vectors%s",
    this->measurement_dir.c_str(), this->output_tag.c_str()
  );
  if (this->save_binned_vectors == "false") {
    this->save_binned_vectors = "";  // transmutation
  } else
  if (this->save_binned_vectors == "true") {
    this->save_binned_vectors = default_bvec_sfilepath;  // transmutation
  } else
  if (this->save_binned_vectors != "") {
    // Check whether path is absolute.
    if (this->save_binned_vectors.rfind("/", 0) != 0 && init) {
      this->save_binned_vectors = this->measurement_dir
        + this->save_binned_vectors;
    }  // transmutation
  }
  trvs::expand_envar_in_path(this->save_binned_vectors);

  // Validate and derive numerical parameters.
  this->volume =
    this->boxsize[0] * this->boxsize[1] * this->boxsize[2];  // derivation
  this->nmesh = static_cast<long long>(this->ngrid[0])
    * this->ngrid[1] * this->ngrid[2];  // derivation

  if (this->volume <= 0.) {
    if (trvs::currTask == 0) {
      trvs::logger.warn(
        "Derived total box volume is non-positive: `volume` = %.6e. "
        "Possible numerical overflow due to large `boxsize`, "
        "or `boxsize` is unset. "
        "In the latter case, the box size will be calculated using "
        "the particle coordinate spans and the box expansion factor.",
        this->volume
      );
    }
    if (this->expand < 1.) {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "The box expansion factor must be >= 1: `expand` = %lg",
          this->expand
        );
      }
      throw trvs::InvalidParameterError(
        "The box expansion factor must be >= 1: `expand` = %lg",
        this->expand
      );
    }
    if (this->expand > 1. and this->catalogue_type == "sim") {
      if (trvs::currTask == 0) {
        trvs::logger.info(
          "The box expansion factor is expected to be unity "
          "for 'sim' catalogue(s): `expand` = %lg",
          this->expand
        );
      }
    }
  }
  if (this->nmesh <= 0) {
    if (trvs::currTask == 0) {
      trvs::logger.warn(
        "Derived total mesh grid number is non-positive: `nmesh` = %lld. "
        "Possible numerical overflow due to large `ngrid`, "
        "or `ngrid` is unset. "
        "In the latter case, the grid cell number will be calculated using "
        "the Nyquist cutoff and the box size.",
        this->nmesh
      );
    }
    if (this->cutoff_nyq < 0.) {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "The Nyquist cutoff must be non-negative: `cutoff_nyq` = %lg",
          this->cutoff_nyq
        );
      }
      throw trvs::InvalidParameterError(
        "The Nyquist cutoff must be non-negative: `cutoff_nyq` = %lg",
        this->cutoff_nyq
      );
    }
    if (this->volume > 0.) {
      set_ngrid_from_cutoff(*this);
    }
  }

  if (this->alignment == "pad") {
    if (this->padfactor < 0.) {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Padding is enabled but the padding factor is negative: "
          "`padfactor` = %lg",
          this->padfactor
        );
      }
      throw trvs::InvalidParameterError(
        "Padding is enabled but the padding factor is negative: "
        "`padfactor` = %lg",
        this->padfactor
      );
    }
    if (this->padscale == "box" && this->padfactor >= 1.) {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Padding is enabled but the %s padding factor is too large "
          "for the box size: `padfactor` = %lg",
          this->padscale.c_str(), this->padfactor
        );
      }
      throw trvs::InvalidParameterError(
        "Padding is enabled but the %s padding factor is too large "
        "for the box size: `padfactor` = %lg",
        this->padscale.c_str(), this->padfactor
      );
    }
    if (this->padscale == "grid" && (
        this->padfactor >= this->ngrid[0]
        || this->padfactor >= this->ngrid[1]
        || this->padfactor >= this->ngrid[2]
    )) {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Padding is enabled but the %s padding factor is too large "
          "for the mesh grid numbers: `padfactor` = %lg",
          this->padscale.c_str(), this->padfactor
        );
      }
      throw trvs::InvalidParameterError(
        "Padding is enabled but the %s padding factor is too large "
        "for the mesh grid numbers: `padfactor` = %lg",
        this->padscale.c_str(), this->padfactor
      );
    }
  }

  if (this->bin_min < 0.) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Lower bin edge must be non-negative.");
    }
    throw trvs::InvalidParameterError(
      "Lower bin edge must be non-negative.\n"
    );
  }
  if (this->bin_min >= this->bin_max) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Lower bin edge must be less than the upper bin edge."
      );
    }
    throw trvs::InvalidParameterError(
      "Lower bin edge must be less than the upper bin edge.\n"
    );
  }
  if (this->space == "fourier") {
    double wavenum_nyquist = M_PI
      * *std::min_element(this->ngrid, this->ngrid + 3)
      / *std::max_element(this->boxsize, this->boxsize + 3);
    if (this->bin_min > wavenum_nyquist) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Lower wavenumber limit exceeds the Nyquist wavenumber %.4f.",
          wavenum_nyquist
        );
      }
    }
  } else
  if (this->space == "config") {
    double separation_nyquist = 2
      * *std::max_element(this->boxsize, this->boxsize + 3)
      / *std::min_element(this->ngrid, this->ngrid + 3);
    if (this->bin_max < separation_nyquist) {
      if (trvs::currTask == 0) {
        trvs::logger.warn(
          "Upper separation limit undershoots the Nyquist scale %.4f.",
          separation_nyquist
        );
      }
    }
  }

  if (this->num_bins < 2) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Number of bins `num_bins` must be >= 2.");
    }
    throw trvs::InvalidParameterError(
      "Number of bins `num_bins` must be >= 2.\n"
    );
  }

  if (
    this->idx_bin < 0
    && this->npoint == "3pt"
    && this->form == "row"
  ) {
    if (trvs::currTask == 0) {
      trvs::logger.error("Fixed row bin index `idx_bin` must be >= 0.");
    }
    throw trvs::InvalidParameterError(
      "Fixed row bin index `idx_bin` must be >= 0.\n"
    );
  }

  // Check for parameter conflicts.
  if (this->binning == "linpad" || this->binning == "logpad") {
    // CAVEAT: See @ref trv::Binning.
    int nbin_pad = 5;

    if (this->num_bins < nbin_pad + 2) {
      if (trvs::currTask == 0) {
        trvs::logger.error(
          "Binning scheme '%s' requires `num_bins` >= %d.",
          this->binning.c_str(), nbin_pad + 2
        );
      }
      throw trvs::InvalidParameterError(
        "Binning scheme '%s' requires `num_bins` >= %d.",
        this->binning.c_str(), nbin_pad + 2
      );
    }
  }

  if (std::abs(this->idx_bin) >= this->num_bins) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Bin index `idx_bin` must be < `num_bins` in absolute value."
      );
    }
    throw trvs::InvalidParameterError(
      "Bin index `idx_bin` must be < `num_bins` in absolute value.\n"
    );
  }

  if (this->npoint == "3pt" && this->interlace == "true") {
    this->interlace = "false";  // transmutation

    if (trvs::currTask == 0) {
      trvs::logger.info(
        "Interlacing is unsupported for three-point measurements. "
        "`interlace` is set to 'false'."
      );
    }
  }

  if (
    (this->statistic_type == "modes" || this->statistic_type == "pairs")
    && this->save_binned_vectors == ""
  ) {
    this->save_binned_vectors = default_bvec_sfilepath;  // transmutation
    if (trvs::currTask == 0) {
      trvs::logger.info(
        "`save_binned_vectors` is overridden, as `statistic_type` is '%s' "
        "so binned vectors are saved as the output to the default path.",
        this->statistic_type.c_str()
      );
    }
  }

  if (trvs::currTask == 0) {
    trvs::logger.stat("Parameters validated.");
  }

  return 0;
}

int ParameterSet::print_to_file(char* out_parameter_filepath) {
  // Create output file.
  std::FILE* ofileptr;
  if (!(ofileptr = std::fopen(out_parameter_filepath, "w"))) {
    if (trvs::currTask == 0) {
      trvs::logger.error(
        "Non-existent or unwritable output directory: %s",
        this->measurement_dir.c_str()
      );
    }
    throw trvs::IOError(
      "Non-existent or unwritable output directory: %s",
      this->measurement_dir.c_str()
    );
  }

  // Define convenience function for printing parameters.
  auto print_par_str = [ofileptr](
    const char* fmt, const std::string& par_val
  ) {
    std::fprintf(ofileptr, fmt, par_val.c_str());
  };
  auto print_par_int = [ofileptr](const char* fmt, int par_val) {
    std::fprintf(ofileptr, fmt, par_val);
  };
  auto print_par_double = [ofileptr](const char* fmt, double par_val) {
    std::fprintf(ofileptr, fmt, par_val);
  };

  // Print parameters to file.
  print_par_str("catalogue_dir = %s\n", this->catalogue_dir);
  print_par_str("measurement_dir = %s\n", this->measurement_dir);
  print_par_str("data_catalogue_file = %s\n", this->data_catalogue_file);
  print_par_str("rand_catalogue_file = %s\n", this->rand_catalogue_file);
  print_par_str("catalogue_columns = %s\n", this->catalogue_columns);
  print_par_str("catalogue_dataset = %s\n", this->catalogue_dataset);
  print_par_str("output_tag = %s\n", this->output_tag);

  print_par_double("boxsize_x = %.3f\n", this->boxsize[0]);
  print_par_double("boxsize_y = %.3f\n", this->boxsize[1]);
  print_par_double("boxsize_z = %.3f\n", this->boxsize[2]);
  print_par_int("ngrid_x = %d\n", this->ngrid[0]);
  print_par_int("ngrid_y = %d\n", this->ngrid[1]);
  print_par_int("ngrid_z = %d\n", this->ngrid[2]);

  print_par_double("volume = %.6e\n", this->volume);
  print_par_int("nmesh = %lld\n", this->nmesh);

  print_par_double("expand = %.4f\n", this->expand);
  print_par_str("alignment = %s\n", this->alignment);
  print_par_str("padscale = %s\n", this->padscale);
  print_par_double("padfactor = %.4f\n", this->padfactor);

  print_par_str("assignment = %s\n", this->assignment);
  print_par_str("interlace = %s\n", this->interlace);
  print_par_int("assignment_order = %d\n", this->assignment_order);

  print_par_str("catalogue_type = %s\n", this->catalogue_type);
  print_par_str("statistic_type = %s\n", this->statistic_type);
  print_par_str("npoint = %s\n", this->npoint);
  print_par_str("space = %s\n", this->space);

  print_par_int("ell1 = %d\n", this->ell1);
  print_par_int("ell2 = %d\n", this->ell2);
  print_par_int("ELL = %d\n", this->ELL);

  print_par_int("i_wa = %d\n", this->i_wa);
  print_par_int("j_wa = %d\n", this->j_wa);

  print_par_str("form = %s\n", this->form);
  print_par_str("norm_convention = %s\n", this->norm_convention);
  print_par_str("binning = %s\n", this->binning);
  print_par_str("shape = %s\n", this->shape);

  print_par_double("bin_min = %.4f\n", this->bin_min);
  print_par_double("bin_max = %.4f\n", this->bin_max);
  print_par_int("num_bins = %d\n", this->num_bins);
  print_par_int("idx_bin = %d\n", this->idx_bin);

  print_par_str("fftw_scheme = %s\n", this->fftw_scheme);
  print_par_str("use_fftw_wisdom = %s\n", this->use_fftw_wisdom.c_str());
  print_par_str("fftw_wisdom_file_f = %s\n", this->fftw_wisdom_file_f.c_str());
  print_par_str("fftw_wisdom_file_b = %s\n", this->fftw_wisdom_file_b.c_str());
  print_par_str("save_binned_vectors = %s\n", this->save_binned_vectors);
  print_par_str("progbar = %s\n", this->progbar);
  print_par_int("verbose = %d\n", this->verbose);
  print_par_int("fftw_planner_flag = %d\n", this->fftw_planner_flag);

  std::fclose(ofileptr);

  if (trvs::currTask == 0) {
    trvs::logger.info(
      "Check used-parameter file for reference: %s", out_parameter_filepath
    );
  }

  return 0;
}

int ParameterSet::print_to_file() {
  // Set output file path to default.
  char ofilepath[1024];
  std::snprintf(
    ofilepath, sizeof(ofilepath), "%sparameters_used%s",
    this->measurement_dir.c_str(), this->output_tag.c_str()
  );

  return ParameterSet::print_to_file(ofilepath);
}

void override_paramset_by_envvars(trv::ParameterSet& params) {
  // Override the output tag.
  char* ev_tag = std::getenv("TRV_OVERRIDE_OUTPUT_TAG");
  if (ev_tag != nullptr) {
    params.output_tag = std::string(ev_tag);
  }

  // Override the FFTW scheme.
  char* ev_fftw_scheme = std::getenv("TRV_OVERRIDE_FFTW_SCHEME");
  if (ev_fftw_scheme != nullptr && !trvs::is_gpu_enabled()) {
    params.fftw_scheme = std::string(ev_fftw_scheme);
  }

  // Override the FFTW wisdom option.
  char* ev_use_fftw_wisdom = std::getenv("TRV_OVERRIDE_USE_FFTW_WISDOM");
  if (ev_use_fftw_wisdom != nullptr && !trvs::is_gpu_enabled()) {
    params.use_fftw_wisdom = std::string(ev_use_fftw_wisdom);
  }

  // Override the verbosity level.
  char* ev_verbose = std::getenv("TRV_OVERRIDE_VERBOSE");
  if (ev_verbose != nullptr) {
    params.verbose = std::stoi(ev_verbose);
  }

  // Override the progress bar.
  char* ev_progbar = std::getenv("TRV_OVERRIDE_PROGBAR");
  if (ev_progbar != nullptr) {
    params.progbar = std::string(ev_progbar);
  }

  params.validate();
}

void set_boxsize_from_expand(
  const double* spans, trv::ParameterSet& params
) {
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    params.boxsize[iaxis] = spans[iaxis] * params.expand;
  }

  params.validate();

  if (trvs::currTask == 0) {
    trvs::logger.info(
      "Box size has been set from particle coordinate spans and "
      "expansion factor: (%.3f, %.3f, %.3f).",
      params.boxsize[0], params.boxsize[1], params.boxsize[2]
    );
  }
}

void set_ngrid_from_cutoff(trv::ParameterSet& params) {
  for (int iaxis = 0; iaxis < 3; iaxis++) {
    int ngrid_;
    if (params.space == "fourier") {
      ngrid_ = static_cast<int>(std::ceil(
        params.boxsize[iaxis] * params.cutoff_nyq / M_PI
      ));
    } else
    if (params.space == "config") {
      ngrid_ = static_cast<int>(std::ceil(
        2. * params.boxsize[iaxis] / params.cutoff_nyq
      ));
    } else {
      throw trvs::InvalidParameterError(
        "Space must be 'fourier' or 'config': `space` = '%s'.",
        params.space.c_str()
      );
    }
    params.ngrid[iaxis] = ngrid_ + (ngrid_ % 2);
  }

  params.validate();

  if (trvs::currTask == 0) {
    trvs::logger.info(
      "Mesh grid numbers have been set from Nyquist cutoff and "
      "box size: (%d, %d, %d).",
      params.ngrid[0], params.ngrid[1], params.ngrid[2]
    );
  }
}

}  // namespace trv
