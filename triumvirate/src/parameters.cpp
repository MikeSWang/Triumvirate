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
 * @file parameters.cpp
 * @author Mike S Wang (https://github.com/MikeSWang)
 *
 */

#include "parameters.hpp"

namespace trv {

int ParameterSet::read_from_file(char* parameter_filepath) {
  /// --------------------------------------------------------------------
  /// Initialisation
  /// --------------------------------------------------------------------

  /// Load parameter file.
  std::string param_filepath = parameter_filepath;

  std::ifstream fin(param_filepath.c_str());

  /// Initialise temporary variables to hold the extracted parameters.
  char catalogue_dir_[1024];
  char measurement_dir_[1024];
  char data_catalogue_file_[1024];
  char rand_catalogue_file_[1024];
  char catalogue_columns_[1024];
  char output_tag_[1024];

  double boxsize_x, boxsize_y, boxsize_z;
  int ngrid_x, ngrid_y, ngrid_z;

  char alignment_[16];
  char padscale_[16];
  char assignment_[16];
  char interlace_[16];

  char catalogue_type_[16];
  char measurement_type_[16];
  char norm_convention_[16];
  char shotnoise_convention_[16];
  char binning_[16];
  char form_[16];

  /// --------------------------------------------------------------------
  /// Extraction
  /// --------------------------------------------------------------------
  std::string line_str;
  char dummy_str[1024];
  while (std::getline(fin, line_str)) {
    /// Check if the line is a parameter assignment.
    if (
      std::sscanf(
        line_str.data(), "%s %s %s", dummy_str, dummy_str, dummy_str
      ) != 3
    ) {
      continue;
    }

    /// Define convenience function for scanning string parameters.
    auto scan_par_str = [line_str, dummy_str](
      const char* par_name, const char* fmt, const char* par_value
    ) {
      if (line_str.find(par_name) != std::string::npos) {
        std::sscanf(
          line_str.data(), fmt, dummy_str, dummy_str, par_value
        );
      }
    };

    /// I/O --------------------------------------------------------------

    scan_par_str("catalogue_dir", "%s %s %s", catalogue_dir_);
    scan_par_str("measurement_dir", "%s %s %s", measurement_dir_);
    scan_par_str("data_catalogue_file", "%s %s %s", data_catalogue_file_);
    scan_par_str("rand_catalogue_file", "%s %s %s", rand_catalogue_file_);
    scan_par_str("catalogue_columns", "%s %s %s", catalogue_columns_);
    scan_par_str("output_tag", "%s %s %s", output_tag_);

    /// Mesh sampling ----------------------------------------------------

    if (line_str.find("boxsize_x") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %lg", dummy_str, dummy_str, &boxsize_x
      );
    }
    if (line_str.find("boxsize_y") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %lg", dummy_str, dummy_str, &boxsize_y
      );
    }
    if (line_str.find("boxsize_z") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %lg", dummy_str, dummy_str, &boxsize_z
      );
    }

    if (line_str.find("ngrid_x") != std::string::npos) {
      std::sscanf(line_str.data(), "%s %s %d", dummy_str, dummy_str, &ngrid_x);
    }
    if (line_str.find("ngrid_y") != std::string::npos) {
      std::sscanf(line_str.data(), "%s %s %d", dummy_str, dummy_str, &ngrid_y);
    }
    if (line_str.find("ngrid_z") != std::string::npos) {
      std::sscanf(line_str.data(), "%s %s %d", dummy_str, dummy_str, &ngrid_z);
    }

    scan_par_str("alignment", "%s %s %s", alignment_);
    scan_par_str("padscale", "%s %s %s", padscale_);

    if (line_str.find("padfactor") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %lg", dummy_str, dummy_str, &this->padfactor
      );
    }

    scan_par_str("assignment", "%s %s %s", assignment_);
    scan_par_str("interlace", "%s %s %s", interlace_);

    /// Measurement ------------------------------------------------------

    scan_par_str("catalogue_type", "%s %s %s", catalogue_type_);
    scan_par_str("measurement_type", "%s %s %s", measurement_type_);
    scan_par_str("norm_convention", "%s %s %s", norm_convention_);
    scan_par_str("shotnoise_convention", "%s %s %s", shotnoise_convention_);
    scan_par_str("binning", "%s %s %s", binning_);
    scan_par_str("form", "%s %s %s", form_);

    if (line_str.find("ell1") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->ell1
      );
    }
    if (line_str.find("ell2") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->ell2
      );
    }
    if (line_str.find("ELL") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->ELL
      );
    }

    if (line_str.find("i_wa") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->i_wa
      );
    }
    if (line_str.find("j_wa") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->j_wa
      );
    }

    if (line_str.find("bin_min") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %lg", dummy_str, dummy_str, &this->bin_min
      );
    }
    if (line_str.find("bin_max") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %lg", dummy_str, dummy_str, &this->bin_max
      );
    }

    if (line_str.find("num_bins") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->num_bins
      );
    }
    if (line_str.find("idx_bin") != std::string::npos) {
      std::sscanf(
        line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->idx_bin
      );
    }
  }

  /// Misc ---------------------------------------------------------------

  if (line_str.find("verbose") != std::string::npos) {
    std::sscanf(
      line_str.data(), "%s %s %d", dummy_str, dummy_str, &this->verbose
    );
  }

  /// --------------------------------------------------------------------
  /// Attribution
  /// --------------------------------------------------------------------

  /// Attribute scalar parameters (directly extracted above).

  /// Attribute string parameters.
  this->catalogue_dir = catalogue_dir_;
  this->measurement_dir = measurement_dir_;
  this->data_catalogue_file = data_catalogue_file_;
  this->rand_catalogue_file = rand_catalogue_file_;
  this->catalogue_columns = catalogue_columns_;
  this->output_tag = output_tag_;

  this->alignment = alignment_;
  this->padscale = padscale_;
  this->assignment = assignment_;
  this->interlace = interlace_;

  this->catalogue_type = catalogue_type_;
  this->measurement_type = measurement_type_;
  this->norm_convention = norm_convention_;
  this->shotnoise_convention = shotnoise_convention_;
  this->binning = binning_;
  this->form = form_;

  /// Attribute derived parameters.
  this->boxsize[0] = boxsize_x;
  this->boxsize[1] = boxsize_y;
  this->boxsize[2] = boxsize_z;

  this->ngrid[0] = ngrid_x;
  this->ngrid[1] = ngrid_y;
  this->ngrid[2] = ngrid_z;

  this->volume = boxsize_x * boxsize_y * boxsize_z;
  this->nmesh = ngrid_x * ngrid_y * ngrid_z;

  /// --------------------------------------------------------------------
  /// Debugging mode
  /// --------------------------------------------------------------------

#ifdef DBG_PARS
  /// Define convenience function for displaying debugged parameters.
  auto debug_par_str = [](std::string name, std::string value) {
    std::cout << name << ": " << value << std::endl;
  };
  auto debug_par_int = [](std::string name, int value) {
    std::cout << name << ": " << value << std::endl;
  };
  auto debug_par_double = [](std::string name, double value) {
    std::cout << name << ": " << value << std::endl;
  };

  /// Display debugged parameters.
  debug_par_str("catalogue_dir", this->catalogue_dir);
  debug_par_str("measurement_dir", this->measurement_dir);
  debug_par_str("data_catalogue_file", this->data_catalogue_file);
  debug_par_str("rand_catalogue_file", this->rand_catalogue_file);
  debug_par_str("catalogue_columns", this->catalogue_columns);
  debug_par_str("output_tag", this->output_tag);

  debug_par_str("alignment", this->alignment);
  debug_par_str("padscale", this->padscale);
  debug_par_str("assignment", this->assignment);
  debug_par_str("interlace", this->interlace);

  debug_par_str("catalogue_type", this->catalogue_type);
  debug_par_str("measurement_type", this->measurement_type);
  debug_par_str("norm_convention", this->norm_convention);
  debug_par_str("shotnoise_convention", this->shotnoise_convention);
  debug_par_str("binning", this->binning);
  debug_par_str("form", this->form);

  debug_par_int("ngrid[0]", this->ngrid[0]);
  debug_par_int("ngrid[1]", this->ngrid[1]);
  debug_par_int("ngrid[2]", this->ngrid[2]);
  debug_par_int("nmesh", this->nmesh);

  debug_par_int("ell1", this->ell1);
  debug_par_int("ell2", this->ell2);
  debug_par_int("ELL", this->ELL);
  debug_par_int("i_wa", this->i_wa);
  debug_par_int("j_wa", this->j_wa);

  debug_par_int("num_bins", this->num_bins);
  debug_par_int("idx_bin", this->idx_bin);

  debug_par_double("boxsize[0]", this->boxsize[0])
  debug_par_double("boxsize[1]", this->boxsize[1])
  debug_par_double("boxsize[2]", this->boxsize[2])
  debug_par_double("volume", this->volume)
  debug_par_double("padfactor", this->padfactor)
  debug_par_double("bin_min", this->bin_min)
  debug_par_double("bin_max", this->bin_max)
#endif  // DBG_PARS

  return this->validate();
}

int ParameterSet::validate() {
  /// Validate and derive string parameters.
  if (this->catalogue_dir != "") {
    this->catalogue_dir += "/";  // transmutation
  }  // any duplicate '/' has no effect
  if (this->catalogue_type == "survey") {
    if (this->data_catalogue_file != "") {
      this->data_catalogue_file = this->catalogue_dir
        + this->data_catalogue_file;  // transmutation
    }
    if (this->rand_catalogue_file != "") {
      this->rand_catalogue_file = this->catalogue_dir
        + this->rand_catalogue_file;  // transmutation
    }
  } else
  if (this->catalogue_type == "random") {
    this->data_catalogue_file = "";  // transmutation
    if (this->rand_catalogue_file != "") {
      this->rand_catalogue_file = this->catalogue_dir
        + this->rand_catalogue_file;  // transmutation
    }
  } else
  if (this->catalogue_type == "sim") {
    if (this->data_catalogue_file != "") {
      this->data_catalogue_file = this->catalogue_dir
        + this->data_catalogue_file;  // transmutation
    }
    this->rand_catalogue_file = "";  // transmutation
  } else {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Catalogue type must be 'survey', 'random' or 'sim': "
        "`catalogue_type` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->catalogue_type.c_str()
      );
    }
  }

  if (!(this->alignment == "centre" || this->alignment == "pad")) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Box alignment must be 'centre' or 'pad': "
        "`alignment` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->alignment.c_str()
      );
    }
  }
  if (!(this->padscale == "box" || this->padscale == "grid")) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Pad scale must be 'box' or 'grid': `padscale` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->padscale.c_str()
      );
    }
  }

  if (!(
    this->assignment == "ngp"
    || this->assignment == "cic"
    || this->assignment == "tsc"
    || this->assignment == "pcs"
  )) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Mesh assignment scheme must be "
        "'ngp', 'cic', 'tsc' or 'pcs': `assignment` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->assignment.c_str()
      );
    }
  }
  if (this->interlace == "true" || this->interlace == "on") {
    this->interlace = "true";  // transmutation
  } else
  if (this->interlace == "false" || this->interlace == "off") {
    this->interlace = "false";  // transmutation
  } else {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Interlacing must be 'true'/'on' or 'false'/'off': "
        "`interlace` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->interlace.c_str()
      );
    }
  }

  if (this->measurement_type == "powspec") {
    this->npoint = "2pt"; this->space = "fourier";  // derivation
  } else
  if (
    this->measurement_type == "2pcf" || this->measurement_type == "2pcf-win"
  ) {
    this->npoint = "2pt"; this->space = "config";  // derivation
  } else
  if (this->measurement_type == "bispec") {
    this->npoint = "3pt"; this->space = "fourier";  // derivation
  } else
  if (
    this->measurement_type == "3pcf"
    || this->measurement_type == "3pcf-win"
    || this->measurement_type == "3pcf-win-wa"
  ) {
    this->npoint = "3pt"; this->space = "config";  // derivation
  } else {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Measurement type is not recognised: "
        "`measurement_type` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->measurement_type.c_str()
      );
    }
  }
  if (!(
    this->norm_convention == "mesh" || this->norm_convention == "particle"
  )) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Normalisation convention must be "
        "'mesh' or 'particle': `norm_convention` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->norm_convention.c_str()
      );
    }
  }
  if (!(
    this->shotnoise_convention == "mesh"
    || this->shotnoise_convention == "particle"
  )) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Shot noise convention convention must be "
        "'mesh' or 'particle': `shotnoise_convention` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->shotnoise_convention.c_str()
      );
    }
  }
  if (!(
    this->binning == "lin"
    || this->binning == "log"
    || this->binning == "linpad"
    || this->binning == "logpad"
    || this->binning == "custom"
  )) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Binning scheme is unrecognised: `binning` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->binning.c_str()
      );
    }
  }
  if (!(this->form == "diag" || this->form == "full")) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] `form` must be either 'full' or 'diag': `form` = '%s'.\n",
        trv::sys::show_timestamp().c_str(),
        this->form.c_str()
      );
    }
  }

  /// Validate numerical parameters.
  if (this->volume < 0.) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Derived total box volume is non-positive: "
        "`volume` = '%d'. Possible numerical overflow due to large `boxsize`, "
        "or `boxsize` is unset.\n",
        trv::sys::show_timestamp().c_str(),
        this->nmesh
      );
    }
  }
  if (this->nmesh <= 0) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Derived total mesh grid number is non-positive: "
        "`nmesh` = '%d'. Possible numerical overflow due to large `ngrid`, "
        "or `ngrid` is unset.\n",
        trv::sys::show_timestamp().c_str(),
        this->nmesh
      );
    }
  }

  if (this->alignment == "pad") {
    if (this->padfactor < 0.) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Padding is enabled but the padding factor is negative: "
        "`padfactor` = '%lg'.\n",
        trv::sys::show_timestamp().c_str(),
        this->padfactor
      );
    }
    if (this->padscale == "box" && this->padfactor >= 1.) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Padding is enabled but the %s padding factor is too large "
        "for the box size: `padfactor` = '%lg'.\n",
        trv::sys::show_timestamp().c_str(),
        this->padscale, this->padfactor
      );
    }
    if (this->padscale == "grid" && (
        this->padfactor >= this->ngrid[0]
        || this->padfactor >= this->ngrid[1]
        || this->padfactor >= this->ngrid[2]
    )) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Padding is enabled but the %s padding factor is too large "
        "for the mesh grid numbers: `padfactor` = '%lg'.\n",
        trv::sys::show_timestamp().c_str(),
        this->padscale, this->padfactor
      );
    }
  }

  if (this->num_bins < 2) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Number of bins `num_bins` must be >= 2.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  if (this->idx_bin < 0 && this->npoint == "3pt" && this->form == "full") {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Fixed bin index `idx_bin` must be >= 0.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  /// Check for parameter conflicts.
  if (this->binning == "linpad" || this->binning == "logpad") {
    /// SEE: See @ref trv::Binning::set_bins().
    int nbin_pad = 5;

    if (this->num_bins < nbin_pad + 2) {
      if (trv::sys::currTask == 0) {
        throw trv::sys::InvalidParameter(
          "[%s ERRO] Binning scheme '%s' requires `num_bins` >= %d.\n",
          trv::sys::show_timestamp().c_str(),
          this->binning.c_str(), nbin_pad + 2
        );
      }
    }
  }

  if (this->idx_bin >= this->num_bins) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::InvalidParameter(
        "[%s ERRO] Bin index `idx_bin` must be < `num_bins`.\n",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  if (this->npoint == "3pt" && this->interlace == "true") {
    this->interlace = "false";  // transmutation

    if (trv::sys::currTask == 0) {
      std::printf(
        "[%s WARN] Interlacing is unsupported for 3-point measurements. "
        "`interlace` is set to 'false'.",
        trv::sys::show_timestamp().c_str()
      );
    }
  }

  if (trv::sys::currTask == 0) {
    std::printf(
      "[%s STAT] Parameters validated.\n",
      trv::sys::show_timestamp().c_str()
    );
  }

  return 0;
}

int ParameterSet::print_to_file(char* out_parameter_filepath) {
  /// Create output file.
  std::FILE* ofileptr;
  if (!(ofileptr = std::fopen(out_parameter_filepath, "w"))) {
    if (trv::sys::currTask == 0) {
      throw trv::sys::IOError(
        "[%s ERRO] Non-existent or unwritable output directory: %s.\n",
        trv::sys::show_timestamp().c_str(),
        this->measurement_dir.c_str()
      );
    }
  }

  /// Define convenience function for printing parameters.
  auto print_par_str = [ofileptr](const char* fmt, std::string par_val) {
    std::fprintf(ofileptr, fmt, par_val.c_str());
  };
  auto print_par_int = [ofileptr](const char* fmt, int par_val) {
    std::fprintf(ofileptr, fmt, par_val);
  };
  auto print_par_double = [ofileptr](const char* fmt, double par_val) {
    std::fprintf(ofileptr, fmt, par_val);
  };

  /// Print parameters to file.
  print_par_str("catalogue_dir = %s\n", this->catalogue_dir);
  print_par_str("measurement_dir = %s\n", this->measurement_dir);
  print_par_str("data_catalogue_file = %s\n", this->data_catalogue_file);
  print_par_str("rand_catalogue_file = %s\n", this->rand_catalogue_file);
  print_par_str("catalogue_columns = %s\n", this->catalogue_columns);
  print_par_str("output_tag = %s\n", this->output_tag);

  print_par_double("boxsize_x = %.2f\n", this->boxsize[0]);
  print_par_double("boxsize_y = %.2f\n", this->boxsize[1]);
  print_par_double("boxsize_z = %.2f\n", this->boxsize[2]);
  print_par_int("ngrid_x = %d\n", this->ngrid[0]);
  print_par_int("ngrid_y = %d\n", this->ngrid[1]);
  print_par_int("ngrid_z = %d\n", this->ngrid[2]);

  print_par_double("volume = %.6e\n", this->volume);
  print_par_int("nmesh = %d\n", this->nmesh);

  print_par_str("alignment = %s\n", this->alignment);
  print_par_str("padscale = %s\n", this->padscale);
  print_par_double("padfactor = %.4f\n", this->padfactor);

  print_par_str("assignment = %s\n", this->assignment);
  print_par_str("interlace = %s\n", this->interlace);

  print_par_str("catalogue_type = %s\n", this->catalogue_type);
  print_par_str("measurement_type = %s\n", this->measurement_type);

  print_par_str("norm_convention = %s\n", this->norm_convention);
  print_par_str("shotnoise_convention = %s\n", this->shotnoise_convention);

  print_par_str("binning = %s\n", this->binning);
  print_par_str("form = %s\n", this->form);

  print_par_str("npoint = %s\n", this->npoint);
  print_par_str("space = %s\n", this->space);

  print_par_int("ell1 = %d\n", this->ell1);
  print_par_int("ell2 = %d\n", this->ell2);
  print_par_int("ELL = %d\n", this->ELL);

  print_par_int("i_wa = %d\n", this->i_wa);
  print_par_int("j_wa = %d\n", this->j_wa);

  print_par_double("bin_min = %.4f\n", this->bin_min);
  print_par_double("bin_max = %.4f\n", this->bin_max);
  print_par_int("num_bins = %d\n", this->num_bins);
  print_par_int("idx_bin = %d\n", this->idx_bin);

  print_par_int("verbose = %d\n", this->verbose);

  std::fclose(ofileptr);

  return 0;
}

int ParameterSet::print_to_file() {
  /// Set output file path to default.
  char ofilepath[1024];
  std::sprintf(
    ofilepath, "%s/parameters_used%s",
    this->measurement_dir.c_str(), this->output_tag.c_str()
  );

  return ParameterSet::print_to_file(ofilepath);
}

}  // namespace trv
