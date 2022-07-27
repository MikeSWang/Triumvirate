/**
 * @file parameters.hpp
 * @brief Program parameter configuration.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_

#include <cstdio>
#include <fstream>
#include <string>

#include <sys/stat.h>

#include "monitor.hpp"

namespace trv {
namespace scheme {

/**
 * Set of parameters.
 *
 * This reads parameters from a file, stores and prints out the extracted
 * parameters, and validates the input.
 *
 */
class ParameterSet {
 public:
  /// I/O.
  std::string catalogue_dir;  ///< catalogue directory
  std::string measurement_dir;  ///< output directory
  std::string data_catalogue_file;  ///< data catalogue file
  std::string rand_catalogue_file;  ///< random catalogue file
  std::string catalogue_header;  ///< catalogue file header (comma-separated)
  std::string output_tag;  ///< output file tag

  /// Sampling.
  double boxsize[3];  ///< boxsize in each dimension
  int ngrid[3];  ///< grid number in each dimension
  std::string alignment = "pad";  /**< box alignment choice:
                                       {'centre', 'pad' (default)} */
  std::string assignment;  ///< mesh assignment scheme: {'ngp', 'cic', 'tsc'}
  std::string norm_convention;  /**< normalisation convention:
                                     {'mesh', 'particle'} */
  std::string shotnoise_convention;  /**< shot noise convention:
                                          {'mesh', 'particle'} */

  /// Measurements.
  std::string catalogue_type;  ///< catalogue type: {'survey', 'mock', 'sim'}
  std::string measurement_type; /**< measurement_type: {
                                       'powspec', '2pcf', '2pcf-win',
                                       'bispec', '3pcf', '3pcf-win',
                                       '3pcf-win-wa'
                                     } */
  int ell1;  ///< spherical degree associated with the first wavevector
  int ell2;  ///< spherical degree associated with the second wavevector
  int ELL;  ///< spherical degree associated with the line-of-sight vector
  int i_wa;  ///< first order of the wide-angle correction term
  int j_wa;  ///< second order of the wide-angle correction term
  std::string binning;  /**< binning interval:
                             {'lin', 'log', 'linpad', 'logpad'} */
  std::string form;  ///< full or diagonal form for the bispectrum
  double kmin;  ///< minimum wavenumber
  double kmax;  ///< maximum wavenumber
  int num_kbin;  ///< number of wavenumber bins
  int ith_kbin;  ///< index of the separation bin measured
  double rmin;  ///< minimum separation
  double rmax;  ///< maximum separation
  int num_rbin;  ///< number of separation bins
  int ith_rbin;  ///< index of the separation bin measured

  /// Derived quantities.
  double volume;  ///< box volume
  int nmesh;  ///< total mesh grid number

  /// Misc.
  int batch_number;  ///< job batch number

  /**
   * Read parameters from a file.
   *
   * @param argv Command-line argument list.
   * @returns Validation exit status.
   */
  int read_from_file(char* argv[]) {
    /// Load the parameter file (assumed to be the second argument,
    /// i.e. after the program file).
    std::string param_file = argv[1];

    std::ifstream fin(param_file.c_str());

    /// Initialise temporary variables to hold the extracted parameters.
    char catalogue_dir_[1024];
    char measurement_dir_[1024];
    char data_catalogue_file_[1024];
    char rand_catalogue_file_[1024];
    char catalogue_header_[1024];
    char output_tag_[1024];
    char catalogue_type_[16];
    char measurement_type_[16];
    char alignment_[16];
    char assignment_[16];
    char norm_convention_[16];
    char shotnoise_convention_[16];
    char binning_[16];
    char form_[16];

    double boxsize_x, boxsize_y, boxsize_z;
    int ngrid_x, ngrid_y, ngrid_z;

    /// Extract parameters from file contents by line parsing.
    std::string str_line;  // string representing the line being parsed
    char str_dummy[1024];  // string placeholder for irrelevant contents
    while (std::getline(fin, str_line)) {
      /// Check if the line is in the correct format for parameter
      /// parsing; process if yes, otherwise move on to the next line.
      if (
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, str_dummy
        ) != 3
      ) {
        continue;
      }

      if (str_line.find("catalogue_dir") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_dir_
        );
      }
      if (str_line.find("measurement_dir") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, measurement_dir_
        );
      }
      if (str_line.find("data_catalogue_file") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s",
          str_dummy, str_dummy, data_catalogue_file_
        );
      }
      if (str_line.find("rand_catalogue_file") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s",
          str_dummy, str_dummy, rand_catalogue_file_
        );
      }
      if (str_line.find("catalogue_header") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_header_
        );
      }
      if (str_line.find("output_tag") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, output_tag_
        );
      }

      if (str_line.find("catalogue_type") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_type_
        );
      }
      if (str_line.find("measurement_type") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, measurement_type_
        );
      }

      if (str_line.find("boxsize_x") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_x
        );
      }
      if (str_line.find("boxsize_y") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_y
        );
      }
      if (str_line.find("boxsize_z") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_z
        );
      }

      if (str_line.find("ngrid_x") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &ngrid_x
        );
      }
      if (str_line.find("ngrid_y") != std::string::npos) {
        std::sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &ngrid_y);
      }
      if (str_line.find("ngrid_z") != std::string::npos) {
        std::sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &ngrid_z);
      }

      if (str_line.find("alignment") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, alignment_
        );
      }
      if (str_line.find("assignment") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, assignment_
        );
      }
      if (str_line.find("norm_convention") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, norm_convention_
        );
      }
      if (str_line.find("shotnoise_convention") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s",
          str_dummy, str_dummy, shotnoise_convention_
        );
      }

      if (str_line.find("ell1") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ell1
        );
      }
      if (str_line.find("ell2") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ell2
        );
      }
      if (str_line.find("ELL") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ELL
        );
      }

      if (str_line.find("i_wa") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->i_wa
        );
      }
      if (str_line.find("j_wa") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->j_wa
        );
      }

      if (str_line.find("binning") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, binning_
        );
      }
      if (str_line.find("form") != std::string::npos) {
        std::sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, form_);
      }

      if (str_line.find("kmin") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->kmin
        );
      }
      if (str_line.find("kmax") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->kmax
        );
      }
      if (str_line.find("num_kbin") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->num_kbin
        );
      }

      if (str_line.find("rmin") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->rmin
        );
      }
      if (str_line.find("rmax") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->rmax
        );
      }
      if (str_line.find("num_rbin") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->num_rbin
        );
      }

      if (str_line.find("ith_kbin") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ith_kbin
        );
      }
      if (str_line.find("ith_rbin") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ith_rbin
        );
      }

      if (str_line.find("batch_number") != std::string::npos) {
        std::sscanf(
          str_line.data(), "%s %s %d",
          str_dummy, str_dummy, &this->batch_number
        );
      }
    }

    /// Store parameters as attributes.
    this->catalogue_dir = catalogue_dir_;
    this->measurement_dir = measurement_dir_;
    this->data_catalogue_file = data_catalogue_file_;
    this->rand_catalogue_file = rand_catalogue_file_;
    this->catalogue_header = catalogue_header_;
    this->output_tag = output_tag_;

    this->catalogue_type = catalogue_type_;
    this->measurement_type = measurement_type_;

    this->alignment = alignment_;
    this->assignment = assignment_;
    this->norm_convention = norm_convention_;
    this->shotnoise_convention = shotnoise_convention_;

    this->binning = binning_;
    this->form = form_;

    this->boxsize[0] = boxsize_x;
    this->boxsize[1] = boxsize_y;
    this->boxsize[2] = boxsize_z;

    this->volume = boxsize_x * boxsize_y * boxsize_z;

    this->ngrid[0] = ngrid_x;
    this->ngrid[1] = ngrid_y;
    this->ngrid[2] = ngrid_z;

    this->nmesh = ngrid_x * ngrid_y * ngrid_z;

    this->set_io_files();

    #ifdef DBG_PARS
    std::cout <<
      "catalogue_dir: " << this->catalogue_dir
    << std::endl;
    std::cout <<
      "measurement_dir: " << this->measurement_dir
    << std::endl;
    std::cout <<
      "data_catalogue_file: " << this->data_catalogue_file
    << std::endl;
    std::cout <<
      "rand_catalogue_file: " << this->rand_catalogue_file
    << std::endl;
    std::cout <<
      "catalogue_header: " << this->catalogue_header
    << std::endl;
    std::cout <<
      "output_tag: " << this->output_tag
    << std::endl;
    std::cout <<
      "boxsize: " << this->boxsize
    << std::endl;
    std::cout <<
      "ngrid: " << this->ngrid
    << std::endl;
    std::cout <<
      "alignment: " << this->alignment
    << std::endl;
    std::cout <<
      "assignment: " << this->assignment
    << std::endl;
    std::cout <<
      "norm_convention: " << this->norm_convention
    << std::endl;
    std::cout <<
      "shotnoise_convention: " << this->shotnoise_convention
    << std::endl;
    std::cout <<
      "catalogue_type: " << this->catalogue_type
    << std::endl;
    std::cout <<
      "measurement_type: " << this->measurement_type
    << std::endl;
    std::cout <<
      "binning: " << this->binning
    << std::endl;
    std::cout <<
      "form: " << this->form
    << std::endl;
    #endif  // DBG_PARS

    return this->validate();
  }

  /**
   * Print out extracted parameters to a file in the measurement output
   * directory.
   *
   * @returns Exit status.
   */
  int printout() {
    /// Set file path.
    char filepath[1024];
    std::sprintf(
      filepath, "%s/parameters_used%s",
      this->measurement_dir.c_str(), this->output_tag.c_str()
    );

    /// Create output file.
    std::FILE* used_param_fileptr;
    if (!(used_param_fileptr = std::fopen(filepath, "w"))) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::IOError(
          "[%s ERRO] Non-existent or unwritable output directory: '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->measurement_dir.c_str()
        );
      }
    }

    /// Print parameters to file.
    std::fprintf(
      used_param_fileptr, "catalogue_dir = %s\n",
      this->catalogue_dir.c_str()
    );
    std::fprintf(
      used_param_fileptr, "measurement_dir = %s\n",
      this->measurement_dir.c_str()
    );
    std::fprintf(
      used_param_fileptr, "data_catalogue_file = %s\n",
      this->data_catalogue_file.c_str()
    );
    std::fprintf(
      used_param_fileptr, "rand_catalogue_file = %s\n",
      this->rand_catalogue_file.c_str()
    );
    std::fprintf(
      used_param_fileptr, "catalogue_header = %s\n",
      this->catalogue_header.c_str()
    );
    std::fprintf(
      used_param_fileptr, "output_tag = %s\n", this->output_tag.c_str()
    );

    std::fprintf(
      used_param_fileptr, "catalogue_type = %s\n",
      this->catalogue_type.c_str()
    );
    std::fprintf(
      used_param_fileptr, "measurement_type = %s\n",
      this->measurement_type.c_str()
    );

    std::fprintf(
      used_param_fileptr, "alignment = %s\n", this->alignment.c_str()
    );
    std::fprintf(
      used_param_fileptr, "assignment = %s\n", this->assignment.c_str()
    );
    std::fprintf(
      used_param_fileptr, "norm_convention = %s\n",
      this->norm_convention.c_str()
    );
    std::fprintf(
      used_param_fileptr, "shotnoise_convention = %s\n",
      this->shotnoise_convention.c_str()
    );

    std::fprintf(used_param_fileptr, "boxsize_x = %.3f\n", this->boxsize[0]);
    std::fprintf(used_param_fileptr, "boxsize_y = %.3f\n", this->boxsize[1]);
    std::fprintf(used_param_fileptr, "boxsize_z = %.3f\n", this->boxsize[2]);

    std::fprintf(used_param_fileptr, "ngrid_x = %d\n", this->ngrid[0]);
    std::fprintf(used_param_fileptr, "ngrid_y = %d\n", this->ngrid[1]);
    std::fprintf(used_param_fileptr, "ngrid_z = %d\n", this->ngrid[2]);

    std::fprintf(used_param_fileptr, "ell1 = %d\n", this->ell1);
    std::fprintf(used_param_fileptr, "ell2 = %d\n", this->ell2);
    std::fprintf(used_param_fileptr, "ELL = %d\n", this->ELL);

    std::fprintf(used_param_fileptr, "i_wa = %d\n", this->i_wa);
    std::fprintf(used_param_fileptr, "j_wa = %d\n", this->j_wa);

    std::fprintf(used_param_fileptr, "binning = %s\n", this->binning.c_str());
    std::fprintf(used_param_fileptr, "form = %s\n", this->form.c_str());

    std::fprintf(used_param_fileptr, "kmin = %.6f\n", this->kmin);
    std::fprintf(used_param_fileptr, "kmax = %.6f\n", this->kmax);
    std::fprintf(used_param_fileptr, "num_kbin = %d\n", this->num_kbin);
    std::fprintf(used_param_fileptr, "ith_kbin = %d\n", this->ith_kbin);

    std::fprintf(used_param_fileptr, "rmin = %.3f\n", this->rmin);
    std::fprintf(used_param_fileptr, "rmax = %.3f\n", this->rmax);
    std::fprintf(used_param_fileptr, "num_rbin = %d\n", this->num_rbin);
    std::fprintf(used_param_fileptr, "ith_rbin = %d\n", this->ith_rbin);

    std::fclose(used_param_fileptr);

    return 0;
  }

  /**
   * Check parameters as valid input.
   *
   * @returns Exit status.
   */
  int validate() {
    /// Validate string parameters.
    if (!(
      this->catalogue_type == "survey"
      || this->catalogue_type == "mock"
      || this->catalogue_type == "sim"
    )) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::InvalidParameter(
          "[%s ERRO] Catalogue type must be 'survey', 'mock' or 'sim': "
          "`catalogue_type` = '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->catalogue_type.c_str()
        );
      }
    }

    if (!(this->alignment == "centre" || this->alignment == "pad")) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::InvalidParameter(
          "[%s ERRO] Box alignment must be 'centre' or 'pad': "
          "`alignment` = '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->alignment.c_str()
        );
      }
    }

    if (!(
      this->assignment == "ngp"
      || this->assignment == "cic"
      || this->assignment == "tsc"
    )) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::InvalidParameter(
          "[%s ERRO] Mesh assignment scheme must be 'ngp', 'cic' or 'tsc': "
          "`assignment` = '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->assignment.c_str()
        );
      }
    }

    if (!(
      this->norm_convention == "mesh" || this->norm_convention == "particle"
    )) {
      if (trv::runtime::currTask == 0) {
        std::printf(
          "[%s WARN] Normalisation convention must be "
          "'mesh' or 'particle': `norm_convention` = '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->norm_convention.c_str()
        );
      }

      /// Set to default convention.
      this->norm_convention == "mesh";

      if (trv::runtime::currTask == 0) {
        std::printf(
          "[%s WARN] Normalisation convention is set to "
          "default value '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->norm_convention.c_str()
        );
      }
    }

    if (!(
      this->shotnoise_convention == "mesh"
      || this->shotnoise_convention == "particle"
    )) {
      if (trv::runtime::currTask == 0) {
        std::printf(
          "[%s WARN] Shot noise convention convention must be "
          "'mesh' or 'particle': `shotnoise_convention` = '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->shotnoise_convention.c_str()
        );
      }

      /// Set to default convention.
      this->shotnoise_convention == "mesh";

      if (trv::runtime::currTask == 0) {
        std::printf(
          "[%s WARN] Shot noise convention is set to "
          "default value '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->shotnoise_convention.c_str()
        );
      }
    }

    if (!(this->form == "diag" || this->form == "full")) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::InvalidParameter(
          "[%s ERRO] `form` must be either 'full' or 'diag': "
          "`form` = '%s'.\n",
          trv::runtime::show_timestamp().c_str(),
          this->form.c_str()
        );
      }
    }

    /// Validate numerical parameters.
    if (this->num_kbin < 2 || this->num_rbin < 2) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::InvalidParameter(
          "[%s ERRO] Number of bins (`num_kbin` or `num_rbin`) "
          "must be >= 2.\n",
          trv::runtime::show_timestamp().c_str()
        );
      }
    }

    if (this->binning == "linpad" || this->binning == "logpad") {
      /// SEE: See `trv::scheme::BinScheme::set_rbin`.
      int nbin_custom = 5;

      if (this->num_rbin < nbin_custom + 2) {
        if (trv::runtime::currTask == 0) {
          throw trv::runtime::InvalidParameter(
            "[%s ERRO] Binning scheme '%s' requires `num_rbin` >= %d.\n",
            trv::runtime::show_timestamp().c_str(),
            this->binning.c_str(), nbin_custom + 2
          );
        }
      }
    }

    if (this->ith_kbin >= this->num_kbin || this->ith_rbin >= this->num_rbin) {
      if (trv::runtime::currTask == 0) {
        throw trv::runtime::InvalidParameter(
          "[%s ERRO] Bin index (`ith_kbin` or `ith_rbin`) must be less than "
          "the number of bins (`num_kbin` or `num_rbin`).\n",
          trv::runtime::show_timestamp().c_str()
        );
      }
    }

    if (trv::runtime::currTask == 0) {
      std::printf(
        "[%s STAT] Parameters validated.\n",
        trv::runtime::show_timestamp().c_str()
      );
    }

    return 0;
  }

 private:
  /**
   * Set I/O paths and directories for data and random catalogue files
   * and processed results.
   */
  void set_io_files() {
    std::string catalogue_dir_ = this->catalogue_dir;
    if (catalogue_dir_ != "") {
      catalogue_dir_ += "/";  // possible duplicate "/" has no effect on path
    }

    /// Set survey data and random catalogue inputs.
    if (this->catalogue_type == "survey") {
      if (this->data_catalogue_file != "") {
        this->data_catalogue_file = catalogue_dir_ + this->data_catalogue_file;
      }
      if (this->rand_catalogue_file != "") {
        this->rand_catalogue_file = catalogue_dir_ + this->rand_catalogue_file;
      }
    }

    if (this->catalogue_type == "sim") {
      this->data_catalogue_file = catalogue_dir_ + this->data_catalogue_file;
      this->rand_catalogue_file = "";
    }

    /// Set mock data and random catalogue inputs.  Make subdirectories
    /// to store the output for each data realisation.
    /// QUEST: Understand the serialisation of catalogue files for
    /// a suite of simulations.
    if (this->catalogue_type == "mock") {
      /// Enumerate the input data catalogue.
      int realisation = trv::runtime::numTasks * this->batch_number
        + trv::runtime::currTask + 1;

      if (realisation > 2048) {
        realisation = 1;
      }

      /// Set input paths.
      char buf_file[2048];

      std::sprintf(
        buf_file, "%s/%s_%04d.dat",
        this->catalogue_dir.c_str(),
        this->data_catalogue_file.c_str(),
        realisation
      );

      this->data_catalogue_file = buf_file;

      std::sprintf(
        buf_file, "%s/%s",
        this->catalogue_dir.c_str(),
        this->rand_catalogue_file.c_str()
      );

      this->rand_catalogue_file = buf_file;

      /// Set the output subdirectory.
      char buf_dir[2048];
      std::sprintf(
        buf_dir, "%s/%04d",
        this->measurement_dir.c_str(), realisation
      );

      struct stat status;  // file permission mode
      if (stat(buf_dir, &status) != 0) {
        /// Make the output subdirectory.
        int ret_status = mkdir(buf_dir, 0777);

        /// Check output subdirectory status and exit upon failure.
        if (ret_status == 0) {
          if (trv::runtime::currTask == 0) {
            std::printf(
              "[%s INFO] Output subdirectory '%s' is successfully made.\n",
              trv::runtime::show_timestamp().c_str(),
              buf_dir
            );
          }
        } else {
          if (trv::runtime::currTask == 0) {
            throw trv::runtime::IOError(
              "[%s ERRO] Failed to make output subdirectory '%s'.\n",
              trv::runtime::show_timestamp().c_str(),
              buf_dir
            );
          }
        }
      }

      this->measurement_dir = buf_dir;
    }
  }
};

}  // trv::scheme::
}  // trv::

#endif  // TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
