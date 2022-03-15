/**
 * @file parameters.hpp
 * @brief Program parameter configuration.
 *
 */

#ifndef TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
#define TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_

#include <cstring>
#include <fstream>
#include <string>

#include <sys/stat.h>

#include "common.hpp"

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
  std::string output_tag;  ///< output file tag

  /// Sampling.
  double boxsize[3];  ///< boxsize in each dimension
  int ngrid[3];  ///< grid number in each dimension
  std::string assignment;  ///< grid assignment scheme: {'NGP', 'CIC', 'TSC'}
  std::string norm_convention;  /**< normalisation convention:
                                     {'grid', 'particle'} */
  /// TODO: Implement shot noise and normalisation convention options.

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
    char output_tag_[1024];
    char catalogue_type_[16];
    char measurement_type_[16];
    char assignment_[16];
    char norm_convention_[16];
    char binning_[16];
    char form_[16];

    double boxsize_x, boxsize_y, boxsize_z;
    int ngrid_x, ngrid_y, ngrid_z;

    /// Extract parameters from file contents by line parsing.
    std::string str_line;  // string representing the line being parsed
    char str_dummy[1024];  // string placeholder for irrelevant contents
    while (getline(fin, str_line)) {
      /// Check if the line is in the correct format for parameter
      /// parsing; process if yes, otherwise move on to the next line.
      if (
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, str_dummy
        ) != 3
      ) {
        continue;
      }

      if (str_line.find("catalogue_dir") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_dir_
        );
      }
      if (str_line.find("measurement_dir") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, measurement_dir_
        );
      }
      if (str_line.find("data_catalogue_file") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s",
          str_dummy, str_dummy, data_catalogue_file_
        );
      }
      if (str_line.find("rand_catalogue_file") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s",
          str_dummy, str_dummy, rand_catalogue_file_
        );
      }
      if (str_line.find("output_tag") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, output_tag_
        );
      }

      if (str_line.find("catalogue_type") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_type_
        );
      }
      if (str_line.find("measurement_type") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, measurement_type_
        );
      }

      if (str_line.find("boxsize_x") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_x
        );
      }
      if (str_line.find("boxsize_y") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_y
        );
      }
      if (str_line.find("boxsize_z") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_z
        );
      }

      if (str_line.find("ngrid_x") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &ngrid_x
        );
      }
      if (str_line.find("ngrid_y") != std::string::npos) {
        sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &ngrid_y);
      }
      if (str_line.find("ngrid_z") != std::string::npos) {
        sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &ngrid_z);
      }

      if (str_line.find("assignment") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, assignment_
        );
      }
      if (str_line.find("norm_convention") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, norm_convention_
        );
      }

      if (str_line.find("ell1") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ell1
        );
      }
      if (str_line.find("ell2") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ell2
        );
      }
      if (str_line.find("ELL") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ELL
        );
      }

      if (str_line.find("i_wa") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->i_wa
        );
      }
      if (str_line.find("j_wa") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->j_wa
        );
      }

      if (str_line.find("binning") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %s", str_dummy, str_dummy, binning_
        );
      }
      if (str_line.find("form") != std::string::npos) {
        sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, form_);
      }

      if (str_line.find("kmin") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->kmin
        );
      }
      if (str_line.find("kmax") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->kmax
        );
      }
      if (str_line.find("num_kbin") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->num_kbin
        );
      }

      if (str_line.find("rmin") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->rmin
        );
      }
      if (str_line.find("rmax") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->rmax
        );
      }
      if (str_line.find("num_rbin") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->num_rbin
        );
      }

      if (str_line.find("ith_kbin") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ith_kbin
        );
      }
      if (str_line.find("ith_rbin") != std::string::npos) {
        sscanf(
          str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ith_rbin
        );
      }

      if (str_line.find("batch_number") != std::string::npos) {
        sscanf(
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
    this->output_tag = output_tag_;

    this->catalogue_type = catalogue_type_;
    this->measurement_type = measurement_type_;

    this->norm_convention = norm_convention_;
    this->assignment = assignment_;

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
    sprintf(
      filepath, "%s/parameters_used%s",
      this->measurement_dir.c_str(), this->output_tag.c_str()
    );

    /// Create output file.
    FILE* used_param_fileptr;
    if (!(used_param_fileptr = fopen(filepath, "w"))) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Output directory '%s' does not exist.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->measurement_dir.c_str()
        );
      }
      return -1;
    }

    /// Print parameters to file.
    fprintf(
      used_param_fileptr, "catalogue_dir = %s\n",
      this->catalogue_dir.c_str()
    );
    fprintf(
      used_param_fileptr, "measurement_dir = %s\n",
      this->measurement_dir.c_str()
    );
    fprintf(
      used_param_fileptr, "data_catalogue_file = %s\n",
      this->data_catalogue_file.c_str()
    );
    fprintf(
      used_param_fileptr, "rand_catalogue_file = %s\n",
      this->rand_catalogue_file.c_str()
    );
    fprintf(
      used_param_fileptr, "output_tag = %s\n", this->output_tag.c_str()
    );

    fprintf(
      used_param_fileptr, "catalogue_type = %s\n",
      this->catalogue_type.c_str()
    );
    fprintf(
      used_param_fileptr, "measurement_type = %s\n",
      this->measurement_type.c_str()
    );

    fprintf(
      used_param_fileptr, "assignment = %s\n", this->assignment.c_str()
    );
    fprintf(
      used_param_fileptr, "norm_convention = %s\n",
      this->norm_convention.c_str()
    );

    fprintf(used_param_fileptr, "boxsize_x = %.2f\n", this->boxsize[0]);
    fprintf(used_param_fileptr, "boxsize_y = %.2f\n", this->boxsize[1]);
    fprintf(used_param_fileptr, "boxsize_z = %.2f\n", this->boxsize[2]);

    fprintf(used_param_fileptr, "ngrid_x = %d\n", this->ngrid[0]);
    fprintf(used_param_fileptr, "ngrid_y = %d\n", this->ngrid[1]);
    fprintf(used_param_fileptr, "ngrid_z = %d\n", this->ngrid[2]);

    fprintf(used_param_fileptr, "ell1 = %d\n", this->ell1);
    fprintf(used_param_fileptr, "ell2 = %d\n", this->ell2);
    fprintf(used_param_fileptr, "ELL = %d\n", this->ELL);

    fprintf(used_param_fileptr, "i_wa = %d\n", this->i_wa);
    fprintf(used_param_fileptr, "j_wa = %d\n", this->j_wa);

    fprintf(used_param_fileptr, "binning = %s\n", this->binning.c_str());
    fprintf(used_param_fileptr, "form = %s\n", this->form.c_str());

    fprintf(used_param_fileptr, "kmin = %.4f\n", this->kmin);
    fprintf(used_param_fileptr, "kmax = %.4f\n", this->kmax);
    fprintf(used_param_fileptr, "num_kbin = %d\n", this->num_kbin);
    fprintf(used_param_fileptr, "ith_kbin = %d\n", this->ith_kbin);

    fprintf(used_param_fileptr, "rmin = %.2f\n", this->rmin);
    fprintf(used_param_fileptr, "rmax = %.2f\n", this->rmax);
    fprintf(used_param_fileptr, "num_rbin = %d\n", this->num_rbin);
    fprintf(used_param_fileptr, "ith_rbin = %d\n", this->ith_rbin);

    fclose(used_param_fileptr);

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
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Catalogue type must be 'survey', 'mock' or 'sim': "
          "`catalogue_type` = %s.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->catalogue_type.c_str()
        );
      }
      return -1;
    }

    if (!(
      this->assignment == "NGP"
      || this->assignment == "CIC"
      || this->assignment == "TSC"
    )) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Grid assignment scheme must be 'NGP', 'CIC' or 'TSC': "
          "`assignment` = %s.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->assignment.c_str()
        );
      }
      return -1;
    }

    if (!(
      this->norm_convention == "grid" || this->norm_convention == "particle"
    )) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[WARN] (+%s) Normalisation convention must be "
          "'grid' or 'particle': `norm_convention` = %s.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->norm_convention.c_str()
        );
      }

      /// Set to default 'grid'.
      this->norm_convention == "grid";

      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[WARN] (+%s) Normalisation convention is set to "
          "default value '%s'.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->norm_convention.c_str()
        );
      }
    }

    if (!(this->form == "diag" || this->form == "full")) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) `form` must be either 'full' or 'diag': "
          "`form` = %s.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
          this->form.c_str()
        );
      }
      return -1;
    }

    /// Validate numerical parameters.
    if (this->num_kbin < 2 || this->num_rbin < 2) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Number of bins (`num_kbin` or `num_rbin`) "
          "must be >= 2.\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      return -1;
    }

    if (this->binning == "linpad" || this->assignment == "logpad") {
      /// HACK: Compare with `BinScheme::set_rbin`.
      int nbin_custom = 5;

      if (this->num_rbin < nbin_custom + 2) {
        if (currTask == 0) {
          clockElapsed = double(clock() - clockStart);
          printf(
            "[ERRO] (+%s) Binning scheme '%s' requires `num_rbin` >= %d.\n",
            calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
            this->binning.c_str(),
            nbin_custom + 2
          );
        }
        return -1;
      }
    }

    if (this->ith_kbin >= this->num_kbin || this->ith_rbin >= this->num_rbin) {
      if (currTask == 0) {
        clockElapsed = double(clock() - clockStart);
        printf(
          "[ERRO] (+%s) Bin index (`ith_kbin` or `ith_rbin`) must be less than "
          "the number of bins (`num_kbin` or `num_rbin`).\n",
          calc_elapsed_time_in_hhmmss(clockElapsed).c_str()
        );
      }
      return -1;
    }

    return 0;
  }

 private:
  /**
   * Set I/O paths and directories for data and random catalogue files
   * and processed results.
   */
  void set_io_files() {
    /// Set survey data and random catalogue inputs.
    if (this->catalogue_type == "survey" || this->catalogue_type == "sim") {
      this->data_catalogue_file =
        this->catalogue_dir + "/" + this->data_catalogue_file;
      this->rand_catalogue_file =
        this->catalogue_dir + "/" + this->rand_catalogue_file;
    }

    /// Set mock data and random catalogue inputs.  Make subdirectories
    /// to store the output for each data realisation.
    /// QUEST: Understand the numbering here.
    if (this->catalogue_type == "mock") {
      /// Enumerate the input data catalogue.
      int realisation = numTasks * this->batch_number + currTask + 1;

      if (realisation > 2048) {
        realisation = 1;
      }

      /// Set input paths.
      char buf_file[2048];

      sprintf(
        buf_file, "%s/%s_%04d.dat",
        this->catalogue_dir.c_str(),
        this->data_catalogue_file.c_str(),
        realisation
      );

      this->data_catalogue_file = buf_file;

      sprintf(
        buf_file, "%s/%s",
        this->catalogue_dir.c_str(),
        this->rand_catalogue_file.c_str()
      );

      this->rand_catalogue_file = buf_file;

      /// Set the output subdirectory.
      char buf_dir[2048];
      sprintf(buf_dir, "%s/%04d", this->measurement_dir.c_str(), realisation);

      struct stat status;  // file permission mode
      if (stat(buf_dir, &status) != 0) {
        /// Make the output subdirectory.
        int ret_status = mkdir(buf_dir, 0777);

        /// Check output subdirectory status and exit upon failure.
        if (ret_status == 0) {
          if (currTask == 0) {
            clockElapsed = double(clock() - clockStart);
            printf(
              "[INFO] (+%s) Output subdirectory '%s' is successfully made.\n",
              calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
              buf_dir
            );
          }
        } else {
          if (currTask == 0) {
            clockElapsed = double(clock() - clockStart);
            printf(
              "[ERRO] (+%s) Failed to make output subdirectory '%s'.\n",
              calc_elapsed_time_in_hhmmss(clockElapsed).c_str(),
              buf_dir
            );
          }
          exit(1);
        }
      }

      this->measurement_dir = buf_dir;
    }
  }
};

#endif  // TRIUMVIRATE_INCLUDE_PARAMETERS_HPP_INCLUDED_
