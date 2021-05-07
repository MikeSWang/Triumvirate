#ifndef TRIUM_PARAMETERS_H_INCLUDED_
#define TRIUM_PARAMETERS_H_INCLUDED_

#ifndef TRIUM_COMMON_H_INCLUDED_
#include "common.hpp"
#endif

/**
 * Set of initial parameters.
 *
 * This reads parameters from a file, stores and prints out the extracted
 * parameters, and validates the input.
 *
 */
class ParameterSet {
 public:
	std::string catalogue_dir;  ///< catalogue directory
	std::string output_dir;  ///< output directory
	std::string data_catalogue_file;  ///< data catalogue file
	std::string rand_catalogue_file;  ///< random catalogue file

	std::string catalogue_type;  ///< type of catalogue: {'survey', 'mock', 'sim'}

	double boxsize[3];  ///< boxsize in each dimension
	double volume;  ///< box volume

	int nmesh[3];  ///< mesh number in each dimension
	int nmesh_tot;  ///< total number of meshes

	std::string assignment;  ///< grid assignment scheme: {'NGP', 'CIC', 'TSC'}

	int ell1;  ///< spherical degree associated with the first wavevector
	int ell2;  ///< spherical degree associated with the second wavevector
	int ELL;  ///< spherical degree associated with the line-of-sight vector
	          // NOTE: standard variable naming convention overriden

	double kmin;  ///< minimum wavenumber boundary
	double kmax;  ///< maximum wavenumber boundary
	int num_kbin;  ///< number of wavenumber bins

	double rmin;  ///< minimum separation-distance boundary
	double rmax;  ///< maximum separation-distance boundary
	int num_rbin;  ///< number of separation-distance bins

	std::string form;  ///< full or diagonal form for the bispectrum

	/// FIXME: The following attributes are irrelevant at the current
	/// development stage.
	int NR;  ///< ???
	int ith_kbin;  ///< ???
	int ith_rbin;  ///< ???

	std::string reconstruction;  ///< reconstruction flag

	double b1_fid;  ///< ???
	double RG;  ///< ???

	/**
	 * Read parameters from a file.
	 *
	 * @param argv Command-line arguments.
	 * @returns Parameter acceptance outcome.
	 */
	int read_parameters(char* argv[]) {
		/// Load the parameter file.
		std::string param_file = argv[1];
		std::ifstream fin(param_file.c_str());

		/// Initialise temporary variables to hold the extracted parameters.
		char catalogue_dir_[1024];
		char output_dir_[1024];
		char data_catalogue_file_[1024];
		char rand_catalogue_file_[1024];
		char catalogue_type_[16];
		char assignment_[16];
		char form_[16];
		char reconstruction_[16];

		double boxsize_x, boxsize_y, boxsize_z;  // boxsize in each dimension
		int nmesh_x, nmesh_y, nmesh_z;  // mesh number in each dimension

		/// Extract parameters from file contents by line parsing.
		std::string str_line;  // string representing the line being parsed
		char str_dummy[1024];  // string placeholder for irrelevant contents

		/// IDEA: Perhaps the following code block could be refactored
		///	using a C++ equivalent to Python's ``zip``.
		while (getline(fin, str_line)) {
				/// Check if the line is in the correct format for parameter
				/// parsing; process if yes, otherwise move on to the next line.
				if (
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, str_dummy) != 3
				) {
					continue;
				}

				if (str_line.find("catalogue_dir") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_dir_);
				}
				if (str_line.find("output_dir") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, output_dir_);
				}
				if (str_line.find("data_catalogue_file") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, data_catalogue_file_);
				}
				if (str_line.find("rand_catalogue_file") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, rand_catalogue_file_);
				}

				if (str_line.find("catalogue_type") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, catalogue_type_);
				}

				if (str_line.find("boxsize_x") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_x);
				}
				if (str_line.find("boxsize_y") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_y);
				}
				if (str_line.find("boxsize_z") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &boxsize_z);
				}

				if (str_line.find("nmesh_x") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &nmesh_x);
				}
				if (str_line.find("nmesh_y") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &nmesh_y);
				}
				if (str_line.find("nmesh_z") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &nmesh_z);
				}

				if (str_line.find("assignment") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, assignment_);
				}

				if (str_line.find("ell1") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ell1);
				}
				if (str_line.find("ell2") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ell2);
				}
				if (str_line.find("ELL") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ELL);
				}

				if (str_line.find("kmin") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->kmin);
				}
				if (str_line.find("kmax") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->kmax);
				}
				if (str_line.find("num_kbin") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->num_kbin);
				}

				if (str_line.find("rmin") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->rmin);
				}
				if (str_line.find("rmax") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->rmax);
				}
				if (str_line.find("num_rbin") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->num_rbin);
				}

				if (str_line.find("form") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, form_);
				}

				/// FIXME: The following parameters are irrelevant at the current
				/// development stage.
				if (str_line.find("ith_kbin") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ith_kbin);
				}
				if (str_line.find("ith_rbin") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->ith_rbin);
				}
				if (str_line.find("NR") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %d", str_dummy, str_dummy, &this->NR);
				}

				if (str_line.find("reconstruction") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %s", str_dummy, str_dummy, reconstruction_);
				}

				if (str_line.find("b1_fid") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->b1_fid);
				}
				if (str_line.find("RG") != std::string::npos) {
					sscanf(str_line.data(), "%s %s %lg", str_dummy, str_dummy, &this->RG);
				}
		}

		/// Store parameters as attributes.
		this->catalogue_dir = catalogue_dir_;
		this->output_dir = output_dir_;
		this->data_catalogue_file = data_catalogue_file_;
		this->rand_catalogue_file = rand_catalogue_file_;

		this->catalogue_type = catalogue_type_;

		this->boxsize[0] = boxsize_x;
		this->boxsize[1] = boxsize_y;
		this->boxsize[2] = boxsize_z;

		this->volume = boxsize_x * boxsize_y * boxsize_z;

		this->nmesh[0] = nmesh_x;
		this->nmesh[1] = nmesh_y;
		this->nmesh[2] = nmesh_z;

		this->nmesh_tot = nmesh_x * nmesh_y * nmesh_z;

		this->assignment = assignment_;
		this->form = form_;
		this->reconstruction = reconstruction_;

		return this->check_parameters();
	}

	/**
	 * Set I/O paths and directories for data and random catalogue files
	 * and processed results.
	 *
	 * @returns Exit status.
	 */
	int set_io_files() {
		/// Set survey data and random catalogue inputs.
		if (this->catalogue_type == "survey") {
			this->data_catalogue_file = this->catalogue_dir + "/" + this->data_catalogue_file;
			this->rand_catalogue_file = this->catalogue_dir + "/" + this->rand_catalogue_file;
		}

		/// Set mock data and random catalogue inputs. Make subdirectories
		/// to store the output for each data realisation.
		if (this->catalogue_type == "mock") {
			/// ???: Enumerate the input data catalogue.
			int realisation = numTasks * this->NR + thisTask + 1;

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

			/// Set output subdirectory.
			char buf_dir[2048];

			sprintf(buf_dir, "%s/%04d", this->output_dir.c_str(), realisation);

			struct stat st;
			if (stat(buf_dir, &st) != 0) {
				/// Make output subdirectory.
				int ret = mkdir(buf_dir, 0777);

				/// Check output subdirectory status and exit upon failure.
				if (ret == 0) {
					if (thisTask == 0) {
						printf("Output subdirectory '%s' is successfully made.\n", buf_dir);
					}
				} else {
					if (thisTask == 0) {
						printf("FATAL: Failed to make output subdirectory '%s'.\n", buf_dir);
					}
					exit(1);
				}
			}

			this->output_dir = buf_dir;
		}

		return 0;
	}

	/**
	 * Print out and check extracted parameters as valid input (i.e. used
	 * parameters).
	 *
	 * @returns Exit status.
	 */
	int check_parameters() {
		/// Print out extracted parameters to a file.
		FILE* used_param_file_ptr;

		char buf[1024];
		sprintf(buf, "%s/used_parameters", this->output_dir.c_str());

		if (!(used_param_file_ptr = fopen(buf, "w"))) {
			if (thisTask == 0) {
				printf(
					"Output directory '%s' does not exist.\n",
					this->output_dir.c_str()
				);
			}
			return -1;
		}

		this->print_parameters(used_param_file_ptr);

		fclose(used_param_file_ptr);

		/// Validate input parameters.
		if (!(this->assignment == "NGP" || this->assignment == "CIC" || this->assignment == "TSC")) {
			if (thisTask == 0) {
				printf(
					"Grid assignment scheme must be 'NGP', 'CIC' or 'TSC': `assignment` = %s.\n",
					this->assignment.c_str()
				);
			}
			return -1;
		}

		if (this->num_kbin < 2 || this->num_rbin < 2) {
			if (thisTask == 0) {
				printf("Number of bins (`num_kbin` or `num_rbin`) must be >= 2.\n");
			}
			return -1;
		}

		if (this->ith_kbin >= this->num_kbin || this->ith_rbin >= this->num_rbin) {
			if (thisTask == 0) {
				printf(
					"Bin index (`ith_kbin` or `ith_rbin`) must be less than "
					"the number of bins (`num_kbin` or `num_rbin`).\n"
				);
			}
			return -1;
		}

		if (!(this->form == "diag" || this->form == "full")) {
			if (thisTask == 0) {
				printf("`form` must be either 'full' or 'diag'.\n");
			}
			return -1;
		}

		return 0;
	}

	/**
	 * Print out extracted parameters which are used in a process.
	 *
	 * @param used_param_file_ptr Pointer to file containing used parameters.
	 * @returns Exit status.
	 */
	int print_parameters(FILE* used_param_file_ptr) {
		/// IDEA: Perhaps the following code block could be refactored
		///	using a C++ equivalent to Python's ``zip``.
		fprintf(used_param_file_ptr, "catalogue_dir = %s\n", this->catalogue_dir.c_str());
		fprintf(used_param_file_ptr, "output_dir = %s\n", this->output_dir.c_str());
		fprintf(used_param_file_ptr, "data_catalogue_file = %s\n", this->data_catalogue_file.c_str());
		fprintf(used_param_file_ptr, "rand_catalogue_file = %s\n", this->rand_catalogue_file.c_str());

		fprintf(used_param_file_ptr, "catalogue_type = %s\n", this->catalogue_type.c_str());

		fprintf(used_param_file_ptr, "boxsize_x = %.2f\n", this->boxsize[0]);
		fprintf(used_param_file_ptr, "boxsize_y = %.2f\n", this->boxsize[1]);
		fprintf(used_param_file_ptr, "boxsize_z = %.2f\n", this->boxsize[2]);

		fprintf(used_param_file_ptr, "nmesh_x = %d\n", this->nmesh[0]);
		fprintf(used_param_file_ptr, "nmesh_y = %d\n", this->nmesh[1]);
		fprintf(used_param_file_ptr, "nmesh_z = %d\n", this->nmesh[2]);

		fprintf(used_param_file_ptr, "assignment = %s\n", this->assignment.c_str());

		fprintf(used_param_file_ptr, "ell1 = %d\n", this->ell1);
		fprintf(used_param_file_ptr, "ell2 = %d\n", this->ell2);
		fprintf(used_param_file_ptr, "ELL = %d\n", this->ELL);

		fprintf(used_param_file_ptr, "kmin = %.4f\n", this->kmin);
		fprintf(used_param_file_ptr, "kmax = %.4f\n", this->kmax);
		fprintf(used_param_file_ptr, "num_kbin = %d\n", this->num_kbin);

		fprintf(used_param_file_ptr, "rmin = %.2f\n", this->rmin);
		fprintf(used_param_file_ptr, "rmax = %.2f\n", this->rmax);
		fprintf(used_param_file_ptr, "num_rbin = %d\n", this->num_rbin);

		fprintf(used_param_file_ptr, "form = %s\n", this->form.c_str());

		fprintf(used_param_file_ptr, "ith_kbin = %d\n", this->ith_kbin);
		fprintf(used_param_file_ptr, "ith_rbin = %d\n", this->ith_rbin);

		fprintf(used_param_file_ptr, "reconstruction = %s\n", this->reconstruction.c_str());
		fprintf(used_param_file_ptr, "b1_fid = %.4f\n", this->b1_fid);
		fprintf(used_param_file_ptr, "RG = %.4f\n", this->RG);

		return 0;
	}
};

#endif
