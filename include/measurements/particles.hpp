#ifndef TRIUM_PARTICLES_H_INCLUDED_
#define TRIUM_PARTICLES_H_INCLUDED_

#ifndef TRIUM_PARAMETERS_H_INCLUDED_
#include "parameters.hpp"
#endif

/**
 * Container of particle data and summary information.
 *
 */
class ParticleBOSSClass {  // FIXME: change class name
 public:
	struct ParticleData {
		double pos[3];  ///< particle position vector
		double w;  ///< particle weight
	}* particles;  ///< particle data

	int n_tot;  ///< total number of particles
	double pos_min[3];  ///< minimum values of particle positions
	double pos_max[3];  ///< maximum values of particle positions

	/**
	 * Return individual particle information.
	 *
	 * @param id Particle ID/index.
	 * @returns Particle information.
	 */
	ParticleData& operator[](int id) {
		return this->particles[id];
	}

	/**
	 * Initialise particle containers.
	 */
	ParticleBOSSClass () {  // FIXME: change class name
		this->particles = NULL;
		this->n_tot = 0;
		this->pos_min[0] = 0.; this->pos_max[0] = 0.;
		this->pos_min[1] = 0.; this->pos_max[1] = 0.;
		this->pos_min[2] = 0.; this->pos_max[2] = 0.;
	}

	/**
	 * Deconstruct particle containers.
	 */
	~ParticleBOSSClass() {  // FIXME: change class name
		finalise_particles();
	}

	/**
	 * Initialise particle data and information.
	 *
	 * @param num Number of particles.
	 */
	void initialise_particles(const int num) {
		/// Check total number of particles.
		if (num <= 0) {
			printf("Number of particles must be >= 0.\n");
			return;
		}
		this->n_tot = num;

		/// Renew particle data.
		delete[] particles; particles = NULL;
		this->particles = new ParticleData[this->n_tot];

		/// Determine memory usage.
		bytes += double(
			this->n_tot * sizeof(struct ParticleData) / 1024. / 1024. / 1024.
		);

		/// Initialise particle data values.
		for (int id = 0; id < this->n_tot; id++) {
			particles[id].pos[0] = 0.;
			particles[id].pos[1] = 0.;
			particles[id].pos[2] = 0.;
			particles[id].w = 0.;
		}
	}

	/**
	 * Finalise particle data and information.
	 */
	void finalise_particles() {
		if (particles != NULL) {  // ???: not this->particles?
			delete[] this->particles; this->particles = NULL;
			bytes -= double(
				this->n_tot * sizeof(struct ParticleData) / 1024. / 1024. / 1024.
			);
		}
	}

	/**
	 * Read in particle data from file.
	 *
	 * @param particles_file Input file path.
	 * @returns Exit code.
	 */
	int read_particles_BOSS(std::string& particles_file) {  // FIXME: change function name

		std::ifstream fin;

		/// Count for the number of lines/particles.
		int num_lines = 0;

		/// Read in data from file.
		fin.open(particles_file.c_str(), std::ios::in);

		if (fin.fail()) {
			printf("Cannot open file '%s'.\n", particles_file.c_str());
			fin.close();
			return -1;
		}

		std::string str_line;
		double x, y, z, w;
		while (getline(fin, str_line)) {
			if (sscanf(str_line.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &w) != 4) {
				continue;
			}
			num_lines++;
		}

		fin.close();

		/// Fill in particle data.
		this->initialise_particles(num_lines);

		fin.open(particles_file.c_str(), std::ios::in);

		num_lines = 0;
		while (getline(fin, str_line)) {
			if (sscanf(str_line.c_str(), "%lf %lf %lf %lf", &x, &y, &z, &w) != 4) {
				continue;
			}

			particles[num_lines].pos[0] = x;
			particles[num_lines].pos[1] = y;
			particles[num_lines].pos[2] = z;
			particles[num_lines].w = w;

			num_lines++;
		}

		fin.close();

		return 0;
	}

	/**
	 * Test the reading in of particle data from file.
	 *
	 * @param particles_file Input file path.
	 * @returns Exit code.
	 */
	int read_particles_test(std::string& particles_file) {

		std::ifstream fin;

		/// Count for the number of lines/particles.
		int num_lines = 0;

		/// Read in data from file.
		fin.open(particles_file.c_str(), std::ios::in);

		if (fin.fail()) {
			printf("Cannot open file '%s'.\n", particles_file.c_str());
			fin.close();
			return -1;
		}

		std::string str_line;
		double x, y, z, vx, vy, vz, mass, dummy;  // FIXME: Redundant?
		int ID, PID;  // FIXME: Redundant?
		while (getline(fin, str_line)) {
			if (sscanf(str_line.c_str(), "%lf %lf %lf", &x, &y, &z) != 3) {
				continue;
			}
			num_lines++;
		}

		fin.close();

		/// Fill in particle data.
		this->initialise_particles(num_lines);

		fin.open(particles_file.c_str(), std::ios::in);

		num_lines = 0;
		while (getline(fin, str_line)) {
			if (sscanf(str_line.c_str(), "%lf %lf %lf", &x, &y, &z) != 3) {
				continue;
			}

			this->particles[num_lines].pos[0] = x;
			this->particles[num_lines].pos[1] = y;
			this->particles[num_lines].pos[2] = z;
			this->particles[num_lines].w = 1.;

			num_lines++;
		}

		fin.close();

		/// Calculate extreme data values.
		this->calc_min_and_max();

		return 0;
	}

	/**
	 * Calculate extreme data values.
	 *
	 * @returns Exit code.
	 */
	int calc_min_and_max() {

		if (particles == NULL) {
			return -1;
		}

		/// Initialise minimum and maximum values with the first
		/// data entry/row.
		double min[3], max[3];

		min[0] = this->particles[0].pos[0]; max[0] = this->particles[0].pos[0];
		min[1] = this->particles[0].pos[1]; max[1] = this->particles[0].pos[1];
		min[2] = this->particles[0].pos[2]; max[2] = this->particles[0].pos[2];

		/// Find/update minimum and maximum values line-by-line.
		for (int id = 0; id < this->n_tot; id++) {
			if (min[0] > particles[id].pos[0]) {
				min[0] = particles[id].pos[0];
			}
			if (min[1] > particles[id].pos[1]) {
				min[1] = particles[id].pos[1];
			}
			if (min[2] > particles[id].pos[2]) {
				min[2] = particles[id].pos[2];
			}

			if (max[0] < particles[id].pos[0]) {
				max[0] = particles[id].pos[0];
			}
			if (max[1] < particles[id].pos[1]) {
				max[1] = particles[id].pos[1];
			}
			if (max[2] < particles[id].pos[2]) {
				max[2] = particles[id].pos[2];
			}
		}

		this->pos_min[0] = min[0]; this->pos_max[0] = max[0];
		this->pos_min[1] = min[1]; this->pos_max[1] = max[1];
		this->pos_min[2] = min[2]; this->pos_max[2] = max[2];

		return 0;
	}

	/**
	 * Calculate the alpha ratio (of weighted number counts or number densities).
	 *
	 * @param particles_data Source-data particle data.
	 * @param particles_rand Random particle data.
	 * @returns alpha Alpha ratio.
	 */
	static double calc_alpha_ratio(ParticleBOSSClass& particles_data, ParticleBOSSClass& particles_rand) {
		double num_wgt_data = 0.;
		for(int id = 0; id < particles_data.n_tot; id++) {
			num_wgt_data += particles_data[id].w;
		}

		double num_wgt_rand = 0.;
		for(int id = 0; id < particles_rand.n_tot; id++) {
			num_wgt_rand += particles_rand[id].w;
		}

		double alpha = num_wgt_data / num_wgt_rand;

		return alpha;
	}

	/**
	 * Calculate power spectrum normalisation.
	 *
	 * @param particles_data Source-data particle data.
	 * @param vol_survey Survey volume.
	 * @returns norm Normalisation constant.
	 */
	static double calc_norm_for_power_spectrum(ParticleBOSSClass& particles_data, double vol_survey) {
		double num_wgt_data = 0.0;
		for(int id = 0; id < particles_data.n_tot; id++) {
			num_wgt_data += particles_data[id].w;
		}

		double norm = vol_survey / num_wgt_data / num_wgt_data;

		return norm;
	}

	/**
	 * Calculate bispectrum normalisation.
	 *
	 * @param particles_data Source-data particle data.
	 * @param vol_survey Survey volume.
	 * @returns norm Normalisation constant.
	 */
	static double calc_norm_for_bispectrum(ParticleBOSSClass& particles_data, double vol_survey) {
		double num_wgt_data = 0.0;
		for(int id = 0; id < particles_data.n_tot; id++) {
			num_wgt_data += particles_data[id].w;
		}

		double norm = vol_survey / num_wgt_data / num_wgt_data;

		norm *= vol_survey / num_wgt_data;

		return norm;
	}

	/**
	 * Offset particle positions.
	 *
	 * @param dpos Position offset vector (subtractive).
	 * @returns Exit code.
	 */
	int offset_particles(const double* dpos) {
		if (particles == NULL) {
			return -1;
		}

		for (int id = 0; id < this->n_tot; id++) {
			this->particles[id].pos[0] -= dpos[0];
			this->particles[id].pos[1] -= dpos[1];
			this->particles[id].pos[2] -= dpos[2];
		}

		return 0;
	}

	/**
	 * Offset particle positions for periodic boundary conditions.
	 *
	 * @param params Input parameter set.
	 * @returns Exit code.
	 */
	int offset_particles_for_periodicity(ParameterSet& params) {
		for (int id = 0; id < this->n_tot; id++) {
			for (int axis = 0; axis < 3; axis++) {
				if (this->particles[id].pos[axis] >= params.boxsize[axis]) {
					this->particles[id].pos[axis] -= params.boxsize[axis];
				} else if (this->particles[id].pos[axis] < 0.) {
					this->particles[id].pos[axis] += params.boxsize[axis];
				}
			}
		}

		/// Recalculate extreme data values.
		this->calc_min_and_max();

		return 0;
	}

	/**
	 * Offset particle positions for FFTs (by grid adjustment).
	 *
	 * @param particles_data Source-data particle data.
	 * @param particles_rand Random particle data.
	 * @param params Input parameter set.
	 * @param factor Offset grid adjustment factor (default is 3.).
	 * @returns Exit code.
	 */
	static int offset_particles_for_fft(
			ParticleBOSSClass& particles_data,
			ParticleBOSSClass& particles_rand,
			ParameterSet& params,
			double factor=3.  /// ???: why 3.?
		) {
		particles_data.calc_min_and_max();
		particles_rand.calc_min_and_max();

		/// ???: Compensate by grid adjustment.
		double dpos[3] = {
			particles_rand.pos_min[0],
			particles_rand.pos_min[1],
			particles_rand.pos_min[2]
		};
		dpos[0] -= factor * params.boxsize[0] / double(params.nmesh[0]);
		dpos[1] -= factor * params.boxsize[1] / double(params.nmesh[1]);
		dpos[2] -= factor * params.boxsize[2] / double(params.nmesh[2]);

		particles_data.offset_particles(dpos);
		particles_rand.offset_particles(dpos);

		/// Recalculate extreme data values.
		particles_data.calc_min_and_max();
		particles_rand.calc_min_and_max();

		return 0;
	}

	/**
	 * Offset particle positions for window function calculations (by
	 * centring inside the pre-determined box).
	 *
	 * @param params Input parameter set.
	 * @returns Exit code.
	 */
	int offset_particles_for_window(ParameterSet& params) {
		double xmid = this->pos_min[0] + (this->pos_max[0] - this->pos_min[0]) / 2.;
		double ymid = this->pos_min[1] + (this->pos_max[1] - this->pos_min[1]) / 2.;
		double zmid = this->pos_min[2] + (this->pos_max[2] - this->pos_min[2]) / 2.;

		double dx[3] = {0., 0., 0.};
		dx[0] = params.boxsize[0]/2. - xmid;
		dx[1] = params.boxsize[1]/2. - ymid;
		dx[2] = params.boxsize[2]/2. - zmid;

		for (int id = 0; id < this->n_tot; id++) {
			for (int axis = 0; axis < 3; axis++) {
				this->particles[id].pos[axis] += dx[axis];
			}
		}

		/// Recalculate extreme data values.
		this->calc_min_and_max();

		return 0;
	}
};

#endif
