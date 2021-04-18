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
	struct particleInfoStruct {
		double pos[3];  ///< particle position vector
		double w;  ///< particle weight
	}* particles;  ///< particle data

	int n_tot;  ///< total number of particles
	double pos_min[3];  ///< minimum values of particle positions
	double pos_max[3];  ///< maximum values of particle positions

	/**
	 * Return individual particle information.
	 *
	 * @param id Particle ID.
	 * @returns Particle information.
	 */
	particleInfoStruct& operator[](int id) {  // FIXME: change function name
		return this->particles[id];
	}

	/**
	 * Initialise particle container.
	 */
	ParticleBOSSClass () {
		this->particles = NULL;
		this->n_tot = 0;
		this->pos_min[0] = 0.; this->pos_max[0] = 0.;
		this->pos_min[1] = 0.; this->pos_max[1] = 0.;
		this->pos_min[2] = 0.; this->pos_max[2] = 0.;
	}

	/**
	 * Deconstruct particle container.
	 */
	~ParticleBOSSClass() {
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
		this->particles = new particleInfoStruct[this->n_tot];

		/// Determine memory usage.
		bytes += double(
			this->n_tot * sizeof(struct particleInfoStruct) / 1024. / 1024. / 1024.
		);

		/// Initialise particle data values.
		for (int i = 0; i < this->n_tot; i++) {
			particles[i].pos[0] = 0.;
			particles[i].pos[1] = 0.;
			particles[i].pos[2] = 0.;
			particles[i].w = 0.;
		}
	}

	/**
	 * Finalise particle data and information.
	 */
	void finalise_particles() {
		if (particles != NULL) {  // ???: not this->particles?
			delete[] this->particles; this->particles = NULL;
			bytes -= double(
				this->n_tot * sizeof(struct particleInfoStruct) / 1024. / 1024. / 1024.
			);
		}
	}

	/**
	 * Read in particle data from file.
	 *
	 * @param fpath_in Input file path.
	 * @returns Exit code.
	 */
	int read_particles_BOSS(std::string& fpath_in) {

		std::ifstream fin;

		/// Count for the number of lines/particles.
		int num_lines = 0;

		/// Read in data from file.
		fin.open(fpath_in.c_str(), std::ios::in);

		if (fin.fail()) {
			printf("Cannot open file '%s'.\n", fpath_in.c_str());
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

		fin.open(fpath_in.c_str(), std::ios::in);

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
	 * @param fpath_in Input file path.
	 * @returns Exit code.
	 */
	int read_particles_test(std::string& fpath_in) {

		std::ifstream fin;

		/// Count for the number of lines/particles.
		int num_lines = 0;

		/// Read in data from file.
		fin.open(fpath_in.c_str(), std::ios::in);

		if (fin.fail()) {
			printf("Cannot open file '%s'.\n", fpath_in.c_str());
			fin.close();
			return -1;
		}

		std::string str_line;
		double x, y, z, vx, vy, vz, mass, dummy;
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

		fin.open(fpath_in.c_str(), std::ios::in);

		num_lines = 0;
		while (getline(fin, str_line)) {
			if (sscanf(str_line.c_str(), "%lf %lf %lf", &x, &y, &z) != 3) {
				continue;
			}

			this->particles[num_lines].pos[0] = x;
			this->particles[num_lines].pos[1] = y;
			this->particles[num_lines].pos[2] = z;
			this->particles[num_lines].w = 1.0;

			num_lines++;
		}

		fin.close();

		/// Calculate data value extremes.
		this->calc_min_and_max();

		return 0;
	}

	/**
	 * Calculate data value extremes.
	 *
	 * @returns Exit code.
	 */
	int calc_min_and_max() {

		if (particles == NULL) {
			return -1;
		}

		/// Initialise minimum and maximum values and set them to the first
		/// data entry/row.
		double min[3], max[3];

		min[0] = this->particles[0].pos[0]; max[0] = this->particles[0].pos[0];
		min[1] = this->particles[0].pos[1]; max[1] = this->particles[0].pos[1];
		min[2] = this->particles[0].pos[2]; max[2] = this->particles[0].pos[2];

		/// Find maximum and minimum line-by-line.
		for (int i = 0; i < this->n_tot; i++) {
			if (min[0] > particles[i].pos[0]) {
				min[0] = particles[i].pos[0];
			}
			if (min[1] > particles[i].pos[1]) {
				min[1] = particles[i].pos[1];
			}
			if (min[2] > particles[i].pos[2]) {
				min[2] = particles[i].pos[2];
			}

			if (max[0] < particles[i].pos[0]) {
				max[0] = particles[i].pos[0];
			}
			if (max[1] < particles[i].pos[1]) {
				max[1] = particles[i].pos[1];
			}
			if (max[2] < particles[i].pos[2]) {
				max[2] = particles[i].pos[2];
			}
		}
		this->pos_min[0] = min[0]; this->pos_max[0] = max[0];
		this->pos_min[1] = min[1]; this->pos_max[1] = max[1];
		this->pos_min[2] = min[2]; this->pos_max[2] = max[2];

		return 0;
	}

	int reset_particles(const double* dP) {

		if (particles == NULL) {
			return -1;
		}

		for (int p = 0; p < this->n_tot; p++) {
			this->particles[p].pos[0] -= dP[0];
			this->particles[p].pos[1] -= dP[1];
			this->particles[p].pos[2] -= dP[2];
		}

		return 0;
	}

	static double calcAlpha(ParticleBOSSClass& P_D, ParticleBOSSClass& P_R) {

		double num_D_weight = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			num_D_weight += P_D[p].w;
		}

		double num_R_weight = 0.0;
		for(int p = 0; p < P_R.n_tot; p++) {
			num_R_weight += P_R[p].w;
		}
		double alpha =  num_D_weight / num_R_weight;

		return alpha;
	}

	static double calcNormalizationForPowerSpectrum(ParticleBOSSClass& P_D, double Vsurvey) {

		double num_D_weight = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			num_D_weight += P_D[p].w;
		}

		double norm = Vsurvey / num_D_weight / num_D_weight;
		return norm;

	}

	static double calcNormalizationForBispectrum(ParticleBOSSClass& P_D, double Vsurvey) {

		double num_D_weight = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			num_D_weight += P_D[p].w;
		}

		double norm = Vsurvey / num_D_weight / num_D_weight;
		norm *= (Vsurvey / num_D_weight);
		return norm;

	}

	static int resetParticleForFFT(ParticleBOSSClass& P_D, ParticleBOSSClass& P_R, ParameterClass& param, double factor=3.0) {

		P_D.calc_min_and_max();
		P_R.calc_min_and_max();

		double dP[3] = {P_R.pos_min[0], P_R.pos_min[1], P_R.pos_min[2]};
		dP[0] -= factor * param.boxsize[0] / double(param.nmesh[0]);
		dP[1] -= factor * param.boxsize[1] / double(param.nmesh[1]);
		dP[2] -= factor * param.boxsize[2] / double(param.nmesh[2]);

		P_D.resetParticle(dP);
		P_R.resetParticle(dP);

		P_D.calc_min_and_max();
		P_R.calc_min_and_max();

		return 0;
	}

	int calcPeriodicBoundary(ParameterClass& param) {

		for(int p = 0; p < this->n_tot; p++) {

			for(int axes = 0; axes < 3; axes++) {
				if (this->particles[p].pos[axes] >= param.boxsize[axes]) {
					this->particles[p].pos[axes] -= param.boxsize[axes];
				} else if (this->particles[p].pos[axes] < 0.0) {
					this->particles[p].pos[axes] += param.boxsize[axes];
				}
			}
		}

		this->calc_min_and_max();

		return 0;
	}

	int resetParticlesForWindowFunction(ParameterClass& param) {

		double x_mid = this->pos_min[0] + (this->pos_max[0] - this->pos_min[0]) / 2.0;
		double y_mid = this->pos_min[1] + (this->pos_max[1] - this->pos_min[1]) / 2.0;
		double z_mid = this->pos_min[2] + (this->pos_max[2] - this->pos_min[2]) / 2.0;

		double dx[3] = {0.0, 0.0, 0.0};

		dx[0] = param.boxsize[0] / 2.0 - x_mid;
		dx[1] = param.boxsize[1] / 2.0 - y_mid;
		dx[2] = param.boxsize[2] / 2.0 - z_mid;

		for(int p = 0; p < this->n_tot; p++) {
			for(int axes = 0; axes < 3; axes++) {
				this->particles[p].pos[axes] += dx[axes];
			}
		}

		this->calc_min_and_max();

		return 0;
	}

 private:
};

#endif
