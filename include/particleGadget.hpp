#ifndef __particleGadget__
#define __particleGadget__

class ParticleBOSSClass {
private:
public:
	struct ParticleInfo {
		double pos[3]; // 3次元の場所。
		double vel[3]; // 3次元の場所。
		double w;
		long long ID;
	} * particles;

	int n_tot;         // 全粒子数
	double pos_max[3]; // 粒子の位置の最小値
	double pos_min[3]; // 粒子の位置の最大値

	double    redshift;
	double    scaleFactor;
	int       numFiles;
	double    boxsize;
	double    omegaM0;
	double    omegaLambda;
	double    hubble;
	double    mass;

	ParticleInfo & operator [] (int id) { return this->particles[id]; } // 粒子の定義

	ParticleBOSSClass() {
		/* initialize */
		this->particles = NULL;
		this->n_tot = 0;
		this->pos_max[0] = 0.0; this->pos_min[0] = 0.0;
		this->pos_max[1] = 0.0; this->pos_min[1] = 0.0;
		this->pos_max[2] = 0.0; this->pos_min[2] = 0.0;

		this->redshift = 0.0;
		this->scaleFactor = 0.0;
		this->numFiles = 0;
		this->boxsize = 0.0;
		this->omegaM0 = 0.0;
		this->omegaLambda = 0.0;
		this->hubble = 0.0;
		this->mass = 0.0;

	}

	~ParticleBOSSClass() {
		this->finalizeParticle();
	}

	struct io_header {
	  int      npart[6];
	  double   mass[6];
	  double   time;
	  double   redshift;
	  int      flag_sfr;
	  int      flag_feedback;
	  int      npartTotal[6];
	  int      flag_cooling;
	  int      num_files;
	  double   BoxSize;
	  double   Omega0;
	  double   OmegaLambda;
	  double   HubbleParam;
	  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 bytess */
	} header;


	void finalizeParticle() {
		if(particles != NULL) {
			delete [] this->particles; this->particles = NULL;
			bytes -= double( sizeof(struct ParticleInfo) * this->n_tot / 1024.0 / 1024.0 / 1024.0);
		}
	}

	int initialise_particles(const int num) {

		if(num <= 0) { printf("Number of particles is <= 0\n"); return -1; }

		/* insert the total number of particles in this->n_tot */
		this->n_tot = num;

		/* allocate particles */
		delete [] particles; particles == NULL;
		this->particles = new ParticleInfo[this->n_tot];

		/* compute memory */
		bytes += double( sizeof(struct ParticleInfo) * this->n_tot / 1024.0 / 1024.0 / 1024.0 );

		/* initialize particles */
		for(int i = 0; i < this->n_tot; i++) {
			particles[i].pos[0] = 0.0;
			particles[i].pos[1] = 0.0;
			particles[i].pos[2] = 0.0;

			particles[i].vel[0] = 0.0;
			particles[i].vel[1] = 0.0;
			particles[i].vel[2] = 0.0;

			particles[i].w = 1.0;

			particles[i].ID = i;
		}

		return 0;
	}


	int read_parameters(std::string & snapshotFileName) {
		FILE *fp;
		char buf[256];
		sprintf(buf, "%s.0", snapshotFileName.c_str());
		if(!(fp = fopen(buf,"rb"))){
			printf("can not open %s\n",buf);
			return -1;
		}

		int dummy = sizeof(header);
		fread(&dummy, sizeof(dummy),1,fp);
		fread(&header, sizeof(header),1,fp);
		fread(&dummy, sizeof(dummy),1,fp);
		fclose(fp);

		int n_tot_sum = 0;
		for(int n = 0; n< 6; n++){
			n_tot_sum += header.npartTotal[n];
		}

		this->n_tot = n_tot_sum;
		this->boxsize = header.BoxSize;
		this->omegaM0 = header.Omega0;
		this->omegaLambda = header.OmegaLambda;
		this->hubble = header.HubbleParam;
		this->redshift = header.redshift;
		this->scaleFactor = header.time;
		this->numFiles = header.num_files;

		double mass_sum = 0.0;
		for(int i = 0; i < 6; i++) {
			mass_sum += this->header.mass[i];
		}
		this->mass = mass_sum;

		return 0;
	}

	int readParticlesFromSnapshotOfGadget(std::string & snapshotFileName) {
		printf("reading snapshot files %s...\n",snapshotFileName.c_str());

		//****** ヘッダーファイルから、パラメータの読み込み。 *****//
		if( this->read_parameters(snapshotFileName)) { printf("fail to read snapshot files\n"); return -1;}

		//****** 粒子の初期化。 *****//
		if ( this->initialise_particles(this->n_tot) ) { printf("fail to read snapshot files\n"); return -1;}

		//****** 粒子の読み込み。 *****//
		int pc = 0, pc_new = 0;
		for(int i = 0; i < this->numFiles; i++, pc = pc_new) {
			FILE *fp;
			char buf[256];
			sprintf(buf,"%s.%d", snapshotFileName.c_str(), i);
			if(!(fp = fopen(buf,"rb"))) {
				printf("can not open %s\n",buf);
				return -1;
			}

			int dummy = sizeof(header);
			fread(&dummy,sizeof(dummy),1,fp);
			fread(&header,sizeof(header),1,fp);
			fread(&dummy,sizeof(dummy),1,fp);

			int NumPart = 0;
			for(int k=0; k< 6; k++){
				NumPart += header.npart[k];
			}

			float *block = new float[3*NumPart];
			long long *blockid = new long long[NumPart];
			dummy = sizeof(float) * 3 * NumPart;
			fread(&dummy,sizeof(dummy),1,fp);
			fread(block,sizeof(float), 3*NumPart, fp);
			fread(&dummy,sizeof(dummy),1,fp);
			pc_new = pc;
			for(int n = 0; n < NumPart; n++) {
			    	for(int k = 0; k < 3; k++) {
				      particles[pc_new].pos[k] =(double) block[3 * n + k];
				}
				pc_new++;
			}

			fread(&dummy,sizeof(dummy),1,fp);
			fread(block,sizeof(float), 3*NumPart, fp);
			fread(&dummy,sizeof(dummy),1,fp);
			pc_new = pc;
			for(int n = 0; n < NumPart; n++) {
			    	for(int k = 0; k < 3; k++) {
				      particles[pc_new].vel[k] = (double) block[3 * n + k] * sqrt(scaleFactor);
				}
				pc_new++;
			}

			dummy = sizeof(long long) * NumPart;
			fread(&dummy,sizeof(dummy),1,fp);
			fread(blockid,sizeof(long long), NumPart, fp);
			fread(&dummy,sizeof(dummy),1,fp);
			pc_new = pc;
			for(int n = 0; n < NumPart; n++) {
				particles[pc_new].ID = (long long) blockid[n];
				pc_new++;
			}
			if(i == numFiles-1) {
			if( pc_new != this->n_tot) {
				printf("fail to read snapshot files...: NumPart == %ld, nTotal = %ld\n",pc_new, this->n_tot);
				return -1;
			}}
			fclose(fp);
		}
		return 0;
	}

	int readParticlesFromRockstar(std::string & fname_in) {

		std::ifstream fin;

		/*****************************************/
		/* count the number of lines (particles) */
		int num_lines = 0;

		fin.open(fname_in.c_str(), std::ios::in);
		if( fin.fail() ) {
			printf("can not open file '%s'...\n", fname_in.c_str());
			fin.close();
			return -1;
		}

		std::string str;
		double x, y, z, vx, vy, vz, mass, dummy;
		int ID, PID;
		while(getline(fin, str)) {
			if( sscanf(str.c_str(), "%d  %lf %lf %lf %lf %lf %lf %lf %lf %lf\
					         %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\
					         %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\
						 %lf %lf %lf %d",
						 &ID,    &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &x,     &y,
						 &z,     &vx,    &vy,    &vz,    &dummy, &dummy, &dummy, &dummy, &dummy, &mass,
						 &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,
						 &dummy, &dummy, &dummy, &PID
						 ) != 34 ) {

				continue;
			}

			if( (12.9 < log10(mass)) && (log10(mass) < 13.1) ) {
				num_lines++;
			}
		}
		fin.close();
		/*****************************************/

		/* initialize particles */
		this->initialise_particles(num_lines);

		/*****************************************/
		/* read particle information */
		num_lines = 0;
		fin.open(fname_in.c_str(), std::ios::in);
		while(getline(fin, str)) {
			if( sscanf(str.c_str(), "%d  %lf %lf %lf %lf %lf %lf %lf %lf %lf\
					         %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\
					         %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\
					         %lf %lf %lf %d",
						 &ID,    &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &x,     &y,
						 &z,     &vx,    &vy,    &vz,    &dummy, &dummy, &dummy, &dummy, &dummy, &mass,
						 &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy,
						 &dummy, &dummy, &dummy, &PID
						 ) != 34 ) {
				continue;
			}
			if( (12.9 < log10(mass)) && (log10(mass) < 13.1) ) {
				particles[num_lines].pos[0] = x;
				particles[num_lines].pos[1] = y;
				particles[num_lines].pos[2] = z;
				particles[num_lines].vel[0] = vx;
				particles[num_lines].vel[1] = vy;
				particles[num_lines].vel[2] = vz;
				particles[num_lines].ID = num_lines;
				num_lines++;

			}
		}
		fin.close();
		/*****************************************/
		return 0;
	}

	int readParticlesRandom(std::string & fname_in) {

		std::ifstream fin;

		/*****************************************/
		/* count the number of lines (particles) */
		int num_lines = 0;

		fin.open(fname_in.c_str(), std::ios::in);
		if( fin.fail() ) {
			printf("can not open file '%s'...\n", fname_in.c_str());
			fin.close();
			return -1;
		}

		std::string str;
		double x, y, z, vx, vy, vz, mass, dummy;
		int ID, PID;
		while(getline(fin, str)) {
			if( sscanf(str.c_str(), "%lf %lf %lf",  &x, &y, &z) != 3 ) {
				continue;
			}

			num_lines++;

		}
		fin.close();
		/*****************************************/

		/* initialize particles */
		this->initialise_particles(num_lines);

		/*****************************************/
		/* read particle information */
		num_lines = 0;
		fin.open(fname_in.c_str(), std::ios::in);
		while(getline(fin, str)) {
			if( sscanf(str.c_str(), "%lf %lf %lf",  &x, &y, &z) != 3 ) {
				continue;
			}
			this->particles[num_lines].pos[0] = x;
			this->particles[num_lines].pos[1] = y;
			this->particles[num_lines].pos[2] = z;
			this->particles[num_lines].w = 1.0;
			num_lines++;
		}
		fin.close();
		/*****************************************/

		this->calcMinAndMax();

		return 0;
	}


	int setSubSamples(ParticleBOSSClass & P_in, int i_box, int j_box, int k_box) {

	    int num_lines = 0;
	    for(int i = 0; i < P_in.n_tot; i++) {
		    if((P_in[i].pos[0] > 500.0*double(i_box)  ) && (P_in[i].pos[1] > 500.0*double(j_box)   ) && (P_in[i].pos[2] > 500.0*double(k_box)  ) &&
		       (P_in[i].pos[0] < 500.0*double(i_box+1)) && (P_in[i].pos[1] < 500.0*double(j_box+1) ) && (P_in[i].pos[2] < 500.0*double(k_box+1))) {
			num_lines++;
		    }
	    }

	    /* initialize particles */
	    this->initialise_particles(num_lines);

	    num_lines = 0;
	    for(int i = 0; i < P_in.n_tot; i++) {
		    if((P_in[i].pos[0] > 500.0*double(i_box)  ) && (P_in[i].pos[1] > 500.0*double(j_box)  ) && (P_in[i].pos[2] > 500.0*double(k_box)  ) &&
		       (P_in[i].pos[0] < 500.0*double(i_box+1)) && (P_in[i].pos[1] < 500.0*double(j_box+1)) && (P_in[i].pos[2] < 500.0*double(k_box+1)) ) {

			this->particles[num_lines].pos[0] = P_in[i].pos[0];
			this->particles[num_lines].pos[1] = P_in[i].pos[1];
			this->particles[num_lines].pos[2] = P_in[i].pos[2];
			this->particles[num_lines].vel[0] = P_in[i].vel[0];
			this->particles[num_lines].vel[1] = P_in[i].vel[1];
			this->particles[num_lines].vel[2] = P_in[i].vel[2];
			this->particles[num_lines].ID = num_lines;
			num_lines++;
		    }
	    }

	}

	void print_parameters() {
		printf("****************************************\n");
		printf("Parameters\n");
		printf("****************************************\n");
		printf("total number of particles =%ld\n", this->n_tot);
		printf("Time = %lf, redshift = %lf\n", this->scaleFactor, this->redshift);
		printf("BoxSize = %lf [Mpc/h]\n", this->boxsize);
		printf("Omega = %lf, Omega_lambda = %lf, Hubble = %lf\n", this->omegaM0, this->omegaLambda, this->hubble);
		printf("mass = %lf [10^10 Msun/h]\n", this->mass);
		printf("****************************************\n");
	}

	int setRSD(int axis) {

		if( (axis > 3) || (axis < 0) ) { return -1; }

		const double a = this->scaleFactor;
		double hubble = 100.0 * sqrt(this->omegaM0 / pow(a,3) + (1.0 - this->omegaM0 - this->omegaLambda) / pow(a,2) + this->omegaLambda);
		double aH = this->scaleFactor * hubble; // aH;
		for (int i = 0; i < this->n_tot; i++) {
			this->particles[i].pos[axis] += particles[i].vel[axis] / aH;
		}

		return 0;
	}

	int calcPeriodicBoundary(ParameterClass & param) {

		for(int p = 0; p < this->n_tot; p++) {

			for(int axes = 0; axes < 3; axes++) {
				if(particles[p].pos[axes] >= param.boxsize[axes]) {
					particles[p].pos[axes] -= param.boxsize[axes];
				} else if(particles[p].pos[axes] < 0.0) {
					particles[p].pos[axes] += param.boxsize[axes];
				}
			}
		}

		return 0;
	}

	int calcMinAndMax() {

		if( particles == NULL) { return -1; }

		double min[3], max[3];

		min[0] = this->particles[0].pos[0]; max[0] = this->particles[0].pos[0];
		min[1] = this->particles[0].pos[1]; max[1] = this->particles[0].pos[1];
		min[2] = this->particles[0].pos[2]; max[2] = this->particles[0].pos[2];

		for (int i = 0; i < this->n_tot; i++) {
			if(min[0] > particles[i].pos[0]) {
				min[0] = particles[i].pos[0];
			}
			if(min[1] > particles[i].pos[1]) {
				min[1] = particles[i].pos[1];
			}
			if(min[2] > particles[i].pos[2]) {
				min[2] = particles[i].pos[2];
			}

			if(max[0] < particles[i].pos[0]) {
				max[0] = particles[i].pos[0];
			}
			if(max[1] < particles[i].pos[1]) {
				max[1] = particles[i].pos[1];
			}
			if(max[2] < particles[i].pos[2]) {
				max[2] = particles[i].pos[2];
			}
		}
		this->pos_min[0] = min[0]; this->pos_max[0] = max[0];
		this->pos_min[1] = min[1]; this->pos_max[1] = max[1];
		this->pos_min[2] = min[2]; this->pos_max[2] = max[2];

		return 0;
	}

	int resetParticle() {

		if( particles == NULL) { return -1; }
		for (int p = 0; p < this->n_tot; p++) {
			this->particles[p].pos[0] -= this->pos_min[0];
			this->particles[p].pos[1] -= this->pos_min[1];
			this->particles[p].pos[2] -= this->pos_min[2];
		}

		return 0;
	}


//	int calcMinAndMax() {
//
//		if( this->particles == NULL) { return -1; }
//
//		double min[3], max[3];
//
//		min[0] = this->particles[0].pos[0]; max[0] = this->particles[0].pos[0];
//		min[1] = this->particles[0].pos[1]; max[1] = this->particles[0].pos[1];
//		min[2] = this->particles[0].pos[2]; max[2] = this->particles[0].pos[2];
//
//		for (int i = 0; i < this->nTotal; i++) {
//			if(min[0] > particles[i].pos[0]) {
//				min[0] = particles[i].pos[0];
//			}
//			if(min[1] > particles[i].pos[1]) {
//				min[1] = particles[i].pos[1];
//			}
//			if(min[2] > particles[i].pos[2]) {
//				min[2] = particles[i].pos[2];
//			}
//
//			if(max[0] < particles[i].pos[0]) {
//				max[0] = particles[i].pos[0];
//			}
//			if(max[1] < particles[i].pos[1]) {
//				max[1] = particles[i].pos[1];
//			}
//			if(max[2] < particles[i].pos[2]) {
//				max[2] = particles[i].pos[2];
//			}
//		}
//		this->posMin[0] = min[0]; this->posMax[0] = max[0];
//		this->posMin[1] = min[1]; this->posMax[1] = max[1];
//		this->posMin[2] = min[2]; this->posMax[2] = max[2];
//
//		printf("min_x = %lf \n", this->posMin[0]);
//		printf("min_y = %lf \n", this->posMin[1]);
//		printf("min_z = %lf \n", this->posMin[2]);
//
//		printf("max_x = %lf \n", this->posMax[0]);
//		printf("max_y = %lf \n", this->posMax[1]);
//		printf("max_z = %lf \n", this->posMax[2]);
//
//		return 0;
//	}
//
//
//
//	std::complex<double> calcSphericalHarmonics(const int ell, const int M, const double mu, const double phi) {
//
//		/* 虚数の定義 */
//		std::complex<double> _I_(0.0, 1.0);
//
//		// Plm は m>=0 しか値を取らない。m<0 の場合は別に考える。
//		//the normalized associated Legendre sqrt((2l+1)/4pi) * sqrt(l-m/l+m) P_lm
//		std::complex<double> Ylm = 0.0;
//		if( M >= 0 ) {
//			Ylm = gsl_sf_legendre_sphPlm(ell, M, mu) * exp( _I_ * double(M) * phi );
//		} else {
//			int s = - M;
//			Ylm = gsl_sf_legendre_sphPlm(ell, s, mu) * exp( _I_ * double(s) * phi );
//			Ylm = pow(-1.0, -s) * std::conj(Ylm);
//		}
//		return Ylm;
//	}
//
//	int setSphericalHarmonics(const int ell, const int m) {
//
//		if(abs(m) > ell) { return -1; }
//
//		for(int i = 0; i < this->nTotal; i++) {
//			double xmag_xyz2 = 0.0;
//			for(int axes = 0; axes < 3; axes++) {
//				xmag_xyz2 += particles[i].pos[axes] * particles[i].pos[axes];
//			}
//			double xmag_xyz = sqrt(xmag_xyz2);
//
//			/* mu の計算。*/
//			double mu = particles[i].pos[2] / xmag_xyz;
//			/* phi の計算 */
//			double phi = atan(particles[i].pos[1]/particles[i].pos[0]);
//			/* atan は -pi/2 < atan < pi/2 までしか範囲を取らないから，o <=phi <= 2pi に変換する。 */
//			if( particles[i].pos[0] > 0 && particles[i].pos[1] < 0 ) {
//				phi = phi + 2.0 * M_PI;
//			}
//			if( particles[i].pos[0] < 0 ) {
//				phi = phi + M_PI;
//			}
//			if(phi < 0.0 || phi > 2.0*M_PI) { return -1; }
//			if(mu < -1.0 || mu > 1.0) { return -1; }
//			/* Ylm の計算 */
//			std::complex<double> Ylm = this->calcSphericalHarmonics(ell, m, mu, phi);
//			/* Y の複素数でweightをかける。*/
//			Ylm = std::conj(Ylm);
//			/* Ylm に sqrt(4pi/(2ell+1)) をかけたものを使う。*/
//			Ylm *= sqrt(4.0 * M_PI/(2.0 * double(ell) + 1.0));
//			/* 粒子に代入。*/
//			particles[i].Ylm = Ylm;
//		}
//		return 0;
//	}
//

	static double calcNormalizationForPowerSpectrum(ParticleBOSSClass & P_D, ParameterClass & param) {

		double num_D_weight = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			num_D_weight += P_D[p].w;
		}

		double norm = param.volume / num_D_weight / num_D_weight;
		return norm;

	}

	static double calcNormalizationForBispectrum(ParticleBOSSClass & P_D, ParameterClass & param) {

		double num_D_weight = 0.0;
		for(int p = 0; p < P_D.n_tot; p++) {
			num_D_weight += P_D[p].w;
		}

		double norm = param.volume / num_D_weight / num_D_weight;
		norm *= (param.volume / num_D_weight);
		return norm;

	}

	static double calcAlpha(ParticleBOSSClass & P_D, ParticleBOSSClass & P_R) {

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



};

template <class TemplateParticle>
class ParticleGadgetClassReduced {
private:
	struct ParticleInfo {
		double pos[3]; // 3次元の場所。
		double vel[3]; // 3次元の速度。
		double Z;      // redshift
		double Weight; //weight
		int    ID;     // ID
		std::complex<double> Ylm; // multipole 計算のための、spherical harmonics
		double T; // kSZ温度
	} * particles;
public:
	const ParticleInfo & operator [] (int id) { return this->particles[id]; } // 粒子の定義。
	int nTotal; // トータル粒子数。
	double posMax[3]; // 粒子の場所の最大値。
	double posMin[3]; // 粒子の場所の最小値。
	double    redshift;
	double    scaleFactor;
	int       numFiles;
	double    boxsize;
	double    omegaM0;
	double    omegaLambda;
	double    hubble;
	double    mass;
/*********************************************/

	ParticleGadgetClassReduced() {
		// 変数を初期化。//
		this->particles = NULL;
		this->nTotal = 0;
		this->posMax[0] = 0.0; this->posMin[0] = 0.0;
		this->posMax[1] = 0.0; this->posMin[1] = 0.0;
		this->posMax[2] = 0.0; this->posMin[2] = 0.0;
		this->redshift = 0.0;
		this->numFiles = 0.0;
		this->boxsize = 0.0;
		this->omegaM0 = 0.0;
		this->omegaLambda = 0.0;
		this->hubble = 0.0;
		this->mass = 0.0;
	}
	~ParticleGadgetClassReduced() {
		//粒子をフリー//
		this->finalizeParticle();
	}
	void finalizeParticle() {
		delete [] this->particles; this->particles = NULL;
	}

	int initialise_particles(const int num) {

		/* n_tot の代入 */
		if(num < 0) { printf("Number of particles is < 0\n"); return -1; }
		this->nTotal = num;

		/* 粒子の定義，メモリ確保，初期化*/
		delete [] particles; particles = NULL;
		this->particles = new ParticleInfo[this->nTotal];
		for(int i = 0; i < this->nTotal; i++) {
			this->particles[i].pos[0] = 0.0;
			this->particles[i].pos[1] = 0.0;
			this->particles[i].pos[2] = 0.0;
			this->particles[i].vel[0] = 0.0;
			this->particles[i].vel[1] = 0.0;
			this->particles[i].vel[2] = 0.0;
			this->particles[i].Z = 0.0;
			this->particles[i].Weight = 0.0;
			this->particles[i].ID = 0;
			this->particles[i].Ylm = 0.0;
			this->particles[i].T = 0.0;
		}

		return 0;
	}

	int checkParticle() {
		if( this->particles == NULL ) {
			return -1;
		} else {
			return 0;
		}
	}

	int getNumberOfParticle() {
		return this->nTotal;
	}

	void print_parameters() {
		printf("****************************************\n");
		printf("Parameters\n");
		printf("****************************************\n");
		printf("total number of particle =%ld\n", this->nTotal);
		printf("Time = %lf, redshift = %lf\n", this->scaleFactor, this->redshift);
		printf("BoxSize = %lf [Mpc/h]\n", this->boxsize);
		printf("Omega = %lf, Omega_lambda = %lf, Hubble = %lf\n", this->omegaM0, this->omegaLambda, this->hubble);
		printf("mass = %lf [10^10 Msun/h]\n", this->mass);
		printf("****************************************\n");
	}

	int calcMinAndMax() {

		if( this->particles == NULL) { return -1; }

		double min[3], max[3];

		min[0] = this->particles[0].pos[0]; max[0] = this->particles[0].pos[0];
		min[1] = this->particles[0].pos[1]; max[1] = this->particles[0].pos[1];
		min[2] = this->particles[0].pos[2]; max[2] = this->particles[0].pos[2];

		for (int i = 0; i < this->nTotal; i++) {
			if(min[0] > particles[i].pos[0]) {
				min[0] = particles[i].pos[0];
			}
			if(min[1] > particles[i].pos[1]) {
				min[1] = particles[i].pos[1];
			}
			if(min[2] > particles[i].pos[2]) {
				min[2] = particles[i].pos[2];
			}

			if(max[0] < particles[i].pos[0]) {
				max[0] = particles[i].pos[0];
			}
			if(max[1] < particles[i].pos[1]) {
				max[1] = particles[i].pos[1];
			}
			if(max[2] < particles[i].pos[2]) {
				max[2] = particles[i].pos[2];
			}
		}
		this->posMin[0] = min[0]; this->posMax[0] = max[0];
		this->posMin[1] = min[1]; this->posMax[1] = max[1];
		this->posMin[2] = min[2]; this->posMax[2] = max[2];

		printf("min_x = %lf \n", this->posMin[0]);
		printf("min_y = %lf \n", this->posMin[1]);
		printf("min_z = %lf \n", this->posMin[2]);

		printf("max_x = %lf \n", this->posMax[0]);
		printf("max_y = %lf \n", this->posMax[1]);
		printf("max_z = %lf \n", this->posMax[2]);

		return 0;
	}







	int setParticle(TemplateParticle & particle) {

		/*************/
		this->nTotal = this->nTotal; // トータル粒子数は変えない。
		/*************/

		this->posMax[0] = particle.posMax[0]; // 粒子の場所の最大値。
		this->posMax[1] = particle.posMax[1]; // 粒子の場所の最大値。
		this->posMax[2] = particle.posMax[2]; // 粒子の場所の最大値。
		this->posMin[0] = particle.posMax[0]; // 粒子の場所の最小値。
		this->posMin[1] = particle.posMax[1]; // 粒子の場所の最小値。
		this->posMin[2] = particle.posMax[2]; // 粒子の場所の最小値。
		this->redshift = particle.redshift;
		this->scaleFactor = particle.scaleFactor;
		this->numFiles = particle.numFiles;
		this->boxsize = particle.boxsize;
		this->omegaM0 = particle.omegaM0;
		this->omegaLambda = particle.omegaLambda;
		this->hubble = particle.hubble;
		this->mass = particle.mass;


		/****************************************************/
		/* random generator の定義。*/
		int seed = 20170112;
		gsl_rng * random_generator = NULL;
		random_generator = gsl_rng_alloc(gsl_rng_mt19937);
		gsl_rng_set(random_generator, seed);
		/* シャッフルするための配列を定義*/
		int shuffle[particle.nTotal];
		/* 値の代入。*/
		for(int p = 0; p < particle.nTotal; p++) { shuffle[p] = p; }
		/* シャッフル */
		gsl_ran_shuffle(random_generator, shuffle, particle.nTotal, sizeof(int));
		/* random generator のフリー*/
		gsl_rng_free(random_generator); random_generator = NULL;
		/****************************************************/

		for(int p = 0; p < this->nTotal; p++) {
			int index = shuffle[p];
			this->particles[p].pos[0] = particle[index].pos[0];
			this->particles[p].pos[1] = particle[index].pos[1];
			this->particles[p].pos[2] = particle[index].pos[2];
			this->particles[p].vel[0] = particle[index].vel[0];
			this->particles[p].vel[1] = particle[index].vel[1];
			this->particles[p].vel[2] = particle[index].vel[2];
			this->particles[p].Z = particle[index].Z;
			this->particles[p].Weight = particle[index].Weight;
			this->particles[p].ID = particle[index].ID;
			this->particles[p].Ylm = particle[index].Ylm;
			this->particles[p].T = particle[index].T;
		}

		return 0;
	}





};

#endif

//
//
//struct ParticleInfo {
//	double pos[3];
//	double vel[3];
////	double acc[3];
//	double dis[3];
//	double ini[3];
//	long long id;
//	double mass;
////	double pot;
//	double rand;
//};
//
//class Particle {
//private:
//	struct io_header {
//	  int      npart[6];
//	  double   mass[6];
//	  double   time;
//	  double   redshift;
//	  int      flag_sfr;
//	  int      flag_feedback;
//	  int      npartTotal[6];
//	  int      flag_cooling;
//	  int      num_files;
//	  double   BoxSize;
//	  double   Omega0;
//	  double   OmegaLambda;
//	  double   HubbleParam;
//	  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 bytess */
//	} header;
//	std::string snapshot_file;
//public:
//
//	ParticleInfo *particles;
//	ParticleInfo & operator [] (long long id) { return particles[id]; }
//	long long n_tot;
//	int       n_sample;
//	double    scale_factor;
//	double    redshift;
//	int       num_files;
//	double    boxsize;
//	double    omega_m0;
//	double    omega_lambda;
//	double    hubble;
//
//	double number_density;
///*********************************************/
//
//	void setParameters() {
//		printf("setParameters\n");
//		long long n_tot_sum = 0;
//		for(int n=0; n< 6; n++){
//			n_tot_sum += header.npartTotal[n];
//		}
//
//		this->n_tot = n_tot_sum;
//		this->n_sample = (int) (pow(n_tot_sum,1.0/3.0) + 1.0e-5);
//		this->boxsize = header.BoxSize;
//		this->omega_m0 = header.Omega0;
//		this->omega_lambda = header.OmegaLambda;
//		this->hubble = header.HubbleParam;
//		this->redshift = header.redshift;
//		this->scale_factor = header.time;
//	}
//
//	Particle() {
//		printf("\n**************************\n");
//		printf("Particle\n");
//		particles = NULL;
//	}
//
//	~Particle() {
//		if(particles != NULL) {
//			printf("freeParticle\n");
//			delete [] particles; particles = NULL;
//		}
//	}
//
//	void setParticleAll(char * snapshot_file, int & bytes) {
//
//		setSnapshotFile(snapshot_file);
//		setHeader();
//		setParameters();
//		printHeader();
//		print_parameters();
//		allocateParticle(bytes);
//		setParticlesFromSnapshotOfGadget();
//		setParticleMass();
//		setLatticeInitialPositions();
//		setParticleDisplacementVector();
//		setPeriodicBoundary();
////		setGaussianRandomDistribution(time(NULL));
//		calcNumberDensity();
//
////		particle.setSnapshotFile(snapshot_file);// 粒子クラスを定義して, 読み込むべきsnapshotファイルをセットする。
////		particle.setHeader(); // ヘッダを読み込む。
////		particle.setParameters(); // 読み込んだヘッダをもとに、粒子クラスにパラメータを入力する。
////		particle.printHeader(); // ヘッダをプリントする。
////		particle.print_parameters(); // 正しくヘッダからパラメータの値を入力できているか確かめるため、パラメータをプリントする。ヘッダと同じ値になるべき。
////		particle.allocateParticle(particle.n_tot, bytes); // 粒子のメモリを割り付け。
////		particle.setParticlesFromSnapshotOfGadget();// スナップショットから、粒子の位置、速度、ID を読み込む。
////	//	particle.setParticlesFromSnapshotOfIC();// スナップショットから、粒子の位置、速度、ID を読み込む。
////	                                                  // 読み込むスナップショットによってここを一番注意して変更すべき。
////	                                                  // Gadget や IC のアウトプットにあわせて、随時変更。(例; id に大して、int <=> long long など)
////		particle.setLatticeInitialPositions();      // 粒子のid を使って、格子状の粒子の初期条件を再現する。
////	                                                  // これは ICs が格子状から始めていなければ意味がない。
////		particle.setParticleDisplacementVector();   // 粒子の位置と、格子初期位置を使って、displacement vector を計算する。線形理論では(ZAでは), 速度と比例する。
////		particle.setPeriodicBoundary(); // 粒子の位置がボックスサイズを超えている場合には、ボックスサイズの中に戻す。
////		particle.setGaussianRandomDistribution(snapshot + 12 * realisation); // ガウシアンランダム分布(正規分布)する、乱数をまく。mu = 0, sigma = 1; seed を引数にもつので、適当に realisation と snapshot を組み合わせる。
////	//	particle.writeParticles("result/particles_g.txt");
//
//	}
//
//	void setSnapshotFile(std::string snapshot_file) {
//		printf("setSnapshotFile\n");
//		this->snapshot_file = snapshot_file;
//	}
//
//	void calcNumberDensity() {
//		printf("calcNumberDensity\n");
//		this->number_density = (double) n_tot / pow(boxsize,3);
//	}
//
//	void freeParticle(int & bytes) {
//		printf("freeParticle\n");
//		bytes -= (int) (sizeof(ParticleInfo) * this->n_tot / 1024 / 1024);
//		if(particles != NULL) {delete [] particles; particles = NULL; }
//		printf("memomry = %ld Mb\n",bytes);
//	}
//
//	void setHeader() {
//		printf("setHeader\n");
//		FILE *fp;
//		char buf[256];
//		sprintf(buf,"%s.0", snapshot_file.c_str());
//		if(!(fp = fopen(buf,"rb"))){
//			fprintf(stderr,"can not open %s\n",buf);
//			endRun();
//		}
//
//		int dummy = sizeof(header);
//		fread(&dummy, sizeof(dummy),1,fp);
//		fread(&header, sizeof(header),1,fp);
//		fread(&dummy, sizeof(dummy),1,fp);
//		fclose(fp);
//
//	}
//
//	void printHeader() {
//		printf("------------------------------------------\n");
//		printf("number of particle = %ld %ld %ld %ld %ld %ld\n", header.npartTotal[0], header.npartTotal[1],
//				                                       header.npartTotal[2],header.npartTotal[3],header.npartTotal[4], header.npartTotal[5]);
//		printf("number of particle in npart= %ld %ld %ld %ld %ld %ld\n", header.npart[0], header.npart[1],
//				                                               header.npart[2],header.npart[3],header.npart[4], header.npart[5]);
//		printf("mass of particle = %lf %lf %lf %lf %lf %lf\n",
//				    header.mass[0],header.mass[1], header.mass[2],header.mass[3],header.mass[4],header.mass[5]);
//		printf("Time = %lf, redshift = %lf\n",header.time,header.redshift);
//		printf("number of files = %ld \n",header.num_files);
//		printf("BoxSize = %lf\n",header.BoxSize);
//		printf("Omega = %lf, Omega_lambda = %lf, Hubble = %lf\n", header.Omega0, header.OmegaLambda, header.HubbleParam);
//		printf("------------------------------------------\n");
//	}
//
//	void print_parameters() {
//		printf("------------------------------------------\n");
//		printf("total number of particle =%ld\n", n_tot);
//		printf("n_sample =%ld\n", n_sample);
//		printf("Time = %lf, redshift = %lf\n",scale_factor, redshift);
//		printf("BoxSize = %lf\n", boxsize);
//		printf("Omega = %lf, Omega_lambda = %lf, Hubble = %lf\n", omega_m0, omega_lambda, hubble);
//		printf("------------------------------------------\n");
//	}
//
//
//	void allocateParticle(int & bytes) {
//		printf("allocateParticle\n");
//		bytes += (int) (sizeof(ParticleInfo) * this->n_tot / 1024 / 1024);
//		printf("memomry = %ld Mb\n",bytes);
//		particles = new ParticleInfo[this->n_tot];
//	}
//
//	void writeParticles(std::string filename) {
//		printf("writeParticles\n");
//		FILE *fp;
//		fp = fopen(filename.c_str(),"w");
//		for(int i = 0; i < n_tot; i++) {
//			fprintf(fp, "%.7e \t %.7e \t %.7e \t %.7e \t %.7e \t %.7e \t %.7e \t %.7e \t %.7e \n",
//					 particles[i].pos[0], particles[i].pos[1] , particles[i].pos[2],
//      			       particles[i].vel[0], particles[i].vel[1] , particles[i].vel[2],
//      			       particles[i].ini[0], particles[i].ini[1] , particles[i].ini[2]);
//		}
//		fclose(fp);
//	}
//
//	long long getNumberOfParticle() { return n_tot; }
//
//	void setParticleDisplacementVector() {
//		printf("setParticleDisplacementVector\n");
//		for(int i = 0; i < n_tot; i++) {
//		for(int axes = 0; axes < 3; axes++) {
//			double disp = particles[i].pos[axes] - particles[i].ini[axes];
//			disp = (disp >   boxsize / 2.0) ? disp - boxsize : disp;
//			disp = (disp < - boxsize / 2.0) ? disp + boxsize : disp;
//			particles[i].dis[axes] = disp;
//		}}
//	}
//
//	void setParticlesFromSnapshotOfGadget() {
//		printf("setParticlesFromSnapshotOfGadget\n");
//		printf("reading snapshot files %s...\n",snapshot_file.c_str());
//		int pc = 0, pc_new = 0;
//		for(int i = 0; i < header.num_files; i++, pc = pc_new) {
//			FILE *fp;
//			char buf[256];
//			sprintf(buf,"%s.%d", snapshot_file.c_str(), i);
//			if(!(fp = fopen(buf,"rb"))) {
//				fprintf(stderr,"can not open %s\n",buf);
//				endRun();
//			}
//
//			int dummy = sizeof(header);
//			fread(&dummy,sizeof(dummy),1,fp);
//			fread(&header,sizeof(header),1,fp);
//			fread(&dummy,sizeof(dummy),1,fp);
//
//			int NumPart = 0;
//			for(int k=0; k< 6; k++){
//				NumPart += header.npart[k];
//			}
//
//
//			float *block = new float[3*NumPart];
//			long long *blockid = new long long[NumPart];
//			dummy = sizeof(float) * 3 * NumPart;
//			fread(&dummy,sizeof(dummy),1,fp);
//			fread(block,sizeof(float), 3*NumPart, fp);
//			fread(&dummy,sizeof(dummy),1,fp);
//			pc_new = pc;
//			for(int n = 0; n < NumPart; n++) {
//			    	for(int k = 0; k < 3; k++) {
//				      particles[pc_new].pos[k] =(double) block[3 * n + k];
//				}
//				pc_new++;
//			}
//
//			fread(&dummy,sizeof(dummy),1,fp);
//			fread(block,sizeof(float), 3*NumPart, fp);
//			fread(&dummy,sizeof(dummy),1,fp);
//			pc_new = pc;
//			for(int n = 0; n < NumPart; n++) {
//			    	for(int k = 0; k < 3; k++) {
//				      particles[pc_new].vel[k] = (double) block[3 * n + k] * sqrt(scale_factor);
//				}
//				pc_new++;
//			}
//
//			dummy = sizeof(long long) * NumPart;
//			fread(&dummy,sizeof(dummy),1,fp);
//			fread(blockid,sizeof(long long), NumPart, fp);
//			fread(&dummy,sizeof(dummy),1,fp);
//			pc_new = pc;
//			for(int n = 0; n < NumPart; n++) {
//				particles[pc_new].id = (long long) blockid[n];
//				pc_new++;
//			}
//			if(i == header.num_files-1) {
//			if( pc_new != n_tot) {
//				printf("fail to read snapshot files...: NumPart == %ld, n_tot = %ld\n",pc_new, n_tot);
//				endRun();
//			}}
//			fclose(fp);
//		}
//	}
//
//	void setParticlesFromSnapshotOfIC() {
//		printf("setParticlesFromSnapshotOfIC\n");
//
//		int pc = 0, pc_new = 0;
//		for(int i = 0; i < header.num_files; i++, pc = pc_new) {
//			if(i == 0) {
//				printf("\n------------------------------------------\n");
//				printf("reading snapshot files...\n\n");
//			}
//
//			FILE *fp;
//			char buf[256];
//			sprintf(buf,"%s.%d", snapshot_file.c_str(), i);
//			if(!(fp = fopen(buf,"rb"))) {
//				fprintf(stderr,"can not open %s\n",buf);
//				endRun();
//			}
//			int dummy = sizeof(header);
//			fread(&dummy,sizeof(dummy),1,fp);
//			fread(&header,sizeof(header),1,fp);
//			fread(&dummy,sizeof(dummy),1,fp);
//
//			int NumPart = 0;
//			for(int k=0; k< 6; k++){
//				NumPart += header.npart[k];
//			}
//
//			dummy = sizeof(float) * 3 * NumPart;
//			fread(&dummy,sizeof(dummy),1,fp);
//			pc_new = pc;
//			for(int k=0; k < 6 ; k++){
//				for(int n=0;n<header.npart[k];n++){
//					float x[3];
//					fread(x,sizeof(float), 3, fp);
//					particles[pc_new].pos[0] = (double) x[0];
//					particles[pc_new].pos[1] = (double) x[1];
//					particles[pc_new].pos[2] = (double) x[2];
//					pc_new++;
//				}
//			}
//			fread(&dummy,sizeof(dummy),1,fp);
//
//			fread(&dummy,sizeof(dummy),1,fp);
//			pc_new = pc;
//			for(int k=0; k < 6 ; k++){
//				for(int n=0;n<header.npart[k];n++){
//					float v[3];
//					fread(v,sizeof(float), 3, fp);
//					particles[pc_new].vel[0] = (double) v[0] * sqrt(scale_factor); // sqrt(time) * comv_vel-> time * comv_vel = physical vel
//					particles[pc_new].vel[1] = (double) v[1] * sqrt(scale_factor);
//					particles[pc_new].vel[2] = (double) v[2] * sqrt(scale_factor);
//					pc_new++;
//				}
//			}
//			fread(&dummy,sizeof(dummy),1,fp);
//
//			dummy = sizeof(long long) * NumPart;
//			fread(&dummy,sizeof(dummy),1,fp);
//			pc_new = pc;
//			for(int k=0; k < 6 ; k++){
//				for(int n=0;n<header.npart[k];n++){
//					long long id;
//					fread(&id,sizeof(long long), 1, fp);
//					particles[pc_new].id = id;
//					pc_new++;
//				}
//			}
//			fread(&dummy,sizeof(dummy),1,fp);
//
//			printf("NumPart = %d, pc_new = %d,TotNumPart = %d\n", NumPart, pc_new, n_tot);
//			fclose(fp);
//		}
//	}
//
//	void setParticleMass() {
//		printf("setParticleMass\n");
//		double mass_sum = 0.0;
//		for(int i = 0; i < 6; i++) {
//			mass_sum += header.mass[i];
//		}
//		for(long long i = 0; i < getNumberOfParticle(); i++) {
//			particles[i].mass = mass_sum;
//		}
//	}
//	void setLatticeInitialPositions() {
//		printf("setLatticeInitialPositions\n");
//		double *temp_x = new double[n_tot];
//		double *temp_y = new double[n_tot];
//		double *temp_z = new double[n_tot];
//		int count = 0;
//		for(int i = 0; i < n_sample; i++){
//	 	for(int j = 0; j < n_sample; j++){
//	      for(int k = 0; k < n_sample; k++){
//			temp_x[count] = (double) i * (boxsize / (double) n_sample);
//		      temp_y[count] = (double) j * (boxsize / (double) n_sample);
//		      temp_z[count] = (double) k * (boxsize / (double) n_sample);
//		      count++;
//		}}}
//
//		for(long long i = 0; i < n_tot; i++){
//			particles[i].ini[0] = temp_x[particles[i].id-1];
//			particles[i].ini[1] = temp_y[particles[i].id-1];
//			particles[i].ini[2] = temp_z[particles[i].id-1];
//		}
//		delete [] temp_x;
//		delete [] temp_y;
//		delete [] temp_z;
//	}
//
//	void setPeriodicBoundary() {
//		printf("setPeriodicBoundary\n");
//		for(long long i = 0; i < n_tot; i++) {
//		for(int axes = 0; axes < 3; axes++) {
//			if ( particles[i].pos[axes] >= boxsize) { particles[i].pos[axes] -= boxsize; }
//			if ( particles[i].pos[axes] <  0.0) { particles[i].pos[axes] += boxsize; }
//		}}
//	}
//
//	double getHubbleParameterTime(double a) {
//		return 100.0 * sqrt(omega_m0 / pow(a,3) + (1.0 - omega_m0 - omega_lambda) / pow(a, 2) + omega_lambda);
//	}
//
//	void setGaussianRandomDistribution(const int seed) {
//		printf("setGaussianRandomDistribution\n");
//		gsl_rng *random_generator;
//		random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
//		gsl_rng_set(random_generator, seed);
//		for(long long i = 0; i < n_tot; i++) {
//			double ampl = gsl_rng_uniform(random_generator);
//			ampl = sqrt(-2.0 * log(ampl));
//			double phase = 2.0 * M_PI * gsl_rng_uniform(random_generator);
//			double rand = ampl * sin(phase);
//			particles[i].rand = rand;
//		}
//		gsl_rng_free(random_generator);
//	}
//
////	void setSmoothedDisplacementVector() {
////
////		fftw_complex *cdata = fftw_alloc_complex(nmesh * nmesh * nmesh);
////		fftw_plan plan = fftw_plan_dft_3d(nmesh, nmesh, nmesh, cdata, cdata, FFTW_FORWARD, FFTW_MEASURE);
////
////		for(int axes1 = 0; axes1 < 3; axes1++) {
////			for(long long i = 0; i < (nmesh * nmesh * nmesh); i++) {
////				cdata[i][0] = 0.0;
////				cdata[i][1] = 0.0;
////			}
////
////			for(long long p = 0; p < particle.getNumberOfParticle(); p++) {
////				double w[3][3];
////				int iw[3][3];
////				for(int axes = 0; axes < 3; axes++) {
////					double xp = particle[p].ini[axes] * (double) nmesh / boxsize;
////					iw[1][axes] = (int) (xp + 0.5);
////					double dx = xp - (double) (iw[1][axes]);
////					w[0][axes] = 0.5 * (0.5 - dx) * (0.5 - dx);
////					w[1][axes] = 0.75 - dx * dx;
////					w[2][axes] = 0.5 * (0.5 + dx) * (0.5 + dx);
////					iw[0][axes] = iw[1][axes] - 1;
////					iw[2][axes] = iw[1][axes] + 1;
////					for(int i = 0; i < 3; i++) {
////						if(iw[i][axes] >= nmesh) { iw[i][axes] -= nmesh; }// >= nmesh is important!
////						if(iw[i][axes] < 0) { iw[i][axes] += nmesh; }
////					}
////				}
////
////				double prefac_dis = particle[p].dis[axes1];
////				for(int i = 0; i < 3; i++) {
////				for(int j = 0; j < 3; j++) {
////				for(int k = 0; k < 3; k++) {
////					int coord = ( iw[i][0] * nmesh + iw[j][1] ) * nmesh + iw[k][2];
////					cdata[coord][0] += prefac_dis * w[i][0] * w[j][1] * w[k][2];
////				}}}
////			} // p
////
////			fftw_execute(plan);
////
////			for(int i = 0; i < nmesh; i++) {
////			for(int j = 0; j < nmesh; j++) {
////			for(int k = 0; k < nmesh; k++) {
////				if((i == 0) && (j ==0) && (k==0)) { continue; }
////				long long coord = ( i * nmesh + j ) * nmesh + k;
////				kvec[0] = (i < nmesh/2) ? (double) i * dk : (double) (i - nmesh) * dk;
////				kvec[1] = (j < nmesh/2) ? (double) j * dk : (double) (j - nmesh) * dk;
////				kvec[2] = (k < nmesh/2) ? (double) k * dk : (double) (k - nmesh) * dk;
////				std::complex<double> rho_temp(cdata[coord][0], cdata[coord][1]);
////				if(space == 0) {
////					rho_r[ND][coord] += kvec[axes1] * kvec[axes2] * rho_temp;
////				} else if (space == 1) {
////					rho_s[ND][coord] += kvec[axes1] * kvec[axes2] * rho_temp;
////				}
////			}}}
////
////
////		}
////
////	}
//
//};
//


