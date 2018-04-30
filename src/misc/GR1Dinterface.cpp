#include "Transport.h"
#include "GridGR1D.h"
#include "GR1DSpectrumArray.h"
#include "Neutrino_GR1D.h"
#include "mpi.h"
#include "omp.h"
#include "global_options.h"
#include <cmath>
#include <iostream>

using namespace std;

double smoothing_timescale;
int GR1D_recalc_every;
const double time_gf = 2.03001708e5;
double GR1D_tau_crit;

extern "C"
void initialize_gr1d_sedonu_(const double *x1i, const int* n_GR1D_zones, const int* M1_imaxradii, const int* ghosts1, Transport** sim){

	// initialize MPI parallelism
	int MPI_myID,MPI_nprocs;
	//MPI_Init(NULL,NULL);
	MPI_Comm_rank( MPI_COMM_WORLD, &MPI_myID );
	MPI_Comm_size( MPI_COMM_WORLD, &MPI_nprocs);
	const int rank0 = (MPI_myID == 0);

	int c_M1_imaxradii = *M1_imaxradii-1;
	if(rank0){
		cout << "# Initializing transport..." << endl;
		cout << "#   Using " << MPI_nprocs << " MPI ranks" << endl;
            #ifdef _OPENMP
		    #pragma omp parallel
		    #pragma omp single
		    cout << "#   Using " << omp_get_num_threads()  << " threads on each MPI rank." << endl;
            #endif
		cout << "#   M1_imaxradii = " << c_M1_imaxradii << endl;
		cout << "#   ghosts1 = " << *ghosts1 << endl;
		cout << "#   n_GR1D_zones = " << *n_GR1D_zones << endl;
		cout << "#   radius range: " << x1i[*ghosts1] << " - " << x1i[c_M1_imaxradii+1] << endl;
		cout << "#   next lowest radius: " << x1i[c_M1_imaxradii] << endl;
	}

	// declare the transport class and save it to the fortran module
	*sim = new Transport;

	// open up the lua parameter file
	Lua lua("param.lua");
	smoothing_timescale = lua.scalar<double>("smoothing_timescale");
	GR1D_recalc_every = lua.scalar<int>("GR1D_recalc_every");

	GR1D_tau_crit = lua.scalar<int>("GR1D_tau_crit");
	if(GR1D_tau_crit < 0) GR1D_tau_crit = INFINITY;
	PRINT_ASSERT(GR1D_tau_crit,>,1);

	// set up the grid
	const int nzones = *M1_imaxradii - *ghosts1;
	GridGR1D* grid = new GridGR1D;
	grid->initialize_grid(x1i,nzones,*ghosts1);
	(*sim)->grid = grid;

	// set up the transport module (includes the grid)
	(*sim)->init(&lua);
	for(unsigned s=0; s<(*sim)->species_list.size(); s++){
		Neutrino_GR1D* tmpSpecies =static_cast<Neutrino_GR1D*>((*sim)->species_list[s]);
		tmpSpecies->ghosts1 = *ghosts1;
		tmpSpecies->n_GR1D_zones = *n_GR1D_zones;
		tmpSpecies->sim = *sim;
	}
	//omp_set_dynamic(true);
	omp_set_num_threads(1);
}

extern "C"
void calculate_mc_closure_(double* q_M1, double* q_M1p, double* q_M1m, double* q_M1_extra,
		double* q_M1_2mom, double* q_M1p_2mom, double* q_M1m_2mom, double* q_M1_extra_2mom,
		const double* eas, const double* rho, const double* T, const double* Ye, const double* v1,
		const double* metricX, const int* iter, const double* dt, const double* rshock, Transport** sim){

	const int nr_GR1D     = static_cast<Neutrino_GR1D*>((*sim)->species_list[0])->n_GR1D_zones;
	const int nghost_GR1D = static_cast<Neutrino_GR1D*>((*sim)->species_list[0])->ghosts1;
	const int nr = (*sim)->grid->rho.size();
	const int ns = (*sim)->species_list.size();
	const int ne = (*sim)->grid->nu_grid_axis.size();
	const double fsmooth = *dt*GR1D_recalc_every/time_gf / smoothing_timescale;
	PRINT_ASSERT(fsmooth,<=,1.0);
	PRINT_ASSERT(fsmooth,>=,0.0);

	// set velocities and opacities
	static_cast<GridGR1D*>((*sim)->grid)->set_fluid(rho, T, Ye, v1);

	bool extract_MC[nr][ns][ne];
	for(unsigned s=0; s<(*sim)->species_list.size(); s++)
		static_cast<Neutrino_GR1D*>((*sim)->species_list[s])->set_eas_external(eas,GR1D_tau_crit,&(extract_MC[0][0][0]),*rshock);

	// do MC calculation
	char* OMP_NUM_THREADS = std::getenv("OMP_NUM_THREADS");
	int nthreads=-1;
	if(OMP_NUM_THREADS!=NULL) nthreads = std::stoi(OMP_NUM_THREADS);
	else nthreads = omp_get_max_threads();
	omp_set_num_threads(nthreads);
	if(*iter%GR1D_recalc_every == 0) (*sim)->step();
	omp_set_num_threads(1);

	// create array for new values so they can be smoothed
	double new_Prr_E[nr][ns][ne];
	double new_Ptt_E[nr][ns][ne];
	double new_Wrrr_Fr[nr][ns][ne];
	double new_Wttr_Fr[nr][ns][ne];


	// set the GR1D quantities
	int indlast = ne*ns*nr_GR1D;
	//#pragma omp for
	unsigned dir_ind[2];
	for(int s=0; s<ns; s++){
		GR1DSpectrumArray* tmpSpectrum = static_cast<GR1DSpectrumArray*>((*sim)->grid->distribution[s]);
		for(int z_ind=0; z_ind<nr; z_ind++){
			dir_ind[0] = z_ind;
			int indr = z_ind+nghost_GR1D;
			const double X = metricX[indr];
			//int inds = indr + s*nr_GR1D;
			for(int ie=0; ie<ne; ie++){
				dir_ind[1] = ie;
				if(extract_MC[z_ind][s][ie]){
					//int inde = inds + ie*ns*nr_GR1D;
					//int indexE    = inde + 0*indlast;
					//int indexFr   = inde + 1*indlast;
					//int indexPrr  = inde + 2*indlast;
					//int indexPtt  = inde + 0*indlast;
					//int indexWrrr = inde + 1*indlast;
					//int indexWttr = inde + 2*indlast;
					//int indexChi  = inde + 3*indlast;
					//int indexPrr_pm = inde + 2*indlast + 0*ne*ns*nr_GR1D*3;
					//int indexChi_pm = inde + 0*indlast + 0*ne*ns*nr_GR1D*1;

					// load up new arrays
					unsigned index = tmpSpectrum->data.direct_index(dir_ind);
					Tuple<double,6> tmp = tmpSpectrum->data[index];
					double Prr  = tmp[2];//q_M1[indexPrr];//
					double Ptt  = tmp[3];//q_M1_extra[indexPtt];//
					double Wrrr = tmp[4];//q_M1_extra[indexWrrr];//
					double Wttr = tmp[5];//q_M1_extra[indexWttr];//

					// enforce local GR consistency
					double P_constraint = Prr/X/X + 2.*Ptt;
					double inv_P_constraint = 1.0 / P_constraint;
					new_Prr_E[z_ind][s][ie] = inv_P_constraint<INFINITY ? Prr * inv_P_constraint : X*X;
					new_Ptt_E[z_ind][s][ie] = inv_P_constraint<INFINITY ? Ptt * inv_P_constraint : 0.0;
					//if(inv_P_constraint == INFINITY && z_ind<150) cout << "warning - inv_P_constraint==INFINITY ir=" << z_ind << " s=" << s << " ie=" << ie << endl;

					double W_constraint = X*X*Wrrr + 2.*Wttr;
					double inv_W_constraint = 1.0 / W_constraint;
					//double Wrrr_Fr = abs(W_constraint)>0 ? Wrrr / W_constraint : 1./X/X;
					//double Wttr_Fr = abs(W_constraint)>0 ? Wttr / W_constraint : 0.0;
					new_Wrrr_Fr[z_ind][s][ie] = abs(inv_W_constraint)<INFINITY ? Wrrr * inv_W_constraint : 1./X/X;
					new_Wttr_Fr[z_ind][s][ie] = abs(inv_W_constraint)<INFINITY ? Wttr * inv_W_constraint : 0.0;
					//if(abs(inv_W_constraint) == INFINITY && z_ind<150) cout << "warning - inv_W_constraint==INFINITY ir=" << z_ind << " s=" << s << " ie=" << ie << endl;
				}
			}
		}
	}

	for(int z_ind=0; z_ind<nr; z_ind++){
		int indr = z_ind+nghost_GR1D;
		const double X = metricX[indr];
		for(int s=0; s<ns; s++){
			int inds = indr + s*nr_GR1D;
			for(int ie=0; ie<ne; ie++){
				int inde = inds + ie*ns*nr_GR1D;
				//int indexE    = inde + 0*indlast;
				int indexFr   = inde + 1*indlast;
				int indexPrr  = inde + 2*indlast;
				int indexPtt  = inde + 0*indlast;
				int indexWrrr = inde + 1*indlast;
				int indexWttr = inde + 2*indlast;
				//int indexChi  = inde + 3*indlast;
				int indexPrr_pm = inde + 2*indlast + 0*ne*ns*nr_GR1D*3;
				//int indexChi_pm = inde + 0*indlast + 0*ne*ns*nr_GR1D*1;

				if(extract_MC[z_ind][s][ie]){
					// spatial smoothing
					double Prr_E   =   new_Prr_E[z_ind][s][ie];
					double Ptt_E   =   new_Ptt_E[z_ind][s][ie];
					double Wrrr_Fr = new_Wrrr_Fr[z_ind][s][ie];
					double Wttr_Fr = new_Wttr_Fr[z_ind][s][ie];

					// temporal smoothing
					Prr_E       = fsmooth*Prr_E                 + (1.-fsmooth)*q_M1[indexPrr];
					Ptt_E       = fsmooth*Ptt_E                 + (1.-fsmooth)*q_M1_extra[indexPtt];
					double Wrrr = fsmooth*Wrrr_Fr*q_M1[indexFr] + (1.-fsmooth)*q_M1_extra[indexWrrr];
					double Wttr = fsmooth*Wttr_Fr*q_M1[indexFr] + (1.-fsmooth)*q_M1_extra[indexWttr];

					// make consistent with GR again
					double P_constraint = Prr_E/X/X + 2.*Ptt_E;
					Prr_E = P_constraint>0 ? Prr_E / P_constraint : X*X;
					Ptt_E = P_constraint>0 ? Ptt_E / P_constraint : 0.0;
					if(Prr_E>X*X and (Prr_E-X*X)/(X*X)<TINY) Prr_E=X*X;

					double W_constraint = X*X*Wrrr + 2.*Wttr;
					Wrrr_Fr = abs(W_constraint)>0 ? Wrrr / W_constraint : 1./X/X;
					Wttr_Fr = abs(W_constraint)>0 ? Wttr / W_constraint : 0.0;
					if(X*X*Wrrr_Fr>1 and (X*X*Wrrr_Fr-1)<TINY) Wrrr_Fr=1./X/X;

					// write out the variables
					q_M1[indexPrr]       = Prr_E;
					q_M1_extra[indexPtt] = Ptt_E;
					q_M1p[indexPrr_pm]   = Prr_E;
					q_M1m[indexPrr_pm]   = Prr_E;
					q_M1_extra[indexWrrr]  = Wrrr;
					q_M1_extra[indexWttr]  = Wttr;
				}
				else{
					q_M1[indexPrr]        = q_M1_2mom[indexPrr];
					q_M1_extra[indexPtt]  = q_M1_extra_2mom[indexPtt];
					q_M1p[indexPrr_pm]    = q_M1p_2mom[indexPrr_pm];
					q_M1m[indexPrr_pm]    = q_M1m_2mom[indexPrr_pm];
					q_M1_extra[indexWrrr] = q_M1_extra_2mom[indexWrrr];
					q_M1_extra[indexWttr] = q_M1_extra_2mom[indexWttr];
				}

				// check that the results are reasonable
				PRINT_ASSERT(q_M1[indexPrr],>=,0);
				PRINT_ASSERT(q_M1[indexPrr]-1,<=,X*X-1);
				PRINT_ASSERT(q_M1_extra[indexPtt],>=,0);
				PRINT_ASSERT(q_M1_extra[indexPtt]-1,<=,0);
				PRINT_ASSERT(q_M1p[indexPrr_pm],>=,0);
				PRINT_ASSERT(q_M1p[indexPrr_pm]-1,<=,X*X-1);
				PRINT_ASSERT(q_M1m[indexPrr_pm],>=,0);
				PRINT_ASSERT(q_M1m[indexPrr_pm]-1,<=,X*X-1);

			}
		}
	}
	cout.flush();
}
