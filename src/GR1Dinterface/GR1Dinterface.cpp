#include "Transport.h"
#include "GridGR1D.h"
#include "GR1DSpectrumArray.h"
#include "Neutrino_GR1D.h"
#include "mpi.h"
#include "omp.h"
#include "global_options.h"
#include <iostream>

using namespace std;

double smoothing_timescale;
int GR1D_recalc_every;
const double time_gf = 2.03001708e5;

extern "C"
void initialize_gr1d_sedonu_(const double *x1i, const int* n_GR1D_zones, const int* M1_imaxradii, const int* ghosts1, Transport** sim){

	// initialize MPI parallelism
	int my_rank,n_procs;
	//MPI_Init(NULL,NULL);
	MPI_Comm_rank( MPI_COMM_WORLD, &my_rank );
	MPI_Comm_size( MPI_COMM_WORLD, &n_procs);
	const int rank0 = (my_rank == 0);


	int c_M1_imaxradii = *M1_imaxradii-1;
	if(rank0){
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

	// set up the grid
	const int nzones = *M1_imaxradii - *ghosts1;
	GridGR1D* grid = new GridGR1D;
	grid->initialize_grid(x1i,nzones,*ghosts1);
	(*sim)->grid = grid;

	// set up the transport module (includes the grid)
	(*sim)->init(&lua);
	for(int s=0; s<(*sim)->species_list.size(); s++){
		Neutrino_GR1D* tmpSpecies =static_cast<Neutrino_GR1D*>((*sim)->species_list[s]);
		tmpSpecies->ghosts1 = *ghosts1;
		tmpSpecies->n_GR1D_zones = *n_GR1D_zones;
	}
	//omp_set_dynamic(true);
	omp_set_num_threads(1);
}

extern "C"
void calculate_mc_closure_(double* q_M1, double* q_M1p, double* q_M1m,
		double* q_M1_extra, double* q_M1_extrap, double* q_M1_extram,
		const double* eas, const double* rho, const double* T, const double* Ye, const double* v1,
		const double* metricX, const int* iter, const double* dt, Transport** sim){

	const int nr_GR1D     = static_cast<Neutrino_GR1D*>((*sim)->species_list[0])->n_GR1D_zones;
	const int nghost_GR1D = static_cast<Neutrino_GR1D*>((*sim)->species_list[0])->ghosts1;
	const int nr = (*sim)->grid->z.size();
	const int ns = (*sim)->species_list.size();
	const int ne = (*sim)->species_list[0]->nu_grid.size();
	const double fsmooth = *dt/time_gf / smoothing_timescale;
	PRINT_ASSERT(fsmooth,<=,1.0);
	PRINT_ASSERT(fsmooth,>=,0.0);

	// set velocities and opacities
	static_cast<GridGR1D*>((*sim)->grid)->set_fluid(rho, T, Ye, v1);

	for(int s=0; s<(*sim)->species_list.size(); s++)
		static_cast<Neutrino_GR1D*>((*sim)->species_list[s])->set_eas_external(eas);

	// do MC calculation
	omp_set_num_threads(std::stoi(std::getenv("OMP_NUM_THREADS")));
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
	for(int z_ind=0; z_ind<nr; z_ind++){
		int indr = z_ind+nghost_GR1D;
		const double X = metricX[indr];
		for(int s=0; s<ns; s++){
			int inds = indr + s*nr_GR1D;
			GR1DSpectrumArray* tmpSpectrum = static_cast<GR1DSpectrumArray*>((*sim)->grid->z[z_ind].distribution[s]);
			for(int ie=0; ie<ne; ie++){
				int inde = inds + ie*ns*nr_GR1D;
				int indexE    = inde + 0*indlast;
				int indexFr   = inde + 1*indlast;
				int indexPrr  = inde + 2*indlast;
				int indexPtt  = inde + 0*indlast;
				int indexWrrr = inde + 1*indlast;
				int indexWttr = inde + 2*indlast;
				//int indexChi  = inde + 3*indlast;
				int indexPrr_pm = inde + 2*indlast + 0*ne*ns*nr_GR1D*3;
				int indexChi_pm = inde + 0*indlast + 0*ne*ns*nr_GR1D*1;

				// load up new arrays
				double Prr  = tmpSpectrum->moments[tmpSpectrum->index(ie,2)];//q_M1[indexPrr];//
				double Ptt  = tmpSpectrum->moments[tmpSpectrum->index(ie,3)];//q_M1_extra[indexPtt];//
				double Wrrr = tmpSpectrum->moments[tmpSpectrum->index(ie,4)];//q_M1_extra[indexWrrr];//
				double Wttr = tmpSpectrum->moments[tmpSpectrum->index(ie,5)];//q_M1_extra[indexWttr];//

				// enforce local GR consistency
				double P_constraint = Prr/X/X + 2.*Ptt;
				double inv_P_constraint = 1.0 / P_constraint;
				new_Prr_E[z_ind][s][ie] = inv_P_constraint<INFINITY ? Prr * inv_P_constraint : X*X;
				new_Ptt_E[z_ind][s][ie] = inv_P_constraint<INFINITY ? Ptt * inv_P_constraint : 0.0;
				//if(inv_P_constraint == INFINITY && z_ind<150) cout << "warning - inv_P_constraint==INFINITY ir=" << z_ind << " s=" << s << " ie=" << ie << endl;

				double W_constraint = X*X*Wrrr + 2.*Wttr;
				double inv_W_constraint = 1.0 / W_constraint;
				double Wrrr_Fr = abs(W_constraint)>0 ? Wrrr / W_constraint : 1./X/X;
				double Wttr_Fr = abs(W_constraint)>0 ? Wttr / W_constraint : 0.0;
				new_Wrrr_Fr[z_ind][s][ie] = abs(inv_W_constraint)<INFINITY ? Wrrr * inv_W_constraint : 1./X/X;
				new_Wttr_Fr[z_ind][s][ie] = abs(inv_W_constraint)<INFINITY ? Wttr * inv_W_constraint : 0.0;
				//if(abs(inv_W_constraint) == INFINITY && z_ind<150) cout << "warning - inv_W_constraint==INFINITY ir=" << z_ind << " s=" << s << " ie=" << ie << endl;
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
				int indexE    = inde + 0*indlast;
				int indexFr   = inde + 1*indlast;
				int indexPrr  = inde + 2*indlast;
				int indexPtt  = inde + 0*indlast;
				int indexWrrr = inde + 1*indlast;
				int indexWttr = inde + 2*indlast;
				//int indexChi  = inde + 3*indlast;
				int indexPrr_pm = inde + 2*indlast + 0*ne*ns*nr_GR1D*3;
				int indexChi_pm = inde + 0*indlast + 0*ne*ns*nr_GR1D*1;

				// spatial smoothing
				double Prr_E   =   new_Prr_E[z_ind][s][ie];
				double Ptt_E   =   new_Ptt_E[z_ind][s][ie];
				double Wrrr_Fr = new_Wrrr_Fr[z_ind][s][ie];
				double Wttr_Fr = new_Wttr_Fr[z_ind][s][ie];

				if(z_ind==0){
					Prr_E   = 0.5*(  new_Prr_E[z_ind][s][ie] +   new_Prr_E[z_ind+1][s][ie]);
					Ptt_E   = 0.5*(  new_Ptt_E[z_ind][s][ie] +   new_Ptt_E[z_ind+1][s][ie]);
					Wrrr_Fr = 0.5*(new_Wrrr_Fr[z_ind][s][ie] + new_Wrrr_Fr[z_ind+1][s][ie]);
					Wttr_Fr = 0.5*(new_Wttr_Fr[z_ind][s][ie] + new_Wttr_Fr[z_ind+1][s][ie]);
				}
				else if(z_ind==nr-1){
					Prr_E   = 0.5*(  new_Prr_E[z_ind][s][ie] +   new_Prr_E[z_ind-1][s][ie]);
					Ptt_E   = 0.5*(  new_Ptt_E[z_ind][s][ie] +   new_Ptt_E[z_ind-1][s][ie]);
					Wrrr_Fr = 0.5*(new_Wrrr_Fr[z_ind][s][ie] + new_Wrrr_Fr[z_ind-1][s][ie]);
					Wttr_Fr = 0.5*(new_Wttr_Fr[z_ind][s][ie] + new_Wttr_Fr[z_ind-1][s][ie]);
				}
				else{
					Prr_E   = 0.25*(  new_Prr_E[z_ind+1][s][ie] +   2.*new_Prr_E[z_ind][s][ie] +   new_Prr_E[z_ind-1][s][ie]);
					Ptt_E   = 0.25*(  new_Ptt_E[z_ind+1][s][ie] +   2.*new_Ptt_E[z_ind][s][ie] +   new_Ptt_E[z_ind-1][s][ie]);
					Wrrr_Fr = 0.25*(new_Wrrr_Fr[z_ind+1][s][ie] + 2.*new_Wrrr_Fr[z_ind][s][ie] + new_Wrrr_Fr[z_ind-1][s][ie]);
					Wttr_Fr = 0.25*(new_Wttr_Fr[z_ind+1][s][ie] + 2.*new_Wttr_Fr[z_ind][s][ie] + new_Wttr_Fr[z_ind-1][s][ie]);
				}

				// temporal smoothing
				Prr_E       = fsmooth*Prr_E                 + (1.-fsmooth)*q_M1[indexPrr];
				Ptt_E       = fsmooth*Ptt_E                 + (1.-fsmooth)*q_M1_extra[indexPtt];
				double Wrrr = fsmooth*Wrrr_Fr*q_M1[indexFr] + (1.-fsmooth)*q_M1_extra[indexWrrr];
				double Wttr = fsmooth*Wttr_Fr*q_M1[indexFr] + (1.-fsmooth)*q_M1_extra[indexWttr];

				// make consistent with GR again
				double P_constraint = Prr_E/X/X + 2.*Ptt_E;
				Prr_E = P_constraint>0 ? Prr_E / P_constraint : X*X;
				Ptt_E = P_constraint>0 ? Ptt_E / P_constraint : 0.0;

				double W_constraint = X*X*Wrrr + 2.*Wttr;
				Wrrr_Fr = abs(W_constraint)>0 ? Wrrr / W_constraint : 1./X/X;
				Wttr_Fr = abs(W_constraint)>0 ? Wttr / W_constraint : 0.0;

				// write out the variables
				q_M1[indexPrr]       = Prr_E;
				q_M1_extra[indexPtt] = Ptt_E;
				q_M1p[indexPrr_pm]   = Prr_E;
				q_M1m[indexPrr_pm]   = Prr_E;
				q_M1_extra[indexWrrr] = Wrrr_Fr * q_M1[indexFr];
				q_M1_extra[indexWttr] = Wttr_Fr * q_M1[indexFr];
			}
		}
	}
}
