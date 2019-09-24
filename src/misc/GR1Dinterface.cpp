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
void initialize_gr1d_sedonu_(const double *x1i, const int* n_GR1D_zones, const int* M1_imaxradii, const int* ghosts1,
			     const double* rho, const double* T, const double* Ye, const double* v1,	const double* metricX, const double* alp, Transport** sim){

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
	static_cast<GridGR1D*>((*sim)->grid)->set_fluid(rho, T, Ye, v1, metricX, alp);

	// set up the transport module (includes the grid)
	(*sim)->init(&lua);
	for(size_t s=0; s<(*sim)->species_list.size(); s++){
		Neutrino_GR1D* tmpSpecies =static_cast<Neutrino_GR1D*>((*sim)->species_list[s]);
		tmpSpecies->ghosts1 = *ghosts1;
		tmpSpecies->n_GR1D_zones = *n_GR1D_zones;
		tmpSpecies->sim = *sim;
	}
	//omp_set_dynamic(true);
#ifdef _OPENMP
	omp_set_num_threads(1);
#endif
}

extern "C"
void calculate_mc_closure_(
		double* q_M1,            // energy density, flux, pressure tensor
		double* q_M1p,           // same as above for plus state
		double* q_M1m,           // same as above for minus state
		double* q_M1_extra,      // P_theta^theta/E=P_phi^phi/E, W^rrr, W_phi^phir, Chi
		double* q_M1_2mom,       // same as above (temporary values previously computed by GR1D closure routine)
		double* q_M1p_2mom,      // ''
		double* q_M1m_2mom,      // ''
		double* q_M1_extra_2mom, // ''
		const double* eas,       // emissivity, absorption opacity, scattering opacity
		const double* rho,       // density
		const double* T,         // temperature
		const double* Ye,        // electron fraction
		const double* v1,        // radial velocity
		const double* metricX,   // g_rr component of metric
		const double* alp,      // lapse
		const int* iter,         // iteration number
		const double* dt,        // timestep in GR1D units
		const double* rshock,    // shock radius (currently unused in set_eas_external)
		Transport** sim){        // pointer to the Sedonu Transport class

	const int nr_GR1D     = static_cast<Neutrino_GR1D*>((*sim)->species_list[0])->n_GR1D_zones;
	const int nghost_GR1D = static_cast<Neutrino_GR1D*>((*sim)->species_list[0])->ghosts1;
	const int nr = (*sim)->grid->rho.size();
	const int ns = (*sim)->species_list.size();
	const int ne = (*sim)->grid->nu_grid_axis.size();
	const double fsmooth = *dt*GR1D_recalc_every/time_gf / smoothing_timescale;
	PRINT_ASSERT(fsmooth,<=,1.0);
	PRINT_ASSERT(fsmooth,>=,0.0);

	// set velocities and opacities
	static_cast<GridGR1D*>((*sim)->grid)->set_fluid(rho, T, Ye, v1, metricX, alp);

	bool extract_MC[nr][ns][ne];
	for(size_t s=0; s<(*sim)->species_list.size(); s++)
		static_cast<Neutrino_GR1D*>((*sim)->species_list[s])->set_eas_external(eas,GR1D_tau_crit,&(extract_MC[0][0][0]),*rshock);

	//===================//
	// do MC calculation //
	//===================//
#ifdef _OPENMP
	char* OMP_NUM_THREADS = std::getenv("OMP_NUM_THREADS");
	int nthreads=-1;
	if(OMP_NUM_THREADS!=NULL) nthreads = std::stoi(OMP_NUM_THREADS);
	else nthreads = omp_get_max_threads();
	omp_set_num_threads(nthreads);
#endif
	if(*iter%GR1D_recalc_every == 0) (*sim)->step();
#ifdef _OPENMP
	omp_set_num_threads(1);
#endif

	// create array for new values so they can be smoothed
	double new_Prr_E[nr][ns][ne];
	double new_Ptt_E[nr][ns][ne];
	double new_Wrrr[nr][ns][ne];
	double new_Wttr[nr][ns][ne];

	// set the GR1D quantities
	const int indlast = ne*ns*nr_GR1D;
	for(int s=0; s<ns; s++){
		GR1DSpectrumArray* tmpSpectrum = static_cast<GR1DSpectrumArray*>((*sim)->grid->distribution[s]);
//#pragma omp parallel for
		for(int z_ind=0; z_ind<nr; z_ind++){
			size_t dir_ind[2];
			dir_ind[0] = z_ind;
			int indr = z_ind+nghost_GR1D;
			const double X = metricX[indr];
			const double vr = v1[indr];
			const double V = vr*X;
			const double WLorentz = 1./sqrt(1.-V*V);
			const double WX = WLorentz*X;

			//int inds = indr + s*nr_GR1D;
			for(int ie=0; ie<ne; ie++){
				dir_ind[1] = ie;
				int inde = s + ie*ns*nr_GR1D;
				int indexE    = inde + 0*indlast;
				int indexFr   = inde + 1*indlast;

				// load the old lab-frame moments
				double E = q_M1[indexE];
				double Fr = q_M1[indexFr]*E;

				if(false){//q_M1[indexFr] > 1e-3 /*extract_MC[z_ind][s][ie]*/){

					// load up new arrays (tetrad-frame moments)
					size_t index = tmpSpectrum->data.direct_index(dir_ind);
					Tuple<double,6> tmp = tmpSpectrum->data[index];
					double J = tmp[0];
					double Hr = tmp[1];
					double Lrr = tmp[2];
					double Ltt = tmp[3];
					double Nrrr = tmp[4];
					double Nttr = tmp[5];

					// get the parameters that make the GR1D and MC moments match
					double denom = (J - Lrr*V*V)*X;
					double alpha = ( E*X*(1.+V)*(1.*V) - 2.*Fr*V     ) / (   denom);
					double beta  = ( Fr*(J+Lrr*V*V)/X  - E*V*(J+Lrr) ) / (Hr*denom);
					J *= alpha;
					Hr *= beta;
					Lrr *= alpha;
					Ltt *= alpha;
					Nrrr *= beta;
					Nttr *= beta;

					// set the lab-frame quantities
					double Prr = WLorentz*WLorentz*X*X * (Lrr + V*(2.*Hr + J*V));
					double Ptt = Ltt;
					double Wrrr = WX*WX*WX * (Nrrr + V*(3.*Lrr + V*(3.*Hr + J*V)));
					double Wrtt = WX * (Nttr + V*Ltt);

					// check for consistency
					PRINT_ASSERT( abs( E  - WLorentz*WLorentz*( J + V*(2.*Hr+Lrr*V) ) ) / abs(E) ,<, 1e-6 );
					PRINT_ASSERT( abs( Fr - WLorentz*WLorentz*X*( Hr*(1.+V*V) + (J+Lrr)*V ) ) ,<, 1e-6 );

					new_Prr_E[z_ind][s][ie] = Prr / E;
					new_Ptt_E[z_ind][s][ie] = Ptt / E;
					new_Wrrr[z_ind][s][ie] = Wrrr;
					new_Wttr[z_ind][s][ie] = Wrtt;
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
				int indexChi  = inde + 3*indlast;
				int indexPrr_pm = inde + 2*indlast + 0*ne*ns*nr_GR1D*3;
				//int indexChi_pm = inde + 0*indlast + 0*ne*ns*nr_GR1D*1;

				// need to set plus/minus chi - used for estimating the characteristic speeds
				// There is no way to compute, e.g., dP/dE to be able to compute characteristics
				// using Monte Carlo.
				q_M1_extra[indexChi] = q_M1_extra_2mom[indexChi];

				// set closure to GR1D value on first iteration and where not using Sedonu
				if( iter==0 or q_M1[indexFr] <= 1e-2 /*(not extract_MC[z_ind][s][ie])*/ ){
					q_M1[indexPrr]        = q_M1_2mom[indexPrr];
					q_M1_extra[indexPtt]  = q_M1_extra_2mom[indexPtt];
					q_M1p[indexPrr_pm]    = q_M1p_2mom[indexPrr_pm];
					q_M1m[indexPrr_pm]    = q_M1m_2mom[indexPrr_pm];
					q_M1_extra[indexWrrr] = q_M1_extra_2mom[indexWrrr];
					q_M1_extra[indexWttr] = q_M1_extra_2mom[indexWttr];
				}

				// Set Monte Carlo closure
				if( q_M1[indexFr] > 1e-2 /*extract_MC[z_ind][s][ie]*/){
					// spatial smoothing
					double Prr_E = new_Prr_E[z_ind][s][ie];
					double Ptt_E = new_Ptt_E[z_ind][s][ie];
					double Wrrr  =  new_Wrrr[z_ind][s][ie];
					double Wttr  =  new_Wttr[z_ind][s][ie];

					// temporal smoothing
					Prr_E = fsmooth*Prr_E + (1.-fsmooth)*q_M1[indexPrr];
					Ptt_E = fsmooth*Ptt_E + (1.-fsmooth)*q_M1_extra[indexPtt];
					Wrrr  = fsmooth*Wrrr  + (1.-fsmooth)*q_M1_extra[indexWrrr];
					Wttr  = fsmooth*Wttr  + (1.-fsmooth)*q_M1_extra[indexWttr];

					// write out the variables
					// do flat interpolation to plus/minus
					// lower order than for other moments, but I don't want to
					// re-implement the TVD and PPM reconstruction algorithms
					q_M1[indexPrr]       = Prr_E;
					q_M1_extra[indexPtt] = Ptt_E;
					q_M1p[indexPrr_pm]   = Prr_E;
					q_M1m[indexPrr_pm]   = Prr_E;
					q_M1_extra[indexWrrr]  = Wrrr;
					q_M1_extra[indexWttr]  = Wttr;

					// check that the results are reasonable
					PRINT_ASSERT(q_M1[indexPrr]         ,>=, 0);
					PRINT_ASSERT(q_M1[indexPrr]-1       ,<=, X*X-1);
					PRINT_ASSERT(q_M1_extra[indexPtt]   ,>=, 0);
					PRINT_ASSERT(q_M1_extra[indexPtt]-1 ,<=, 0);
					PRINT_ASSERT(q_M1p[indexPrr_pm]     ,>=, 0);
					PRINT_ASSERT(q_M1p[indexPrr_pm]-1   ,<=, X*X-1);
					PRINT_ASSERT(q_M1m[indexPrr_pm]     ,>=, 0);
					PRINT_ASSERT(q_M1m[indexPrr_pm]-1   ,<=, X*X-1);

				}


			} // energy
		} // species
	} // z_ind
	cout.flush();
}
