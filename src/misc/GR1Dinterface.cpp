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

double GR1D_smoothing_timescale;
int GR1D_recalc_every;
const double time_gf = 2.03001708e5;
ScalarMultiDArray<double, 3> Ptt_E_tet, Wrtt_Fr_tet;
double tolerance = 1e-6;

struct Moments{
	double E, Fr, Prr, Ptt, Wrrr, Wrtt;
};

Moments toLab(const Moments tet, double vr, double X){
	Moments lab;
	const double X2 = X*X;
	const double X3 = X2*X;
	const double V = vr*X;
	const double V2 = V*V;
	const double WLorentz = 1./sqrt(1.-V2);
	const double W2 = WLorentz*WLorentz;
	const double W3 = W2*WLorentz;

	lab.E = W2*( tet.E + V*(2.*tet.Fr+tet.Prr*V) );
	lab.Fr = W2*X*( tet.Fr*(1.+V*V) + (tet.E+tet.Prr)*V );
	lab.Prr = W2*X2 * (tet.Prr + V*(2.*tet.Fr + tet.E*V));
	lab.Ptt = tet.Ptt;
	lab.Wrrr = W3/X3 * (tet.Wrrr + V*(3.*tet.Prr + V*(3.*tet.Fr + tet.E*V)));
	lab.Wrtt = WLorentz/X * (tet.Wrtt + V*tet.Ptt);

	PRINT_ASSERT(lab.E,>=,0);
	PRINT_ASSERT(lab.Prr,>=,0);
	PRINT_ASSERT(lab.Ptt,>=,0);

	return lab;
}
Moments toTet(const Moments lab, double vr, double X){
	Moments tet;

	// account for floor in GR1D
	const double X2 = X*X;
	const double X3 = X*X*X;
	const double V = vr*X;
	const double V2 = V*V;
	const double V4 = V2*V2;
	const double WLorentz = 1./sqrt(1.-V2);
	const double W2 = WLorentz*WLorentz;
	const double W3 = W2*WLorentz;

	tet.E = W2 * (lab.E - 2.*vr*lab.Fr + vr*vr*lab.Prr );
	tet.Fr = W2/X * ( lab.Fr*(1.+V2) - vr*(lab.Prr+lab.E*X2) );
	tet.Prr = W2/X2 * (lab.Prr - 2.*lab.Fr*vr*X2 + lab.E*V2*X2 );
	tet.Ptt = lab.Ptt;
	tet.Wrtt = lab.Wrtt*X/WLorentz - lab.Ptt*V;
	tet.Wrrr = lab.Wrrr*X3/W3 - 1./X2*( lab.Prr*V*(3.+W2*V4) - V2*X*lab.Fr*(1.+2.*W2) + W2*V2*X2*V*lab.E);

	// We know the GR1D closure isn't great, so we will fix it.
	if(tet.Fr!=0){
	  double Wmag = tet.Wrrr+2.*tet.Wrtt;
	  tet.Wrrr *= tet.Fr/Wmag;
	  tet.Wrtt *= tet.Fr/Wmag;
	}
	else{
	  tet.Wrrr = 0;
	  tet.Wrtt = 0;
	}

	PRINT_ASSERT(tet.E,>,0);
	PRINT_ASSERT(tet.Prr,>=,0);
	PRINT_ASSERT(tet.Ptt,>=,0);
	PRINT_ASSERT(abs(tet.Fr),<=,tet.E);
	PRINT_ASSERT(abs(tet.E  - (tet.Prr +2.*tet.Ptt ))/tet.E ,<, tolerance);
	PRINT_ASSERT(abs(tet.Fr - (tet.Wrrr+2.*tet.Wrtt))/tet.E ,<, tolerance);

	return tet;
}
void applyClosure(Moments& tet, Moments& lab, double Ptt_E_tet, double Wrtt_Fr_tet, double vr, double X){

	if(tet.E==0) return;

	// apply the closure in the tetrad frame
	tet.Ptt  = tet.E  *            Ptt_E_tet ;
	tet.Prr  = tet.E  * (1. - 2.*  Ptt_E_tet);
	tet.Wrtt = tet.Fr *          Wrtt_Fr_tet ;
	tet.Wrrr = tet.Fr * (1. - 2.*Wrtt_Fr_tet);

	// get the parameters that make the lab and tetrad moments match
	const double V = vr*X;
	double denom = ( tet.E - tet.Prr*V*V)*X;
	double alpha = ( lab.E*X*(1.+V*V) - 2.*lab.Fr*V ) / denom;
	double beta  = ( lab.Fr*(tet.E+tet.Prr*V*V) - lab.E*V*X*(tet.E+tet.Prr) ) / (tet.Fr*denom);
	PRINT_ASSERT(denom,>,0);
	PRINT_ASSERT(alpha,>,0);

	// modify the tetrad such that it follows the closure relation AND matches lab-frame
	tet.E    *= alpha;
	tet.Prr  *= alpha;
	tet.Ptt  *= alpha;
	tet.Fr   *= beta;
	tet.Wrrr *= beta;
	tet.Wrtt *= beta;

	// get lab moments. First two moments should match by construction of the closure.
	double Eold = lab.E;
	double Frold = lab.Fr;
	lab = toLab(tet, vr, X);
	PRINT_ASSERT(abs(lab.E  - Eold )/Eold ,<, tolerance);
	PRINT_ASSERT(abs(lab.Fr - Frold)/Eold ,<, tolerance);
}

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
	GR1D_smoothing_timescale = lua.scalar<double>("GR1D_smoothing_timescale");
	GR1D_recalc_every = lua.scalar<int>("GR1D_recalc_every");

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

	// resize arrays to be used for time-averaging
	const size_t ns = (*sim)->species_list.size();
	const size_t nr = (*sim)->grid->rho.size();
	const size_t ne = (*sim)->grid->nu_grid_axis.size();
	vector<Axis> axes(3);
	axes[0] = Axis(0, ns, ns);
	axes[1] = Axis(0, nr, nr);
	axes[2] = Axis(0, ne, ne);
	Ptt_E_tet.set_axes(axes);
	Wrtt_Fr_tet.set_axes(axes);

	// set initial value to 1/3
	for(size_t i=0; i<ns*nr*ne; i++){
	  Ptt_E_tet[i] = 1./3.;
	  Wrtt_Fr_tet[i] = 1./3.;
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
	const size_t nr = (*sim)->grid->rho.size();
	const size_t ns = (*sim)->species_list.size();
	const size_t ne = (*sim)->grid->nu_grid_axis.size();
	const double fsmooth = *dt*GR1D_recalc_every/time_gf / GR1D_smoothing_timescale;
	PRINT_ASSERT(fsmooth,<=,1.0);
	PRINT_ASSERT(fsmooth,>=,0.0);

	// set velocities and opacities
	static_cast<GridGR1D*>((*sim)->grid)->set_fluid(rho, T, Ye, v1, metricX, alp);

	bool extract_MC[nr][ns][ne];
	for(size_t s=0; s<(*sim)->species_list.size(); s++)
		static_cast<Neutrino_GR1D*>((*sim)->species_list[s])->set_eas_external(eas,&(extract_MC[0][0][0]),*rshock);

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

	// set the GR1D quantities
	for(size_t s=0; s<ns; s++){
		GR1DSpectrumArray* tmpSpectrum = static_cast<GR1DSpectrumArray*>((*sim)->grid->distribution[s]);
        #pragma omp parallel for
		for(size_t z_ind=0; z_ind<nr; z_ind++){
			size_t dir_ind[2];
			dir_ind[0] = z_ind;
			int indr = z_ind+nghost_GR1D;
			const double X = metricX[indr];
			const double lapse = alp[indr];
			const double vr = v1[indr]/pc::c;

			for(size_t ie=0; ie<ne; ie++){
				dir_ind[1] = ie;
				int ind_sre = indr + s*nr_GR1D + ie*ns*nr_GR1D;
				int indexE      = ind_sre + 0*ne*ns*nr_GR1D;
				int indexFr     = ind_sre + 1*ne*ns*nr_GR1D;
				int indexPrr    = ind_sre + 2*ne*ns*nr_GR1D;
				int indexPtt    = ind_sre + 0*ne*ns*nr_GR1D;
				int indexWrrr   = ind_sre + 1*ne*ns*nr_GR1D;
				int indexWrtt   = ind_sre + 2*ne*ns*nr_GR1D;
				int indexChi    = ind_sre + 3*ne*ns*nr_GR1D;
				int indexPrr_pm = ind_sre + 2*ne*ns*nr_GR1D + 0*ne*ns*nr_GR1D*3;

				// load the old lab-frame moments
				Moments lab;
				lab.E = q_M1[indexE];
				lab.Fr = q_M1[indexFr];
				lab.Prr = q_M1_2mom[indexPrr] * lab.E;
				lab.Ptt  = q_M1_extra_2mom[indexPtt] * lab.E;
				lab.Wrrr = q_M1_extra_2mom[indexWrrr];
				lab.Wrtt = q_M1_extra_2mom[indexWrtt];

				const double crit_fluxfac = 1e-2;

				// load up tetrad-frame moments (just make zero if not much stuff around)
				Moments tet;
				if(lab.E < 1e-98){
					tet.E = tet.Fr = tet.Prr = tet.Ptt = tet.Wrrr = tet.Wrtt = 0;
				}
				else{
					if(q_M1[indexFr] > crit_fluxfac){
						size_t index = tmpSpectrum->data.direct_index(dir_ind);
						Tuple<double,6> tmp = tmpSpectrum->data[index];
						tet.E    = tmp[0];
						tet.Fr   = tmp[1];
						tet.Prr  = tmp[2];
						tet.Ptt  = tmp[3];
						tet.Wrrr = tmp[4];
						tet.Wrtt = tmp[5];
					}
					else tet = toTet(lab, vr, X);
					PRINT_ASSERT(tet.E,>,0);
					PRINT_ASSERT(tet.Prr,>=,0);
					PRINT_ASSERT(tet.Ptt,>=,0);
					PRINT_ASSERT(abs(tet.Fr),<=,tet.E);
					PRINT_ASSERT(abs(tet.E  - (tet.Prr +2.*tet.Ptt ))/tet.E ,<, tolerance);
					PRINT_ASSERT(abs(tet.Fr - (tet.Wrrr+2.*tet.Wrtt))/tet.E ,<, tolerance);
				}

				// get tetrad-frame closure relations
				double Ptt_E_tetnew   = (tet.E ==0 ? 1./3. : tet.Ptt  / tet.E);
				double Wrtt_Fr_tetnew = (tet.Fr==0 ? 1./3. : tet.Wrtt / tet.Fr);
				PRINT_ASSERT(Ptt_E_tetnew ,<=, 0.5);
				PRINT_ASSERT(Ptt_E_tetnew ,>=, 0.0);
				PRINT_ASSERT(Wrtt_Fr_tetnew ,<=, 0.5);
				PRINT_ASSERT(Wrtt_Fr_tetnew ,>=, 0.0);

				// temporal smoothing of closure relations
				size_t indices[3] = {s, z_ind, ie};
				size_t index = Ptt_E_tet.direct_index(indices);
				Ptt_E_tet[index]   = fsmooth*Ptt_E_tetnew   + (1.-fsmooth)*Ptt_E_tet[index];
				Wrtt_Fr_tet[index] = fsmooth*Wrtt_Fr_tetnew + (1.-fsmooth)*Wrtt_Fr_tet[index];
				PRINT_ASSERT(Ptt_E_tet[index] ,<=, 0.5);
				PRINT_ASSERT(Ptt_E_tet[index] ,>=, 0.0);
				PRINT_ASSERT(Wrtt_Fr_tet[index] ,<=, 0.5);
				PRINT_ASSERT(Wrtt_Fr_tet[index] ,>=, 0.0);

				// get the closed lab state
				applyClosure(tet, lab, Ptt_E_tet[index], Wrtt_Fr_tet[index], vr, X);

				// need to set plus/minus chi - used for estimating the characteristic speeds
				// There is no way to compute, e.g., dP/dE to be able to compute characteristics
				// using Monte Carlo.
				q_M1_extra[indexChi] = q_M1_extra_2mom[indexChi];

				// write out the variables
				// do flat interpolation to plus/minus
				// lower order than for other moments, but I don't want to
				// re-implement the TVD and PPM reconstruction algorithms
				q_M1[indexPrr]        = lab.Prr / lab.E;
				q_M1_extra[indexPtt]  = lab.Ptt / lab.E;
				q_M1_extra[indexWrrr] = lab.Wrrr;
				q_M1_extra[indexWrtt] = lab.Wrtt;
				q_M1p[indexPrr_pm]    = lab.Prr / lab.E;
				q_M1m[indexPrr_pm]    = lab.Prr / lab.E;

			} // energy
		} // species
	} // z_ind



	cout.flush();
}

extern "C"
void write_sedonu_closure_(const hid_t* file_id){
  H5::H5File file(*file_id);
  Ptt_E_tet.write_HDF5(file, "SedonuGR_Ptt_E_tet");
  Wrtt_Fr_tet.write_HDF5(file, "SedonuGR_Wrtt_Fr_tet");
}

extern "C"
void read_sedonu_closure_(const hid_t* file_id){
  H5::H5File file(*file_id);
  Ptt_E_tet.read_HDF5(file, "SedonuGR_Ptt_E_tet", Ptt_E_tet.axes);
  Wrtt_Fr_tet.read_HDF5(file, "SedonuGR_Wrtt_Fr_tet", Wrtt_Fr_tet.axes);
}
