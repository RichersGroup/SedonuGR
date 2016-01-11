-- what are we simulating?

verbose      = 1
radiative_eq = 0
do_visc      = 1
do_annihilation = 0
solve_T      = 0
solve_Ye     = 0
do_relativity = 0
reflect_outer= 0
output_hdf5 = 0
max_n_iter = 1

importance_bias = 1
min_importance = 1e-4
bias_path_length = 1
max_path_length_boost = 2
min_packet_energy = 10
max_packet_energy = 1e500


-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "neutron_star.mod"
nulib_table = "/data/tables/NuLib/MCNuTrans/Helmholtz_full/NuLib.h5"
nulib_eos_filename = ""
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_distribution = 0

-- spectrum parameters

nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- distribution parameters
distribution_nmu = 1
distribution_nphi = 1

-- particle creation parameters

n_emit_core    = 0
n_emit_therm   = 2e6
max_particles  = 9e6
do_emit_by_bin = 0
nut_cdf_cutoff = 0
cdf_interpolation_order = 0
visc_specific_heat_rate = 1e22
emissions_per_timestep = 10

-- particle propagation parameters

step_size = 0.4

-- core parameters

core_emit_method = 1
r_core = 0
T_core = 0
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity    = -1
nut_grey_abs_frac   = -1
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.1
brent_itmax = 100
brent_tolerance = 1e-5