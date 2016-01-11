-- what are we simulating?

verbose       = 1
radiative_eq  = 0
solve_T       = 0
solve_Ye      = 0
do_visc       = 0
do_annihilation = 0
do_relativity = 0
reflect_outer = 0
output_hdf5 = 0
max_n_iter = 1

-- biasing
importance_bias = 0
bias_path_length = 1
max_path_length_boost = 10
min_packet_energy = 0
max_packet_energy = 1e100

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "empty_sphere.mod"
nulib_table = "/data/tables/NuLib/MCNuTrans/LS220_full/NuLib.h5"
nulib_eos_filename = "/data/tables/EOS/HShen.h5"
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_distribution = 0

-- spectrum parameters

nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- distribution parameters

distribution_nmu = 2
distribution_nphi = 2

-- particle creation parameters

n_emit_core_per_bin    = 1000
n_emit_therm_per_bin   = 0
max_particles  = 2e7
emissions_per_timestep = 1
cdf_interpolation_order = 0
nut_cdf_cutoff    = 0
do_emit_by_bin = 1

-- particle propagation parameters

step_size = 0.4

-- inner source

core_emit_method = 1
r_core = 9.99999e5
T_core = 10
core_nue_chem_pot = 20
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity  = -1
nut_grey_abs_frac = -1
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.5
brent_itmax = 100
brent_tolerance = 0.01