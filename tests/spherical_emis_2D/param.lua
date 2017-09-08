-- Included Physics

do_visc       = 0
do_relativity = 0
do_annihilation = 0
radiative_eq  = 1
reflect_outer = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

-- Opacity and Emissivity

neutrino_type = "NuLib"
nulib_table = "NuLib.h5"
nulib_eos = "LS220.h5"
opac_interp_method = 0
cdf_interpolation_order = 1
cdf_cutoff    = 0

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

distribution_type = "Polar"
distribution_polar_basis = 0
distribution_nmu = 2
distribution_nphi = 2

-- Grid and Model

grid_type = "Grid2DSphere"
Grid2DSphere_ignore_theta_min_dist = 0
model_type = "custom"
model_file = "empty_sphere.mod"

-- Output

write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_zones_distribution = 0
output_hdf5 = 0

-- Particle Creation

max_particles  = 2e7
n_subcycles = 1
do_emit_by_bin = 1
n_emit_core_per_bin    = 1000
n_emit_therm_per_bin   = 0

-- Inner Source

r_core = 9.99999e5
core_emit_method = 1
T_core = 10
core_nue_chem_pot = 10
core_lum_multiplier = 1.0

-- General Controls

verbose       = 1
max_n_iter =  1
step_size = 0.4
max_time_hours = -1

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_energy = 1e24
max_packet_energy = 1e99
exponential_decay = 0

-- Random Walk

randomwalk_sphere_size = 0