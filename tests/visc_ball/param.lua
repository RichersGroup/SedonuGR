-- Included Physics

do_visc       = 1
do_relativity = 0
do_annihilation = 0
radiative_eq  = 0
reflect_outer = 0
visc_specific_heat_rate = 1e22

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

model_type = "custom"
grid_type = "Grid1DSphere"
model_file = "neutron_star.mod"

-- Output

write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_zones_distribution = 0
output_hdf5 = 0

-- Particle Creation

max_particles  = 9e6
n_subcycles = 1
do_emit_by_bin = 0
n_emit_core    = 0
n_emit_therm   = 2e6

-- Inner Source

r_core = 0
core_emit_method = 1
T_core = 0
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- General Controls

verbose       = 1
max_n_iter =  1
step_size = 0.4
max_time_hours = -1

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_energy = 1e20
max_packet_energy = 1e99
max_path_length_boost = 2
min_importance = 1e-4
exponential_decay = 0

-- Random Walk

randomwalk_sphere_size = 0