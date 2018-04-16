-- Included Physics

do_visc       = 0
do_relativity = 0
do_annihilation = 0
radiative_eq  = 1
reflect_outer = 0
use_scattering_kernels = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 1
equilibrium_damping = 0.5
equilibrium_itmax = 100
equilibrium_tolerance = 1e-5

-- Opacity and Emissivity

neutrino_type = "NuLib"
nulib_table = "NuLib.h5"
nulib_eos = "SFHo.h5"
opac_interp_method = 0
cdf_interpolation_order = 1
cdf_cutoff    = 0

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

distribution_type = "Polar"
distribution_polar_basis = 1
distribution_nmu = 2
distribution_nphi = 2

-- Grid and Model

grid_type = "Grid1DSphere"
model_type = "custom"
model_file = "neutron_star.mod"

-- Output

write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_zones_distribution = 0
output_hdf5 = 0

-- Particle Creation

max_particles  = 2e7
n_subcycles = 1
n_emit_core_per_bin    = 1e3
n_emit_therm_per_bin   = 0

-- Inner Source

r_core = 9e5
core_emit_method = 1
T_core = {10,10,10}
core_chem_pot = {0,0,0}
core_lum_multiplier = {1.0,1.0,1.0}

-- General Controls

verbose       = 1
max_n_iter =  20
min_step_size = 0.01
max_step_size = 0.1
max_time_hours = -1

-- Biasing

bias_path_length = 0
min_packet_weight = 0.1

-- Random Walk

do_randomwalk = 1
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 200
randomwalk_min_optical_depth = 100
randomwalk_interpolation_order = 1
