-- Included Physics

do_visc       = 0
do_relativity = 0
do_randomwalk = 0
do_annihilation = 0
radiative_eq  = 1
reflect_outer = 1
use_scattering_kernels = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 1
equilibrium_damping = 0.0
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
write_rays_every    = 0
write_spectra_every = 0
output_zones_distribution = 0
output_hdf5 = 0

-- Particle Creation

max_particles  = 0
n_subcycles = 1
n_emit_therm_per_bin  = 0
n_emit_core_per_bin = 0

-- Inner Source

r_core = 0
core_emit_method = 1
T_core = 0
core_nue_chem_pot = 0
core_lum_multiplier = 0.0

-- General Controls

verbose       = 1
max_n_iter = 1
min_step_size = 0.01
max_step_size = 0.1
max_time_hours = -1

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_weight = 0.01

-- Random Walk

randomwalk_sphere_size = 0
