-- Included Physics

do_visc       = 0
do_relativity = 0
do_GR = 1
do_annihilation = 0
radiative_eq  = 1
reflect_outer = 0
use_scattering_kernels = 1

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

-- Opacity and Emissivity

neutrino_type = "grey"
grey_abs_frac = 0
grey_opac = 0
opac_interp_method=0
cdf_cutoff = 0
cdf_interpolation_order=0
nugrid_n=1
nugrid_start=0
nugrid_stop=1e99

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

distribution_type = "Polar"
distribution_nmu = 2
distribution_nphi = 2
distribution_polar_basis = 0

-- Grid and Model

grid_type = "Grid1DSphere"
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
n_emit_core_per_bin    = 0
n_emit_therm_per_bin   = 0
max_time_hours = -1

-- Inner Source

r_core = 0
core_emit_method = 1
T_core = 10
core_nue_chem_pot = 10
core_lum_multiplier = 1.0

-- General Controls

verbose       = 0
max_n_iter =  1
min_step_size = 0.01
max_step_size = 0.01

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_number = 0
max_packet_number = 1e99
exponential_decay = 0

-- Random Walk

randomwalk_sphere_size = 0

Grid1DSchwarzschild_r_sch = 1.0
initial_xup = {1.5,0,0,0}
initial_kup = {0,1,0,1}
Grid1DSphere_radial_interpolation_method = 1