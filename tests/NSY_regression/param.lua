-- Included Physics

do_visc       = 0
do_relativity = 1
do_annihilation = 0
radiative_eq  = 0
use_scattering_kernels = 0
reflect_outer = 0
do_randomwalk = 1

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

-- Opacity and Emissivity

opacity_dir = "NSY/opacities"
opac_interp_method = 0
cdf_interpolation_order = 3
cdf_cutoff    = 0

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

neutrino_type="Nagakura"
distribution_polar_basis = 1
distribution_type = "RadialMoments"
distribution_moment_order = 3
nugrid_filename               = "NSY/meshdata/Emesh_ascii_edge.dat"
nugrid_n = -1

-- Grid and Model

grid_type = "Grid1DSphere"
model_type="Nagakura"
Grid1DSphere_Nagakura_rgrid_file="NSY/meshdata/rmesh_ascii_edge.dat"
model_file="NSY/data_Sherwood_format/douSherw.00120"

-- Output

write_zones_every   = 1
write_rays_every    = 0
write_spectra_every = 1
output_zones_distribution = 1
output_hdf5 = 1

-- Particle Creation

max_particles  = 2e7
n_subcycles = 1
do_emit_by_bin = 1
n_emit_core_per_bin    = 0
n_emit_therm_per_bin = 1
min_packet_weight = 0.1

-- Inner Source

r_core = 0

-- General Controls

verbose       = 1
max_n_iter =  1
min_step_size = 0.01
max_step_size = 0.4
max_time_hours = -1

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_energy = 1e21
max_packet_energy = 1e99
max_path_length_boost = 0
min_importance = 1e-4
exponential_decay = 0

-- Random Walk

randomwalk_sphere_size = 0.01
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 100
randomwalk_min_optical_depth = 12
randomwalk_absorption_depth_limit = 1.0
randomwalk_interpolation_order = 1
randomwalk_n_isotropic = 20
