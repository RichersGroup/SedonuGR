-- Included Physics

do_GR = 0
do_visc       = 0
do_randomwalk = 1
do_relativity = 1
do_annihilation = 0
use_scattering_kernels = 0
radiative_eq  = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

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
distribution_nmu = 2
distribution_nphi = 2
distribution_polar_basis = 0

-- Grid and Model

grid_type = "Grid0DIsotropic"
Grid0DIsotropic_rho = 0
Grid0DIsotropic_T   = 0
Grid0DIsotropic_Ye  = 0
Grid0DIsotropic_vx  = 0
Grid0DIsotropic_vy  = 0
Grid0DIsotropic_vz  = 0

-- Output

write_zones_every   = 1
write_rays_every    = 0
write_spectra_every = 1
output_zones_distribution = 0
output_hdf5 = 0

-- Particle Creation

max_particles  = 2e6
n_subcycles = 1
do_emit_by_bin = 0
n_emit_core_per_bin    = 0
n_emit_therm_per_bin   = 500

-- Inner Source

r_core = 0

-- General Controls

verbose       = 1
min_step_size = 0.01
max_step_size = 0.4

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_weight = 0.01

-- Random Walk

randomwalk_sphere_size = 0.1
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 200
randomwalk_min_optical_depth = 100
randomwalk_interpolation_order = 1
randomwalk_n_isotropic=0