-- Included Physics

do_visc       = 0
do_relativity = 0
do_annihilation = 0
radiative_eq  = 0
reflect_outer = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

-- Opacity and Emissivity

grey_opacity  = -1
grey_abs_frac = -1
nulib_table = "/home/srichers/software/NuLib/NuLib_noscat.h5"
nulib_eos = "/home/srichers/software/sedonu-devel/external/tables/EOS/SFHo.h5"
opac_interp_method = 0
cdf_interpolation_order = 1
cdf_cutoff    = 0

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

distribution_nmu = 2
distribution_nphi = 2

-- Grid and Model

grid_type = "Grid0DIsotropic"
model_file = ""
Grid0DIsotropic_rho = 0
Grid0DIsotropic_T   = 0
Grid0DIsotropic_Ye  = 0
Grid0DIsotropic_vx  = 0
Grid0DIsotropic_vy  = 0
Grid0DIsotropic_vz  = 0

-- Output

write_zones_every   = 0
write_rays_every    = 0
write_spectra_every = 1

-- Particle Creation

max_particles  = 2e6
n_subcycles = 10
do_emit_by_bin = 1
n_emit_core_per_bin    = 0
n_emit_therm_per_bin   = 1000

-- Inner Source

r_core = 0
core_emit_method = 1
T_core = 0
core_nue_chem_pot = 0
core_lum_multiplier = 0.0

-- General Controls

verbose       = 1
max_n_iter = 1
step_size = 0.4

-- Biasing

importance_bias = 0
bias_path_length = 0
min_packet_energy = 1e24
max_packet_energy = 1e99
max_path_length_boost = 2
min_importance = 1e-4
exponential_decay = 0

-- Random Walk

randomwalk_sphere_size = 0 --0.4
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 200
randomwalk_min_optical_depth = 100
randomwalk_interpolation_order = 1