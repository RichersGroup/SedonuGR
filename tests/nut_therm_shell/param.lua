-- what are we simulating?

verbose      = 1
iterative    = 1
radiative_eq = 1
solve_T      = 0
solve_Ye     = 1
do_visc      = 0
do_relativity = 0
reflect_outer= 1
do_annihilation = 0

-- distribution parameters

distribution_nmu = 4
distribution_nphi = 8

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "shell.mod"
nulib_table = "/data/tables/NuLib/Helmholtz_full/NuLib.h5"
nulib_eos_filename = ""
write_zones_every   = 5
write_rays_every    = -1
write_spectra_every = -1
output_distribution = 0

-- spectrum parameters

nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- particle creation parameters

n_initial      = 0
n_emit_core    = 2e5
n_emit_therm   = 0
max_particles  = 2e5
emissions_per_timestep = 1
cdf_interpolation_order = 0
nut_cdf_cutoff = 0
ratio_emit_by_bin = 0

-- particle propagation parameters

max_n_steps = 20
dt = -1
step_size = 0.4

-- inner source

core_emit_method = 1
r_core = 9.9999e5
T_core = 1e-10
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity    = -2e-16
nut_grey_abs_frac   = -1
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.7
brent_itmax = 100
brent_tolerance = 1e-10
