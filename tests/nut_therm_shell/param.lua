-- what are we simulating?

verbose      = 0
do_photons   = 0
do_neutrinos = 1
iterative    = 1
radiative_eq = 1
solve_T      = 1
solve_Ye     = 1
do_visc      = 0
reflect_outer= 1
do_distribution = 0

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "shell.mod"
nulib_table = "/home/srichers/tables/NuLib_HShen_noscat_highnures.h5"
write_zones_every   = 5
write_rays_every    = -1
write_spectra_every = -1

-- spectrum parameters

nut_spec_time_grid  = {1,1,1}
nut_spec_nu_grid    = {0,2e22,1e20}
nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- particle creation parameters

n_initial      = 0
n_emit_core    = 1e5
n_emit_therm   = 0
max_particles  = 1e5
emissions_per_timestep = 1
cdf_interpolation_order = 0
nut_cdf_cutoff = 0
ratio_emit_by_bin = 1

-- particle propagation parameters

max_n_steps = 10
dt = -1
step_size = 0.4

-- inner source

r_core = 9.9999e5
T_core = 1e-10
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity    = -2e-16
nut_grey_abs_frac   = -1
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.1
brent_itmax = 100
brent_tolerance = 1e-10
