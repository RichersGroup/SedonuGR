-- what are we simulating?

verbose      = 1
do_photons   = 0
do_neutrinos = 1
iterative    = 1
radiative_eq = 1
solve_T      = SOLVE_T_HERE
solve_Ye     = SOLVE_YE_HERE
do_visc      = 0
reflect_outer= 1

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "shell.mod"
nulib_table = "/data/tables/NuLib/NuLib_HShen_noscat.h5"
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
n_emit_core    = 1e6
n_emit_therm   = 0
max_particles  = 1e6
emissions_per_timestep = 1
cdf_interpolation_order = 1
nut_cdf_cutoff = 1e-2

-- particle propagation parameters

max_n_steps = 5
dt = -1
step_size = 0.4

-- inner source

r_core = 9.9999e5
T_core = TEMP_HERE
core_nue_chem_pot = MUNU_HERE
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity    = -2e-16
nut_grey_abs_frac   = -1
nut_nugrid_start    = 1.0
nut_nugrid_stop     = 6.3e22
nut_nugrid_n        = 48

-- equilibrium solver parameters

damping = 0.1
brent_itmax = 100
brent_tolerance = 1e-5
