-- what are we simulating?

verbose       = 1
do_photons    = 0
do_neutrinos  = 1
radiative_eq  = 0
iterative     = 1
solve_T       = 0
solve_Ye      = 0
do_visc       = 0
reflect_outer = 0

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "empty_sphere.mod"
nulib_table = "/data/tables/NuLib/NuLib_HShen_noscat.h5"
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1

-- spectrum parameters

nut_spec_time_grid  = {1,1,1}
nut_spec_nu_grid    = {0,5e22,2e20}
nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- particle creation parameters

n_initial      = 0
n_emit_core    = 1e6
n_emit_therm   = 0
max_particles  = 1e6
emissions_per_timestep = 1

-- particle propagation parameters

max_n_steps = 1
dt = -1
step_size = 0.4

-- inner source

r_core = 1
T_core = 10
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity  = -1
nut_grey_abs_frac = -1
nut_nugrid_start  = 1
nut_nugrid_stop   = 2e10
nut_nugrid_n      = 10
nut_cdf_cutoff    = 1e-2

-- equilibrium solver parameters

damping = 0.5
brent_itmax = 100
brent_tolerance = 0.01