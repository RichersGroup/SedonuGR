-- what are we simulating?

verbose      = 1
do_photons   = 0
do_neutrinos = 1
iterative    = 1
radiative_eq = 0
solve_T      = 0
solve_Ye     = 0
do_visc=0
reflect_outer= 0
do_distribution = 0

-- input/output files

grid_type = "grid_0D_isotropic"
model_file = ""
rho = RHO_HERE
T   = TEMP_HERE
Ye  = YE_HERE
nulib_table = NULIB_HERE
write_zones_every   = 1
write_rays_every    = -1
write_spectra_every = -1

-- spectrum parameters

nut_spec_time_grid  = {1,1,1}
nut_spec_nu_grid    = {1,1,1}
nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- particle creation parameters

n_initial      = 0
n_emit_core    = 0
n_emit_therm   = NPARTICLES_HERE
n_emit_decay   = 0
max_particles  = NPARTICLES_HERE
ratio_emit_by_zone = 0
emissions_per_timestep = 1
cdf_interpolation_order = 1

-- particle propagation parameters

max_n_steps = 1
dt = -1
step_size = 1.0

-- inner source

r_core = 0
core_nue_chem_pot = 0
T_core = 0

-- opacity parameters

nut_grey_opacity  =  -1
nut_grey_abs_frac =  -1
nut_cdf_cutoff    = 1e-2
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.1
brent_itmax = 100
brent_tolerance = 0.01