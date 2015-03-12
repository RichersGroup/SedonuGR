-- what are we simulating?

verbose      = 1
iterative    = 1
radiative_eq = 0
do_visc      = 1
do_annihilation = 0
solve_T      = 1
solve_Ye     = 0
do_relativity = 0
reflect_outer= 0

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "neutron_star.mod"
nulib_table = "/data/tables/NuLib/Helmholtz_full/NuLib.h5"
nulib_eos_filename = ""
write_zones_every   = 10
write_rays_every    = 10
write_spectra_every = 10

-- spectrum parameters

nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- distribution parameters
distribution_nmu = 1
distribution_nphi = 1

-- particle creation parameters

n_initial      = 0
n_emit_core    = 0
n_emit_therm   = 1e6
n_emit_decay   = 0
max_particles  = 1e6
ratio_emit_by_bin = 0
nut_cdf_cutoff = 1e-2
cdf_interpolation_order = 1
visc_specific_heat_rate = 1e22
emissions_per_timestep = 1
ratio_emit_by_zone = 0

-- particle propagation parameters

max_n_steps = 10
dt = -1
step_size = 0.4

-- core parameters

core_emit_method = 1
r_core = 0
T_core = 0
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity    =  -1
nut_grey_abs_frac   =  -1
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.0
brent_itmax = 100
brent_tolerance = 0.01