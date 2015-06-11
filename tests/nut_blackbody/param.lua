-- what are we simulating?

verbose       = 1
radiative_eq  = 0
iterative     = 1
solve_T       = 0
solve_Ye      = 0
do_visc       = 0
do_annihilation = 0
do_relativity = 0
reflect_outer = 0

-- input/output files

grid_type = "grid_0D_isotropic"
model_file = ""
rho = 0
T   = 0
Ye  = 0
nulib_table = "/data/tables/NuLib/Helmholtz_noscat/NuLib.h5"
nulib_eos_filename = "/data/tables/EOS/HShen.h5"
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_distribution = 0

-- spectrum parameters

nut_spec_n_mu       = 1
nut_spec_n_phi      = 1

-- distribution parameters

distribution_nmu = 2
distribution_nphi = 2

-- particle creation parameters

n_initial      = 0
n_emit_core    = 0
n_emit_therm   = 1e6
max_particles  = 1e6
emissions_per_timestep = 1
cdf_interpolation_order = 1
nut_cdf_cutoff    = 0
ratio_emit_by_bin = 0

-- particle propagation parameters

max_n_steps = 1
dt = -1
step_size = 0.4

-- inner source

core_emit_method = 1
r_core = 0
T_core = 0
core_nue_chem_pot = 0
core_lum_multiplier = 0.0

-- opacity parameters

nut_grey_opacity  = -1
nut_grey_abs_frac = -1
opac_interp_method = 0

-- equilibrium solver parameters

damping = 0.5
brent_itmax = 100
brent_tolerance = 0.01