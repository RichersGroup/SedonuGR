
-- Included Physics

do_GR         = 0
do_visc       = 0
do_relativity = 0
do_annihilation = 0
do_randomwalk = 0
use_scattering_kernels = 1
radiative_eq  = 1
reflect_outer = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

-- Opacity and Emissivity

neutrino_type = "grey"
nugrid_n = 20
nugrid_start = 0
nugrid_stop = 150
grey_opac = 0
grey_abs_frac = 0
--nulib_table = "NuLib.h5"
--nulib_eos = "SFHo.h5"
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

grid_type = "Grid1DSphere"
model_type = "custom"
model_file = "empty_sphere.mod"

-- Output

write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1
output_zones_distribution = 0
output_hdf5 = 0

-- Particle Creation

max_particles  = 2e7
n_subcycles = 1
do_emit_by_bin = 1
n_emit_core_per_bin    = 1000
n_emit_therm_per_bin   = 0
max_time_hours = -1

-- Inner Source

r_core = 1.5
core_emit_method = 1
T_core = {10,10}
core_chem_pot = {10,-10}
core_lum_multiplier = {1.0,1.0}

-- General Controls

verbose       = 1
max_n_iter =  1
min_step_size = 0.1
max_step_size = 0.4

-- Biasing

min_packet_weight = 0

-- Random Walk

randomwalk_sphere_size = 0