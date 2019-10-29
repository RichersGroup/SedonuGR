
-- Included Physics

do_relativity = 0
do_annihilation = 0
do_randomwalk = 0
reflect_outer = 0

-- Equilibrium Solving

equilibrium_T  = 0
equilibrium_Ye = 0

-- Opacity and Emissivity

neutrino_type = "grey"
nugrid_n = 50
nugrid_start = 0
nugrid_stop = 100
grey_opac = 0
grey_abs_frac = 0
grey_chempot = 0

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

distribution_type = "Polar"
distribution_nmu = 2
distribution_nphi = 2

-- Grid and Model

grid_type = "Grid1DSphere"
model_type = "custom"
model_file = "empty_sphere.mod"

-- Output

write_zones_every   = 1

-- Particle Creation

n_subcycles = 1
n_emit_core_per_bin    = 100
n_emit_therm_per_bin   = 0
max_time_hours = -1

-- Inner Source

r_core = 1.5
T_core = {10}
core_chem_pot = {10}
core_lum_multiplier = {1.0}

-- General Controls

verbose       = 1
max_n_iter =  1
min_step_size = 0.1
max_step_size = 0.1
absorption_depth_limiter = 1.0

-- Biasing

min_packet_weight = 0

-- Random Walk

