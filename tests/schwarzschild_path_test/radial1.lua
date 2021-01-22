-- Included Physics

do_randomwalk = 0
do_annihilation = 0
reflect_outer = 0

-- Opacity and Emissivity

neutrino_type = "grey"
Neutrino_grey_abs_frac = 0
Neutrino_grey_opac = 0
Neutrino_grey_chempot = 0
nugrid_n=1
nugrid_start=0
nugrid_stop=1e99

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
n_emit_core_per_bin    = 0
n_emit_therm_per_bin   = 0
max_time_hours = -1

-- Inner Source

r_core = 0
T_core = 10
core_lum_multiplier = 1.0

-- General Controls

verbose       = 0
max_n_iter =  1
min_step_size = 0.01
max_step_size = 0.1

-- Biasing

min_packet_weight = 0

-- Random Walk

Schwarzschild_initial_xup = {1.5,0,0,0}
Schwarzschild_initial_kup = {1,0,0,1}
