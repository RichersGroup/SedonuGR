
-- Included Physics

do_annihilation = 0
do_randomwalk = 0
reflect_outer = 0

-- Opacity and Emissivity

neutrino_type = "grey"
Neutrino_grey_opac  = 4
Neutrino_grey_abs_frac = 1
Neutrino_grey_chempot = 0.
nugrid_start = 10
nugrid_stop = 10.001
nugrid_n = 1

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

distribution_type = "Moments"

-- Grid and Model

grid_type = "Grid1DSphere"
model_type = "custom"
model_file = "uniform_sphere.mod"

-- Output

write_zones_every   = 1

-- Particle Creation

n_subcycles = 10
n_emit_core_per_bin    = 0 --100
n_emit_therm_per_bin   = 10
max_time_hours = -1

-- Inner Source

r_core = 0 --1.5e5
T_core = {10}
core_chem_pot = {0}
core_lum_multiplier = {1.0}

-- General Controls

verbose       = 1
max_n_iter =  1
min_step_size = .4 --0.01
max_step_size = 0.4

-- Biasing

min_packet_weight = 0.01

-- Random Walk

randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 200
randomwalk_min_optical_depth = 5
