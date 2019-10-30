-- Included Physics

do_annihilation = 0
reflect_outer = 0
do_randomwalk = 1

-- Opacity and Emissivity

opacity_dir = "NSY/opacities"

-- Escape Spectra

spec_n_mu       = 1
spec_n_phi      = 1

-- Distribution Function

neutrino_type="Nagakura"
distribution_type = "RadialMoments"
nugrid_filename               = "NSY/meshdata/Emesh_ascii_edge.dat"
nugrid_n = -1

-- Grid and Model

grid_type = "Grid1DSphere"
model_type="Nagakura"
Grid1DSphere_Nagakura_rgrid_file="NSY/meshdata/rmesh_ascii_edge.dat"
model_file="NSY/data_Sherwood_format/douSherw.00120"

-- Output

write_zones_every   = 1

-- Particle Creation

n_subcycles = 1
n_emit_core_per_bin    = 0
n_emit_therm_per_bin = 1
min_packet_weight = 0.5

-- Inner Source

r_core = 0

-- General Controls

verbose       = 1
max_n_iter =  1
min_step_size = 0.01
max_step_size = 0.4
max_time_hours = -1
absorption_depth_limiter = 1.0

-- Random Walk

randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 100
randomwalk_min_optical_depth = 12
randomwalk_interpolation_order = 1
