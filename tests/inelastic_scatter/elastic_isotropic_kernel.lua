verbose = 1
reflect_outer = 0
do_annihilation = 0
neutrino_type = "grey"
nugrid_n = 10
nugrid_start = 0
nugrid_stop = 10
grey_abs_frac = 0
grey_opac = 1
grey_chempot = 0

-- output parameters
write_zones_every = 1

-- bias parameters
min_packet_weight = 0.707106781 --0.707106781 -- 1/sqrt(2)

-- distribution parameters
distribution_type = "Moments"

-- input/output files
grid_type = "Grid0DIsotropic"
Grid0DIsotropic_rho = 0
Grid0DIsotropic_T   = 0
Grid0DIsotropic_Ye  = 0

-- spectrum parameters
spec_n_mu = 1
spec_n_phi = 1

-- particle creation parameters
n_emit_core_per_bin = 0
n_emit_therm_per_bin = 1
n_subcycles = 16
r_core = 0 --7e5
max_n_iter = 1
max_time_hours = -1

-- particle propagation parameters
min_step_size = 0.05
max_step_size = 0.5
absorption_depth_limiter = 1.

-- randomwalk
do_randomwalk = 1
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 100
randomwalk_min_optical_depth = 12
randomwalk_interpolation_order = 1
