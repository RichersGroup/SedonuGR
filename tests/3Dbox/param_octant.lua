do_relativity = 1
verbose = 1
equilibrium_T = 0
equilibrium_Ye = 0
do_visc = 0
reflect_outer = 0
do_annihilation = 0
use_scattering_kernels = 0

-- opacity stuff
neutrino_type = "grey"
nugrid_n = 20
grey_chempot = 10
nugrid_start = 0
nugrid_stop = 200
nugrid_n = 50
grey_opac = 1
grey_abs_frac = 1

-- output parameters
write_zones_every = 1

-- bias parameters
min_packet_weight = 0.707106781 --0.707106781 -- 1/sqrt(2)

-- distribution parameters
distribution_type = "Moments"

-- input/output files
grid_type = "Grid3DCart"
model_type = "THC"
Grid3DCart_THC_reflevel = 0
Grid3DCart_reflect_x=1
Grid3DCart_reflect_y=1
Grid3DCart_reflect_z=1
Grid3DCart_rotate_quadrant = 0
Grid3DCart_rotate_hemisphere_x = 0
Grid3DCart_rotate_hemisphere_y = 0
model_file = "stationary.h5"

-- spectrum parameters
spec_n_mu = 1
spec_n_phi = 1

-- particle creation parameters
n_emit_core_per_bin = 0
n_emit_therm_per_bin = 100
max_particles = 2e8
n_subcycles = 1
core_emit_method = 0
r_core = 0 --7e5
max_n_iter = 1
max_time_hours = -1

-- particle propagation parameters
min_step_size = 0.05
max_step_size = 0.5
absorption_depth_limiter = 1000000.0

-- randomwalk
do_randomwalk = 1
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 100
randomwalk_min_optical_depth = 6
randomwalk_interpolation_order = 1
