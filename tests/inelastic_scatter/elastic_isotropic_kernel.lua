do_relativity = 1
verbose = 1
equilibrium_T = 0
equilibrium_Ye = 0
do_visc = 0
reflect_outer = 0
do_annihilation = 0
use_scattering_kernels = 1
neutrino_type = "NuLib"
nulib_table = "/panfs/ds08/sxs/scurtis/Tables/BHBlp_tables/simple_table/NuLib_rho82_temp65_ye51_ng15_ns3_version1.0_20181125.h5"

-- output parameters
write_zones_every = 1

-- bias parameters
min_packet_weight = 0.707106781 --0.707106781 -- 1/sqrt(2)

-- distribution parameters
distribution_type = "Moments"

-- input/output files
nulib_eos = "/panfs/ds08/sxs/scurtis/EOS/BHB_lpEOS_rho234_temp180_ye60_version_1.02_20140422.h5"
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
max_particles = 2e8
n_subcycles = 16
core_emit_method = 0
r_core = 0 --7e5
max_n_iter = 1
max_time_hours = -1

-- particle propagation parameters
min_step_size = 0.05
max_step_size = 0.5

-- randomwalk
do_randomwalk = 1
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 100
randomwalk_min_optical_depth = 12
randomwalk_interpolation_order = 1
