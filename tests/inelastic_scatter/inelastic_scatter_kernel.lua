verbose = 1
reflect_outer = 0
do_annihilation = 0
neutrino_type = "NuLib"
nulib_eos = "/panfs/ds08/sxs/scurtis/EOS/Hempel_SFHoEOS_rho222_temp180_ye60_version_1.1_20120817.h5"
nulib_table = "/panfs/ds08/sxs/scurtis/Tables/SFHo_tables/simple_with_Ielastic/NuLib_rho82_temp65_ye51_ng15_ns3_Itemp65_Ieta61_version1.0_20181209.h5"

-- output parameters
write_zones_every = 1

-- bias parameters
min_packet_weight = 0.707106781 --0.707106781 -- 1/sqrt(2)

-- distribution parameters
distribution_type = "Moments"

-- input/output files
grid_type = "Grid0DIsotropic"
Grid0DIsotropic_rho = 1e14
Grid0DIsotropic_T   = 14.3
Grid0DIsotropic_Ye  = 0.11

-- spectrum parameters
spec_n_mu = 1
spec_n_phi = 1

-- particle creation parameters
n_emit_core_per_bin = 0
n_emit_therm_per_bin = 10
n_subcycles = 1 
r_core = 0 --7e5
max_n_iter = 1
max_time_hours = -1

-- particle propagation parameters
min_step_size = 0.05
max_step_size = 0.5
absorption_depth_limiter = 1000000.

-- randomwalk
do_randomwalk = 1
randomwalk_max_x = 2
randomwalk_sumN = 1000
randomwalk_npoints = 100
randomwalk_min_optical_depth = 12
randomwalk_interpolation_order = 1
