-- what are we simulating?
do_photons   = 0
do_neutrinos = 1

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "neutrino_test.mod"
output_file =  "my_lightcurve.dat"          -- output light curve file
nulib_table = "../../external/NuLib/NuLib_LS220.h5"

-- time stepping (all times in days)
day = 3600.0*24

tstep_start =   1.0
tstep_min   =   0.01*day          -- smallest allowed time step
tstep_max   =   0.1*day           -- largest allowed time step
tstep_del   =   0.1               -- largest allowed step fraction of time
n_times     =  1000               -- maximum nuber of steps taken
t_stop      =  5*day             -- time to stop calculation, in days
t_delta     =   1.0               -- time spacing in ouput light curve
write_out   =  1.0*day

n_photons   =  1e5           -- total number of photon packets to use
init_particles = 0
step_size          = 0.4

radiative_eq   = 1;   -- set 
iterate        = 2;   -- set to do an iterative (time independent) calc

-- opacity parameters
nut_grey_opacity = -1  -- grey opacity - set to negative to turn off
nut_epsilon      = -1  -- absorption fraction - set to negative to turn off

-- inner source
r_core = 12e3
L_core = 1e43
T_core = 3.5e11  -- 30 MeV
n_inject = 1e3
rho_core = 3e14
Ye_core = 0.35

-- output spectrum
nut_spec_time_grid   = {1,1,1}
nut_spec_log_nu_grid =   {3e26,4e28,1e27}
nut_n_mu             =    1           -- number of cos(theta) bins in output spectrum 
nut_n_phi            =    1           -- number of phi bins in output spectrum


