-- what are we simulating?
do_photons   = 0
do_neutrinos = 1

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "neutron_star.mod"
output_file =  "my_lightcurve.dat"          -- output light curve file
nulib_table = "../../external/tables/NuLib.h5"

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

init_particles = 0
step_size      = 0.4

radiative_eq   = 1;   -- set 
iterate        = 10;   -- set to do an iterative (time independent) calc
damping        = 0.5;

solve_T  = 1
solve_Ye = 1

-- opacity parameters
nut_grey_opacity = -1  -- grey opacity - set to negative to turn off
nut_epsilon      = -1  -- absorption fraction - set to negative to turn off

-- inner source
r_core = 11.5e5
L_core = 1e46
T_core = 3.5e11  -- 30 MeV
n_inject = 6e3

-- output spectrum
nut_spec_time_grid   = {1,1,1}
nut_spec_log_nu_grid =   {19,22,0.1}
nut_n_mu             =    1           -- number of cos(theta) bins in output spectrum 
nut_n_phi            =    1           -- number of phi bins in output spectrum


