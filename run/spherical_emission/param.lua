-- what are we simulating?
do_photons = 1

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "lucy.mod"
output_file =  "my_lightcurve.dat"          -- output light curve file

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
n_mu        =    1           -- number of cos(theta) bins in output spectrum 
n_phi       =    1           -- number of phi bins in output spectrum
init_particles = 0

step_size          = 0.4

radiative_eq   = 1;   -- set 
iterate        = 2;   -- set to do an iterative (time independent) calc

-- inner source
r_core = 2e15
L_core = 1e43
T_core = 10000
n_inject = 1e3

-- opacities
gray_abs_opacity    =  0.1           -- optical grey opacity in cm^2/g
gray_scat_opacity = -1
epsilon         =  1.0

-- opacity grid
nu_start        =  2e13
nu_stop         =  1e15
n_nu            =  500


-- output spectrum
spec_time_grid = {1,1,1}
spec_nu_grid =   {1e14,1e15,1e13}
use_transport = 0