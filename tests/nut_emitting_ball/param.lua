-- what are we simulating?
do_photons     = 0
do_neutrinos   = 1
solve_T        = 0
solve_Ye       = 0
radiative_eq   = 1;   -- set to enforce radiative equilibrium. 
iterate        = 2;   -- set to do an iterative (time independent) calc

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "neutron_star.mod"
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

-- spectrum parameters
nut_spec_time_grid   = {1,1,1}
nut_spec_log_nu_grid = {20,22.5,0.1}
nut_n_mu             = 1           -- number of cos(theta) bins in output spectrum 
nut_n_phi            = 1           -- number of phi bins in output spectrum


-- particle creation parameters
init_particles = 0                 -- # particles spawned from the pre-existing 'radiation energy' in each zone
n_emit_core    = 0                 -- # particles to emit from core each timestep
n_emit_heat    = 0                 -- # particles to emit from zones each timestep ("actual" emission, ignored if radiative_eq)
n_emit_visc    = 1e4               -- # particles to emit from zones each timestep (from viscosity, ignored if !radiative_eq)
n_emit_decay   = 0                 -- # particles to emit from zones each timestep (from non-thermal processes)
visc_specific_heat_rate = 1.0

-- particle propagation parameters
step_size = 0.4                    -- move at most step_size*min_grid_length at a time


-- opacity parameters
nut_grey_opacity    =  -1        -- optical grey opacity (cm^2/g)
nut_epsilon         =  -1        -- absorption fraction


-- equilibrium solver parameters
damping = 0.5                    -- changes in values between iterations are decreased by this factor