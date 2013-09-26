-- what are we simulating?

do_photons     = 1
do_neutrinos   = 0
solve_T        = 0
solve_Ye       = 0
radiative_eq   = 1;   -- set to enforce radiative equilibrium. 
iterate        = 1;   -- set to do an iterative (time independent) calc

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "lucy.mod"

-- spectrum parameters
spec_n_mu      =    1              -- number of cos(theta) bins in output spectrum 
spec_n_phi     =    1              -- number of phi bins in output spectrum
spec_time_grid = {1,1,1}           -- {start, stop, bin width}
spec_nu_grid   = {0,2e15,2e13}     -- {start, stop, bin width}


-- particle creation parameters
init_particles = 0                 -- # particles spawned from the pre-existing 'radiation energy' in each zone
n_emit_core    = 0                 -- # particles to emit from core each timestep
n_emit_heat    = 0                 -- # particles to emit from zones each timestep ("actual" emission, ignored if radiative_eq)
n_emit_visc    = 1e4               -- # particles to emit from zones each timestep (from viscosity, ignored if !radiative_eq)
n_emit_decay   = 0                 -- # particles to emit from zones each timestep (from non-thermal processes)
max_particles  = 1e6
visc_specific_heat_rate = 1.0      -- specific heating rate due to viscosity (erg/s/g)

-- particle propagation parameters
step_size = 0.4                    -- move at most step_size*min_grid_length at a time


-- inner source
r_core = 0                       -- core radius (cm)
L_core = 0                       -- core luminosity (erg/s)
T_core = 0                       -- core temperature (K)

-- opacity parameters
grey_opacity    =  0.1           -- optical grey opacity (cm^2/g)
epsilon         =  1.0           -- absorption fraction
nu_start        =  0             -- photon opacity grid start
nu_stop         =  2e15          -- photon opacity grid stop
n_nu            =  500           -- photon opacity grid number of points along frequency


-- equilibrium solver parameters
damping = 0.5                    -- changes in values between iterations are decreased by this factor
brent_itmax = 100
brent_tolerance = 0.01