-- what are we simulating?

do_photons     = 0                 -- simulate photons?
do_neutrinos   = 1                 -- simulate neutrinos?
radiative_eq   = 1                 -- set to enforce radiative equilibrium. 
iterate        = 1                -- set to do an iterative (time independent) calc
solve_T        = 0                 -- (if radiative_eq) solves each zone's temperature based on its absorbed energy
solve_Ye       = 0                 -- (if radiative_eq) solves each zone's Ye based on its absorbed lepton number

-- input/output files

grid_type = "grid_1D_sphere"       -- grid geometry. Must match grid geometry in model file if used  
model_file  =  "neutron_star.mod"  -- model file. "custom" --> use hard coded model
nulib_table = "../../external/tables/NuLib.h5" -- NuLib opacity/emissivity table

-- spectrum parameters

spec_time_grid   = {1,1,1}         -- spectrum time grid parameters {start,stop,delta} (start==stop --> single bin catch-all)
spec_nu_grid = {0,3e22,1e20}    -- spectrum frequency grid parameters {start,stop,delta} (start==stop --> single bin catch-all)
spec_n_mu             = 1          -- number of cos(theta) bins in output spectrum 
spec_n_phi            = 1          -- number of phi bins in output spectrum

-- particle creation parameters

init_particles = 0                 -- # particles spawned from the pre-existing 'radiation energy' in each zone
n_emit_core    = 0                 -- # particles to emit from core each timestep
n_emit_heat    = 0                 -- # particles to emit from zones each timestep ("actual" emission, ignored if radiative_eq)
n_emit_visc    = 1e4               -- # particles to emit from zones each timestep (from viscosity, ignored if !radiative_eq)
n_emit_decay   = 0                 -- # particles to emit from zones each timestep (from non-thermal processes)
max_particles  = 1e6
visc_specific_heat_rate = 1.0

-- opacity parameters

nut_grey_opacity    =  -1          -- optical grey opacity (cm^2/g)
nut_epsilon         =  -1          -- absorption fraction

-- particle propagation parameters

step_size = 0.4                    -- move at most step_size*min_grid_length at a time

-- equilibrium solver parameters

damping = 0.5                      -- changes in values between iterations are decreased by this factor
