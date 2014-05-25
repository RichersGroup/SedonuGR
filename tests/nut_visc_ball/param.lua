-- what are we simulating?

do_photons   = 0                 -- simulate photons?
do_neutrinos = 1                 -- simulate neutrinos?
steady_state = 1                 -- iterative calculation (solve for steady-state configuration)? 
radiative_eq = 0
do_visc      = 1
solve_T      = 1                 -- (if iterative) solves each zone's temperature based on its absorbed energy
solve_Ye     = 1                 -- (if iterative) solves each zone's Ye based on its absorbed lepton number
reflect_outer= 0

-- input/output files

grid_type = "grid_1D_sphere"       -- grid geometry. Must match grid geometry in model file if used  
model_file  =  "neutron_star.mod"  -- model file. "custom" --> use hard coded model
nulib_table = "../../external/tables/NuLib_LS220_rho150_temp90_ye60_ng24_ns3_version1.0_20140327.h5" -- NuLib opacity/emissivity table
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1

-- spectrum parameters

nut_spec_time_grid  = {1,1,1}          -- spectrum time grid parameters {start,stop,delta} (start==stop --> single bin catch-all)
nut_spec_nu_grid    = {0,2e22,1e20}    -- spectrum frequency grid parameters {start,stop,delta} (start==stop --> single bin catch-all)
nut_spec_n_mu       = 1                -- number of cos(theta) bins in output spectrum 
nut_spec_n_phi      = 1                -- number of phi bins in output spectrum

-- particle creation parameters

n_emit_core    = 0                 -- # particles to emit from core each timestep
n_emit_therm   = 1e6                 -- # particles to emit from zones each timestep
n_emit_decay   = 0                 -- # particles to emit from zones each timestep (from non-thermal processes)
max_particles  = 1e6
visc_specific_heat_rate = 1e22

-- particle propagation parameters

max_n_steps = 10
dt = -1
step_size = 0.4                    -- move at most step_size*min_grid_length at a time

-- opacity parameters

nut_grey_opacity    =  -1          -- optical grey opacity (cm^2/g)
nut_epsilon         =  -1          -- absorption fraction

-- equilibrium solver parameters

damping = 0.0                      -- changes in values between iterations are decreased by this factor
brent_itmax = 100
brent_tolerance = 0.01