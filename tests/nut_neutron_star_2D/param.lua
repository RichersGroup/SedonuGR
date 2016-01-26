-- what are we simulating?

verbose      = 1
do_photons   = 0                 -- simulate photons?
do_neutrinos = 1                 -- simulate neutrinos?
radiative_eq = 1
steady_state = 1                 -- iterative calculation (solve for steady-state configuration)? 
solve_T      = 1                 -- (if iterative) solves each zone's temperature based on its absorbed energy
solve_Ye     = 0                 -- (if iterative) solves each zone's Ye based on its absorbed lepton number
do_visc=0
reflect_outer= 0

-- input/output files

grid_type = "grid_2D_sphere"       -- grid geometry. Must match grid geometry in model file if used  
model_file  =  "custom"  -- model file. "custom" --> use hard coded model
nulib_table = "../../external/tables/NuLib_LS220_rho100_temp60_ye40_ng24_ns3_Itemp10_Ieta10_version1.0_20140629.h5" -- NuLib opacity/emissivity table
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1

-- spectrum parameters

nut_spec_time_grid   = {1,1,1}         -- spectrum time grid parameters {start,stop,delta} (start==stop --> single bin catch-all)
nut_spec_nu_grid = {0,1e22,1e20}    -- spectrum frequency grid parameters {start,stop,delta} (start==stop --> single bin catch-all)
nut_spec_n_mu             = 32          -- number of cos(theta) bins in output spectrum 
nut_spec_n_phi            = 32          -- number of phi bins in output spectrum

-- particle creation parameters

n_emit_core    = 1e7               -- # particles to emit from core each timestep
n_emit_therm   = 0                 -- # particles to emit from zones each timestep
n_emit_decay   = 0                 -- # particles to emit from zones each timestep (from non-thermal processes)
max_particles  = 1e7

-- particle propagation parameters

max_n_iter = 20
dt = -1
step_size = 0.4                    -- move at most step_size*min_grid_length at a time

-- inner source

r_core = 9e5                       -- core radius (cm)
T_core = 10 -- 10MeV           -- core temperature (K)
core_nue_chem_pot = 0
core_lum_multiplier = 1.0

-- opacity parameters

nut_grey_opacity    =  -1          -- optical grey opacity (cm^2/g)
nut_grey_abs_frac   =  -1          -- absorption fraction

-- equilibrium solver parameters

damping         = 0.1              -- changes in values between iterations are decreased by this factor
brent_tolerance = 0.01             -- how close must the measured emission be to the predicted?
brent_itmax     = 100              -- how many iterations of brent cycle (both inner and outer loop) are allowed
