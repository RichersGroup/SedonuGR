-- what are we simulating?

do_photons   = 1                 -- simulate photons?
do_neutrinos = 0                 -- simulate neutrinos?
iterate      = 1                 -- iterative calculation (solve for steady-state configuration)? 
solve_T      = 1                 -- (if iterative) solves each zone's temperature based on its absorbed energy
solve_Ye     = 0                 -- (if iterative) solves each zone's Ye based on its absorbed lepton number
do_visc      = 0

-- input/output files

grid_type = "grid_1D_sphere"
model_file  =  "lucy.mod"
write_zones_every   = 1
write_rays_every    = 1
write_spectra_every = 1

-- spectrum parameters

phot_spec_n_mu           =    1         -- number of cos(theta) bins in output spectrum 
phot_spec_n_phi          =    1         -- number of phi bins in output spectrum
phot_spec_time_grid = {1,1,1}           -- {start, stop, bin width}
phot_spec_nu_grid   = {0,2e15,2e13}     -- {start, stop, bin width}

-- particle creation parameters

n_emit_core    = 1e4               -- # particles to emit from core each timestep
n_emit_therm   = 0                 -- # particles to emit from zones each timestep ("actual" emission, ignored if radiative_eq)
n_emit_decay   = 0                 -- # particles to emit from zones each timestep (from non-thermal processes)
max_particles = 1e6

-- particle propagation parameters

max_n_steps = 2
dt = 1
step_size = 0.4                    -- move at most step_size*min_grid_length at a time

-- inner source

r_core = 2e15                       -- core radius (cm)
L_core = 1e43                       -- core luminosity (erg/s)
T_core = 10000                      -- core temperature (K)

-- opacity parameters

phot_grey_opacity    =  0.0 --1           -- optical grey opacity (cm^2/g)
phot_epsilon         =  1.0           -- absorption fraction
phot_nu_start        =  0             -- photon opacity grid start
phot_nu_stop         =  2e15          -- photon opacity grid stop
phot_n_nu            =  100           -- photon opacity grid number of points along frequency

-- equilibrium solver parameters

damping = 0.5                    -- changes in values between iterations are decreased by this factor
brent_itmax = 100
brent_tolerance = 0.01