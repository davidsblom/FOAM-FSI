-- Configuration file for Ateles --


-- This is a configuration file for the Finite Volume / Discontinuous Galerkin Solver ATELES.
-- It provides a testcase for the simulation of Euler equations in a homogenous media. The simulation domain
-- is a periodic cube with edge length 2.0. Therefore this is a very good way to verify your algorithmic implementations,
-- since this testcase does not involve any boundary conditions.
-- The testcase simulates the temporal development of Gaussian pulse in density. Since we
-- are considering a very simple domain, an analytic solution is well known and given as Lua functions in this script.
-- Therefore we suggest this testcase to verify one of the following topics
-- ... algorihtmic correctness
-- ... spatial and temporal order of your scheme
-- ... diffusion and dispersion relations of your scheme
-- ... and many more.
-- This testcase can be run in serial (only one execution unit) or in parallel (with multiple mpi ranks).
-- To specify the number of ranks please modify nprocs variable. To calculate a grid convergence behavior please modify the
-- level variable. An increment of one will half the radius of your elements.

timestep_info = 1

-- Check for Nans and unphysical values
check =  {
           interval = 1,
         }

logging = { level = 5 }
cubeLength = 2.0
-- global simulation options
simulation_name='ateles_lineareuler'
fin_time = 0.012
degree = 4
sim_control = {
             time_control = {
                  min = 0,
                  max = {iter = 910}, -- final simulation time
                  interval = {iter=1}
                }
}

-- table for preCICE
precice = {
           accessor = 'Ateles_acoustic',
           configFile ='precice.xml',
          }

-- Mesh definitions --
mesh = 'mesh_right/'

---- Restart settings
restart = {
--            -- file to restart from
--            read = './restart/lineareuler/ateles_lineareuler_lastHeader.lua',
--            -- folder to write restart data to
            write = './restart/',
            -- temporal definition of restart write
            time_control = { min = 0,
                             max = {iter = 900},
                             interval = {iter = 100}
                           }
          }

-- timing settings (i.e. output for performance measurements, this table is otional)
timing = {
          folder = './',                  -- the folder for the timing results
          filename = 'timing_lineareuler.res'         -- the filename of the timing results
         }

-- Equation definitions --
velocityX = 0.0
press = 100000
dens = 1.225
equation = {
    name   = 'linearEuler',
    therm_cond = 2.555e-02,
    isen_coef = 1.4,
    r      = 296.0,
    background = {
                      density=dens,
                      velocityX = 0.0,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                      pressure = press,
                 }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
    -- the spatial discretization scheme
    spatial =  {
               name = 'modg',           -- we use the modal discontinuous Galerkin scheme
               modg_space = 'Q',        -- the polynomial space Q or P
               m = degree,                   -- the maximal polynomial degree for each spatial direction
               },
    -- the temporal discretization scheme
    temporal = {
                name = 'explicitSSPRungeKutta',  -- we use ssp explicit runge kutta in time
                steps = 2,
               -- how to control the timestep
               control = {
                          name = 'cfl',   -- the name of the timestep control mechanism
                          cfl  = 0.8,     -- Courant�Friedrichs�Lewy number
                         },
               },
}

function dens(x,y,z)
  return x
end

projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
              lobattoPoints = false,  -- if lobatto points should be used, default = false
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
           -- blocksize = 32,        -- for fpt, default -1
           -- fftMultiThread = false -- for fpt, logical, default false
             }

-- This is a very simple example to define constant boundary condtions.
-- Transport velocity of the pulse in x direction.
initial_condition = { density = 0.0,
                      velocityX = 0.0,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                      pressure = 0.0,
                    }


 -- Boundary definitions
 boundary_condition = {
                         {
                           label = 'wall_1',
                           kind = 'primitives',
                           density = 0.0,
                           v_x = 0.0,
                           v_y = 0.0,
                           v_z = 0.0,
                           pressure = 0.0
                         }
                         ,
                         {
                           label = 'wall_2',
                           kind = 'primitives',
                           density = 0.0,
                           v_x = 0.0,
                           v_y = 0.0,
                           v_z = 0.0,
                           pressure = 0.0
                         }
                         ,
                         {
                           label = 'precice_rightmesh',
                           kind = 'precice',
                           precice_mesh = 'AcousticSurface_Ateles',

                           exchange_data_read =  {'Acoustics_Density',
                                                  'Acoustics_Velocity_X',
                                                  'Acoustics_Velocity_Y',
                                                  'Acoustics_Velocity_Z',
                                                  'Acoustics_Pressure'
                                                 }
                         }
                         ,
                         {
                           label = 'wall_4',
                           kind = 'primitives',
                           density = 0.0,
                           v_x = 0.0,
                           v_y = 0.0,
                           v_z = 0.0,
                           pressure = 0.0
                         }
                         ,
                         {
                           label = 'wall_5',
                           kind = 'primitives',
                           density = 0.0,
                           v_x = 0.0,
                           v_y = 0.0,
                           v_z = 0.0,
                           pressure = 0.0
                         }
                         ,
                         {
                           label = 'wall_6',
                           kind = 'primitives',
                           density = 0.0,
                           v_x = 0.0,
                           v_y = 0.0,
                           v_z = 0.0,
                           pressure = 0.0
                         }
                       }
