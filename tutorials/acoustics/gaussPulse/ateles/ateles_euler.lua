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
degree = 4
cubeLength = 2.0
-- global simulation options
simulation_name='ateles_euler'
fin_time = 0.012
sim_control = {
             time_control = {
                              min = 0,
                              max = {iter=910}, -- final simulation ti
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
--            read = './restart_lob_nearest/left/ateles_left_lastHeader.lua',
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
          filename = 'timing_euler.res'         -- the filename of the timing results
         }

-- Equation definitions --
equation = {
    name   = 'euler',
    therm_cond = 2.555e-02,
    isen_coef = 1.4,
    r      = 296.0,
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
    -- the spatial discretization scheme
    spatial =  {
               name = 'modg',           -- we use the modal discontinuous Galerkin scheme
               modg_space = 'Q',        -- the polynomial space Q or P
               m = degree,                    -- the maximal polynomial degree for each spatial direction
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
projection = {
              kind = 'fpt',  -- 'fpt' or 'l2p', default 'l2p'
              -- for fpt the  nodes are automatically 'chebyshev'
              -- for lep the  nodes are automatically 'gauss-legendre'
              factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
              lobattoPoints = true,  -- if lobatto points should be used, default = false
           -- blocksize = 32,        -- for fpt, default -1
           -- fftMultiThread = false -- for fpt, logical, default false
             }

-- This is a very simple example to define constant boundary condtions.
-- Transport velocity of the pulse in x direction.
velocityX = 0.0
press = 100000
dens = 1.225

function ic_pressure_gauss (x,y,z)
    d= (x+1*x+1)+y*y+z*z
    return( press + 100* math.exp(-d/0.2*math.log(2)) )
end

initial_condition = { density = dens,
                      pressure = ic_pressure_gauss,
                      --pressure = press,
                      velocityX = velocityX,
                      velocityY = 0.0,
                      velocityZ = 0.0,
                    }

 -- Boundary definitions
boundary_condition = {
                         {
                           label = 'wall_1',
                           kind = 'outflow',
                           pressure = press,
                         }
                         ,
                         {
                           label = 'wall_2',
                           kind = 'outflow',
                           pressure = press,
                         }
                         ,
                         {
                           label = 'wall_4',
                           kind = 'outflow',
                           pressure = press,
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
                                                  },
                          exchange_data_write =  {'Acoustics_Density',
                                                 'Acoustics_Velocity_X',
                                                 'Acoustics_Velocity_Y',
                                                 'Acoustics_Velocity_Z',
                                                 'Acoustics_Pressure'
                                                 }
                         }
                         ,
                         {
                           label = 'wall_5',
                           kind = 'outflow',
                           pressure = press,
                         }
                         ,
                         {
                           label = 'wall_6',
                           kind = 'outflow',
                           pressure = press,
                         }
                       }
