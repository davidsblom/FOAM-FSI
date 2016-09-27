timestep_info = 1

logging = {level=10}

-- Check for Nans and unphysical values
check =  { interval = 1 }

-- global simulation options
simulation_name='ateles'
-- general fluid's parameter
isen_coef  = 1.4
r          = 296.0
-- flow state
dens = 1.225
press = 100000
velocityX = 0.0
velocityY = 0.0
velocityZ = 0.0
-- numerical parameter
degree = 1
dt = 1e-5
tmax = 1.0

function ic_pressure_gauss (x,y,z)
    d= (x+1*x+1)+y*y+z*z
    return( press + 50* math.exp(-d/0.1*math.log(2)) )
end

sim_control = {
             time_control = {
                  min = 0,
                  max = {iter=100}, --tmax, -- final simulation time
                  interval = {iter = 10}, -- final simulation time
                }
}


tracking = {
  {
    label = 'center',
    folder = './harvester/',
    shape = {kind = 'all'},
    variable = {'density','pressure','temperature'},
    time_control = {min = 0, max = tmax, interval = {iter=1}},
    output = { format = 'vtk'},
  }
}

-- table for preCICE
precice = {
           accessor = 'Ateles',
           configFile ='precice.xml',
          }

---- Restart settings
restart = {
  write = './restart/',
  time_control = {
    min = 0,
    max = tmax,
    interval = {iter=1}
  }
}

-- Variable system definintion--
characteristic = 0.0
function relax_velocity(x,y,z,t)
  return {0.0, 0.0, 0.0}
end
relax_temperature = 0.0

variable = {
  {
     name = 'characteristic',
     ncomponents = 1,
     vartype = "st_fun",
     st_fun = characteristic
  },
  {
     name = 'relax_velocity',
     ncomponents = 3,
     vartype = "st_fun",
     st_fun = relax_velocity
  },
  {
     name = 'relax_temperature',
     ncomponents = 1,
     vartype = "st_fun",
     st_fun = relax_temperature
  },
  {
     name = 'grad_pressure',
     ncomponents = 3,
     vartype = 'operation',
     operation = {
       kind = 'gradient',
       input_varname = 'pressure',
     }
  },
  {
     name = 'grad_velocity',
     ncomponents = 9,
     vartype = 'operation',
     operation = {
       kind = 'gradient',
       input_varname = 'velocity',
     }
  },
  -- {
  --   name = 'Temp',
  --   ncomponents = 1,
  --   vartype = 'operation',
  --   operation = {
  --     kind = 'extract',
  --     input_varname = 'temperature',
  --     input_varindex = {1}
  --   }
  -- },
  {
    name = 'grad_temp',
    ncomponents = 3,
    vartype = 'operation',
    operation = {
      kind = 'gradient',
      --input_varname = 'pressure',
      input_varname = 'temperature',
    }
  },
 {
    name = 'coupling',
    ncomponents = 5,
    vartype = 'st_fun',
    st_fun = {
      predefined = 'precice',
      precice_mesh = 'Ateles_Acoustics',
      write_varname = {
                       'Acoustics_Temperature_Gradient',
                       'Acoustics_Velocity_X_Gradient',
                       'Acoustics_Velocity_Y_Gradient',
                       'Acoustics_Velocity_Z_Gradient',
                       'Acoustics_Pressure_Gradient',
                       },
      read_varname = {'Acoustics_Density',
                      'Acoustics_Velocity_X',
                      'Acoustics_Velocity_Y',
                      'Acoustics_Velocity_Z',
                      'Acoustics_Pressure'}
    }
  },
  -- write to precice
  {
     name = 'Acoustics_Temperature_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'grad_temp',
       input_varindex = {1}
     }
  },
  {
     name = 'Acoustics_Velocity_X_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'grad_velocity',
       input_varindex = {1}
     }
  },
  {
     name = 'Acoustics_Velocity_Y_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'grad_velocity',
       input_varindex = {4}
     }
  },
  {
     name = 'Acoustics_Velocity_Z_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'grad_velocity',
       input_varindex = {8}
     }
  },
  {
     name = 'Acoustics_Pressure_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'grad_pressure',
       input_varindex = {1}
     }
  },
--  -- read from precice
 {
    name = 'Acoustics_Density',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'coupling',
      input_varindex = {1}
    }
 },
 {
    name = 'Acoustics_Velocity_X',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'coupling',
      input_varindex = {2}
    }
 },
 {
    name = 'Acoustics_Velocity_Y',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'coupling',
      input_varindex = {3}
    }
 },
 {
    name = 'Acoustics_Velocity_Z',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'coupling',
      input_varindex = {4}
    }
 },
 {
    name = 'Acoustics_Pressure',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'coupling',
      input_varindex = {5}
    }
 },
}

-- Mesh definitions --
mesh = 'mesh/'

-- timing settings (i.e. output for performance measurements, this table is otional)
timing_file = 'timing_right.res'         -- the filename of the timing results

-- Equation definitions --
equation = {
  penalization = {
    global = {
      kind = 'const',
      const = {0.0, 0.0, 0.0, 0.0, 0.0},
    },
  },
  name       = 'euler',
  isen_coef  = isen_coef,
  r          = r,
  material = {
    characteristic = 'characteristic',
    relax_velocity = 'relax_velocity',
    relax_temperature = 'relax_temperature'
  }
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',        -- we use the modal discontinuous Galerkin scheme
    modg_space = 'Q',        -- the polynomial space Q or P
    m = degree,                   -- the maximal polynomial degree for each spatial direction
  },
--  stabilization = {
--    {
--      name = 'spectral_viscosity',
--      alpha = 36,
--      order = 50
--    },
--  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitSSPRungeKutta',  --'explicitEuler',
    steps = 2,
    -- how to control the timestep
    control = {
--      name = 'cfl',   -- the name of the timestep control mechanism
--      cfl  = 0.2,     -- Courant�Friedrichs�Lewy number
--      cfl_visc  = 0.2,     -- Courant�Friedrichs�Lewy number
      name = 'fixed',
      dt = dt
    },
  },
}

-- ...the general projection table
projection = {
  kind = 'l2p',  -- 'fpt' or 'l2p', default 'l2p'
  factor = 1.0,          -- dealising factpr for fpt, oversampling factor for l2p, float, default 1.0
}

-- This is a very simple example to define constant boundary condtions.
-- Transport velocity of the pulse in x direction.
initial_condition = {
  density = dens,
 pressure = press,
  -- pressure = ic_pressure_gauss,
  velocityX = velocityX,
  velocityY = velocityY,
  velocityZ = velocityZ,
}

 -- Boundary definitions
boundary_condition = {
  {
    label = 'precice_mesh',
--    kind = 'outflow',
--    pressure = press,
    kind = 'primitives',
    density = 'Acoustics_Density',
    velocityX = 'Acoustics_Velocity_X',
    velocityY = 'Acoustics_Velocity_Y',
    velocityZ = 'Acoustics_Velocity_Z',
    pressure = 'Acoustics_Pressure'
  }
  ,
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
