---- testcase setup------------------------------------------------------------
gamma      = 1.4
r          = 296.0
therm_cond = 0
mu         = 0
ip_param   = 4.0

function ic_pressure_gauss (x,y,z)
    d= (x+1*x+1)+y*y+z*z
    return( press + 100* math.exp(-d/0.2*math.log(2)) )
end

dens = 1.225
velX = 0.0
velY = 0.0
velZ = 0.0
press = 100000
c = math.sqrt(gamma* press / dens)
gp_center = {0.0, 0.0, 0.0}
gp_halfwidth = 3.0
gp_amplitude = 0.001
gp_background = 0.0
gp_c = c

-----simulation setup-----------------------------------------------------------
simulation_name='ateles'
logging = {level=10}
timestep_info = 1
check =  { interval = 1 }
tmax = 0.012
-- order of scheme
degree = 2
-- use fix timestep
dt = 1.25e-02
sim_control = {
  time_control = {
    min = 0,
    max = tmax,
    interval = {iter = 100}, -- final simulation time
  }
}

-- Mesh definitions --
mesh = 'mesh/'

-- timing setting
timing_file = 'timing.res'

-- table for preCICE
precice = {
           accessor = 'Ateles',
           configFile ='precice.xml',
          }

restart = {
  write = 'restart/',
  time_control = {
    min = 0,
    max = tmax,
    interval = {iter = 100}, -- final simulation time
  }
}
-- Variable system definintion--
characteristic = 0.0
function relax_velocity(x,y,z,t)
  return {0.0,0.0, 0.0}
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
     name = 'coupling',
     ncomponents = 10,
     vartype = 'st_fun',
     st_fun = {
       predefined = 'precice',
       precice_mesh = 'Ateles_Acoustics',
       read_varname = { 'Acoustics_Density_Gradient',
                        'Acoustics_Velocity_X_Gradient',
                        'Acoustics_Velocity_Y_Gradient',
                        'Acoustics_Velocity_Z_Gradient',
                        'Acoustics_Pressure_Gradient'},
       write_varname = {'Acoustics_Density',
                        'Acoustics_Velocity_X',
                        'Acoustics_Velocity_Y',
                        'Acoustics_Velocity_Z',
                        'Acoustics_Pressure'}
    }
  },
  -- write to precice
  {
    name = 'Acoustics_Density',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'density',
      input_varindex = {1}
    }
  },
  {
    name = 'Acoustics_Velocity_X',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'velocity',
      input_varindex = {1}
    }
  },
  {
    name = 'Acoustics_Velocity_Y',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'velocity',
      input_varindex = {2}
    }
  },
  {
    name = 'Acoustics_Velocity_Z',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'velocity',
      input_varindex = {3}
    }
  },
  {
    name = 'Acoustics_Pressure',
    ncomponents = 1,
    vartype = 'operation',
    operation = {
      kind = 'extract',
      input_varname = 'pressure',
      input_varindex = {1}
    }
  },
  -- read from precice
  {
     name = 'Acoustics_Density_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'coupling',
       input_varindex = {1}
     }
  },
  {
     name = 'Acoustics_Velocity_X_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'coupling',
       input_varindex = {2}
     }
  },
  {
     name = 'Acoustics_Velocity_Y_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'coupling',
       input_varindex = {3}
     }
  },
  {
     name = 'Acoustics_Velocity_Z_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'coupling',
       input_varindex = {4}
     }
  },
  {
     name = 'Acoustics_Pressure_Gradient',
     ncomponents = 1,
     vartype = 'operation',
     operation = {
       kind = 'extract',
       input_varname = 'coupling',
       input_varindex = {5}
     }
  }
}
-- Scheme definitions --
scheme = {
  -- the spatial discretization scheme
  spatial =  {
    name = 'modg',
    modg_space = 'Q',
    m = degree,
  },
  -- the temporal discretization scheme
  temporal = {
    name = 'explicitSSPRungeKutta',
    steps = 2,
    control = {
      name = 'fixed',
      dt = dt
    },
  },
}
-- the general projection table --
projection = {
  kind = 'l2p',
  factor = 1.0,
}

-- Equation definitions --
penalization_eps = 8.0/(degree+1)
penalization_alpha = 1.0
equation = {
  penalization ={
    global = {
      kind = 'const',
      const = {0.0, 0.0, 0.0,0.0, 0.0},
    },
  },
  name       = 'navier_stokes',
  isen_coef  = gamma,
  r          = r,
  therm_cond = therm_cond,
  mu         = mu,
  ip_param   = ip_param,
  -- Parameters of the penalization
  porosity             = penalization_eps,
  viscous_permeability = penalization_alpha*penalization_eps,
  thermal_permeability = penalization_alpha*penalization_eps,
  material = {
    characteristic = 'characteristic',
    relax_velocity = 'relax_velocity',
    relax_temperature = 'relax_temperature'
  },
}
-- (cv) heat capacity and (r) ideal gas constant
equation["cv"] = equation["r"] / (equation["isen_coef"] - 1.0)

-- initial condition --
initial_condition = {
  density = dens,
  velocityX = velX,
  velocityY = velY,
  velocityZ = velZ,
  pressure = ic_pressure_gauss,
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
    label = 'precice_mesh',
    kind = 'gradients',
    grad_density   = 'Acoustics_Density_Gradient',
    grad_velocityX = 'Acoustics_Velocity_X_Gradient',
    grad_velocityY = 'Acoustics_Velocity_Y_Gradient',
    grad_velocityZ = 'Acoustics_Velocity_Z_Gradient',
    grad_pressure  = 'Acoustics_Pressure_Gradient'
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
