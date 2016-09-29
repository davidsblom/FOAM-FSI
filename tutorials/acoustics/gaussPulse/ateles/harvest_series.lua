-- Use the configuration of the original simulation run.
require 'ateles_harvester'
logging = {level = 10}

-- Set the restart data to harvest.
restart.read = 'restart/ateles_header_20.000E-06.lua'

-- Subsampling
ply_sampling = { nlevels = 1 }

-- Example tracking to generate vtk files:
tracking = {
  { label = 'visu',
    variable = {'density', 'velocity', 'pressure', 'momentum', 'energy', 'mach_number', 'temperature'},
    shape = {kind='global'},
    folder = 'harvester/',
    output = {format = 'vtk', write_pvd=false}
  }
}
