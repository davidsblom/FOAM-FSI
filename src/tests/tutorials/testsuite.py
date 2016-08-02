#!/usr/bin/env python
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import os, shutil, subprocess
import xml.etree.ElementTree

# Run the tutorials for one time step

os.chdir( "../../../tutorials/fsi" )
mainDir = os.getcwd()
for tutorial in os.listdir("."):
    print tutorial

    os.chdir( mainDir + "/" + tutorial )

    if not os.path.isfile( "fluid/system/controlDict" ): continue 
    
    controlDict = ParsedParameterFile( "fluid/system/controlDict" )
    controlDict['endTime'] = controlDict['deltaT']
    controlDict['writeInterval'] = 1
    controlDict['writeControl'] = "timeStep"
    controlDict['startFrom'] = "startTime"
    controlDict.writeFile()

    status = subprocess.call( "./Allrun", shell = True )
    assert status == 0
    
    simulationCompleted = False

    with open( "fluid/log.fsiFoam" ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True

    assert simulationCompleted == True
    
    if not os.path.isfile( "Allrun_precice" ): continue

    # Edit precice configuration file
    f = open( "fluid/constant/preCICE.xml", 'r' )
    precice = f.readlines()
    f.close()
    i = 0
    for line in precice:
         if "<max-timesteps value=" in line:
             line = "<max-timesteps value=\"1\" />\n"
             precice[i] = line
         i += 1
    f = open( "fluid/constant/preCICE.xml", 'w' )
    f.write( ''.join(precice) )
    f.close()

    controlDict = ParsedParameterFile( "solid/system/controlDict" )
    controlDict['endTime'] = controlDict['deltaT']
    controlDict['writeInterval'] = 1
    controlDict['writeControl'] = "timeStep"
    controlDict['startFrom'] = "startTime"
    controlDict.writeFile()
    
    status = subprocess.call( "./Allrun_precice", shell = True )
    assert status == 0
    
    simulationCompleted = False

    with open( "fluid/log.fsiFluidFoam" ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True
    
    assert simulationCompleted == True
    
    simulationCompleted = False

    with open( "solid/log.fsiSolidFoam" ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True

    assert simulationCompleted == True
