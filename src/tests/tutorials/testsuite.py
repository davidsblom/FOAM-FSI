#!/usr/bin/env python
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import os, shutil, subprocess

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
    
    status = subprocess.call( "./Allrun_precice", shell = True )
    assert status == 0
    
    simulationCompleted = False

    with open( "fluid/log.fsiFluidFoam" ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True
    
    assert simulationCompleted == True
