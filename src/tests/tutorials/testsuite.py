#!/usr/bin/env python

from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
import os, shutil, subprocess, yaml
import xml.etree.ElementTree

# Run the tutorials for one time step

os.chdir( "../../../tutorials/fsi" )
mainDir = os.getcwd()
for tutorial in sorted( os.listdir(".") ):

    if tutorial != "cylinderFlap_FSI3_manifoldMapping": continue

    if not os.path.isdir( mainDir + "/" + tutorial ): continue
    os.chdir( mainDir + "/" + tutorial )

    if not os.path.isfile( "fluid/system/controlDict" ) and not os.path.isfile( "fluid-level-1/system/controlDict" ):
        continue

    try:
        controlDict = ParsedParameterFile( "fluid/system/controlDict" )
    except:
        controlDict = ParsedParameterFile( "fluid-level-1/system/controlDict" )
    controlDict['endTime'] = controlDict['deltaT']
    controlDict['writeInterval'] = 1
    controlDict['writeControl'] = "timeStep"
    controlDict['startFrom'] = "startTime"
    controlDict.writeFile()

    try:
        stream = open( "fluid/constant/fsi.yaml", "r" )
        fsi = yaml.load( stream )
        stream.close()
        fsi["coupling-scheme-implicit"]["max-iterations"] = 3
        stream = open( "fluid/constant/fsi.yaml", "w" )
        yaml.dump( fsi, stream, default_flow_style=False )
        stream.close()
    except:
        stream = open( "fluid-level-1/constant/fsi.yaml", "r" )
        fsi = yaml.load( stream )
        stream.close()
        fsi["multi-level-acceleration"]["levels"][0]["max-iterations"] = 3
        fsi["multi-level-acceleration"]["levels"][1]["max-iterations"] = 3
        stream = open( "fluid-level-1/constant/fsi.yaml", "w" )
        yaml.dump( fsi, stream, default_flow_style=False )
        stream.close()

    status = subprocess.call( "./Allrun", shell = True )
    if status != 0: subprocess.call( "cat fluid/log.fsiFoam", shell = True )
    if status != 0: subprocess.call( "cat fluid-level-1/log.fsiFoam", shell = True )
    assert status == 0

    simulationCompleted = False

    fileName = "fluid/log.fsiFoam"

    if not os.path.isfile( fileName ):
        fileName = "fluid-level-1/log.fsiFoam"

    with open( fileName ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True
            if "End" in line: simulationCompleted = True

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
         if "<max-iterations value=" in line:
             line = "<max-iterations value=\"3\" />\n"
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
    if status != 0:
        subprocess.call( "cat fluid/log.fsiFluidFoam", shell = True )
        subprocess.call( "cat solid/log.fsiSolidFoam", shell = True )
        subprocess.call( "cat fluid/log.reconstructPar", shell = True )
    assert status == 0

    simulationCompleted = False

    with open( "fluid/log.fsiFluidFoam" ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True
            if "End" in line: simulationCompleted = True

    assert simulationCompleted == True

    simulationCompleted = False

    with open( "solid/log.fsiSolidFoam" ) as f:
        for line in f:
            assert "assert" not in line
            if "Finalising parallel run" in line: simulationCompleted = True
            if "End" in line: simulationCompleted = True

    assert simulationCompleted == True
