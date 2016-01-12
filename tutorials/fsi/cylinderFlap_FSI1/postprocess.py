#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import os
import numpy as np
import matplotlib.pyplot as plt
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile

fileName = 'fluid/forces/0/forces.dat'

liftList = []
dragList = []

controlDict = ParsedParameterFile( 'fluid/system/controlDict' )
deltaT = controlDict['deltaT']

index = 0

with open ( fileName ) as FILE:

    for line in FILE:

        if index > 0:

            subString = line[ line.find( '(((' ) + 3 :-1 ]
            pressureForce = float( subString[ : subString.find( ' ' ) ] )

            pressureForce = np.array( subString[ : subString.find( ')' ) ].split(' '), dtype = np.float64 )

            subString = subString[ subString.find( '(' ) + 1 :  ]
            viscousForce = np.array( subString[ : subString.find( ')' ) ].split(' '), dtype = np.float64 )

            lift = pressureForce[1] + viscousForce[1]
            drag = pressureForce[0] + viscousForce[0]

            liftList.append( lift )
            dragList.append( drag )

        index += 1

fileName = "fluid/probesSolid/solid/0/U"

U = open( fileName, 'r' ).readlines()

Ux = []
Uy = []

for line in U:

    if "#" in line:
        continue

    line = line.strip()
    time = float( line[ : line.find(" ") ] )

    disp = line[ line.find(" ") : ].strip(" ()").split( " " )
    index = 0
    for u in disp:
        disp[index] = float( u )
        index += 1

    Ux.append( disp[0] )
    Uy.append( disp[1] )

timeList = np.arange( 1, len(liftList)+1 ) * deltaT

plt.figure( figsize=(20, 10) )

plt.subplot(221)
plt.grid()
plt.plot( timeList, liftList )
plt.xlabel( "Time [s]" )
plt.ylabel( "Lift [N]" )
(x1, x2, y1, y2) = plt.axis()
plt.axis( (timeList[0], x2, y1, y2) )

plt.subplot(222)
plt.grid()
plt.plot( timeList, dragList )
plt.xlabel( "Time [s]" )
plt.ylabel( "Drag [N]" )
(x1, x2, y1, y2) = plt.axis()
plt.axis( (timeList[0], x2, y1, y2) )

plt.subplot(223)
plt.plot( timeList, Ux )
plt.grid()
plt.xlabel( "Time [s]" )
plt.ylabel( "Displacement x [m]" )
(x1, x2, y1, y2) = plt.axis()
plt.axis( (timeList[0], x2, y1, y2) )

plt.subplot(224)
plt.plot( timeList, Uy )
plt.grid()
plt.xlabel( "Time [s]" )
plt.ylabel( "Displacement y [m]" )
(x1, x2, y1, y2) = plt.axis()
plt.axis( (timeList[0], x2, y1, y2) )

plt.savefig( "forces-displacements.pdf", bbox_inches = 'tight' )
plt.close()
