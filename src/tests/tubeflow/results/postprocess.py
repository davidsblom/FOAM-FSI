#!/usr/bin/env python

import os
import numpy as np
from tabulate import tabulate

table = []

referenceTime = 1e50

for resultFile in sorted( os.listdir( "." ) ):

    if not ".log" in resultFile: continue

    with open( resultFile ) as file:

        computationalTime = 1e50

        for line in file:

            if "elapsed time" in line:
                computationalTime = float( line[15:-1] )

        if computationalTime < referenceTime:
          referenceTime = computationalTime

for resultFile in sorted( os.listdir( "." ) ):

    if not ".log" in resultFile: continue

    with open( resultFile ) as file:

        solver = "solver"
        parallel = False
        nbReuse = 0
        computationalTime = 0
        nbLevels = 0
        nbTimeSteps = 0
        noIterLevel0 = 0
        noIterLevel1 = 0
        noIterLevel2 = 0
        speedup = 0

        for line in file:

            if "solver" in line:
                solver = line[9:-1]

            if "parallel" in line:
                if line[11:-1] == "0": parallel = False
                if line[11:-1] == "1": parallel = True

            if "nbIter" in line and "level" not in line:
                noIterLevel1 = int( line[9:-1] )

            if "nbIter" in line and "level" in line:
                if "level 0" in line:
                    noIterLevel2 = int( line[16:-1] )
                if "level 1" in line:
                    noIterLevel1 = int( line[16:-1] )
                if "level 2" in line:
                    noIterLevel0 = int( line[16:-1] )

            if "nbLevels" in line:
                nbLevels = int( line[10:-1] )

            if "nbReuse" in line:
                nbReuse = int( line[10:-1] )

            if "elapsed time" in line:
                computationalTime = float( line[15:-1] )

            if "nbTimeSteps" in line:
                nbTimeSteps = int( line[13:-1] )

        if solver == "IQN-ILS" and nbLevels == 3: continue

        solver += "(" + str(nbReuse) + ")"

        if nbLevels == 2:
            noIterLevel0 = noIterLevel1
            noIterLevel1 = 0

        duration = computationalTime / referenceTime

        noIterLevel0 = float( noIterLevel0 ) / float( nbTimeSteps )
        noIterLevel1 = float( noIterLevel1 ) / float( nbTimeSteps )
        noIterLevel2 = float( noIterLevel2 ) / float( nbTimeSteps )

        assert nbTimeSteps == 100

        row = [ solver, noIterLevel0, noIterLevel2, duration ]

        if not parallel:
          table.append( row )

header = ["Algorithm", "Iterations fine model", "Iterations coarse model", "Duration" ]
print tabulate( table, headers = header, floatfmt = ".1f" )

# print tabulate( table, headers = header, floatfmt = ".1f", tablefmt = "latex" )
