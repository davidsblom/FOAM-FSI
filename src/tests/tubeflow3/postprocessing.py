#!/usr/bin/env python

import numpy as np
import os, tabulate

os.chdir( 'data3' )

# --------------------------------

def readReuse( content ):

    for line in content:
        if 'nbReuse =' in line:
            return int( line[10:-1] )

    return 0

# --------------------------------

def readLevels( content ):

    nbLevels = 0
    for line in content:
        if 'nbLevels =' in line:
            return int( line[11:-1] )

# --------------------------------

def readLabel( content ):

    nbLevels = readLevels( content )
    label = str( nbLevels ) + " ("

    for i in np.arange( nbLevels ):
        for line in content:

            if 'level ' + str( i ) + ' N' in line:
                label += line[12:-1]

        if i < nbLevels - 1:
            label += ", "

    label += ")"

    return label

# --------------------------------

def getLabel( file ):

    parallel = False
    if 'parallel_1' in file:
        parallel = True

    label = "S-"
    if parallel:
        label = "P-"
    if "Anderson" in file:
        label += "Anderson"
    if "Aitken" in file:
        label += "Aitken"

    return label

# --------------------------------

def readAvgIter( content, level ):

    nbLevels = readLevels( content )

    for line in content:

        if 'level ' + str( level ) + ' avgIter' in line:
            return float( line[18:-1] )

# --------------------------------

def avgIter( content ):

    nbLevels = readLevels( content )

    for line in content:

        if 'avgIter' in line:
            return float( line[10:-1] )

# --------------------------------

def readData( content, row ):

    nbLevels = readLevels( content )
    nbReuse = readReuse( content )
    label = readLabel( content )

    for i in np.arange( nbLevels ):

        nbIter = readAvgIter( content, i )

        if nbReuse == 0:
            row[1+i+3-nbLevels] = nbIter
        if nbReuse == 1:
            row[4+i+3-nbLevels] = nbIter
        if nbReuse == 4:
            # row[7+i+3-nbLevels] = nbIter
            row[4+i+3-nbLevels] = nbIter

    row[0] = label

# --------------------------------

for dt in ["0.01", "0.1"]:

    table = []

    for file in sorted( os.listdir( '.' ) ):

        row = ["","","","",""]

        if 'MM_' in file: continue
        if 'complete_jacobian' in file: continue
        if 'dt_' + dt in file: continue

        f = open( file, "r" )
        content = f.readlines()

        nbIter = avgIter( content )
        nbReuse = readReuse( content )

        if nbReuse != 0: continue

        label = getLabel( file )

        for file in sorted( os.listdir( '.' ) ):

            if 'MM_' in file: continue
            if 'dt_' + dt in file: continue
            if getLabel( file ) != label: continue

            f = open( file, "r" )
            content = f.readlines()

            nbIter = avgIter( content )

            if 'complete_jacobian' in file:

                row[4] = nbIter

            if 'complete_jacobian' not in file:

                nbReuse = readReuse( content )

                row[0] = label
                if nbReuse == 0:
                    row[1] = nbIter
                    if 'Aitken' in file:
                        table.append( row )
                if nbReuse == 1:
                    row[2] = nbIter
                if nbReuse == 4:
                    row[3] = nbIter
                    table.append( row )

    headers = ["Histories",0,1,4,"all"]

    if dt == "0.1":
        caption = 'dt = 0.01'
    if dt == "0.01":
        caption = 'dt = 0.1'

    print "\n\\begin{table}\n\centering"
    print "\caption{"+caption+"}"
    print tabulate.tabulate( table, headers, floatfmt=".1f", tablefmt="latex_booktabs" )
    print "\end{table}"

    # print tabulate.tabulate( table, headers, floatfmt=".1f" )

# --------------------------------

for dt in ["0.1", "0.01"]:

    for parallel in [1,0]:

        for solver in ["fluid_linearized", "fluid_non-linear"]:

            table = []

            for iLevel in [2,3]:

                for file in sorted( os.listdir( '.' ) ):

                    row = ["","","","","","","","","",""]

                    if file[0:3] != 'MM_': continue
                    if 'complete_jacobian' in file: continue
                    if 'parallel_' + str( parallel ) in file: continue
                    if 'dt_' + dt in file: continue
                    if solver in file: continue

                    f = open( file, "r" )
                    content = f.readlines()

                    if readLevels( content ) != iLevel: continue
                    if readReuse( content ) != 0: continue

                    label = readLabel( content )

                    for file in sorted( os.listdir( '.' ) ):

                        if solver in file: continue
                        if 'dt_' + dt in file: continue
                        if 'parallel_' + str( parallel ) in file: continue

                        f = open( file, "r" )
                        content = f.readlines()

                        if readLevels( content ) != iLevel: continue
                        if readLabel( content ) != label: continue
                        if readReuse( content ) == 1: continue

                        if 'complete_jacobian' not in file:

                            readData( content, row )

                            if readReuse( content ) == 4:
                                table.append( row )
                        else:
                            nbLevels = readLevels( content )

                            for i in np.arange( nbLevels ):
                                # row[10+i+3-nbLevels] = readAvgIter( content, i )
                                row[7+i+3-nbLevels] = readAvgIter( content, i )


            # headers = [ "Levels / Histories", "0--1", "0--2", "0--3", "1--1", "1--2", "1--3", "4--1", "4--2", "4--3", "all--1", "all--2", "all--3" ]
            headers = [ "Levels / Histories", "0--1", "0--2", "0--3", "4--1", "4--2", "4--3", "all--1", "all--2", "all--3" ]

            if dt == "0.1":
                caption = 'dt = 0.01'
            if dt == "0.01":
                caption = 'dt = 0.1'

            if parallel == 1:
                caption += ", serial coupling"

            if parallel == 0:
                caption += ", parallel coupling"

            if solver == "fluid_linearized":
                caption += ", low fidelity model: non-linear fluid model"

            if solver == "fluid_non-linear":
                caption +=  ", low fidelity model: linearized fluid model"

            print "\n\\begin{table}\n\centering"
            print "\caption{"+caption+"}"
            print tabulate.tabulate( table, headers, floatfmt=".1f", tablefmt="latex_booktabs" )
            print "\end{table}"

            # print "\n", caption
            # print tabulate.tabulate( table, headers, floatfmt=".1f" )
