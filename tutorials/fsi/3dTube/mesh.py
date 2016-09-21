#!/usr/bin/env python

import argparse, mako.template

parser = argparse.ArgumentParser( description='Generate blockMesh configuration file for fluid and solid domains' )
parser.add_argument('level', type=int, help='Mesh refinement level (0,1,2,...,N)' )

args = parser.parse_args()

template = mako.template.Template( filename = 'fluid/constant/polyMesh/blockMeshDict.template' )

f = open( 'fluid/constant/polyMesh/blockMeshDict', 'w' )
output = template.render( level = args.level )
f.write( output )
f.close()

template = mako.template.Template( filename = 'solid/constant/polyMesh/blockMeshDict.template' )

f = open( 'solid/constant/polyMesh/blockMeshDict', 'w' )
output = template.render( level = args.level )
f.write( output )
f.close()
