#!/usr/bin/env python

from operator import itemgetter

fileName = 'tests_0.log'

tests = []

with open( fileName ) as f:
    for line in f:
        if ' OK ]' in line:
            substring = line[ line.find( '(' ) + 1 : ]
            timing = int( substring[ : substring.find( 'ms' ) ] )
            name = line[ line.find( 'OK ]' ) + 5 : line.find( '(' )-1 ]
            tests.append( (timing,name) )

tests = sorted( tests, key=itemgetter(0), reverse=True )

for i in range( 20 ):
    print tests[i][0], 'ms', tests[i][1]
