#!/usr/bin/env python

import os, subprocess, multiprocessing, time, sys, argparse

nbCores = int( os.environ['WM_NCOMPPROCS'] )

parser = argparse.ArgumentParser( description='Run the test suite' )
parser.add_argument('testsuite', help='which testsuite: testsuite-dealii, testsuite-rbf, testsuite-spacemapping, testsuite-fsi, or testsuite-sdc' )
args = parser.parse_args()

runs = []

for i in range( nbCores ):
    run = subprocess.Popen( "GTEST_TOTAL_SHARDS=" + str(nbCores) + " GTEST_SHARD_INDEX=" + str(i) + " " + args.testsuite + " --gtest_throw_on_failure > tests_" + str(i) + ".log 2>&1", shell = True )
    runs.append( run )

i = 0
returnCode = 0
for run in runs:
    code = None
    while code is None:
        print i
        i += 1
        sys.stdout.flush()
        code = run.poll()
        if code > returnCode:
            returnCode = code
        if code is not None:
            break
        time.sleep( 1 )
    print 'Run finished'
print 'All runs finished'

for i in range( nbCores ):
    subprocess.call("tail -n 20 tests_" + str(i) + ".log 2>&1", shell=True)

if returnCode == 0:
    print "Finished successfully"
else:
    print "Tests failed"
exit( returnCode )
