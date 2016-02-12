#!/usr/bin/env python

import os, subprocess, multiprocessing, time

nbCores = multiprocessing.cpu_count()

runs = []

for i in range( nbCores ):
    run = subprocess.Popen( "GTEST_TOTAL_SHARDS=" + str(nbCores) + " GTEST_SHARD_INDEX=" + str(i) + " tests > tests_" + str(i) + ".log 2>&1", shell = True )
    runs.append( run )

returnCode = 0
for run in runs:
    while ( not run.poll() ):
        time.sleep( 5 )
        print '.'
        code = run.poll()
        if code > returnCode:
            returnCode = code

for i in range( nbCores ):
    subprocess.call("tail -n 50 tests_" + str(i) + ".log 2>&1", shell=True)

exit( returnCode )
