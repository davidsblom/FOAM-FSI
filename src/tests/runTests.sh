#!/bin/bash

rm -rf *.log

NB_CORES=`grep -c ^processor /proc/cpuinfo`

NB_CORES=8

CORE_COUNT=`expr $NB_CORES - 1`

for i in `seq 0 $CORE_COUNT`
  do
    GTEST_TOTAL_SHARDS=$NB_CORES GTEST_SHARD_INDEX=${i} tests > tests_${i}.log 2>&1 &
  done

wait

for i in `seq 0 $CORE_COUNT`
  do
    tail tests_${i}.log
  done
