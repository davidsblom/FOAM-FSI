#!/bin/bash

set -e
set -x

git submodule update --init --recursive

# foam-extend
cd /home
tar xf foam-extend-3.2.tar.gz && rm foam-extend-3.2.tar.gz
source foam-extend-3.2/etc/prefs.sh && source foam-extend-3.2/etc/bashrc

# boost
cd /pipeline/source/src/thirdParty
ln -s /home/boost_1_61_0
ln -s /home/boost_1_61_0 boost
cp boost/stage/lib/* $FOAM_LIBBIN/

# petsc
ln -s /home/petsc-3.6.4
ln -s /home/petsc-3.6.4 petsc
cp petsc/x86_64/lib/lib* ${FOAM_LIBBIN}/

# deal.II
ln -s /home/dealii-8.4.1
ln -s /home/dealii-8.4.1 dealii
cp dealii/bin/lib/libdeal* ${FOAM_LIBBIN}/

# googletest
ln -s /home/googletest-release-1.8.0
ln -s /home/googletest-release-1.8.0 gtest
cp gtest/install/lib/lib* ${FOAM_LIBBIN}/

# yaml-cpp
ln -s /home/yaml-cpp-release-0.5.3
ln -s /home/yaml-cpp-release-0.5.3 yaml-cpp
cp yaml-cpp/install/lib/libyaml-cpp* ${FOAM_LIBBIN}/
cd ../../

./Allwmake.firstInstall
