#!/bin/bash
# script to clean up output given by c++ solver

echo "loading modules..."
module load cmake
module load intel/18.0

echo "removing old cmake data..."
rm -R CMakeFiles
rm CMakeCache.txt
echo "making make with cmake"
cmake -DDEAL_II_DIR=~/testbuild/dealii .
#cmake -D CMAKE_CXX_COMPILER=/afs/crc.nd.edu/x86_64_linux/intel/15.0/bin/icc -DDEAL_II_DIR=~/dealii .

