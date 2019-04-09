#!/bin/sh

#  build various solvers

echo "===== Cleaning old bins ====="
test -d bin || mkdir -p bin
rm -f bin/akmaxsat*
rm -f bin/CCLS*

echo "===== Building akmaxsat ====="
cd src/akmaxsat_1.1
make cleanup
make
mv akmaxsat ../../bin/
cd ../../

echo "===== Building CCLS - an incomplete solver (2014 version) ====="
cd src/CCLS2014
make cleanup
make
mv CCLS2014 ../../bin/
cd ../../

echo "===== Building CCLS_to_akmaxsat ====="
cd src/CCLS_to_akmaxsat
make cleanup
make
mv CCLS_to_akmaxsat ../../bin/
cd ../../
