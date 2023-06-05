#!/bin/bash

echo "This is only for bootstrapping from scratch faster, ... not a build system. Use at your own risk."

sleep 1
echo -n "."
sleep 1
echo -n "."
sleep 1
echo -n "."
sleep 1

compile_pkg () {
cd $1
mkdir -p build
cd build
cmake -DCMAKE_CXX_STANDARD=17 ..
make install
cd ../..
}

cd source
compile_pkg CompareTrueMEProcessor
compile_pkg ZHHPostRecoMEProcessor
cd ..
