#!/bin/bash

compile_pkg () {
cd $1
mkdir -p build
cd build
cmake -DCMAKE_CXX_STANDARD=17 ..
make install
cd ../..
}

cd source
compile_pkg CompareMEProcessor
compile_pkg ZHHPostRecoMEProcessor
cd ..
