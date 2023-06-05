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
compile_pkg AddNeutralPFOCovMat
compile_pkg CheatedMCOverlayRemoval
compile_pkg HdecayMode
compile_pkg JetErrorAnalysis
compile_pkg LeptonErrorAnalysis
compile_pkg LeptonPairing
compile_pkg PreSelection
compile_pkg ZHHKinfitProcessors
compile_pkg ZHHPostRecoMEProcessor
cd ..
