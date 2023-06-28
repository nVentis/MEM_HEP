#!/bin/bash

cd /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/TestPhyssimK4H/Physsim/build
cmake -DCMAKE_CXX_STANDARD=17 ..
make install

cd /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/source/CompareMEProcessor/build
cmake -DCMAKE_CXX_STANDARD=17 ..
make install
cd ../../../
rm compare_out_mcparticle.root
Marlin marlin_compare_true_me.xml
mv compare_out.root compare_out_mcparticle.root

# Assumes we're using the correct conda environment (py311)
# python 