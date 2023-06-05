#!/bin/bash

#source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
#source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
source /cvmfs/ilc.desy.de/key4hep/setup.sh

echo "Using current directory as relative path for libraries, which is ${PWD}"

export MARLIN_DLL=$MARLIN_DLL:$PWD/source/ZHHPostRecoMEProcessor/lib/libZHHPostRecoMEProcessor.so:$PWD/source/CompareTrueMEProcessor/lib/libCompareTrueMEProcessor.so

export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim:$CMAKE_PREFIX_PATH
export LD_LIBRARY=$LCIO/lib64:$LD_LIBRARY_PATH
