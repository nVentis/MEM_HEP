#!/bin/bash

#source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
#source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
source /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/TestPhyssimK4H/only_physsim.sh

echo "Using current directory as relative path for libraries, which is ${PWD}"

export MARLIN_DLL=$MARLIN_DLL:$PWD/source/AddNeutralPFOCovMat/lib/libAddNeutralPFOCovMat.so:$PWD/source/LeptonErrorAnalysis/lib/libLeptonErrorAnalysis.so:$PWD/source/CheatedMCOverlayRemoval/lib/libCheatedMCOverlayRemoval.so:$PWD/source/LeptonPairing/lib/libLeptonPairing.so:$PWD/source/HdecayMode/lib/libHdecayMode.so:$PWD/source/PreSelection/lib/libPreSelection.so:$PWD/source/JetErrorAnalysis/lib/libJetErrorAnalysis.so:$PWD/source/ZHHKinfitProcessors/lib/libZHHKinfitProcessors.so:$PWD/source/Misclustering/lib/libMisclustering.so

export LD_LIBRARY=$LCIO/lib64:$LD_LIBRARY_PATH
