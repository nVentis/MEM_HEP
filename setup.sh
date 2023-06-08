#!/bin/bash

#source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
#source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
source /cvmfs/ilc.desy.de/key4hep/setup.sh

echo "Using current directory as relative path for libraries, which is ${PWD}"

# Modules of ZHH
export ILD_ANASOFT_ZHH=/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/ZHH

export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/build/lib/libPhyssim.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/AddNeutralPFOCovMat/lib/libAddNeutralPFOCovMat.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/LeptonErrorAnalysis/lib/libLeptonErrorAnalysis.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/CheatedMCOverlayRemoval/lib/libCheatedMCOverlayRemoval.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/LeptonPairing/lib/libLeptonPairing.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/HdecayMode/lib/libHdecayMode.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/PreSelection/lib/libPreSelection.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/JetErrorAnalysis/lib/libJetErrorAnalysis.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/ZHHKinfitProcessors/lib/libZHHKinfitProcessors.so
export MARLIN_DLL=$MARLIN_DLL:$ILD_ANASOFT_ZHH/source/Misclustering/lib/libMisclustering.so

# Modules of MEM_HEP
export MARLIN_DLL=$MARLIN_DLL:$PWD/source/ZHHPostRecoMEProcessor/lib/libZHHPostRecoMEProcessor.so
export MARLIN_DLL=$MARLIN_DLL:$PWD/source/CompareTrueMEProcessor/lib/libCompareTrueMEProcessor.so

export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim:$CMAKE_PREFIX_PATH
export LD_LIBRARY=$LCIO/lib64:$LD_LIBRARY_PATH
