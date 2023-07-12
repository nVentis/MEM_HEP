#!/bin/bash

#source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
#source /afs/desy.de/project/ilcsoft/sw/x86_64_gcc82_centos7/v02-02-03/init_ilcsoft.sh
source /cvmfs/ilc.desy.de/key4hep/setup.sh

echo "Using current directory as relative path for libraries, which is ${PWD}"

# Modules of ZHH
export ILD_ANASOFT_ZHH=/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/ZHH
export MEM_HEP=/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP

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
export MARLIN_DLL=$MARLIN_DLL:$MEM_HEP/source/ZHHPostRecoMEProcessor/lib/libZHHPostRecoMEProcessor.so
export MARLIN_DLL=$MARLIN_DLL:$MEM_HEP/source/CompareMEProcessor/lib/libCompareMEProcessor.so
export MARLIN_DLL=$MARLIN_DLL:$MEM_HEP/source/SaveFinalKinematics/lib/libSaveFinalKinematics.so

# Other modules
export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/b/bliewert/public/yradkhorrami/SLDecayCorrection/build/lib/libSLDecayCorrection.so

export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim:$CMAKE_PREFIX_PATH
export LD_LIBRARY=$LCIO/lib64:$LD_LIBRARY_PATH

# CPATH header files; marlin and streamlog
export CPATH=$CPATH:/cvmfs/ilc.desy.de/key4hep/spackages/marlin/1.19/x86_64-centos7-gcc11.2.0-opt/aj5n37vac5zzw4eil7raxolernp6o4vi/include
export CPATH=$CPATH:/cvmfs/ilc.desy.de/key4hep/spackages/ilcutil/1.7/x86_64-centos7-gcc11.2.0-opt/nnyawmpqncnl6tu3nfhufy4seyzdlfj2/include
export CPATH=$CPATH:/cvmfs/ilc.desy.de/key4hep/spackages/root/6.26.06/x86_64-centos7-gcc11.2.0-opt/dctcyvzmo7xg4dehiooyfl24oevtaids/include

# MadGraph and madjax
# export PATH=/nfs/dust/ilc/user/bliewert/MG5_aMC_v3_5_0/bin:$PATH
# export PATH=/nfs/dust/ilc/user/bliewert/MG5_aMC_v2_9_15/bin:$PATH
export PATH=/nfs/dust/ilc/user/bliewert/MG5_aMC_v2_8_1/bin:$PATH
export PATH=/afs/desy.de/user/b/bliewert/.local/bin:$PATH

# This is incompatible with py311 environment!