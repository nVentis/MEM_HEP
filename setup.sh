#!/bin/bash

source /cvmfs/ilc.desy.de/key4hep/setup.sh
#source /cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/key4hep-stack/2023-05-24/x86_64-centos7-gcc12.3.0-opt/7emhu/setup.sh
#export JUPYTER_PATH=${ROOTSYS}/etc/notebook:${JUPYTER_PATH}

export Python_INCLUDE_DIR=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/python/3.10.10/x86_64-centos7-gcc12.3.0-opt/3hqba/include/python3.10
export Python_LIBRARY=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/python/3.10.10/x86_64-centos7-gcc12.3.0-opt/3hqba/lib/libpython3.10.so
export Python_EXECUTABLE=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/python/3.10.10/x86_64-centos7-gcc12.3.0-opt/3hqba/bin/python
export Python_FIND_VIRTUALENV=STANDARD
export Python_ROOT_DIR=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/python/3.10.10/x86_64-centos7-gcc12.3.0-opt/3hqba
export Python_FIND_STRATEGY=LOCATION
export PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/python/3.10.10/x86_64-centos7-gcc12.3.0-opt/3hqba/bin/python:$PATH

# Modules of ZHH
export ILD_ANASOFT_ZHH=/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/ZHH
export MEM_HEP=/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP

export MARLIN_DLL=/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim/build/lib/libPhyssim.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/AddNeutralPFOCovMat/lib/libAddNeutralPFOCovMat.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/LeptonErrorAnalysis/lib/libLeptonErrorAnalysis.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/CheatedMCOverlayRemoval/lib/libCheatedMCOverlayRemoval.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/LeptonPairing/lib/libLeptonPairing.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/HdecayMode/lib/libHdecayMode.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/PreSelection/lib/libPreSelection.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/JetErrorAnalysis/lib/libJetErrorAnalysis.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/ZHHKinfitProcessors/lib/libZHHKinfitProcessors.so:$MARLIN_DLL
export MARLIN_DLL=$ILD_ANASOFT_ZHH/source/Misclustering/lib/libMisclustering.so:$MARLIN_DLL

# Modules of MEM_HEP
#export MARLIN_DLL=$MARLIN_DLL:$MEM_HEP/source/ZHHPostRecoMEProcessor/lib/libZHHPostRecoMEProcessor.so
export MARLIN_DLL=$MEM_HEP/source/CompareMEProcessor/lib/libCompareMEProcessor.so:$MARLIN_DLL
export MARLIN_DLL=$MEM_HEP/source/SaveFinalKinematics/lib/libSaveFinalKinematics.so:$MARLIN_DLL
export MARLIN_DLL=$MEM_HEP/source/JetConvProcessor/build/lib/libJetConvProcessor.so:$MARLIN_DLL
#export MARLIN_DLL=$MARLIN_DLL:$MEM_HEP/source/Use_TrueJet/build/lib/libUse_TrueJet.so

# Other modules
#export MARLIN_DLL=$MARLIN_DLL:/afs/desy.de/user/b/bliewert/public/yradkhorrami/SLDecayCorrection/build/lib/libSLDecayCorrection.so

export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/ILCSoft/Physsim:$CMAKE_PREFIX_PATH
export LD_LIBRARY_PATH=$LCIO/lib64:$LD_LIBRARY_PATH

#export GSL_ROOT_DIR=/afs/desy.de/user/b/bliewert/public/DevRepositories/gsl

# CPATH header files; marlin and streamlog
#export CPATH=/afs/desy.de/user/b/bliewert/public/ILCSoft/MarlinUtil/source/include:$CPATH
export CPATH=/cvmfs/ilc.desy.de/key4hep/releases/089d775cf2/marlin/1.19/x86_64-centos7-gcc12.3.0-opt/d6jkp/include:$CPATH
export CPATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/ilcutil/1.7/x86_64-centos7-gcc12.3.0-opt/b3vqf/include:$CPATH
export CPATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/root/6.28.04/x86_64-centos7-gcc12.3.0-opt/owni5/include:$CPATH
export CPATH=/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/analysis/cffi/include:$CPATH

# MadGraph and madjax
export PATH=/nfs/dust/ilc/user/bliewert/MG5_aMC_v3_5_1/bin:$PATH
# export PATH=/nfs/dust/ilc/user/bliewert/MG5_aMC_v2_9_15/bin:$PATH
#export PATH=/nfs/dust/ilc/user/bliewert/MG5_aMC_v2_8_1/bin:$PATH
#export PATH=/afs/desy.de/user/b/bliewert/.local/bin:$PATH

# This is incompatible with py311 environment!

# PyTorch
export PYTHONPATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages:$PYTHONPATH
export LD_LIBRARY_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages/torch/lib:$LD_LIBRARY_PATH
export CMAKE_PREFIX_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages/torch/share/cmake:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/py-torch/1.13.1/x86_64-centos7-gcc12.3.0-opt/i7cg4/lib/python3.10/site-packages/torch/include/torch/csrc/api:$CMAKE_PREFIX_PATH

# PyTorch Geometric dependencies
export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pytorch_scatter/share/cmake:$CMAKE_PREFIX_PATH
export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pytorch_sparse/share/cmake:$CMAKE_PREFIX_PATH

# pybind11
export CMAKE_PREFIX_PATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pybind11/share/cmake:$CMAKE_PREFIX_PATH
export CPATH=/afs/desy.de/user/b/bliewert/public/DevLocal/pybind11/include:$CPATH

# lcfiplus
#export LD_LIBRARY_PATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/lcfiplus/0.10.1/x86_64-centos7-gcc12.3.0-opt/7tnbl/lib:$LD_LIBRARY_PATH
#export CPATH=/cvmfs/ilc.desy.de/key4hep/releases/2023-05-23/lcfiplus/0.10.1/x86_64-centos7-gcc12.3.0-opt/7tnbl/include:$CPATH


