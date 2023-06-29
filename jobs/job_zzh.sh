#!/bin/bash
# source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-02/init_ilcsoft.sh &&
source /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/setup.sh
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zzh
cd /nfs/dust/ilc/user/bliewert/fullflow_v3/zzh

rm -rf job${2}
mkdir job${2}
cd ./job${2}

#ls -la
#cp ${1} .
#sleep 2
#ls -la
#filename=`ls *.slcio`
filename=${1}
bname=`basename -s .slcio $filename`

Marlin "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/scripts/e2e2hh.xml" --global.LCIOInputFiles="${filename}" --constant.OutputDirectory="." --constant.OutputBaseName="${bname}" --global.MaxRecordNumber=0
mv *.root ../root
mv *.slcio ../slcio

cd ..
rm -rf job${2}
