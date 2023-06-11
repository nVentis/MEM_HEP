#!/bin/bash
# source /cvmfs/ilc.desy.de/sw/x86_64_gcc82_centos7/v02-02-02/init_ilcsoft.sh &&
source /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/setup.sh

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

Marlin "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/jobs/2023-06-08_a/marlin_file.xml" --global.LCIOInputFiles="${filename}" --constant.OutputBaseName="${bname}"
mv *.root ../
mv *.slcio ../

cd ..
rm -rf job${2}
