#!/bin/bash

source /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/setup.sh

cd /nfs/dust/ilc/user/bliewert/fullflow_v3/
rm -rf comparison
mkdir comparison
cd ./comparison

mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/root/prod

#filename=${1}
#bname=`basename -s .slcio $filename`

Marlin "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/scripts/comparison_prod.xml"
# --global.LCIOInputFiles="${filename}"

cd ..
