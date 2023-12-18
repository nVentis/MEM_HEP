#!/bin/bash

source /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/setup.sh

filename=${1}
bname=`basename -s .slcio $filename`

rm -rf /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/${bname}
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/${bname}/root/prod
cd /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/${bname}

Marlin "/afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/scripts/comparison.xml" --global.LCIOInputFiles="${filename}" --constant.Verbosity="MESSAGE" --constant.RunMode="prod" --constant.RunAll="false"