#!/bin/bash

source /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/setup.sh

rm -rf /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/root/prod
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison/root/prod
cd /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison

#filename=${1}
#bname=`basename -s .slcio $filename`

stf="comparison_prod.xml"

rm $stf && echo "$stf deleted" || echo "$stf not deleted; skipping"

# Read in comparison.xml and adapt for prod
steering_file=$(cat /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/scripts/comparison.xml)
files=$(cat /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/jobs/input_files/comparison_prod.txt)
echo "${steering_file/\$\{InputFiles\}/$files}" > $stf

[[ -f $stf ]] && echo "Created steering file $stf; starting Marlin" || echo "Could not create steering file $stf; aborting" && return 1

Marlin $stf --constant.Verbosity="DEBUG" --constant.RunMode="prod"