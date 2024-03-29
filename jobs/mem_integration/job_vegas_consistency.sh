#!/bin/bash

#export PATH="/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/bin/"
source /afs/desy.de/user/b/bliewert/.zshrc
conda activate /nfs/dust/ilc/user/bliewert/miniconda3/envs/py37

mkdir -p /nfs/dust/ilc/user/bliewert/vegas_consistency/results
mkdir -p /nfs/dust/ilc/user/bliewert/vegas_consistency/log

cd /nfs/dust/ilc/user/bliewert/vegas_consistency/results

event=${1}
job=${2}

mkdir -p event_$event
cd event_$event

python /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/cli mem integrate $event ./result.txt --me_type=1 --with_perms=0 --job=${job}