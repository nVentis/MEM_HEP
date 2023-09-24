#!/bin/bash

#export PATH="/nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311/bin/"
source /afs/desy.de/user/b/bliewert/.zshrc
conda activate /nfs/dust/ilc/user/bliewert/.mambaforge/envs/py311

cd /nfs/dust/ilc/user/bliewert/mem_integrate/results

event=${1}

mkdir -p event_$event
cd event_$event

python /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/cli mem integrate $event ./result.txt