rm -rf /nfs/dust/ilc/user/bliewert/fullflow_v3/zhh
rm -rf /nfs/dust/ilc/user/bliewert/fullflow_v3/zzh
#rm -rf /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison

mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zhh/slcio
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zhh/root
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zhh/log

mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zzh/slcio
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zzh/root
mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/zzh/log

#mkdir -p /nfs/dust/ilc/user/bliewert/fullflow_v3/comparison

condor_submit /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/jobs/zhh.sub
condor_submit /afs/desy.de/user/b/bliewert/public/MarlinWorkdirs/MEM_HEP/jobs/zzh.sub