#!/bin/bash
echo -n  "Min number of processors: "
read min_np
echo -n "Max number of processors: "
read max_np

for (( procs=min_np; procs<=max_np; ++procs))
do
    echo $procs
    mkdir -p results/speedup/sim_${procs}
    sed "s/NUMPROCS/$procs/" sim_template.sub > ./results/speedup/sim_${procs}/script.sub
    qsub ./results/speedup/sim_${procs}/script.sub
done
exit 0
