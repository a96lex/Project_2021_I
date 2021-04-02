#!/bin/bash
echo -n  "Min number of processors (power): "
read min_np
echo -n "Max number of processors (power): "
read max_np

for (( i=min_np; i<=max_np; ++i))
do
    procs=$((2**i))
    echo $procs
    mkdir -p results/speedup/sim_${procs}
    sed "s/NUMPROCS/$procs/" sim_template.sub > ./results/speedup/sim_${procs}/script.sub
    qsub ./results/speedup/sim_${procs}/script.sub
done
exit 0
