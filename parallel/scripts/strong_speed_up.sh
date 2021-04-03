#!/bin/bash
#$ -S /bin/bash
#$ -q master.q
#$ -N test
#$ -cwd

echo -n  "Min number of processors: "
read min_np
echo -n "Max number of processors: "
read max_np
for (( i=min_np; i <= max_np; ++i ))
do
    start=`date +%s.%N`
    mpirun -np $i main.x input_template.txt 1>/dev/null 2>&1 || exit 1
    end=`date +%s.%N`
    runtime=$( echo "$end - $start" | bc -l )
    echo "processors: $i, time: $runtime seconds" 
done
