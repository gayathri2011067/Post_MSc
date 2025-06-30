#!/bin/bash
mkdir -p logs
times=(0 1 2 3 4 5 6 7)
i=0
for param in 0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0; do
    level=$((700 + i))
    echo "bash run.sh $param $level > logs/run_${param}_${level}.out 2>&1" | at now + ${times[$i]} minute
    ((i++))
done
