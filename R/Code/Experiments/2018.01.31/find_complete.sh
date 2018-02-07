#!/bin/bash

touch incomplete.txt

for i in $(seq 1 23); do
    if grep -q "Execution halted" slurm-21750412_${i}.out; then
	echo -n ","$i >> incomplete.txt;
    fi
done

for i in $(seq 25 106); do
    if grep -q "Execution halted" slurm-21750412_${i}.out; then
	echo -n ","$i >> incomplete.txt;
    fi
done

