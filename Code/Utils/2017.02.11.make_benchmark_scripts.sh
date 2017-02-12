#!/bin/sh

while read l; do
    if [ "$l" != "abalone" ]; then
	sed "s/abalone/$l/g" ../Experiments/2017.02.11.benchmarks/run_abalone.m > ../Experiments/2017.02.11.benchmarks/run_$l.m
    fi
done < ../../Data/Benchmarks/Benchmark_names.txt