#!/bin/bash

while read d
do
if [ $d != abalone ]
    then
    sed -e "s/abalone/$d/g" run_abalone_2017_05_27.m > run_${d}_2017_05_27.m
fi
done < ~/Benchmarks/Data/uci/Names.txt