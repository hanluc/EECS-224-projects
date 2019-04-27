#!/bin/bash
#$ -N qsort
#$ -q class 
#$ -pe openmp 8-64
#$ -m beas

make quicksort
perf stat ./quicksort 1000
