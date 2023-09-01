#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --tasks-per-node=8
#SBATCH --mem=16gb
#SBATCH --time=3:00:00
#SBATCH --partition=medium
#SBATCH --account=proj-fs0002

# params
bampath=$1
prefix=$2
mosdepth_bin_path=$3

threads=8

${mosdepth_bin_path} --no-per-base --fast-mode --threads ${threads} ${prefix} ${bampath}
