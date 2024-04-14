#!/bin/sh
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --job-name=REMDcalc

python main.py