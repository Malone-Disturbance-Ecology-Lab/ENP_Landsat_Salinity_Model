#!/bin/bash
#SBATCH --job-name=enp_longr
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=logs/longreduce_%j.out
module load miniconda; conda activate sal30
python long_masked_reduce.py $HOME/enp_salinity
