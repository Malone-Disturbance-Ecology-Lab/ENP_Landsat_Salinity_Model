#!/bin/bash
#SBATCH --job-name=enp_m250
#SBATCH --partition=day
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/m250_%j.out
module load miniconda; conda activate sal30
python reduce_monthly250.py $HOME/enp_salinity
