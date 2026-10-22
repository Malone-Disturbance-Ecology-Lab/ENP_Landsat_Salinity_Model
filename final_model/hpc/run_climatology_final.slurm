#!/bin/bash
#SBATCH --job-name=enp_clim
#SBATCH --partition=day
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=logs/clim_%j.out
module load miniconda; conda activate sal30
python climatology_final.py $HOME/enp_salinity
