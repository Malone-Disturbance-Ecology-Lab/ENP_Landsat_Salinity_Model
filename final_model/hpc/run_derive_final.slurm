#!/bin/bash
#SBATCH --job-name=enp_derive
#SBATCH --partition=day
#SBATCH --array=0-308
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/derive_%A_%a.out
module load miniconda; conda activate sal30
python derive_final.py $SLURM_ARRAY_TASK_ID $HOME/enp_salinity
