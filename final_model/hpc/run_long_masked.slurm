#!/bin/bash
#SBATCH --job-name=enp_long
#SBATCH --partition=day
#SBATCH --array=0-74
#SBATCH --cpus-per-task=2
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=logs/long_%A_%a.out
module load miniconda; conda activate sal30
python long_masked_hpc.py $SLURM_ARRAY_TASK_ID $HOME/enp_salinity
