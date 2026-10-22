#!/bin/bash
#SBATCH --job-name=enp_wdf
#SBATCH --partition=day
#SBATCH --array=153-308
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --output=logs/wdf_%A_%a.out
module load miniconda; conda activate sal30
python wetdays_final.py $SLURM_ARRAY_TASK_ID $HOME/enp_salinity
