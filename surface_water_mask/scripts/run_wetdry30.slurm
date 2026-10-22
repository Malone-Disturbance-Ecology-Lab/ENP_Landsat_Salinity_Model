#!/bin/bash
# Daily wet/dry mask + wet-day monthly salinity, one task per month (0 = 2000-04 ... 308 = 2025-12).
#SBATCH --job-name=enp_wd30
#SBATCH --partition=day
#SBATCH --array=0-308
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=01:00:00
#SBATCH --output=logs/wd30_%A_%a.out
module load miniconda; conda activate sal30
python wetdry30_hpc.py $SLURM_ARRAY_TASK_ID $HOME/enp_salinity/wetdry/in $HOME/enp_salinity
