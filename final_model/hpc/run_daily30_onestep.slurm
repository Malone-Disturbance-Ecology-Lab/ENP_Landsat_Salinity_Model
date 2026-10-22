#!/bin/bash
# Daily 30 m maps from the final one-step (no-spectra, sea level) model, one array task per month.
# Index 0 = 2000-04 ... 308 = 2025-12 (309 months).
#SBATCH --job-name=enp_os30
#SBATCH --partition=day
#SBATCH --array=0-308
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/os30_%A_%a.out
module load miniconda; conda activate sal30
IN=$HOME/enp_salinity/hpc_inputs
OUT=$HOME/enp_salinity/daily30m_onestep_final
export EDEN_DIR=$HOME/enp_salinity/data/EDEN_GRIDS_Q NTHREAD=$SLURM_CPUS_PER_TASK
mkdir -p logs
python daily30_onestep_hpc.py $SLURM_ARRAY_TASK_ID $IN $OUT
