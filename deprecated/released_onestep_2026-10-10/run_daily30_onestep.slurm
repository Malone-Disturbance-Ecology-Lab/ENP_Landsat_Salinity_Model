#!/bin/bash
# Daily 30 m one-step (no-spectra) hindcast maps, one array task per month.
# Index 0 = 2000-04 ... 152 = 2012-12 (153 months). Indices 153-308 = 2013-01..2025-12 if a one-step daily record is also wanted then.
#SBATCH --job-name=enp_os30
#SBATCH --partition=day
#SBATCH --array=0-152
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=04:00:00
#SBATCH --output=logs/os30_%A_%a.out
module load miniconda; conda activate sal30
IN=${IN:-$HOME/enp_salinity/hpc_inputs}
OUT=${OUT:-$HOME/enp_salinity/daily30m_onestep}
export EDEN_DIR=${EDEN_DIR:-$HOME/enp_salinity/data/EDEN_GRIDS_Q} NTHREAD=$SLURM_CPUS_PER_TASK
mkdir -p logs
python daily30_onestep_hpc.py $SLURM_ARRAY_TASK_ID $IN $OUT
