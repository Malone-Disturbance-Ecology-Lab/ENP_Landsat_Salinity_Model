#!/bin/bash
#SBATCH --job-name=enp_sal30
#SBATCH --partition=day
#SBATCH --array=0-155
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sal30_%A_%a.out
# CHANGE AS NEEDED: environment and paths
module load miniconda; conda activate sal30
IN=${IN:-$HOME/enp_salinity/hpc_inputs}          # prepared inputs (see README_HPC.md)
OUT=${OUT:-$HOME/enp_salinity/daily30m}
export NTHREAD=$SLURM_CPUS_PER_TASK TILE_ROWS=64
mkdir -p logs
python daily30_hpc.py $SLURM_ARRAY_TASK_ID $IN $OUT
