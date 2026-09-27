#!/bin/bash
#SBATCH --job-name=enp_sal30
#SBATCH --array=0-155
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=06:00:00
#SBATCH --output=logs/sal30_%A_%a.out
# CHANGE AS NEEDED: environment and paths
module load miniconda 2>/dev/null; conda activate enp_salinity 2>/dev/null
IN=${IN:-$HOME/project/enp_salinity/hpc_inputs}          # prepared inputs (see README_HPC.md)
OUT=${OUT:-/path/to/MaloneLab/Research/ENP_Landsat_Salinity_Model/best_two_step/daily30m}
export NTHREAD=$SLURM_CPUS_PER_TASK TILE_ROWS=64
mkdir -p logs
python daily30_hpc.py $SLURM_ARRAY_TASK_ID $IN $OUT
