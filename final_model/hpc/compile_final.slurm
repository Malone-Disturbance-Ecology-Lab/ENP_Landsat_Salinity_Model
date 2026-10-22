#!/bin/bash
#SBATCH --job-name=enp_os30_compile
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/os30_compile_%j.out
module load miniconda; conda activate sal30
mkdir -p logs
python compile_final.py $HOME/enp_salinity/hpc_inputs
