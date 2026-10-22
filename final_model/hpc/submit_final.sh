#!/bin/bash
# Run on Bouchet from ~/enp_salinity/hpc: compile the model, then the 309 daily 30 m tasks once the compile succeeds.
cd ~/enp_salinity/hpc
j=$(sbatch --parsable compile_final.slurm) && echo compile $j
a=$(sbatch --parsable --dependency=afterok:$j run_daily30_onestep.slurm) && echo array $a
