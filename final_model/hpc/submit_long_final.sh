#!/bin/bash
# Run on Bouchet from ~/enp_salinity/hpc after the daily 30 m array (job id as $1): trends + events per tile, then the reduce step.
cd ~/enp_salinity/hpc
t=$(sbatch --parsable --dependency=afterok:$1 run_long_masked.slurm) && echo tiles $t
r=$(sbatch --parsable --dependency=afterok:$t run_long_reduce.slurm) && echo reduce $r
