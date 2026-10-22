#!/bin/bash
cd ~/enp_salinity/hpc
a=$(sbatch --parsable run_derive_final.slurm) && echo derive $a
r=$(sbatch --parsable --dependency=afterok:$a run_reduce_monthly250.slurm) && echo monthly250 $r
