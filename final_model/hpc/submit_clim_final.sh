#!/bin/bash
cd ~/enp_salinity/hpc
a=$(sbatch --parsable run_wetdays_final.slurm) && echo wetdays $a
c=$(sbatch --parsable --dependency=afterok:$a run_climatology_final.slurm) && echo climatology $c
