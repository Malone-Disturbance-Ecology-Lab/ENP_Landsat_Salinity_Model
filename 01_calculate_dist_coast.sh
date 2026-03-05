#!/bin/bash

### Run on day partition ###
#SBATCH -p day

### Request 2 nodes ###
#SBATCH --nodes=2

### Request 2 cores per node ###
#SBATCH --cpus-per-task=2

### Request 1 day of run time ###
#SBATCH -t 1-00:00:00

### Receive email notifications when this job starts/completes ###
#SBATCH --mail-type=all

### Remove any loaded modules ###
module reset 

### Load R module ###
module load R/4.4.1-foss-2022b

### Execute the R script ###
Rscript 01_calculate_dist_coast.R