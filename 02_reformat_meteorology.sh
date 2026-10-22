#!/bin/bash

### Run on day partition ###
#SBATCH -p day

### Request 1 node ###
#SBATCH --nodes=1

### Request 32 GB of mem/cpu ###
#SBATCH --mem-per-cpu=32G

### Request 1 day of run time ###
#SBATCH -t 1-00:00:00

### Receive email notifications when this job starts/completes ###
#SBATCH --mail-type=all

### Remove any loaded modules ###
module reset 

### Load R module ###
module load R/4.4.1-foss-2022b

### Execute the R script ###
Rscript 02_reformat_meteorology.R