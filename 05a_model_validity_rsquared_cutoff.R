## --------------------------------------------- ##
#         Model Validity: R-squared cutoff
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks for the ideal R-squared cutoff value to be
## considered a "good" value in order to find out which 
## stations the model consistently performs well for.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)

## --------------------------------------------- ##
#                 Harmonizing -----
## --------------------------------------------- ##

path <- file.path("model_validity_results_timeframe_diff", "random_forest")

# List files 
files_to_harmonize <- list.files(path, pattern = "obs_vs_pred_", full.names = T)

files_to_harmonize

results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) %>%
  # Round numbers
  dplyr::mutate(coeff_det = round(coeff_det, digits = 4)) %>%
  # Remove infinite values
  dplyr::filter(!if_any(everything(), is.infinite)) 


## --------------------------------------------- ##
#   Exploration: Label R^2 > 0.5 as "good" -----
## --------------------------------------------- ##

results_0.5_cutoff <- results_harmonized %>%
  # Group by prediction year and station
  dplyr::group_by(pred_year, station_name) %>%
  # Get the average R^2
  dplyr::summarize(avg_coeff_det = mean(coeff_det)) %>%
  # Round numbers
  dplyr::mutate(avg_coeff_det = round(avg_coeff_det, digits = 4)) %>%
  # Denote "good" R^2 values as > 0.5
  # otherwise "bad" R^2 value as <= 0.5
  dplyr::mutate(does_well = case_when(
    avg_coeff_det > 0.5 ~ 1,
    avg_coeff_det <= 0.5 ~ 0,
    T ~ NA
  ))

results_0.5_cutoff_v2 <- results_0.5_cutoff %>%
  # Group by station
  dplyr::group_by(station_name) %>%
  # Find how many years did model predict well for each station
  # (Max number is 13 since there are 13 prediction years from 2013-2025)
  dplyr::summarize(count_does_well = sum(does_well))

# When R^2 > 0.5 is labelled as "good",
# SEVENPALM, ENPWP, TAYLORS3, ENPCW, FLAB44 are the consistent "good" stations
# (count_does_well > = 3)

## --------------------------------------------- ##
#     Exploration: Is > 0.5 too high? -----
## --------------------------------------------- ##

results_ranked <- results_harmonized %>%
  # Group by station
  dplyr::group_by(station_name) %>%
  # Order rows by descending R^2
  dplyr::arrange(desc(coeff_det), .by_group = T) %>%
  # Rank the R^2
  dplyr::mutate(rank = row_number()) %>%
  # Find the middle rank 
  dplyr::mutate(middle_rank = round(max(rank)/2)) %>%
  # Ungroup
  dplyr::ungroup() %>%
  # Filter to find the R^2 value at the middle rank
  dplyr::filter(rank == middle_rank)

results_ranked_v2 <- results_ranked %>%
  # Filter to only positive R^2 values
  dplyr::filter(coeff_det > 0 ) %>%
  dplyr::arrange(desc(coeff_det))

# Plot to check
ggplot() + 
  geom_density(aes(x = coeff_det), data = results_ranked_v2)

# Find the value at the peak of the density plot
dens <- density(results_ranked_v2$coeff_det)
peak_x <- dens$x[which.max(dens$y)]
# 0.33 

## --------------------------------------------- ##
#   Exploration: Label R^2 > 0.33 as "good" -----
## --------------------------------------------- ##

results_0.33_cutoff <- results_harmonized %>%
  # Group by prediction year and station
  dplyr::group_by(pred_year, station_name) %>%
  # Get the average R^2
  dplyr::summarize(avg_coeff_det = mean(coeff_det)) %>%
  # Round numbers
  dplyr::mutate(avg_coeff_det = round(avg_coeff_det, digits = 4)) %>%
  # Denote "good" R^2 values as > 0.33
  # otherwise "bad" R^2 value as <= 0.33
  dplyr::mutate(does_well = case_when(
    avg_coeff_det > 0.33 ~ 1,
    avg_coeff_det <= 0.33 ~ 0,
    T ~ NA
  )) %>%
  # Group by station
  dplyr::group_by(station_name) %>%
  # Find how many years did model predict well for each station
  # (Max number is 13 since there are 13 prediction years from 2013-2025)
  dplyr::summarize(count_does_well = sum(does_well))

# More stations getting labelled as consistently "good"
# 21 stations where count_does_well >= 3

write_csv(results_0.33_cutoff, "results_0.33_cutoff.csv")
