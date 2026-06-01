## --------------------------------------------- ##
#           Model Validity: Random Forest
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks for the validity of random forest models
## by calculating R-squared for observed vs predicted for every station.
## This tests how well does the model capture temporal patterns in the last year of the training set
## across different timeframe lengths.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(caret)

# Create new folders to store results
dir.create(path = file.path("model_validity_results_timeframe"), showWarnings = F)
dir.create(path = file.path("model_validity_results_timeframe", "random_forest"), showWarnings = F)

DBSAL_orig <- readr::read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) 

DBSAL <- DBSAL_orig %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  dplyr::filter(!if_any(everything(), is.infinite)) 

Fmask_lookup <- readr::read_csv("HLSL30-020-Fmask-lookup.csv") %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL %>%
  # Filter out cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) 

# Create a version of the data with filled out band values
# (Landsat and Sentinel values taken from the nearest previous day)
DBSAL_filled <- DBSAL_orig %>%
  dplyr::arrange(formatted_date) %>%
  dplyr::group_by(station) %>%
  tidyr::fill(starts_with("B0"), .direction = "down") %>%
  tidyr::fill(Fmask, .direction = "down") %>%
  tidyr::fill(instrument, .direction = "down") %>%
  dplyr::ungroup() %>%
  # Calculate indices
  dplyr::mutate(NDVI = (B05 - B04) / (B05 + B04),
                SI = (B03*B04)^0.5,
                NLI = (B05^2 - B04)/(B05^2 + B04),
                SRSI = ((NDVI - 1)^2 + SI^2)^0.5,
                S7 = (B06 - B07)/(B06 + B07),
                CRSI = ((B05*B04-B03*B02)/((B05*B04+B03*B02)))^0.5,
                NDSI = (B05 - B06) / (B05 + B06)) 

## --------------------------------------------- ##
#                 Modeling -----
## --------------------------------------------- ##

set.seed(77)

calc_timeframe <- function(start_year, interval, training_data, filled_data){
  
  # Subset to interval
  DBSAL_some_years <- training_data %>%
    dplyr::filter(formatted_date >= as.Date(paste0(start_year,"-01-01")) & formatted_date < as.Date(paste0(start_year+interval,"-01-01")))
  
  # Split data set into training and test sets
  sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
  train_some_years <- DBSAL_some_years[sample1,] 
  test_some_years <- DBSAL_some_years[-sample1,]
  
  # Random forest modeling
  # Picked variables based on VSURF results
  # Chosen variables: distCoast+slope+B03+srad+CRSI+SI+NDSI+B06+NLI+tavg
  rf_model <- randomForest(salinity ~ distCoast+slope+B03+srad+CRSI+SI+NDSI+B06+NLI+tavg,
                           data = train_some_years,
                           importance = TRUE)
  
  # Sanity check
  my_pred <- predict(rf_model, newdata = test_some_years)
  mse <- mean((my_pred - test_some_years$salinity)^2)
  r_squared <- 1 - sum((test_some_years$salinity - my_pred)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)
  
  # Predict for the last year in the interval
  DBSAL_end_year <- filled_data %>%
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+interval-1,"-01-01")) & formatted_date < as.Date(paste0(start_year+interval,"-01-01"))) %>% 
    # Remove missing values
    na.omit() %>%
    # Remove infinite values
    dplyr::filter(!if_any(everything(), is.infinite)) 
  
  # Calculate coefficient of determination for all stations altogether
  my_pred2 <- predict(rf_model, newdata = DBSAL_end_year)
  mse2 <- mean((my_pred2 - DBSAL_end_year$salinity)^2)
  r_squared2 <- 1 - sum((DBSAL_end_year$salinity - my_pred2)^2) / sum((DBSAL_end_year$salinity - mean(DBSAL_end_year$salinity))^2)
  
  # Add a new column for predicted values
  DBSAL_end_year_v2 <- DBSAL_end_year %>%
    dplyr::mutate(pred = my_pred2)
  
  # Create lists to store our results
  r_sq1_list <- list()
  r_sq2_list <- list()
  stations <- list()
  n_rows_list <- list()
  
  # For each station...
  for (i in seq_along(unique(DBSAL_end_year_v2$station))){
    
    # Save station name
    stations[[i]] <- unique(DBSAL_end_year_v2$station)[i]
    
    message(paste("station:", unique(DBSAL_end_year_v2$station)[i]))
    
    station_sub <- DBSAL_end_year_v2 %>%
      # Filter to one station
      dplyr::filter(station == unique(DBSAL_end_year_v2$station)[i]) 
    
    # Save number of rows for that station
    n_rows_list[[i]] <- nrow(station_sub) 
    
    # Grab Pearson's correlation squared
    r_sq1 <- caret::postResample(obs = station_sub$salinity, pred = station_sub$pred)
    r_sq1_list[[i]] <- r_sq1[2]
    message(paste("r_sq1:", r_sq1[2]))
    
    # Calculate coefficient of determination
    r_sq2 <- 1 - sum((station_sub$salinity - station_sub$pred)^2) / sum((station_sub$salinity - mean(station_sub$salinity))^2)
    r_sq2_list[[i]] <- r_sq2
    message(paste("r_sq2:", r_sq2))
    
  }
  
  # Save results in data frame
  results <- data.frame(interval_length = interval,
                        training_years = paste0(start_year, "-", start_year+interval-1),
                        pred_year = paste0(start_year+interval-1),
                        variables = paste(rownames(rf_model$importance), collapse = ", "),
                        station_name = unlist(stations),
                        n_rows = unlist(n_rows_list),
                        pearsons_sq = round(unlist(r_sq1_list), digits = 4),
                        coeff_det = round(unlist(r_sq2_list), digits = 4))
  
  # Add an additional row to save results for all stations altogether
  results <- results %>%
    dplyr::add_row(interval_length = interval,
                   training_years = paste0(start_year, "_", start_year+interval-1),
                   pred_year = paste0(start_year+interval-1),
                   variables = paste(rownames(rf_model$importance), collapse = ", "),
                   station_name = "all",
                   n_rows = nrow(DBSAL_end_year_v2),
                   pearsons_sq = caret::postResample(obs = DBSAL_end_year_v2$salinity, pred = DBSAL_end_year_v2$pred)[2],
                   coeff_det = r_squared2)
  
  return(results)
  
}

# Set the range of years we're interested in
# For example, setting my_years <- 2013:2025 will cover every timeframe
# from 2013-2025, 2014-2025, 2015-2025, ..., 2025-2025

# Therefore, the best way to cover every possible timeframe is to run like so:
# my_years <- 2013:2025
# my_years <- 2013:2024
# my_years <- 2013:2023
# ...
# my_years <- 2013:2014
# my_years <- 2013

my_years <- 2013:2025
# Set target year to be the last year in the timeframe
my_target <- my_years[length(my_years)]

# For every year in the timeframe range...
for (i in my_years){
  message("on: ", i)
  
  # Set it as the start year
  my_start_year <- i
  # Set the target as the last year in the timeframe
  target <- my_target
  # Calculate the interval length from start year to last year
  my_interval <- target+1-my_start_year
  
  # Run modelling function for this specified timeframe
  df <- calc_timeframe(start_year = my_start_year,
                       interval = my_interval,
                       training_data = DBSAL_v2,
                       filled_data = DBSAL_filled)
  
  # Export to CSV
  readr::write_csv(df, file.path("model_validity_results_timeframe",
                                 "random_forest",
                                 paste0("obs_vs_pred_", my_start_year, "_", my_start_year+my_interval-1, ".csv")))
  
}

## --------------------------------------------- ##
#                 Harmonizing -----
## --------------------------------------------- ##

path <- file.path("model_validity_results_timeframe", "random_forest")

# List files 
files_to_harmonize <- list.files(path, pattern = "obs_vs_pred_", full.names = T)

files_to_harmonize

results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) %>%
  # Drop Pearson's correlation squared
  dplyr::select(-pearsons_sq) %>%
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