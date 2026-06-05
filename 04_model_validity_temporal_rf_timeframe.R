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

# Create new folders to store results
dir.create(path = file.path("model_validity_results_timeframe"), showWarnings = F)
# Folder for first attempt
dir.create(path = file.path("model_validity_results_timeframe", "random_forest"), showWarnings = F)
# Folder for second attempt that included calculating difference betw actual values and predictions
dir.create(path = file.path("model_validity_results_timeframe_diff", "random_forest"), showWarnings = F)

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
  #sanity_check_pred <- predict(rf_model, newdata = test_some_years)
  #sanity_check_mse <- mean((sanity_check_pred - test_some_years$salinity)^2)
  #sanity_check_r_squared <- 1 - sum((test_some_years$salinity - sanity_check_pred)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)
  
  # Predict for the last year in the interval
  DBSAL_end_year <- filled_data %>%
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+interval-1,"-01-01")) & formatted_date < as.Date(paste0(start_year+interval,"-01-01"))) %>% 
    # Remove missing values
    na.omit() %>%
    # Remove infinite values
    dplyr::filter(!if_any(everything(), is.infinite)) 
  
  # Calculate coefficient of determination for all stations altogether
  my_pred <- predict(rf_model, newdata = DBSAL_end_year)
  #my_mse <- mean((my_pred - DBSAL_end_year$salinity)^2)
  r_squared <- 1 - sum((DBSAL_end_year$salinity - my_pred)^2) / sum((DBSAL_end_year$salinity - mean(DBSAL_end_year$salinity))^2)
  
  DBSAL_end_year_v2 <- DBSAL_end_year %>%
    # Add a new column for predicted values
    # Add a new column for the difference between actual value and prediction
    dplyr::mutate(pred = my_pred,
                  diff = salinity - my_pred)
  
  # Create lists to store our results
  r_sq_list <- list()
  stations <- list()
  n_rows_list <- list()
  avg_diff_list <- list()
  min_sal_list <- list()
  mean_sal_list <- list()
  max_sal_list <- list()
  min_pred_list <- list()
  mean_pred_list <- list()
  max_pred_list <- list()
  
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
    
    # Find average difference
    avg_diff_list[[i]] <- mean(station_sub$diff)
    
    # Find min, mean, max salinity
    min_sal_list[[i]] <- min(station_sub$salinity)
    mean_sal_list[[i]] <- mean(station_sub$salinity)
    max_sal_list[[i]] <- max(station_sub$salinity)
    
    # Find min, mean, max prediction
    min_pred_list[[i]] <- min(station_sub$pred)
    mean_pred_list[[i]] <- mean(station_sub$pred)
    max_pred_list[[i]] <- max(station_sub$pred)    
    
    # Calculate coefficient of determination
    station_r_sq <- 1 - sum((station_sub$salinity - station_sub$pred)^2) / sum((station_sub$salinity - mean(station_sub$salinity))^2)
    r_sq_list[[i]] <- station_r_sq
    message(paste("station_r_sq:", station_r_sq))
    
  }
  
  # Save results in data frame
  results <- data.frame(interval_length = interval,
                        training_years = paste0(start_year, "-", start_year+interval-1),
                        pred_year = paste0(start_year+interval-1),
                        variables = paste(rownames(rf_model$importance), collapse = ", "),
                        station_name = unlist(stations),
                        n_rows = unlist(n_rows_list),
                        coeff_det = round(unlist(r_sq_list), digits = 4),
                        avg_diff = round(unlist(avg_diff_list), digits = 2),
                        min_sal = unlist(min_sal_list),
                        mean_sal = round(unlist(mean_sal_list), digits = 2),
                        max_sal = unlist(max_sal_list),
                        min_pred = round(unlist(min_pred_list), digits = 2),
                        mean_pred = round(unlist(mean_pred_list), digits = 2),
                        max_pred = round(unlist(max_pred_list), digits = 2))
  
  # Add an additional row to save results for all stations altogether
  results <- results %>%
    dplyr::add_row(interval_length = interval,
                   training_years = paste0(start_year, "_", start_year+interval-1),
                   pred_year = paste0(start_year+interval-1),
                   variables = paste(rownames(rf_model$importance), collapse = ", "),
                   station_name = "all",
                   n_rows = nrow(DBSAL_end_year_v2),
                   coeff_det = round(r_squared, digits = 4),
                   avg_diff = round(mean(DBSAL_end_year_v2$diff), digits = 2),
                   min_sal = min(DBSAL_end_year_v2$salinity),
                   mean_sal = round(mean(DBSAL_end_year_v2$salinity), digits = 2),
                   max_sal = max(DBSAL_end_year_v2$salinity),
                   min_pred = round(min(DBSAL_end_year_v2$pred), digits = 2),
                   mean_pred = round(mean(DBSAL_end_year_v2$pred), digits = 2),
                   max_pred = round(max(DBSAL_end_year_v2$pred), digits = 2))
  
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

my_years <- 2013
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
  readr::write_csv(df, file.path("model_validity_results_timeframe_diff",
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
  # Round numbers
  dplyr::mutate(coeff_det = round(coeff_det, digits = 4)) %>%
  # Remove infinite values
  dplyr::filter(!if_any(everything(), is.infinite)) 
