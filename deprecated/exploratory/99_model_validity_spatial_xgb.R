## --------------------------------------------- ##
#           Model Validity: XGBoost
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks for the validity of XGBoost models
## by calculating R-squared between each pair of stations for both observed and predicted data.
## This tests how well does the model capture spatial patterns.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(xgboost)
library(parallel)
library(caret)

# Create new folders to store results
dir.create(path = file.path("model_validity_results"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "xgboost"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "xgboost", "2016_2020"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "xgboost", "2017_2021"), showWarnings = F)

DBSAL <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

## --------------------------------------------------------------- ##
#    Observed Data: Getting R-Squared for Pairs of Stations -----
## --------------------------------------------------------------- ##

set.seed(77)
start_year <- 2017
my_interval <- 5

# Specify the next years to compare our future predictions to
next_years <- 5

# Get unique combinations of stations where order doesn't matter
grid <- expand.grid(unique(DBSAL$station), unique(DBSAL$station))
unique_pairs <- grid[!duplicated(t(apply(grid, 1, sort))), ] %>%
  # Remove rows with the same station twice
  dplyr::filter(Var1 != Var2) %>%
  # Rename columns
  dplyr::rename(station1 = Var1) %>%
  dplyr::rename(station2 = Var2) %>%
  # Create new column to store our R-squared values
  dplyr::mutate(pearsons_sq = NA) %>%
  # Convert columns to character type
  dplyr::mutate(station1 = as.character(station1),
                station2 = as.character(station2)) %>%
  # Create new data type column
  dplyr::mutate(type = "obs", .before = station1)

for (i in 1:nrow(unique_pairs)){
  station_sub1 <- DBSAL %>%
    # Filter to one station
    dplyr::filter(station == unique_pairs[i,]$station1) %>%
    # Filter to specified next years to compare our predictions to 
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval+next_years,"-01-01"))) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, salinity)
  
  station_sub2 <- DBSAL %>%
    # Filter to another station
    dplyr::filter(station == unique_pairs[i,]$station2) %>%
    # Filter to specified next years to compare our predictions to 
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval+next_years,"-01-01"))) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, salinity)
  
  if(nrow(station_sub1) == 0){
    unique_pairs[i,]$pearsons_sq <- NA
    next
  } else if (nrow(station_sub2) == 0){
    unique_pairs[i,]$pearsons_sq <- NA
    next
  } else {
    
    # Find dates that both stations have in common
    common_dates <- intersect(as.character(station_sub1$formatted_date), as.character(station_sub2$formatted_date))
    
    if (length(common_dates) == 0){
      unique_pairs[i,]$pearsons_sq <- NA
      next
    }
    
    station_sub1_v2 <- station_sub1 %>%
      # Filter to rows that have common dates
      dplyr::filter(formatted_date %in% common_dates) %>%
      # Get distinct rows (there are some duplicates because one day can have both Landsat+Sentinel)
      dplyr::distinct() %>%
      dplyr::group_by(formatted_date) %>%
      # Calculate mean daily salinity
      dplyr::mutate(mean_daily_sal = mean(salinity)) %>%
      # Drop salinity column
      dplyr::select(-salinity) %>%
      # Get distinct rows again (this will drop rows with same dates)
      dplyr::distinct()
    
    station_sub2_v2 <- station_sub2 %>%
      # Filter to rows that have common dates
      dplyr::filter(formatted_date %in% common_dates) %>%
      # Get distinct rows (there are some duplicates because one day can have both Landsat+Sentinel)
      dplyr::distinct() %>%
      dplyr::group_by(formatted_date) %>%
      # Calculate mean daily salinity
      dplyr::mutate(mean_daily_sal = mean(salinity)) %>%
      # Drop salinity column
      dplyr::select(-salinity) %>%
      # Get distinct rows again (this will drop rows with same dates)
      dplyr::distinct()
    
    # One way to get R-squared
    model <- lm(station_sub1_v2$mean_daily_sal ~ station_sub2_v2$mean_daily_sal)
    unique_pairs[i,]$pearsons_sq <- summary(model)$r.squared
  }
}

unique_pairs <- unique_pairs %>%
  dplyr::mutate(next_n_years = next_years)

# Export to CSV
readr::write_csv(unique_pairs, file.path("model_validity_results", 
                                         "xgboost",
                                         paste0(start_year, "_", start_year+my_interval-1),
                                         paste0("obs_station_pairs_rsquared_", start_year, "_", start_year+my_interval-1, "_next_", next_years, ".csv")))

## --------------------------------------------------------------- ##
#    Predicted Data: Getting R-Squared for Pairs of Stations -----
## --------------------------------------------------------------- ##

Fmask_lookup <- read_csv(file.path("appeears_landsat_lookup", "HLSL30-020-Fmask-lookup.csv")) %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL %>%
  # Filter out cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) %>%
  # Grab only continuous values
  dplyr::filter(grab == 0)

# Subset to interval
DBSAL_some_years <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date(paste0(start_year,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval,"-01-01")))

# Split data set into training and test sets
sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
train_some_years <- DBSAL_some_years[sample1,] 
test_some_years <- DBSAL_some_years[-sample1,]

# XGBoost modeling
# Picked variables based on VSURF results
# 2016: selected_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "tavg", "B06")
# 2017: selected_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "NLI", "B06", "tavg")
selected_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "NLI", "B06", "tavg")

if (start_year == 2016){
  my_nrounds = 20
} else if (start_year == 2017){
  my_nrounds = 30
}

xgb_model <- xgboost(train_some_years[, selected_vars], train_some_years$salinity, 
                     objective = "reg:squarederror",
                     nrounds = my_nrounds)

# Sanity check
my_pred <- predict(xgb_model, newdata = test_some_years)
mse <- mean((my_pred - test_some_years$salinity)^2)
r_squared <- 1 - sum((test_some_years$salinity - my_pred)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)
# 0.72 

# Predict for next X years
DBSAL_next_years <- DBSAL %>%
  dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval+next_years,"-01-01")))

my_pred2 <- predict(xgb_model, newdata = DBSAL_next_years)
mse2 <- mean((my_pred2 - DBSAL_next_years$salinity)^2)
r_squared2 <- 1 - sum((DBSAL_next_years$salinity - my_pred2)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
# 0.56  

# Add a new column for predicted values
DBSAL_next_years_v2 <- DBSAL_next_years %>%
  dplyr::mutate(pred = my_pred2)

# Get unique combinations of stations where order doesn't matter
grid2 <- expand.grid(unique(DBSAL$station), unique(DBSAL$station))
unique_pairs2 <- grid2[!duplicated(t(apply(grid2, 1, sort))), ] %>%
  # Remove rows with the same station twice
  dplyr::filter(Var1 != Var2) %>%
  # Rename columns
  dplyr::rename(station1 = Var1) %>%
  dplyr::rename(station2 = Var2) %>%
  # Create new column to store our R-squared values
  dplyr::mutate(pearsons_sq = NA) %>%
  # Convert columns to character type
  dplyr::mutate(station1 = as.character(station1),
                station2 = as.character(station2)) %>%
  # Create new data type column
  dplyr::mutate(type = "pred", .before = station1)

for (i in 1:nrow(unique_pairs2)){
  station_sub1 <- DBSAL_next_years_v2 %>%
    # Filter to one station
    dplyr::filter(station == unique_pairs2[i,]$station1) %>%
    # Filter to specified next years to compare our predictions to 
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval+next_years,"-01-01"))) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, pred)
  
  station_sub2 <- DBSAL_next_years_v2 %>%
    # Filter to another station
    dplyr::filter(station == unique_pairs2[i,]$station2) %>%
    # Filter to specified next years to compare our predictions to 
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval+next_years,"-01-01"))) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, pred)
  
  if(nrow(station_sub1) == 0){
    unique_pairs2[i,]$pearsons_sq <- NA
    next
  } else if (nrow(station_sub2) == 0){
    unique_pairs2[i,]$pearsons_sq <- NA
    next
  } else {
    
    # Find dates that both stations have in common
    common_dates <- intersect(as.character(station_sub1$formatted_date), as.character(station_sub2$formatted_date))
    
    if (length(common_dates) == 0){
      unique_pairs2[i,]$pearsons_sq <- NA
      next
    }
    
    station_sub1_v2 <- station_sub1 %>%
      # Filter to rows that have common dates
      dplyr::filter(formatted_date %in% common_dates) %>%
      # Get distinct rows (there are some duplicates because one day can have both Landsat+Sentinel)
      dplyr::distinct() %>%
      dplyr::group_by(formatted_date) %>%
      # Calculate mean daily salinity
      dplyr::mutate(mean_daily_sal = mean(pred)) %>%
      # Drop pred column
      dplyr::select(-pred) %>%
      # Get distinct rows again (this will drop rows with same dates)
      dplyr::distinct()
    
    station_sub2_v2 <- station_sub2 %>%
      # Filter to rows that have common dates
      dplyr::filter(formatted_date %in% common_dates) %>%
      # Get distinct rows (there are some duplicates because one day can have both Landsat+Sentinel)
      dplyr::distinct() %>%
      dplyr::group_by(formatted_date) %>%
      # Calculate mean daily salinity
      dplyr::mutate(mean_daily_sal = mean(pred)) %>%
      # Drop pred column
      dplyr::select(-pred) %>%
      # Get distinct rows again (this will drop rows with same dates)
      dplyr::distinct()
    
    # One way to get R-squared
    model <- lm(station_sub1_v2$mean_daily_sal ~ station_sub2_v2$mean_daily_sal)
    unique_pairs2[i,]$pearsons_sq <- summary(model)$r.squared
  }
}

unique_pairs2 <- unique_pairs2 %>%
  dplyr::mutate(next_n_years = next_years)

# Export to CSV
readr::write_csv(unique_pairs2, file.path("model_validity_results", 
                                          "xgboost",
                                          paste0(start_year, "_", start_year+my_interval-1),
                                          paste0("pred_station_pairs_rsquared_", start_year, "_", start_year+my_interval-1, "_next_", next_years, ".csv")))

## --------------------------------------------- ##
#                 Harmonizing -----
## --------------------------------------------- ##

start_year <- 2017
my_interval <- 5

selected_interval <- paste0(start_year, "_", start_year+my_interval-1)

next_years <- 5

distance_station_pairs <- readr::read_csv("distance_station_pairs.csv")

path <- file.path("model_validity_results", "xgboost", selected_interval)

# List files 
files_to_harmonize <- c(file.path(path, paste0("obs_station_pairs_rsquared_", selected_interval, "_next_", next_years, ".csv")),
                        (file.path(path, paste0("pred_station_pairs_rsquared_", selected_interval, "_next_", next_years, ".csv"))))

files_to_harmonize

results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) %>%
  # Create new column to indicate model type
  dplyr::mutate(model_name = "xgb", .before = type) %>%
  # Pivot wider to reorganize
  tidyr::pivot_wider(names_from = type,
                     values_from = pearsons_sq) %>%
  # Rename pivoted columns
  dplyr::rename(rsq_for_obs_station_pairs = obs) %>%
  dplyr::rename(rsq_for_pred_station_pairs = pred) %>%
  # Round values
  dplyr::mutate(rsq_for_obs_station_pairs = round(rsq_for_obs_station_pairs, digits = 4),
                rsq_for_pred_station_pairs = round(rsq_for_pred_station_pairs, digits = 4)) %>%
  # Combine with distance between stations info
  dplyr::left_join(distance_station_pairs) %>%
  # Drop station coordinates
  dplyr::select(-station1_coord, -station2_coord) %>%
  # Create new column with difference between R-squared values for observed and predicted data
  dplyr::mutate(rsq_abs_diff = abs(rsq_for_obs_station_pairs - rsq_for_pred_station_pairs), .before = dist_m)

# Export to CSV
readr::write_csv(results_harmonized, file.path("model_validity_results", 
                                               "xgboost",
                                               selected_interval,
                                               paste0("xgb_station_pairs_", selected_interval, "_next_", next_years, "_harmonized.csv")))
