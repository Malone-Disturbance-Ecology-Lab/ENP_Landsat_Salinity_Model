## --------------------------------------------- ##
#           Model Validity: Random Forest
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks for the validity of random forest models
## by calculating R-squared between each pair of stations for both observed and predicted data.
## This tests how well does the model capture spatial patterns in the last year of the training set.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(caret)

# Create new folders to store results
dir.create(path = file.path("model_validity_results_present"), showWarnings = F)
dir.create(path = file.path("model_validity_results_present", "random_forest"), showWarnings = F)
dir.create(path = file.path("model_validity_results_present", "random_forest", "2016_2020"), showWarnings = F)
dir.create(path = file.path("model_validity_results_present", "random_forest", "2017_2021"), showWarnings = F)

DBSAL_orig <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) 

DBSAL <- DBSAL_orig %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

Fmask_lookup <- read_csv(file.path("appeears_landsat_lookup", "HLSL30-020-Fmask-lookup.csv")) %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL %>%
  # Filter out cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) #%>%
# # Grab only continuous values
# dplyr::filter(grab == 0)

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

## --------------------------------------------------------------- ##
#    Observed Data: Getting R-Squared for Pairs of Stations -----
## --------------------------------------------------------------- ##

set.seed(77)
start_year <- 2017
my_interval <- 6

# Find the stations that have data in the last year of the specified interval
DBSAL_end_year <- DBSAL_filled %>%
  dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval-1,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval,"-01-01"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 
  
# Get unique combinations of stations where order doesn't matter
grid <- expand.grid(unique(DBSAL_end_year$station), unique(DBSAL_end_year$station))
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

n_rows_list <- list() 

for (i in 1:nrow(unique_pairs)){
  station_sub1 <- DBSAL_filled %>%
    # Filter to one station
    dplyr::filter(station == unique_pairs[i,]$station1) %>%
    # Filter to specified years to compare our predictions to 
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval-1,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval,"-01-01"))) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, salinity)
  
  station_sub2 <- DBSAL_filled %>%
    # Filter to another station
    dplyr::filter(station == unique_pairs[i,]$station2) %>%
    # Filter to specified years to compare our predictions to 
    dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval-1,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval,"-01-01"))) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, salinity)
  
  if(nrow(station_sub1) == 0){
    unique_pairs[i,]$pearsons_sq <- NA
    n_rows_list[[i]] <- NA
    next
  } else if (nrow(station_sub2) == 0){
    unique_pairs[i,]$pearsons_sq <- NA
    n_rows_list[[i]] <- NA
    next
  } else {
    
    # Find dates that both stations have in common
    common_dates <- intersect(as.character(station_sub1$formatted_date), as.character(station_sub2$formatted_date))
    
    if (length(common_dates) == 0){
      unique_pairs[i,]$pearsons_sq <- NA
      n_rows_list[[i]] <- NA
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
    
    n_rows_list[[i]] <- nrow(station_sub2_v2)
    
    if (all(is.na(station_sub1_v2$mean_daily_sal)) | all(is.na(station_sub2_v2$mean_daily_sal))){
      unique_pairs[i,]$pearsons_sq <- NA
      next
    }
    
    # One way to get R-squared
    model <- tryCatch({
      lm(station_sub1_v2$mean_daily_sal ~ station_sub2_v2$mean_daily_sal)
    }, error = function(e){
      return(NA)
    })
    
    if (all(is.na(model))){
      unique_pairs[i,]$pearsons_sq <- NA
      next
    } else {
      unique_pairs[i,]$pearsons_sq <- summary(model)$r.squared
    }
  }
}

unique_pairs <- unique_pairs %>%
  dplyr::mutate(n_rows_obs = unlist(n_rows_list), .before = pearsons_sq) %>%
  dplyr::rename(rsq_for_obs_station_pairs = pearsons_sq)


# Export to CSV
readr::write_csv(unique_pairs, file.path("model_validity_results_present", 
                                         "random_forest",
                                         paste0(start_year, "_", start_year+my_interval-1),
                                         paste0("obs_station_pairs_rsquared_", start_year, "_", start_year+my_interval-1, ".csv")))

## --------------------------------------------------------------- ##
#    Predicted Data: Getting R-Squared for Pairs of Stations -----
## --------------------------------------------------------------- ##

# Subset to interval
DBSAL_some_years <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date(paste0(start_year,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval,"-01-01")))

# Split data set into training and test sets
sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
train_some_years <- DBSAL_some_years[sample1,] 
test_some_years <- DBSAL_some_years[-sample1,]

# Random forest modeling
# Picked variables based on VSURF results
# 2016 - 2020: distCoast+slope+B03+srad+CRSI+SI+NDSI+B06+tavg
# 2016 - 2020 with grab: distCoast+slope+B03+srad+SI+CRSI+NDSI+tavg
# 2016 - 2021 with grab: distCoast+slope+B03+srad+SI+CRSI+NDSI+NLI+B05+tavg
# 2017 - 2021 with grab: distCoast+slope+B03+srad+CRSI+SI+NDSI+B06+NLI+tavg
# 2017 - 2022 with grab: distCoast+slope+B03+srad+CRSI+SI+NDSI+B06+NLI+tavg
rf_model <- randomForest(salinity ~ distCoast+slope+B03+srad+CRSI+SI+NDSI+B06+NLI+tavg,
                         data = train_some_years,
                         importance = TRUE)

# Sanity check
my_pred <- predict(rf_model, newdata = test_some_years)
mse <- mean((my_pred - test_some_years$salinity)^2)
r_squared <- 1 - sum((test_some_years$salinity - my_pred)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)

# Predict for the last year in the interval
my_pred2 <- predict(rf_model, newdata = DBSAL_end_year)
mse2 <- mean((my_pred2 - DBSAL_end_year$salinity)^2)
r_squared2 <- 1 - sum((DBSAL_end_year$salinity - my_pred2)^2) / sum((DBSAL_end_year$salinity - mean(DBSAL_end_year$salinity))^2)

# Add a new column for predicted values
DBSAL_end_year_v2 <- DBSAL_end_year %>%
  dplyr::mutate(pred = my_pred2)

# Get unique combinations of stations where order doesn't matter
grid2 <- expand.grid(unique(DBSAL_end_year$station), unique(DBSAL_end_year$station))
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

n_rows_list2 <- list()

for (i in 1:nrow(unique_pairs2)){
  station_sub1 <- DBSAL_end_year_v2 %>%
    # Filter to one station
    dplyr::filter(station == unique_pairs2[i,]$station1) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, pred)
  
  station_sub2 <- DBSAL_end_year_v2 %>%
    # Filter to another station
    dplyr::filter(station == unique_pairs2[i,]$station2) %>%
    # Select relevant columns
    dplyr::select(station, formatted_date, pred)
  
  if(nrow(station_sub1) == 0){
    unique_pairs2[i,]$pearsons_sq <- NA
    n_rows_list2[[i]] <- NA
    next
  } else if (nrow(station_sub2) == 0){
    unique_pairs2[i,]$pearsons_sq <- NA
    n_rows_list2[[i]] <- NA
    next
  } else {
    
    # Find dates that both stations have in common
    common_dates <- intersect(as.character(station_sub1$formatted_date), as.character(station_sub2$formatted_date))
    
    if (length(common_dates) == 0){
      unique_pairs2[i,]$pearsons_sq <- NA
      n_rows_list2[[i]] <- NA
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
    
    n_rows_list2[[i]] <- nrow(station_sub2_v2) 
    
    if (all(is.na(station_sub1_v2$mean_daily_sal)) | all(is.na(station_sub2_v2$mean_daily_sal))){
      unique_pairs2[i,]$pearsons_sq <- NA
      next
    }
    
    # One way to get R-squared
    model <- tryCatch({
      lm(station_sub1_v2$mean_daily_sal ~ station_sub2_v2$mean_daily_sal)
    }, error = function(e){
      return(NA)
    })
    
    if (all(is.na(model))){
      unique_pairs2[i,]$pearsons_sq <- NA
      next
    } else {
      unique_pairs2[i,]$pearsons_sq <- summary(model)$r.squared
    }
  }
}

unique_pairs2 <- unique_pairs2 %>%
  dplyr::mutate(n_rows_pred = unlist(n_rows_list2), .before = pearsons_sq)%>%
  dplyr::rename(rsq_for_pred_station_pairs = pearsons_sq)

# Export to CSV
readr::write_csv(unique_pairs2, file.path("model_validity_results_present", 
                                          "random_forest",
                                          paste0(start_year, "_", start_year+my_interval-1),
                                          paste0("pred_station_pairs_rsquared_", start_year, "_", start_year+my_interval-1, "_grab.csv")))

## --------------------------------------------- ##
#                 Harmonizing -----
## --------------------------------------------- ##

start_year <- 2017
my_interval <- 6

selected_interval <- paste0(start_year, "_", start_year+my_interval-1)

distance_station_pairs <- readr::read_csv("distance_station_pairs.csv")

path <- file.path("model_validity_results_present", "random_forest", selected_interval)

# List files 
files_to_harmonize <- c(file.path(path, paste0("obs_station_pairs_rsquared_", selected_interval, ".csv")),
                        (file.path(path, paste0("pred_station_pairs_rsquared_", selected_interval, "_grab.csv"))))

files_to_harmonize

results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  purrr::map(.f = ~select(.x, -type)) %>%
  # Combine them together
  purrr::reduce(dplyr::left_join) %>%
  # Create new column to indicate model type
  dplyr::mutate(model_name = "rf", .before = station1) %>%
  # Round values
  dplyr::mutate(rsq_for_obs_station_pairs = round(rsq_for_obs_station_pairs, digits = 4),
                rsq_for_pred_station_pairs = round(rsq_for_pred_station_pairs, digits = 4)) %>%
  dplyr::relocate(n_rows_pred, .before = rsq_for_obs_station_pairs) %>%
  # Combine with distance between stations info
  dplyr::left_join(distance_station_pairs) %>%
  # Drop station coordinates
  dplyr::select(-station1_coord, -station2_coord) %>%
  # Combine again, but with the opposite permutation of station names
  dplyr::left_join(distance_station_pairs, by = c("station1" = "station2",
                                                  "station2" = "station1")) %>%
  # Drop station coordinates
  dplyr::select(-station1_coord, -station2_coord) %>%
  # Coalesce distance columns into one
  dplyr::mutate(dist_m = case_when(
    !is.na(dist_m.x) ~ dist_m.x, 
    !is.na(dist_m.y) ~ dist_m.y,
    T ~ NA
  )) %>%
  dplyr::select(-dist_m.x, -dist_m.y) %>%
  # Create new column with difference between R-squared values for observed and predicted data
  dplyr::mutate(rsq_abs_diff = abs(rsq_for_obs_station_pairs - rsq_for_pred_station_pairs), .before = dist_m) 

# Export to CSV
readr::write_csv(results_harmonized, file.path("model_validity_results_present", 
                                               "random_forest",
                                               selected_interval,
                                               paste0("rf_station_pairs_", selected_interval, "_harmonized.csv")))
