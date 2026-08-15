## --------------------------------------------- ##
#                 XGBoost Models
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script uses XGBoost to predict salinity in the Everglades
## by training on chunks of 4 and 5-year data intervals

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(caret)
library(xgboost)
library(parallel)

# Create new folders to store results
dir.create(path = file.path("xgb_5_years_results"), showWarnings = F)

DBSAL <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

Fmask_lookup <- read_csv(file.path("appeears_landsat_lookup", "HLSL30-020-Fmask-lookup.csv")) %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL %>%
  # Filter out cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) %>%
  # Grab only continuous values
  dplyr::filter(grab == 0)

## --------------------------------------------- ##
#                 Modeling -----
## --------------------------------------------- ##

set.seed(77)

# Got these predictors from VSURF results
#selected_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "NLI", "B06", "tavg")
selected_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "tavg", "B06")
start_years <- c(2016)
my_interval <- 5

for (i in start_years) {
  
  selected_interval <- paste0(i, "-", i+my_interval-1)
  
  message(paste0("interval from: ", selected_interval))
  
  # Subset to interval
  DBSAL_some_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval,"-01-01")))
  
  # Split data set into training and test sets
  sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
  train_some_years <- DBSAL_some_years[sample1,] 
  test_some_years <- DBSAL_some_years[-sample1,]
  
  # nrounds = 20 best for 2016-2020
  # nrounds = 30 best for 2017-2021
  xgb_model <- xgboost(train_some_years[, selected_vars], train_some_years$salinity, 
                       objective = "reg:squarederror",
                       nrounds = 20)
  
  # Sanity check
  my_pred1 <- predict(xgb_model, test_some_years)
  mse1 <- mean((my_pred1 - test_some_years$salinity)^2)
  r_squared1 <- 1 - sum((test_some_years$salinity - my_pred1)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)
  
  # Predict for next 2 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+2,"-01-01")))
  
  my_pred2 <- predict(xgb_model, newdata = DBSAL_next_years)
  mse2 <- mean((my_pred2 - DBSAL_next_years$salinity)^2)
  r_squared2 <- 1 - sum((DBSAL_next_years$salinity - my_pred2)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  # Predict for next 3 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+3,"-01-01")))
  
  my_pred3 <- predict(xgb_model, newdata = DBSAL_next_years)
  mse3 <- mean((my_pred3 - DBSAL_next_years$salinity)^2)
  r_squared3 <- 1 - sum((DBSAL_next_years$salinity - my_pred3)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  # Predict for next 4 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+4,"-01-01")))
  
  my_pred4 <- predict(xgb_model, newdata = DBSAL_next_years)
  mse4 <- mean((my_pred4 - DBSAL_next_years$salinity)^2)
  r_squared4 <- 1 - sum((DBSAL_next_years$salinity - my_pred4)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  # Predict for next 5 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+5,"-01-01")))
  
  my_pred5 <- predict(xgb_model, newdata = DBSAL_next_years)
  mse5 <- mean((my_pred5 - DBSAL_next_years$salinity)^2)
  r_squared5 <- 1 - sum((DBSAL_next_years$salinity - my_pred5)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  # Save results
  result <-data.frame(interval_length = paste0(my_interval, " years"),
                      years = selected_interval,
                      nrow_for_specified_years = nrow(DBSAL_some_years),
                      xgb_vars = paste(selected_vars, collapse = ", "),
                      R_squared_test_set = r_squared1,
                      R_squared_next_2_years = r_squared2,
                      R_squared_next_3_years = r_squared3,
                      R_squared_next_4_years = r_squared4,
                      R_squared_next_5_years = r_squared5)
  
  # Export as CSV
  readr::write_csv(result, 
                   file.path(paste0("xgb_", my_interval, "_years_results"), paste0("xgb_", selected_interval, ".csv")))
}
