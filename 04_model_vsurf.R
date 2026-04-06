## --------------------------------------------- ##
#                 VSURF Models
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script uses VSURF to choose the best predictors to predict salinity in the Everglades
## by training on chunks of 4 and 5-year data intervals

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(VSURF)
library(caret)
library(xgboost)
library(parallel)

# Create new folders to store results
dir.create(path = file.path("vsurf_4_years_results"), showWarnings = F)
dir.create(path = file.path("vsurf_5_years_results"), showWarnings = F)

DBSAL <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

Fmask_lookup <- read_csv("HLSL30-020-Fmask-lookup.csv") %>%
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

#start_years <- c(2014, 2015, 2016, 2017)
start_years <- c(2020)
my_interval <- 5

for (i in start_years) {
  
  selected_interval <- paste0(i, "-", i+my_interval-1)
  
  message(paste0("interval from: ", selected_interval))
  
  # Subset to interval
  DBSAL_some_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval,"-01-01")))
  
  sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
  train_some_years <- DBSAL_some_years[sample1,] 
  test_some_years <- DBSAL_some_years[-sample1,]
  
  rf_vsurf <- VSURF(train_some_years[, 8:27], 
                    train_some_years$salinity,
                    ntree = 1000,
                    RFimplem = "randomForest", 
                    clusterType = "PSOCK", 
                    verbose = TRUE,
                    ncores = detectCores() - 1,
                    parallel = TRUE)
  
  # variables selected by VSURF
  selected_vars <- paste(names(train_some_years[, 8:27][rf_vsurf$varselect.pred]), collapse = ", ")
  
  predictions1 <- predict(rf_vsurf, newdata = test_some_years, step = "pred")
  mse1 <- mean((predictions1 - test_some_years$salinity)^2)
  r_squared1 <- 1 - sum((test_some_years$salinity - predictions1)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)
  
  # Predict for next 2 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+2,"-01-01")))
  
  
  predictions2 <- predict(rf_vsurf, newdata = DBSAL_next_years, step = "pred")
  mse2 <- mean((predictions2 - DBSAL_next_years$salinity)^2)
  r_squared2 <- 1 - sum((DBSAL_next_years$salinity - predictions2)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  
  # Predict for next 3 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+3,"-01-01")))
  
  
  predictions3 <- predict(rf_vsurf, newdata = DBSAL_next_years, step = "pred")
  mse3 <- mean((predictions3 - DBSAL_next_years$salinity)^2)
  r_squared3 <- 1 - sum((DBSAL_next_years$salinity - predictions3)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  
  # Predict for next 4 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+4,"-01-01")))
  
  
  predictions4 <- predict(rf_vsurf, newdata = DBSAL_next_years, step = "pred")
  mse4 <- mean((predictions4 - DBSAL_next_years$salinity)^2)
  r_squared4 <- 1 - sum((DBSAL_next_years$salinity - predictions4)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  
  # Predict for next 5 years
  DBSAL_next_years <- DBSAL_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i+my_interval,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval+5,"-01-01")))
  
  
  predictions5 <- predict(rf_vsurf, newdata = DBSAL_next_years, step = "pred")
  mse5 <- mean((predictions5 - DBSAL_next_years$salinity)^2)
  r_squared5 <- 1 - sum((DBSAL_next_years$salinity - predictions5)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
  
  # Save results
  result <-data.frame(interval_length = paste0(my_interval, " years"),
                years = selected_interval,
                nrow_for_specified_years = nrow(DBSAL_some_years),
                VSURF_vars = selected_vars,
                R_squared_test_set = r_squared1,
                R_squared_next_2_years = r_squared2,
                R_squared_next_3_years = r_squared3,
                R_squared_next_4_years = r_squared4,
                R_squared_next_5_years = r_squared5)
  
  # Export as CSV
  readr::write_csv(result, 
                   file.path(paste0("vsurf_", my_interval, "_years_results"), paste0("vsurf_", selected_interval, ".csv")))
}
