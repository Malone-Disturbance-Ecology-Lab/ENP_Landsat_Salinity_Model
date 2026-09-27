## --------------------------------------------- ##
#           Model Validity: Random Forest
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks for the validity of random forest models
## by calculating R-squared for observed vs predicted for every station.
## This tests how well does the model capture temporal patterns.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(caret)

# Create new folders to store results
dir.create(path = file.path("model_validity_results"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "random_forest"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "random_forest", "2016_2020"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "random_forest", "2017_2021"), showWarnings = F)
dir.create(path = file.path("model_validity_results", "random_forest", "2013_2020"), showWarnings = F)

DBSAL <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

Fmask_lookup <- read_csv(file.path("appeears_landsat_lookup", "HLSL30-020-Fmask-lookup.csv")) %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL #%>%
  # # Filter out cloudy days
  # dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) #%>%
  # # Grab only continuous values
  # dplyr::filter(grab == 0)

## --------------------------------------------- ##
#                 Modeling -----
## --------------------------------------------- ##

set.seed(77)
start_year <- 2013
my_interval <- 8

# Specify the next years to predict for (need to do 2 through 5)
next_years <- 3

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
# 2017 - 2021: distCoast+slope+B03+srad+CRSI+SI+NDSI+NLI+B06+tavg
# 2013 - 2020 with grab samples: distCoast+slope+B03+srad+SI+CRSI+NDSI+B06+B07+tavg
# 2013 - 2020 with grab and cloudy rows: distCoast+slope+srad+B03+SI+tavg
rf_model <- randomForest(salinity ~ distCoast+slope+srad+B03+SI+tavg,
                         data = train_some_years,
                         importance = TRUE)

# Sanity check
my_pred <- predict(rf_model, newdata = test_some_years)
mse <- mean((my_pred - test_some_years$salinity)^2)
r_squared <- 1 - sum((test_some_years$salinity - my_pred)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)
# 0.73

# Predict for next X years
DBSAL_next_years <- DBSAL %>%
  dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval+next_years,"-01-01")))

my_pred2 <- predict(rf_model, newdata = DBSAL_next_years)
mse2 <- mean((my_pred2 - DBSAL_next_years$salinity)^2)
r_squared2 <- 1 - sum((DBSAL_next_years$salinity - my_pred2)^2) / sum((DBSAL_next_years$salinity - mean(DBSAL_next_years$salinity))^2)
# 0.67 

# Add a new column for predicted values
DBSAL_next_years_v2 <- DBSAL_next_years %>%
  dplyr::mutate(pred = my_pred2)

# Create lists to store our results
r_sq1_list <- list()
r_sq2_list <- list()
stations <- list()

for (i in seq_along(unique(DBSAL$station))){
  
  stations[[i]] <- unique(DBSAL$station)[i]
  
  message(paste("station:", unique(DBSAL$station)[i]))
  
  station_sub <- DBSAL_next_years_v2 %>%
    # Filter to one station
    dplyr::filter(station == unique(DBSAL$station)[i]) 
  
  # Plot to check
  # ggplot() + geom_point(aes(x=station_sub$salinity, y=station_sub$pred)) +
  #   geom_smooth(method='lm')
  
  # One way to get R-squared
  # model <- lm(station_sub$salinity ~ station_sub$pred)
  # summary(model)$r.squared
  
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
results <- data.frame(training_years = paste0(start_year, "-", start_year+my_interval-1),
                   variables = paste(rownames(rf_model$importance), collapse = ", "),
                   station_name = unlist(stations),
                   # Change column names here ----------------
                   pearsons_sq_next_3_years = round(unlist(r_sq1_list), digits = 4),
                   coeff_det_next_3_years = round(unlist(r_sq2_list), digits = 4))

# Export to CSV
readr::write_csv(results, file.path("model_validity_results", 
                                    "random_forest",
                                 paste0(start_year, "_", start_year+my_interval-1),
                                 paste0("obs_vs_pred_next_",next_years,"_years_cloudy.csv")))

## --------------------------------------------- ##
#                 Harmonizing -----
## --------------------------------------------- ##

# # List files 
# files_to_harmonize <- list.files(path = file.path("model_validity_results",
#                                                   "random_forest",
#                                                   "2017_2021"),
#                                  pattern = "obs_vs_pred_next_", full.names = T)
# 
# results_harmonized <- files_to_harmonize %>%
#   # Read in files as CSVs
#   purrr::map(read.csv) %>%
#   # Combine them together
#   purrr::reduce(dplyr::left_join, by = c("training_years", "variables", "station_name")) %>%
#   # Create new column to indicate model type
#   dplyr::mutate(model_name = "rf", .before = training_years)
# 
# # Export to CSV
# readr::write_csv(results_harmonized, file.path("model_validity_results", 
#                                                "random_forest",
#                                                "2017_2021",
#                                                "rf_obs_vs_pred_2017_2021_harmonized.csv"))
