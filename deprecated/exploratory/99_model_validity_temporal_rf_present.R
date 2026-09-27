## --------------------------------------------- ##
#           Model Validity: Random Forest
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks for the validity of random forest models
## by calculating R-squared for observed vs predicted for every station.
## This tests how well does the model capture temporal patterns in the last year of the training set.

## Notes:
## Including grab samples in training data seem to improve predictions for grab stations for current years,
## but will provide overall worse predictions for any future years.

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
dir.create(path = file.path("model_validity_results_present", "random_forest", "2016_2021"), showWarnings = F)
dir.create(path = file.path("model_validity_results_present", "random_forest", "2017_2021"), showWarnings = F)
dir.create(path = file.path("model_validity_results_present", "random_forest", "2017_2022"), showWarnings = F)

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

## --------------------------------------------- ##
#                 Modeling -----
## --------------------------------------------- ##

set.seed(77)
start_year <- 2017
my_interval <- 6

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
DBSAL_end_year <- DBSAL_filled %>%
  dplyr::filter(formatted_date >= as.Date(paste0(start_year+my_interval-1,"-01-01")) & formatted_date < as.Date(paste0(start_year+my_interval,"-01-01"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

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

for (i in seq_along(unique(DBSAL_end_year_v2$station))){
  
  stations[[i]] <- unique(DBSAL_end_year_v2$station)[i]
  
  message(paste("station:", unique(DBSAL_end_year_v2$station)[i]))
  
  station_sub <- DBSAL_end_year_v2 %>%
    # Filter to one station
    dplyr::filter(station == unique(DBSAL_end_year_v2$station)[i]) 
  
  n_rows_list[[i]] <- nrow(station_sub) 
  
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

# ggplot() + geom_point(aes(x=DBSAL_end_year_v2$salinity, y=DBSAL_end_year_v2$pred, color = DBSAL_end_year_v2$station)) +
#   geom_smooth(method='lm')

# Save results in data frame
results <- data.frame(training_years = paste0(start_year, "-", start_year+my_interval-1),
                   variables = paste(rownames(rf_model$importance), collapse = ", "),
                   station_name = unlist(stations),
                   n_rows = unlist(n_rows_list),
                   pearsons_sq = round(unlist(r_sq1_list), digits = 4),
                   coeff_det = round(unlist(r_sq2_list), digits = 4))

results <- results %>%
  dplyr::add_row(training_years = paste0(start_year, "_", start_year+my_interval-1),
                 variables = paste(rownames(rf_model$importance), collapse = ", "),
                 station_name = "all",
                 n_rows = nrow(DBSAL_end_year_v2),
                 pearsons_sq = caret::postResample(obs = DBSAL_end_year_v2$salinity, pred = DBSAL_end_year_v2$pred)[2],
                 coeff_det = r_squared2)

# Export to CSV
readr::write_csv(results, file.path("model_validity_results_present", 
                                    "random_forest",
                                 paste0(start_year, "_", start_year+my_interval-1),
                                 paste0("obs_vs_pred_", start_year, "_", start_year+my_interval-1, "_grab.csv")))

