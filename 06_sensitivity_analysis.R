## --------------------------------------------- ##
#             Sensitivity Analysis
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script is for performing sensitivity analysis on our model.
## It investigates how the model predicts when all other variables 
## are held constant (aside from our variable of interest).

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(rlang)

# CHANGE AS NEEDED (also see towards bottom of script) ----------

# Use the modelling data on MaloneLab Server? 0 for no, 1 for yes 
use_data_on_server <- 1

# ---------------------------------------------------------------

# Create new folders to store results
dir.create(path = file.path("sensitivity_analysis"), showWarnings = F)
dir.create(path = file.path("sensitivity_analysis", "random_forest"), showWarnings = F)
# Folders that includes full results (specific predictions) and summary results
dir.create(path = file.path("sensitivity_analysis", "random_forest", "full_results"), showWarnings = F)
dir.create(path = file.path("sensitivity_analysis", "random_forest", "summary_results"), showWarnings = F)

# Variable-specific folders
vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")

for (a_var in vars){
  dir.create(path = file.path("sensitivity_analysis", "random_forest", "full_results", a_var), showWarnings = F)
  dir.create(path = file.path("sensitivity_analysis", "random_forest", "summary_results", a_var), showWarnings = F)
}

if (use_data_on_server == 1){
  # Point to the MaloneLab Server project folder
  project_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Landsat_Salinity_Model") 
} else {
  # Point to our current working directory
  # (It should be your local project folder)
  project_folder <- getwd()
}

# Find the non-grab stations
non_grab_stations <- readr::read_csv(file.path(project_folder, "DBHydro_lonlat.csv")) %>%
  dplyr::filter(grab == 0) %>%
  dplyr::pull(station)

# Read in modelling data
DBSAL_orig <- readr::read_csv(file.path(project_folder, "DBSAL.csv"), col_types = cols(formatted_date = col_date(format = "%Y-%m-%d")))

DBSAL <- DBSAL_orig %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  dplyr::filter(!if_any(everything(), is.infinite)) %>%
  # Remove grab stations
  dplyr::filter(station %in% non_grab_stations)

Fmask_lookup <- readr::read_csv(file.path("appeears_landsat_lookup", "HLSL30-020-Fmask-lookup.csv")) %>%
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
                NDSI = (B05 - B06) / (B05 + B06)) %>%
  # Remove grab stations
  dplyr::filter(station %in% non_grab_stations)

# Select only columns of interest to avoid clutter
DBSAL_filled_v2 <- DBSAL_filled %>%
  dplyr::select(station, formatted_date, salinity, distCoast, slope, B03, srad, CRSI, SI, NDSI, B06, NLI, tavg)

## --------------------------------------------- ##
#             Modeling Function -----
## --------------------------------------------- ##

calc_timeframe <- function(start_year, interval, training_data, filled_data){
  
  set.seed(77)
  
  # Subset to interval
  DBSAL_some_years <- training_data %>%
    dplyr::filter(formatted_date >= as.Date(paste0(start_year,"-01-01")) & formatted_date < as.Date(paste0(start_year+interval,"-01-01")))
  
  # Split data set into training and test sets
  sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
  train_some_years <- DBSAL_some_years[sample1,] 
  test_some_years <- DBSAL_some_years[-sample1,]
  
  set.seed(77)
  
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
  
  full_results <- DBSAL_end_year_v2 %>%
    # Bin the actual salinity values to save in our full results
    dplyr::mutate(salinity_bin = dplyr::case_when(
      salinity >= 0 & salinity <= 5 ~ "0 to 5",
      salinity > 5 & salinity <= 10 ~ "05 to 10",
      salinity > 10 & salinity <= 15 ~ "10 to 15",
      salinity > 15 & salinity <= 20 ~ "15 to 20",
      salinity > 20 & salinity <= 25 ~ "20 to 25",
      salinity > 25 & salinity <= 30 ~ "25 to 30",
      salinity > 30 & salinity <= 35 ~ "30 to 35",
      salinity > 35 & salinity <= 40 ~ "35 to 40",
      salinity > 40 & salinity <= 45 ~ "40 to 45",
      salinity > 45 & salinity <= 50 ~ "45 to 50",
      salinity > 50 ~ "50 and above"
    ), .after = salinity) %>%
    dplyr::mutate(salinity_bin = as.factor(salinity_bin))
  
  # Create lists to store our summary results
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
  for (i in seq_along(unique(full_results$station))){
    
    # Save station name
    stations[[i]] <- unique(full_results$station)[i]
    
    #message(paste("station:", unique(full_results$station)[i]))
    
    station_sub <- full_results %>%
      # Filter to one station
      dplyr::filter(station == unique(full_results$station)[i])
    
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
    #message(paste("station_r_sq:", station_r_sq))
    
  }
  
  # Save summary results in data frame
  summary_results <- data.frame(interval_length = interval,
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
  summary_results <- summary_results %>%
    dplyr::add_row(interval_length = interval,
                   training_years = paste0(start_year, "_", start_year+interval-1),
                   pred_year = paste0(start_year+interval-1),
                   variables = paste(rownames(rf_model$importance), collapse = ", "),
                   station_name = "all",
                   n_rows = nrow(full_results),
                   coeff_det = round(r_squared, digits = 4),
                   avg_diff = round(mean(full_results$diff), digits = 2),
                   min_sal = min(full_results$salinity),
                   mean_sal = round(mean(full_results$salinity), digits = 2),
                   max_sal = max(full_results$salinity),
                   min_pred = round(min(full_results$pred), digits = 2),
                   mean_pred = round(mean(full_results$pred), digits = 2),
                   max_pred = round(max(full_results$pred), digits = 2))
  
  # Return our full and summary results
  return(list(full = full_results,
              summary = summary_results))
  
}

## --------------------------------------------- ##
#                 Execution -----
## --------------------------------------------- ##

# Set the range of years we're interested in
# For example, setting my_years_series <- list(2013:2025) will train the model on 
# the training years 2013-2025, then 2014-2025, then 2015-2025, ..., and finally 2025-2025

# Therefore, the best way to cover every possible timeframe range is to run like so:
# my_years_series <- list(2013:2025, 2013:2024, 2013:2023, 2013:2022, 2013:2021,
#                         2013:2020, 2013:2019, 2013:2018, 2013:2017, 2013:2016, 
#                         2013:2015, 2013:2014, 2013)

# Then run every possible timeframe range for each variable of interest

# CHANGE AS NEEDED ----------------------------

my_years_series <- list(2013:2025)
#my_years_series <- list(2013:2024, 2013:2023, 2013:2022, 2013:2021)
#my_years_series <- list(2013:2020, 2013:2019, 2013:2018, 2013:2017)
#my_years_series <- list(2013:2016, 2013:2015, 2013:2014, 2013)

# Set variable of interest (we let this vary)
my_var <- "distCoast"
# ---------------------------------------------

# For every range in the series...
for (my_years in my_years_series){
  # Set target year to be the last year in the timeframe
  my_target <- my_years[length(my_years)]
  message("target: ", my_target)
  
  # List all model variables
  all_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")
  
  # List variables we hold constant
  const_vars <- setdiff(all_vars, my_var)
  
  # List percentile levels
  all_levels <- c("low", "med", "high")
  
  # For every year in the timeframe range...
  for (i in my_years){
    message("on: ", i)
    
    # Set it as the start year
    my_start_year <- i
    
    # Calculate the interval length from start year to last year
    my_interval <- my_target+1-my_start_year
    
    # For every level...
    for (j in all_levels){
      
      if (j == "low"){
        # Set low percentile
        my_percentile <- 0.2
        
        # Create a fake dataframe with variables constant at that low percentile
        fake_data <- DBSAL_filled_v2 %>% 
          dplyr::mutate(!!sym(const_vars[1]) := quantile(!!sym(const_vars[1]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[2]) := quantile(!!sym(const_vars[2]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[3]) := quantile(!!sym(const_vars[3]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[4]) := quantile(!!sym(const_vars[4]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[5]) := quantile(!!sym(const_vars[5]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[6]) := quantile(!!sym(const_vars[6]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[7]) := quantile(!!sym(const_vars[7]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[8]) := quantile(!!sym(const_vars[8]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[9]) := quantile(!!sym(const_vars[9]), my_percentile, na.rm = T)[[1]])
        
      } else if (j == "med"){
        # Set median percentile
        my_percentile <- 0.5
        
        # Create a fake dataframe with variables constant at that median percentile
        fake_data <- DBSAL_filled_v2 %>% 
          dplyr::mutate(!!sym(const_vars[1]) := quantile(!!sym(const_vars[1]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[2]) := quantile(!!sym(const_vars[2]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[3]) := quantile(!!sym(const_vars[3]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[4]) := quantile(!!sym(const_vars[4]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[5]) := quantile(!!sym(const_vars[5]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[6]) := quantile(!!sym(const_vars[6]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[7]) := quantile(!!sym(const_vars[7]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[8]) := quantile(!!sym(const_vars[8]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[9]) := quantile(!!sym(const_vars[9]), my_percentile, na.rm = T)[[1]])
        
      } else if (j == "high"){
        # Set high percentile
        my_percentile <- 0.8
        
        # Create a fake dataframe with variables constant at that high percentile
        fake_data <- DBSAL_filled_v2 %>% 
          dplyr::mutate(!!sym(const_vars[1]) := quantile(!!sym(const_vars[1]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[2]) := quantile(!!sym(const_vars[2]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[3]) := quantile(!!sym(const_vars[3]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[4]) := quantile(!!sym(const_vars[4]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[5]) := quantile(!!sym(const_vars[5]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[6]) := quantile(!!sym(const_vars[6]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[7]) := quantile(!!sym(const_vars[7]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[8]) := quantile(!!sym(const_vars[8]), my_percentile, na.rm = T)[[1]],
                        !!sym(const_vars[9]) := quantile(!!sym(const_vars[9]), my_percentile, na.rm = T)[[1]])
      }
      
      # Run modelling function for this specified timeframe
      df <- calc_timeframe(start_year = my_start_year,
                           interval = my_interval,
                           training_data = DBSAL_v2,
                           filled_data = fake_data)
      
      # Export the full results that has the specific predictions
      readr::write_csv(df$full, file.path("sensitivity_analysis",
                                          "random_forest",
                                          "full_results",
                                          my_var,
                                          paste0(my_var, "_", j, "_obs_vs_pred_", my_start_year, "_", my_start_year+my_interval-1, "_full_results.csv")))
      
      # Export the summary results by station
      readr::write_csv(df$summary, file.path("sensitivity_analysis",
                                             "random_forest",
                                             "summary_results",
                                             my_var,
                                             paste0(my_var, "_", j, "_obs_vs_pred_", my_start_year, "_", my_start_year+my_interval-1, "_summary.csv")))
    }
  }
}

