## --------------------------------------------- ##
#      Sensitivity Analysis: Visualizations
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script summarizes and visualizes the results from
## the sensitivity analysis.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(plotly)

# Create new folders to store results
dir.create(path = file.path("sensitivity_analysis_visualizations"), showWarnings = F)

## --------------------------------------------- ##
#     Harmonizing and Visualizing Function -----
## --------------------------------------------- ##

harmonize_visualize <- function(var_interest, path_to_files){
  
  message("Variable: ", var_interest)
  
  #  ------------ Harmonizing ------------
  
  message("Now harmonizing for level: high")
  
  # List high files 
  high_files_to_harmonize <- list.files(path_to_files, pattern = "high_obs_vs_pred_", full.names = T)
  
  high_results_harmonized <- high_files_to_harmonize %>%
    # Read in files as CSVs
    purrr::map(read.csv) %>%
    # Combine them together
    purrr::list_rbind(x = .) %>%
    # Format date column
    dplyr::mutate(formatted_date = as.Date(formatted_date))
  
  message("Now harmonizing for level: medium")
  
  # List medium files 
  med_files_to_harmonize <- list.files(path_to_files, pattern = "med_obs_vs_pred_", full.names = T)
  
  med_results_harmonized <- med_files_to_harmonize %>%
    # Read in files as CSVs
    purrr::map(read.csv) %>%
    # Combine them together
    purrr::list_rbind(x = .) %>%
    # Format date column
    dplyr::mutate(formatted_date = as.Date(formatted_date))
  
  message("Now harmonizing for level: low")
  
  # List low files 
  low_files_to_harmonize <- list.files(path_to_files, pattern = "low_obs_vs_pred_", full.names = T)
  
  low_results_harmonized <- low_files_to_harmonize %>%
    # Read in files as CSVs
    purrr::map(read.csv) %>%
    # Combine them together
    purrr::list_rbind(x = .) %>%
    # Format date column
    dplyr::mutate(formatted_date = as.Date(formatted_date))
  
  message("Now summarizing")
  
  # Summarize high results
  high_avg <- high_results_harmonized %>%
    dplyr::group_by(formatted_date, station, !!sym(var_interest)) %>%
    # Create a column for average daily prediction
    dplyr::summarize(high = mean(pred)) %>%
    dplyr::ungroup()
  
  # Summarize medium results
  med_avg <- med_results_harmonized %>%
    dplyr::group_by(formatted_date, station, !!sym(var_interest)) %>%
    # Create a column for average daily prediction
    dplyr::summarize(med = mean(pred)) %>%
    dplyr::ungroup()
  
  # Summarize low results
  low_avg <- low_results_harmonized %>%
    dplyr::group_by(formatted_date, station, !!sym(var_interest)) %>%
    # Create a column for average daily prediction
    dplyr::summarize(low = mean(pred)) %>%
    dplyr::ungroup()
  
  # Combine altogether
  out_df <- list(high_avg, med_avg, low_avg) %>%
    purrr::reduce(dplyr::left_join) 
  
  # Pivot longer for ggplot
  out_df_v2 <- out_df %>%
    tidyr::pivot_longer(cols = c(high, med, low),
                        names_to = "level",
                        values_to = "avg_daily_pred") 
  
  #  ------------ Visualizing ------------
  
  message("Now creating plot")
  
  # Plot
  p <- ggplot() +
    geom_point(aes(x = !!sym(var_interest), y = avg_daily_pred, color = level), data = out_df_v2) +
    #theme(axis.text.x = element_text(angle = 45, size = 10)) +
    labs(title = paste("Average daily prediction vs.", var_interest))
  
  # Export
  ggsave(file.path("sensitivity_analysis_visualizations", paste0("avgpred_vs_", var_interest, ".png")), 
         p, height = 7, width = 10)
  
}

## --------------------------------------------- ##
#                 Execution -----
## --------------------------------------------- ##

# CHANGE AS NEEDED ----------------------------

# List all variables
all_vars <- c("distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")

# For all variables...
for (i in all_vars){
  
  # Select variable of interest
  my_var <- i
  
  # Point to folder with model results
  my_path <- file.path("sensitivity_analysis", "random_forest", "full_results", my_var)
  
  # Visualize
  harmonize_visualize(var_interest = my_var, path_to_files = my_path)
}

# ---------------------------------------------
