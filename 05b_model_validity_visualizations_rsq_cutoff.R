## ---------------------------------------------------- ##
#  Model Validity: Visualizations for R-squared Cutoffs
## ---------------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script plots the selected model variables for each station
## to see if there are any patterns for "good" or "bad" stations.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(caret)
library(rlang)

# Create new folders to store results
dir.create(path = file.path("model_validity_visualizations"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "cutoff_0.5"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "cutoff_0.33"), showWarnings = F)

DBSAL_orig <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>%
  # Removing outliers to see general patterns more easily
  # Removing 4 outlier points where salinity had values of 3400, 1380, 157
  dplyr::filter(salinity < 100 | is.na(salinity)) %>%
  # Removing 1 outlier where CRSI is 8792
  dplyr::filter(CRSI < 8000 | is.na(CRSI))

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
#               Visualization -----
## --------------------------------------------- ##

# Create plots when the R^2 cutoff to be considered "good" is 0.5 ---------

check_0.5 <- DBSAL_filled %>%
  # Select relevant variables
  dplyr::select(station, formatted_date, grab, salinity, distCoast, slope, B03, srad, CRSI, SI, NDSI, B06, NLI, tavg) %>%
  # Denote "good" stations with a "1", "bad" stations with a "0"
  dplyr::mutate(good = dplyr::case_when(
    station == "SEVENPALM" | station == "ENPWP" | station == "TAYLORS3" | station == "ENPCW" | station == "FLAB44" ~ 1,
    T ~ 0
  ))

# List relevant model variables
vars <- c("salinity", "distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")

# For each variable...
for (a_var in vars){
  
  check_0.5_v2 <- check_0.5 %>%
    dplyr::select(formatted_date, station, !!sym(a_var), good) %>%
    na.omit()
  
  # Plot variable, facet by station
  p <- ggplot(check_0.5_v2, aes(x = formatted_date, y = !!sym(a_var), color = good)) +
    geom_line() +
    facet_wrap(~station)
  
  # Save plot
  ggsave(file.path("model_validity_visualizations", "cutoff_0.5", paste0(a_var, "_plot.png")), p, height = 12, width = 18)
}

# Create plots when the R^2 cutoff to be considered "good" is 0.33 ---------

good_stations <- c("2290930", "ENPCW", "ENPLO", "SEVENPALM", "TAYLORUPS",
                   "ENPWP", "TAYLORS3", "TROUT CR_B", "FLAB37", "JOEBAY2E",
                   "MCCORMICK", "MUD_CRKM", "ENPTC", "ENPTR", "FLAB29",
                   "HIGHWAY_CR", "TTI57", "ENPGI", "FLAB30", "FLAB44",
                   "TTI51B")

check_0.33 <- DBSAL_filled %>%
  # Select relevant variables
  dplyr::select(station, formatted_date, grab, salinity, distCoast, slope, B03, srad, CRSI, SI, NDSI, B06, NLI, tavg) %>%
  # Denote "good" stations with a "1", "bad" stations with a "0"
  dplyr::mutate(good = dplyr::case_when(
    station %in% good_stations ~ 1,
    T ~ 0
  ))

# List relevant model variables
vars <- c("salinity", "distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")

# For each variable...
for (a_var in vars){
  
  check_0.33_v2 <- check_0.33 %>%
    dplyr::select(formatted_date, station, !!sym(a_var), good) %>%
    na.omit()
  
  # Plot variable, facet by station
  p <- ggplot(check_0.33_v2, aes(x = formatted_date, y = !!sym(a_var), color = good)) +
    geom_line() +
    facet_wrap(~station)
  
  # Save plot
  ggsave(file.path("model_validity_visualizations", "cutoff_0.33", paste0(a_var, "_plot.png")), p, height = 12, width = 18)
}
