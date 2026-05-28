## --------------------------------------------- ##
#         Model Validity: Visualizations
## --------------------------------------------- ##
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

dir.create(path = file.path("model_validity_visualizations"), showWarnings = F)

DBSAL_orig <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d")))  %>%
  # Removing outliers to see general patterns more easily
  # Removing 4 outlier points where salinity had values of 3400, 1380, 157
  filter(salinity < 100) %>%
  # Removing 1 outlier where CRSI is 8792
  filter(CRSI < 8000)

DBSAL <- DBSAL_orig %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite))

Fmask_lookup <- read_csv("HLSL30-020-Fmask-lookup.csv") %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL %>%
  # Filter out cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) 

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

check <- DBSAL_filled %>%
  select(station, formatted_date, grab, salinity, distCoast, slope, B03, srad, CRSI, SI, NDSI, B06, NLI, tavg) %>%
  mutate(good = case_when(
    station == "SEVENPALM" | station == "ENPWP" | station == "TAYLORS3" | station == "ENPCW" | station == "FLAB44" ~ 1,
    T ~ 0
  ))

vars <- c("salinity", "distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")

for (a_var in vars){
  
  p <- ggplot(check, aes(x = formatted_date, y = !!sym(a_var), color = good)) +
    geom_line() +
    facet_wrap(~station)
  
  ggsave(file.path("model_validity_visualizations", paste0(a_var, "_plot.png")), p, height = 12, width = 18)
}
