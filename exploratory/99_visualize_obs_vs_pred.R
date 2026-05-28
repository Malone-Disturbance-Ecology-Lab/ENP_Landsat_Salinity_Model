## --------------------------------------------- ##
#     Visualize Obs vs. Pred R-Squared Values
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script plots observed vs. predicted R-squared values for each station.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(ggrepel)
library(sf)

rf_obs_vs_pred <- read_csv(file.path("model_validity_results", "random_forest", "2016_2020", "rf_obs_vs_pred_2016_2020_harmonized.csv"))

ENP <- sf::st_read(file.path("Everglades_NP_4326", "Everglades_NP_4326.shp"))

DBHydro_lonlat <- read_csv("DBHydro_lonlat.csv") %>%
  dplyr::select(-grab) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  dplyr::left_join(rf_obs_vs_pred, by = c("station" = "station_name"))

ggplot() +
  geom_sf(mapping = aes(geometry = geometry), data = ENP) +
  geom_sf(mapping = aes(geometry = geometry), data = DBHydro_lonlat) +
  # geom_sf_label(mapping = aes(label = coeff_det_next_3_years), data = DBHydro_lonlat, size = 2) +
  geom_label_repel(mapping = aes(label = coeff_det_next_3_years, geometry = geometry), data = DBHydro_lonlat,
                size = 2, stat = "sf_coordinates")

# Zoom in on plot 
                