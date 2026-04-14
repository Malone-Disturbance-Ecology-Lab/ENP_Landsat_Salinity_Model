## --------------------------------------------- ##
#       Calculate Distance Between Stations
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script calculates the distance between each unique pair of stations.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(sf)
library(tidyverse)

DBSAL <- readr::read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

## --------------------------------------------- ##
#           Calculating Distance -----
## --------------------------------------------- ##

# Get unique combinations of stations where order doesn't matter
grid <- expand.grid(unique(DBSAL$station), unique(DBSAL$station))
station_pairs <- grid[!duplicated(t(apply(grid, 1, sort))), ] %>%
  # Remove rows with the same station twice
  dplyr::filter(Var1 != Var2) %>%
  # Rename columns
  dplyr::rename(station1 = Var1) %>%
  dplyr::rename(station2 = Var2) %>%
  # Convert columns to character type
  dplyr::mutate(station1 = as.character(station1),
                station2 = as.character(station2))

# Read in station coordinates
DBHydro_lonlat <- readr::read_csv("DBHydro_lonlat.csv") %>%
  dplyr::select(-grab) %>%
  sf::st_as_sf(coords = c("longitude", "latitude"), crs = 4326)

# Combine unique pairs of stations with coordinate info
distance_station_pairs <- station_pairs %>%
  dplyr::left_join(DBHydro_lonlat, by = c("station1" = "station")) %>%
  dplyr::left_join(DBHydro_lonlat, by = c("station2" = "station")) %>%
  dplyr::rename(station1_coord = geometry.x) %>%
  dplyr::rename(station2_coord = geometry.y) %>%
  # Calculate distance between each unique pair of stations
  dplyr::mutate(dist_m = sf::st_distance(station1_coord, station2_coord, by_element = TRUE))
  
# Export
readr::write_csv(distance_station_pairs, "distance_station_pairs.csv")
