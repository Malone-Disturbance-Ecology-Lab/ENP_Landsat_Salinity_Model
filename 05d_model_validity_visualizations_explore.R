## ------------------------------------------------------- ##
#  Model Validity: Visualizations for Exploring Good/Bad 
## ------------------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script plots:
## - A map of good stations
## - Salinity bins
## - Average predictions vs prediction year, colored by good or bad stations
## - Good vs. bad, years and stations
## - Salinity vs. average daily prediction for all prediction years

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(plotly)
library(sf)

# Create new folders to store results
dir.create(path = file.path("model_validity_visualizations"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "explore"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "explore", "sal_vs_pred"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "explore", "diff_sal_bins"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "explore", "model_variables_by_year"), showWarnings = F)

## --------------------------------------------- ##
#         Looking at Map of Stations -----
## --------------------------------------------- ##

DBHydro_sf <- sf::st_read(file.path("ENP_DBHydro_sf", "ENP_DBHydro_sf.shp"))

ENP <- sf::st_read(file.path("Everglades_NP_4326", "Everglades_NP_4326.shp"))

results_0.33 <- read_csv("results_0.33_cutoff.csv")

good_stations <- results_0.33 %>%
  dplyr::filter(count_does_well >= 5) %>%
  dplyr::filter(station_name != "all") %>%
  dplyr::pull(station_name)

bad_stations <- results_0.33 %>%
  dplyr::filter(!(station_name %in% good_stations)) %>%
  dplyr::pull(station_name)

DBHydro_good_stations <- DBHydro_sf %>%
  filter(station %in% good_stations)

DBHydro_bad_stations <- DBHydro_sf %>%
  filter(station %in% bad_stations)

ggplot() +
  geom_sf(aes(geometry = geometry), data = ENP) +
  geom_sf(aes(geometry = geometry), data = DBHydro_good_stations, color = "green") +
  geom_sf(aes(geometry = geometry), data = DBHydro_bad_stations, color = "red")

## --------------------------------------------- ##
#           Looking at Full Results -----
## --------------------------------------------- ##

# Point to folder with model results
full_path <- file.path("model_validity_results_timeframe_diff", "random_forest", "full_results")

# List files 
files_to_harmonize <- list.files(full_path, pattern = "obs_vs_pred_", full.names = T)

files_to_harmonize

full_results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) 

# Remove FLAB29 entirely because it had an outlier where salinity is 1380
# And it greatly affects the scale of the visualizations
full_results_harmonized_v2 <- full_results_harmonized %>%
  dplyr::filter(station != "FLAB29") 

cor.test(full_results_harmonized_v2$salinity, full_results_harmonized_v2$diff)

cor(full_results_harmonized_v2$B03, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$srad, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$CRSI, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$SI, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$NDSI, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$B06, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$NLI, full_results_harmonized_v2$diff)
cor(full_results_harmonized_v2$tavg, full_results_harmonized_v2$diff)

p <- ggplot() +
  geom_density(aes(x = diff), data = full_results_harmonized_v2) +
  facet_wrap(~salinity_bin) +
  coord_cartesian(xlim = c(-35, 35))

ggsave(file.path("model_validity_visualizations", "explore", "diff_sal_bins", paste0("diff_sal_bins_all_years.png")), 
       p, height = 6, width = 8)

# ggplot() +
#   geom_point(aes(x = salinity, y = diff), data = full_results_harmonized_v2)

my_years <- 2013:2025

for (i in my_years){
  sub_results <- full_results_harmonized_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(i,"-01-01")) & formatted_date < as.Date(paste0(i+1,"-01-01")))
  
  p <- ggplot() +
    geom_density(aes(x = diff), data = sub_results) +
    facet_wrap(~salinity_bin) +
    labs(title = paste("Difference betw. salinity and prediction, faceted by salinity bins,", i))
  
  ggsave(file.path("model_validity_visualizations", "explore", "diff_sal_bins", paste0("diff_sal_bins_", i, ".png")), 
         p, height = 6, width = 8)
}

## --------------------------------------------- ##
#         Looking at Summary Results -----
## --------------------------------------------- ##

results_0.33 <- read_csv("results_0.33_cutoff.csv")

good_stations <- results_0.33 %>%
  dplyr::filter(count_does_well >= 5) %>%
  dplyr::filter(station_name != "all") %>%
  dplyr::pull(station_name)

# Point to folder with model results
summary_path <- file.path("model_validity_results_timeframe_diff", "random_forest", "summary_results")

# List files 
files_to_harmonize2 <- list.files(summary_path, pattern = "obs_vs_pred_", full.names = T)

files_to_harmonize2

summary_results_harmonized <- files_to_harmonize2 %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) %>%
  # Round numbers
  dplyr::mutate(coeff_det = round(coeff_det, digits = 4)) %>%
  # Remove infinite values
  dplyr::filter(!if_any(everything(), is.infinite)) %>%
  mutate(pred_year = as.factor(pred_year)) %>%
  # Remove FLAB29 entirely because it had an outlier where salinity is 1380
  # And it greatly affects the scale of the visualizations
  dplyr::filter(station_name != "FLAB29") 

sal_summary <- summary_results_harmonized %>%
  # Group by prediction year and station
  dplyr::group_by(pred_year, station_name) %>%
  # Summarize stats further
  dplyr::summarize(avg_coeff_det = round(mean(coeff_det), digits = 4),
                   avg_sal = mean(mean_sal),
                   avg_pred = round(mean(mean_pred), digits = 2)) %>%
  dplyr::mutate(good_or_bad = case_when(
    station_name %in% good_stations ~ "good",
    station_name == "all" ~ "good",
    T ~ "bad"
  ))

p <- ggplot(aes(x = pred_year, y = avg_pred, group = 1, color = good_or_bad), data = sal_summary) + 
  geom_line()+
  geom_point() +
  facet_wrap(~station_name) +
  theme(axis.text.x = element_text(angle = 45, size = 8)) +
  labs(title = "Average predictions vs. prediction years, colored by good or bad stations")

ggsave(file.path("model_validity_visualizations", "explore", paste0("pred_vs_pred_year.png")), 
       p, height = 16, width = 20)

summary_results_harmonized_v2 <- summary_results_harmonized %>%
  # Bin the R^2 values 
  dplyr::mutate(coeff_det_bin = dplyr::case_when(
    coeff_det < 0 ~ "less than 0",
    coeff_det > 0 & coeff_det <= .10 ~ "0 to .10",
    coeff_det > .10 & coeff_det <= .20 ~ ".10 to .20",
    coeff_det > .20 & coeff_det <= .30 ~ ".20 to .30",
    coeff_det > .30 & coeff_det <= .40 ~ ".30 to .40",
    coeff_det > .40 & coeff_det <= .50 ~ ".40 to .50",
    coeff_det > .50 & coeff_det <= .60 ~ ".50 to .60",
    coeff_det > .60 & coeff_det <= .70 ~ ".60 to .70",
    coeff_det > .70 & coeff_det <= .80 ~ ".70 to .80",
    coeff_det > .80 & coeff_det <= .90 ~ ".80 to .90",
    coeff_det > .90 ~ ".90 to 1"
  ), .after = coeff_det) %>%
  dplyr::filter(avg_diff < 300)


ggplot() +
  geom_density(aes(x = avg_diff), data = summary_results_harmonized_v2) +
  facet_wrap(~coeff_det_bin)


## --------------------------------------------- ##
# Looking at Good vs. Bad Years and Stations -----
## --------------------------------------------- ##

DBSAL_orig <- readr::read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d")))

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

good_year <- 2022
bad_year <- 2024
a_var <- "salinity"

filled_good <- DBSAL_filled %>%
  dplyr::filter(formatted_date >= as.Date(paste0(good_year,"-01-01")) & formatted_date < as.Date(paste0(good_year+1,"-01-01")))  %>%
  dplyr::select(formatted_date, station, !!sym(a_var)) %>%
  na.omit() %>%
  dplyr::mutate(year = as.factor(good_year)) %>%
  dplyr::mutate(date = as.Date(stringr::str_sub(formatted_date, start = 6), "%m-%d"))

filled_bad <- DBSAL_filled %>%
  dplyr::filter(formatted_date >= as.Date(paste0(bad_year,"-01-01")) & formatted_date < as.Date(paste0(bad_year+1,"-01-01")))  %>%
  dplyr::select(formatted_date, station, !!sym(a_var)) %>%
  na.omit() %>%
  dplyr::mutate(year = as.factor(bad_year)) %>%
  dplyr::mutate(date = as.Date(stringr::str_sub(formatted_date, start = 6), "%m-%d"))

t <- rbind(filled_good, filled_bad) %>%
  select(-formatted_date) %>%
  filter(station == "EASTSIDE" | station == "ENPBK" | station == "ENPWB"|
           station == "SEVENPALM" | station == "ENPWP" | station == "TAYLORS3")


p <- ggplot(t, aes(x = date, y = !!sym(a_var), group = year, color = year)) +
  geom_line() +
  facet_wrap(~station) +
  theme(axis.text.x = element_text(angle = 45, size = 6)) +
  scale_x_date(date_labels = "%b") 

ggplotly(p)

## --------------------------------------------- ##
#    Looking at salinity vs. prediction -----
## --------------------------------------------- ##

# Point to folder with model results
full_path <- file.path("model_validity_results_timeframe_diff", "random_forest", "full_results")

# List files 
files_to_harmonize <- list.files(full_path, pattern = "obs_vs_pred_", full.names = T)

files_to_harmonize

full_results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) 

# Remove FLAB29 entirely because it had an outlier where salinity is 1380
# And it greatly affects the scale of the visualizations
full_results_harmonized_v2 <- full_results_harmonized %>%
  # Remove outliers
  dplyr::filter(station != "FLAB29") %>%
  dplyr::filter(salinity < 150) %>%
  dplyr::select(formatted_date, station, salinity, pred) %>%
  dplyr::mutate(formatted_date = as.Date(formatted_date)) %>%
  na.omit() %>%
  dplyr::group_by(formatted_date, station, salinity) %>%
  dplyr::summarize(avg_daily_pred = mean(pred)) %>%
  ungroup()


years <- 2013:2025

for (i in years){
  start_year <- i
  full_results_harmonized_v3 <- full_results_harmonized_v2 %>%
    dplyr::filter(formatted_date >= as.Date(paste0(start_year,"-01-01")) & formatted_date < as.Date(paste0(start_year+1,"-01-01")))
  
  p <- ggplot() +
    geom_line(aes(x = formatted_date, y = salinity, group = 1),color = "darkgreen", data=full_results_harmonized_v3) +
    geom_line(aes(x = formatted_date, y = avg_daily_pred, group = 2), color = "darkred",data=full_results_harmonized_v3) +
    facet_wrap(~station) +
    theme(axis.text.x = element_text(angle = 45, size = 6)) +
    labs(title = paste("salinity (green) vs. avg daily prediction (red),", start_year))
  
  ggsave(file.path("model_validity_visualizations", "explore", "sal_vs_pred", paste0("sal_vs_avgpred_", start_year, ".png")), p, height = 12, width = 18)
}

## --------------------------------------------- ##
#     Looking at Model Variables by Year -----
## --------------------------------------------- ##

# DBSAL_orig <- readr::read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) 
# 
# # Create a version of the data with filled out band values
# # (Landsat and Sentinel values taken from the nearest previous day)
# DBSAL_filled <- DBSAL_orig %>%
#   dplyr::arrange(formatted_date) %>%
#   dplyr::group_by(station) %>%
#   tidyr::fill(starts_with("B0"), .direction = "down") %>%
#   tidyr::fill(Fmask, .direction = "down") %>%
#   tidyr::fill(instrument, .direction = "down") %>%
#   dplyr::ungroup() %>%
#   # Calculate indices
#   dplyr::mutate(NDVI = (B05 - B04) / (B05 + B04),
#                 SI = (B03*B04)^0.5,
#                 NLI = (B05^2 - B04)/(B05^2 + B04),
#                 SRSI = ((NDVI - 1)^2 + SI^2)^0.5,
#                 S7 = (B06 - B07)/(B06 + B07),
#                 CRSI = ((B05*B04-B03*B02)/((B05*B04+B03*B02)))^0.5,
#                 NDSI = (B05 - B06) / (B05 + B06))
# 
# # List relevant model variables
# vars <- c("salinity", "distCoast", "slope", "B03", "srad", "CRSI", "SI", "NDSI", "B06", "NLI", "tavg")
# 
# a_year <- 2022
# 
# for (a_var in vars){
#   
#   filled_year <- DBSAL_filled %>%
#     dplyr::filter(formatted_date >= as.Date(paste0(a_year,"-01-01")) & formatted_date < as.Date(paste0(a_year+1,"-01-01")))  %>%
#     dplyr::select(formatted_date, station, !!sym(a_var)) %>%
#     na.omit()
#   
#   # Plot variable, facet by station
#   p <- ggplot(filled_year, aes(x = formatted_date, y = !!sym(a_var), group = 1)) +
#     geom_line() +
#     facet_wrap(~station) +
#     theme(axis.text.x = element_text(angle = 45, size = 6))
#   
#   # Save plot
#   ggsave(file.path("model_validity_visualizations", "explore", "model_variables_by_year", paste0(a_year), paste0(a_var, "_", a_year, "_plot.png")), p, height = 12, width = 18)
# }
