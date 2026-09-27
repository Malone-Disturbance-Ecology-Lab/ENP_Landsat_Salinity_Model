## --------------------------------------------- ##
#        Clean DBSAL and split train / test
## --------------------------------------------- ##
# Purpose:
## Cleans the harmonized salinity + surface reflectance data (DBSAL.csv) and
## splits it into a training set (<= 2024) and a test set (2025).
##
## Cleaning steps:
##  1. Remove implausible salinity values (> max_plausible_sal psu).
##     Florida Bay reached ~71 psu during the 2015 drought (ENPGB), so the
##     cutoff is set above that to keep real hypersaline events.
##  2. Mask surface reflectance for pixels flagged in Fmask as cloud, cirrus,
##     cloud shadow, adjacent to cloud/shadow, or snow/ice. Small negative
##     reflectance (>= -0.01; dark-water noise, ~94% water pixels) is set to 0;
##     reflectance < -0.01 or > 1 is masked.
##  3. Collapse to one row per station per day:
##       - salinity = mean of all readings that day (n_sal_readings, sal_sd kept)
##       - reflectance = mean of the clear Landsat/Sentinel images that day,
##         indices recomputed from the averaged bands
##  4. Split station-days with both salinity and clear reflectance into
##     train (<= train_end_year) and test (test_year).
##
## Inputs:  DBSAL.csv
## Outputs (in cleaned_data/):
##   - DBSAL_daily_clean.csv        all station-days after cleaning
##   - DBSAL_train_2013_2024.csv    matched salinity + clear imagery, <= 2024
##   - DBSAL_test_2025.csv          matched salinity + clear imagery, 2025
##   - removed_salinity_values.csv  log of salinity values removed in step 1
##   - cleaning_summary.csv         row counts at each step

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)

# CHANGE AS NEEDED ----------------------------------------------

# Use the data on MaloneLab Server? 0 for no, 1 for yes
use_data_on_server <- 1

# Salinity values above this (psu) are treated as errors
max_plausible_sal <- 80

# Last training year and the test year
train_end_year <- 2024
test_year <- 2025

# ---------------------------------------------------------------

if (use_data_on_server == 1){
  project_folder <- file.path("/", "Volumes", "malonelab", "Research", "ENP_Landsat_Salinity_Model")
} else {
  project_folder <- getwd()
}

out_folder <- file.path(project_folder, "cleaned_data")
dir.create(out_folder, showWarnings = F)

DBSAL <- readr::read_csv(file.path(project_folder, "DBSAL.csv"),
                         col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"),
                                          instrument = col_character()))

band_cols <- c("B01", "B02", "B03", "B04", "B05", "B06", "B07")
index_cols <- c("NDVI", "SI", "NLI", "SRSI", "S7", "CRSI", "NDSI")
static_cols <- c("elevation", "slope", "distCoast")
met_cols <- c("precip", "srad", "tavg")

summary_log <- tibble(step = "raw DBSAL rows", n = nrow(DBSAL))

## --------------------------------------------- ##
#        1. Remove implausible salinity -----
## --------------------------------------------- ##

removed_sal <- DBSAL %>%
  dplyr::filter(!is.na(salinity) & (salinity > max_plausible_sal | salinity < 0)) %>%
  dplyr::distinct(station, formatted_date, grab, salinity)

readr::write_csv(removed_sal, file.path(out_folder, "removed_salinity_values.csv"))

DBSAL_v1 <- DBSAL %>%
  # Set bad values to NA rather than dropping the row (the image is still valid)
  dplyr::mutate(salinity = ifelse(!is.na(salinity) & (salinity > max_plausible_sal | salinity < 0),
                                  NA_real_, salinity))

summary_log <- summary_log %>%
  dplyr::add_row(step = paste0("salinity values removed (> ", max_plausible_sal, " psu)"),
                 n = nrow(removed_sal))

## --------------------------------------------- ##
#        2. Mask flagged pixels -----
## --------------------------------------------- ##

# HLS Fmask bits: 0 cirrus, 1 cloud, 2 adjacent to cloud/shadow,
#                 3 cloud shadow, 4 snow/ice, 5 water, 6-7 aerosol
bad_fmask <- function(fmask){
  f <- as.integer(fmask)
  bitwAnd(f, 1L) > 0 | bitwAnd(f, 2L) > 0 | bitwAnd(f, 4L) > 0 |
    bitwAnd(f, 8L) > 0 | bitwAnd(f, 16L) > 0
}

# Image-level table: one row per station x date x instrument
images <- DBSAL_v1 %>%
  dplyr::filter(!is.na(instrument)) %>%
  dplyr::distinct(station, formatted_date, instrument, Fmask, dplyr::across(dplyr::all_of(band_cols))) %>%
  # Set tiny negative reflectance (dark-water noise) to 0
  dplyr::mutate(dplyr::across(dplyr::all_of(band_cols), ~ ifelse(!is.na(.x) & .x < 0 & .x >= -0.01, 0, .x))) %>%
  dplyr::mutate(fmask_flagged = is.na(Fmask) | bad_fmask(Fmask),
                refl_out_of_range = dplyr::if_any(dplyr::all_of(band_cols), ~ is.na(.x) | .x < -0.01 | .x > 1),
                clear = !fmask_flagged & !refl_out_of_range)

summary_log <- summary_log %>%
  dplyr::add_row(step = "station-date-instrument images", n = nrow(images)) %>%
  dplyr::add_row(step = "images removed: Fmask cloud/cirrus/shadow/adjacent/snow", n = sum(images$fmask_flagged)) %>%
  dplyr::add_row(step = "images removed: reflectance < -0.01 or > 1 or missing (not already flagged)",
                 n = sum(!images$fmask_flagged & images$refl_out_of_range)) %>%
  dplyr::add_row(step = "clear images kept", n = sum(images$clear))

## --------------------------------------------- ##
#        3. One row per station per day -----
## --------------------------------------------- ##

# Reflectance: average the clear images that day, then recompute indices
daily_refl <- images %>%
  dplyr::filter(clear) %>%
  dplyr::group_by(station, formatted_date) %>%
  dplyr::summarise(instrument = dplyr::case_when(
                     dplyr::n_distinct(instrument) > 1 ~ "Both",
                     T ~ dplyr::first(instrument)),
                   n_images = dplyr::n(),
                   dplyr::across(dplyr::all_of(band_cols), ~ mean(.x)),
                   .groups = "drop") %>%
  dplyr::mutate(NDVI = (B05 - B04) / (B05 + B04),
                SI = (B03*B04)^0.5,
                NLI = (B05^2 - B04)/(B05^2 + B04),
                SRSI = ((NDVI - 1)^2 + SI^2)^0.5,
                S7 = (B06 - B07)/(B06 + B07),
                CRSI = ((B05*B04-B03*B02)/((B05*B04+B03*B02)))^0.5,
                NDSI = (B05 - B06) / (B05 + B06))

# Salinity: the full join in 03 repeats each reading once per image that day,
# so divide the row count by the number of images to get the true reading count
daily_sal <- DBSAL_v1 %>%
  dplyr::filter(!is.na(salinity)) %>%
  dplyr::group_by(station, formatted_date) %>%
  dplyr::summarise(n_images_raw = max(1, dplyr::n_distinct(instrument, na.rm = T)),
                   n_rows = dplyr::n(),
                   salinity_mean = mean(salinity),
                   sal_sd = sd(salinity),
                   sal_min = min(salinity),
                   sal_max = max(salinity),
                   any_grab = any(grab == 1, na.rm = T),
                   .groups = "drop") %>%
  dplyr::mutate(n_sal_readings = n_rows / n_images_raw) %>%
  dplyr::select(-n_images_raw, -n_rows) %>%
  dplyr::rename(salinity = salinity_mean)

# Meteorology (daily) and static site variables
daily_met <- DBSAL_v1 %>%
  dplyr::distinct(station, formatted_date, dplyr::across(dplyr::all_of(met_cols))) %>%
  dplyr::filter(!dplyr::if_all(dplyr::all_of(met_cols), is.na)) %>%
  dplyr::group_by(station, formatted_date) %>%
  dplyr::summarise(dplyr::across(dplyr::all_of(met_cols), ~ mean(.x, na.rm = T)), .groups = "drop")

site_static <- DBSAL_v1 %>%
  dplyr::filter(!dplyr::if_all(dplyr::all_of(static_cols), is.na)) %>%
  dplyr::distinct(station, dplyr::across(dplyr::all_of(static_cols))) %>%
  dplyr::group_by(station) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup()

daily_clean <- dplyr::full_join(daily_sal, daily_refl, by = c("station", "formatted_date")) %>%
  dplyr::left_join(daily_met, by = c("station", "formatted_date")) %>%
  dplyr::left_join(site_static, by = "station") %>%
  dplyr::mutate(year = lubridate::year(formatted_date),
                has_salinity = !is.na(salinity),
                has_clear_image = !is.na(B01)) %>%
  dplyr::relocate(year, has_salinity, has_clear_image, .after = formatted_date) %>%
  dplyr::arrange(formatted_date, station)

readr::write_csv(daily_clean, file.path(out_folder, "DBSAL_daily_clean.csv"))

summary_log <- summary_log %>%
  dplyr::add_row(step = "station-days (all)", n = nrow(daily_clean)) %>%
  dplyr::add_row(step = "station-days with salinity", n = sum(daily_clean$has_salinity)) %>%
  dplyr::add_row(step = "station-days with >1 salinity reading (averaged)",
                 n = sum(daily_clean$n_sal_readings > 1, na.rm = T)) %>%
  dplyr::add_row(step = "station-days with a clear image", n = sum(daily_clean$has_clear_image)) %>%
  dplyr::add_row(step = "station-days with both Landsat & Sentinel (averaged)",
                 n = sum(daily_clean$instrument == "Both", na.rm = T))

## --------------------------------------------- ##
#        4. Train / test split -----
## --------------------------------------------- ##

matched <- daily_clean %>%
  dplyr::filter(has_salinity, has_clear_image) %>%
  # Require bands and site variables only. Indices and meteorology are left as-is:
  #  - CRSI is undefined (sqrt of a negative) for ~1/3 of matched days, mostly
  #    saltier open-water sites; requiring it would bias the data toward fresh sites
  #  - gridMET is missing at a few open-water stations (ENPJK, ENPLR, ENPPK, ENPBA)
  tidyr::drop_na(dplyr::all_of(c(band_cols, static_cols))) %>%
  # Infinite index values (zero denominators) become NA
  dplyr::mutate(dplyr::across(dplyr::all_of(index_cols), ~ ifelse(is.infinite(.x), NA_real_, .x)))

train <- matched %>% dplyr::filter(year <= train_end_year)
test  <- matched %>% dplyr::filter(year == test_year)

readr::write_csv(train, file.path(out_folder, paste0("DBSAL_train_", min(train$year), "_", train_end_year, ".csv")))
readr::write_csv(test,  file.path(out_folder, paste0("DBSAL_test_", test_year, ".csv")))

summary_log <- summary_log %>%
  dplyr::add_row(step = "matched station-days (salinity + clear image + site vars)", n = nrow(matched)) %>%
  dplyr::add_row(step = paste0("train rows (<= ", train_end_year, ")"), n = nrow(train)) %>%
  dplyr::add_row(step = "train stations", n = dplyr::n_distinct(train$station)) %>%
  dplyr::add_row(step = paste0("test rows (", test_year, ")"), n = nrow(test)) %>%
  dplyr::add_row(step = "test stations", n = dplyr::n_distinct(test$station)) %>%
  dplyr::add_row(step = "test stations not in train",
                 n = length(setdiff(unique(test$station), unique(train$station)))) %>%
  dplyr::add_row(step = paste0("matched rows after ", test_year, " (not used)"),
                 n = sum(matched$year > test_year)) %>%
  dplyr::add_row(step = "matched rows with CRSI = NA", n = sum(is.na(matched$CRSI))) %>%
  dplyr::add_row(step = "matched rows with meteorology = NA", n = sum(is.na(matched$precip)))

readr::write_csv(summary_log, file.path(out_folder, "cleaning_summary.csv"))
print(summary_log, n = Inf)
