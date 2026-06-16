## ------------------------------------------------------- ##
#  Model Validity: Visualizations for Checking Predictions
## ------------------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script checks if my predictions are over or underestimating
## and plots the difference between predicted & actual salinity value.

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(RColorBrewer)

# Create new folders to store results
dir.create(path = file.path("model_validity_visualizations"), showWarnings = F)
dir.create(path = file.path("model_validity_visualizations", "check_pred"), showWarnings = F)

# Point to folder with model results
path <- file.path("model_validity_results_timeframe_diff", "random_forest", "summary_results")

# List files 
files_to_harmonize <- list.files(path, pattern = "obs_vs_pred_", full.names = T)

files_to_harmonize

results_harmonized <- files_to_harmonize %>%
  # Read in files as CSVs
  purrr::map(read.csv) %>%
  # Combine them together
  purrr::list_rbind(x = .) %>%
  # Round numbers
  dplyr::mutate(coeff_det = round(coeff_det, digits = 4)) %>%
  # Remove infinite values
  dplyr::filter(!if_any(everything(), is.infinite)) %>%
  mutate(pred_year = as.factor(pred_year))

## --------------------------------------------- ##
#               Plotting -----
## --------------------------------------------- ##

# Checking if I over or underestimated more ----------------------

over_under_count <- results_harmonized %>%
  # If avg_diff is negative, mark it as overestimating
  dplyr::mutate(over = case_when(
    avg_diff < 0 ~ 1,
    T ~ 0
  )) %>%
  # If avg_diff is positive, mark it as underestimating
  dplyr::mutate(under = case_when(
    avg_diff > 0 ~ 1,
    T ~ 0
  )) %>%
  # Sum the count
  dplyr::mutate(under_count = sum(under),
                over_count = sum(over))

message("Overestimating count: ", unique(over_under_count$over_count), "\n",
        "Underestimating count: ", unique(over_under_count$under_count))

# Remove FLAB29 entirely because it had an outlier where salinity is 1380
# And it greatly affects the scale of the visualizations
results_harmonized_v2 <- results_harmonized %>%
  dplyr::filter(station_name != "FLAB29") 

# Looks like there are more data points that overestimated
ggplot()+
  geom_point(aes(x = coeff_det, y = avg_diff), data = results_harmonized_v2)

# Could see there are more overestimated points
ggplot() + 
  geom_density(aes(x = avg_diff), data = results_harmonized_v2)

# Facet by station to check whether I over or underestimate for each station
ggplot() + 
  geom_density(aes(x = avg_diff), data = results_harmonized_v2)+
  facet_wrap(~station_name)

# Plotting prediction year vs salinity ----------------------

sal_summary <- results_harmonized_v2 %>%
  # Group by prediction year and station
  dplyr::group_by(pred_year, station_name) %>%
  # Summarize stats further
  dplyr::summarize(avg_coeff_det = round(mean(coeff_det), digits = 4),
            avg_sal = mean(mean_sal),
            avg_pred = round(mean(mean_pred), digits = 2))

# Facet by station
p1 <- ggplot(data = sal_summary) + 
  geom_segment(aes(x = pred_year, xend = pred_year, y = avg_sal, yend = avg_pred), color = "ivory4") +
  geom_point(aes(x = pred_year, y = avg_sal), color = "darkgreen") +
  geom_point(aes(x = pred_year, y = avg_pred), color = "red3") +
  coord_flip() +
  facet_wrap(~station_name) +
  ylab("actual salinity value (green) or prediction (red)")

ggsave(file.path("model_validity_visualizations", "check_pred", "pred_year_vs_sal.png"), 
       p1, height = 16, width = 20)

# Plotting prediction year vs R-squared ----------------------

# Filter to all stations to check overall distribution first
all <- results_harmonized_v2 %>%
  filter(station_name == "all")

ggplot() + 
  geom_point(aes(x=pred_year, y = coeff_det, color = avg_diff), data = all) +
  facet_wrap(~station_name, scales = "free") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_color_binned(type = "viridis")

# Facet by station
# Include the negative R-squared values just to see
p2 <- ggplot() + 
  geom_point(aes(x=pred_year, y = coeff_det, color = avg_diff), data = results_harmonized_v2) +
  facet_wrap(~station_name, scales = "free") +
  labs(title = "Positive and negative R-squared values included") +
  theme(axis.text.x = element_text(angle = 45, size = 6)) +
  scale_color_binned(type = "viridis") 

# Remove negative R-squared values because they bloat the plots
results_harmonized_v3 <- results_harmonized_v2 %>%
  dplyr::filter(coeff_det > 0)

# Facet by station
# Positive R-squared values only
p3 <- ggplot() + 
  geom_point(aes(x=pred_year, y = coeff_det, color = avg_diff), data = results_harmonized_v3) +
  facet_wrap(~station_name, scales = "free") +
  labs(title = "Positive R-squared values only") +
  theme(axis.text.x = element_text(angle = 45, size = 8)) +
  scale_color_binned(type = "viridis") 

ggsave(file.path("model_validity_visualizations", "check_pred", "pred_year_vs_rsq.png"), 
       p2, height = 16, width = 20)

ggsave(file.path("model_validity_visualizations", "check_pred", "pred_year_vs_rsq_positive.png"), 
       p3, height = 16, width = 20)

