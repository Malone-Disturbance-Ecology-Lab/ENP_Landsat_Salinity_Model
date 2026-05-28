## --------------------------------------------- ##
#             Neural Network Models
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script uses neuralnet to predict salinity in the Everglades
## by training on chunks of 4 and 5-year data intervals

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(neuralnet)
library(caret)

DBSAL <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

Fmask_lookup <- read_csv("HLSL30-020-Fmask-lookup.csv") %>%
  # Find the Fmask values for cloudy days
  dplyr::filter(Cloud == "Yes")

DBSAL_v2 <- DBSAL %>%
  # Filter out cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) %>%
  # Grab only continuous values
  dplyr::filter(grab == 0) 

## --------------------------------------------- ##
#                 Modeling -----
## --------------------------------------------- ##

set.seed(77)

i <- 2016
my_interval <- 5

# Subset to interval
DBSAL_some_years <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date(paste0(i,"-01-01")) & formatted_date < as.Date(paste0(i+my_interval,"-01-01")))

sample1 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
train_some_years <- DBSAL_some_years[sample1,] 
test_some_years <- DBSAL_some_years[-sample1,]

# Normalize data
preProc <- preProcess(DBSAL_some_years, method = c('center', 'scale'))
DBSAL_some_years_preProc <- predict(preProc, DBSAL_some_years)

# Grab normalized training and test sets
DBSAL_some_years_preProc_train <- DBSAL_some_years_preProc[sample1,]
DBSAL_some_years_preProc_test <- DBSAL_some_years_preProc[-sample1,]


nn_model <- neuralnet(salinity ~ distCoast+ slope+B03+srad+CRSI+SI+NDSI+tavg+B06, 
                      data = DBSAL_some_years_preProc_train, 
                      hidden = c(6, 4), # Two hidden layers with 6 and 4 neurons
                      threshold = 0.2,
                      linear.output = TRUE,) # TRUE for regression, FALSE for classification

# Predict using the compute function
predictions <- compute(nn_model, DBSAL_some_years_preProc_test)
r_squared <- 1 - sum((DBSAL_some_years_preProc_test$salinity - predictions$net.result)^2) / sum((DBSAL_some_years_preProc_test$salinity - mean(DBSAL_some_years_preProc_test$salinity))^2)

# 0.54

nn_model2 <- neuralnet(salinity ~ distCoast+ slope+B03+srad+CRSI+SI+NDSI+tavg+B06, 
                       data = DBSAL_some_years_preProc_train, 
                       hidden = c(7, 5),
                       threshold = 0.2,
                       linear.output = TRUE,) 

predictions2 <- compute(nn_model2, DBSAL_some_years_preProc_test)
r_squared2 <- 1 - sum((DBSAL_some_years_preProc_test$salinity - predictions2$net.result)^2) / sum((DBSAL_some_years_preProc_test$salinity - mean(DBSAL_some_years_preProc_test$salinity))^2)

# 0.56

nn_model3 <- neuralnet(salinity ~ distCoast+ slope+B03+srad+CRSI+SI+NDSI+tavg+B06, 
                       data = DBSAL_some_years_preProc_train, 
                       hidden = c(7, 5, 3),
                       threshold = 0.2,
                       linear.output = TRUE,) 

predictions3 <- compute(nn_model3, DBSAL_some_years_preProc_test)
r_squared3 <- 1 - sum((DBSAL_some_years_preProc_test$salinity - predictions3$net.result)^2) / sum((DBSAL_some_years_preProc_test$salinity - mean(DBSAL_some_years_preProc_test$salinity))^2)

# didn't converge