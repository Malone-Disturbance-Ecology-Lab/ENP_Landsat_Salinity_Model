## --------------------------------------------- ##
#                 Explore Models
## --------------------------------------------- ##
# Script author(s): Angel Chen

# Purpose:
## This script investigates tree-based models for predicting salinity in the Everglades
## with predictors precip, solar radiation, avg temperature, elevation, slope,
## distance to coast, Landsat+Sentinel bands, NDVI, SI, NLI, SRSI, S7, CRSI, NDSI

## --------------------------------------------- ##
#               Housekeeping -----
## --------------------------------------------- ##

library(tidyverse)
library(randomForest)
library(VSURF)
library(caret)
library(xgboost)

DBSAL <- read_csv("DBSAL.csv", col_types = cols(formatted_date = col_date(format = "%Y-%m-%d"))) %>% 
  # Remove missing values
  na.omit() %>%
  # Remove infinite values
  filter(!if_any(everything(), is.infinite)) 

# Check multicollinearity 
DBSAL_cor <- cor(DBSAL[,c(8:27)], use = "pairwise.complete.obs")
corrplot::corrplot(DBSAL_cor)

Fmask_lookup <- read_csv(file.path("appeears_landsat_lookup", "HLSL30-020-Fmask-lookup.csv")) %>%
  # Find the Fmask values for high aerosol, cloudy days
  dplyr::filter((`Aerosol level` == "High aerosol" & Cloud == "Yes") |
                (`Aerosol level` == "Moderate aerosol" & Cloud == "Yes"))

DBSAL_v2 <- DBSAL %>%
  # Filter out high aerosol & cloudy days
  dplyr::filter(!(Fmask %in% Fmask_lookup$Value)) %>%
  # Getting only the cont salinity values increased r^2 from .71-.72 to .77 ?
  dplyr::filter(grab == 0)
  

## --------------------------------------------- ##
#           Modeling: For All Years -----
## --------------------------------------------- ##

set.seed(1)
sample1 <- sample(1:nrow(DBSAL_v2), 0.8*nrow(DBSAL_v2)) 
train_all_years <- DBSAL_v2[sample1,] 
test_all_years <- DBSAL_v2[-sample1,]


rf_model <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                           B01 + B02 + B03 + B04 + B05 + B06 + B07 + 
                           NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI,
                         data = train_all_years,
                         importance = TRUE)

varImpPlot(rf_model)

my_pred <- predict(rf_model, newdata = test_all_years)
mse <- mean((my_pred - test_all_years$salinity)^2)
r_squared <- 1 - sum((test_all_years$salinity - my_pred)^2) / sum((test_all_years$salinity - mean(test_all_years$salinity))^2)


rf_model2 <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                           NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI,
                         data = train_all_years,
                         importance = TRUE)

varImpPlot(rf_model2)

my_pred2 <- predict(rf_model2, newdata = test_all_years)
mse2 <- mean((my_pred2 - test_all_years$salinity)^2)
r_squared2 <- 1 - sum((test_all_years$salinity - my_pred2)^2) / sum((test_all_years$salinity - mean(test_all_years$salinity))^2)

rf_model3 <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                            B01 + B02 + B03 + B04 + B05 + B06 + B07,
                          data = train_all_years,
                          importance = TRUE)

varImpPlot(rf_model3)

my_pred3 <- predict(rf_model3, newdata = test_all_years)
mse3 <- mean((my_pred3 - test_all_years$salinity)^2)
r_squared3 <- 1 - sum((test_all_years$salinity - my_pred3)^2) / sum((test_all_years$salinity - mean(test_all_years$salinity))^2)

# All did around .78

## --------------------------------------------- ##
#       Modeling: For Subset of Years -----
## --------------------------------------------- ##

# Subset to 2016 and make predictions for 2017
DBSAL_one_year <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date("2016-01-01") & formatted_date < as.Date("2017-01-01"))

sample2 <- sample(1:nrow(DBSAL_one_year), 0.8*nrow(DBSAL_one_year)) 
train_one_year <- DBSAL_one_year[sample2,] 
test_one_year <- DBSAL_one_year[-sample2,]


rf_model4 <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                            NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI,
                          data = train_one_year,
                          importance = TRUE)

varImpPlot(rf_model4)

my_pred4 <- predict(rf_model4, newdata = test_one_year)
mse4 <- mean((my_pred4 - test_one_year$salinity)^2)
r_squared4 <- 1 - sum((test_one_year$salinity - my_pred4)^2) / sum((test_one_year$salinity - mean(test_one_year$salinity))^2)


DBSAL_2017 <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date("2017-01-01") & formatted_date < as.Date("2018-01-01"))

my_pred5 <- predict(rf_model4, newdata = DBSAL_2017)
mse5 <- mean((my_pred5 - DBSAL_2017$salinity)^2)
r_squared5 <- 1 - sum((DBSAL_2017$salinity - my_pred5)^2) / sum((DBSAL_2017$salinity - mean(DBSAL_2017$salinity))^2)

# Did horribly

# What about including all predictors?
rf_model5 <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                            B01 + B02 + B03 + B04 + B05 + B06 + B07 +
                            NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI,
                          data = train_one_year,
                          importance = TRUE)

varImpPlot(rf_model5)

my_pred6 <- predict(rf_model5, newdata = test_one_year)
mse6 <- mean((my_pred6 - test_one_year$salinity)^2)
r_squared6 <- 1 - sum((test_one_year$salinity - my_pred6)^2) / sum((test_one_year$salinity - mean(test_one_year$salinity))^2)

my_pred7 <- predict(rf_model5, newdata = DBSAL_2017)
mse7 <- mean((my_pred7 - DBSAL_2017$salinity)^2)
r_squared7 <- 1 - sum((DBSAL_2017$salinity - my_pred7)^2) / sum((DBSAL_2017$salinity - mean(DBSAL_2017$salinity))^2)

# Still pretty bad

# Subset to 2013-2019
DBSAL_some_years <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date("2013-01-01") & formatted_date < as.Date("2020-01-01"))

sample3 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
train_some_years <- DBSAL_some_years[sample3,] 
test_some_years <- DBSAL_some_years[-sample3,]


rf_model6 <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                            B01 + B02 + B03 + B04 + B05 + B06 + B07 +
                            NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI,
                          data = train_some_years,
                          importance = TRUE)

my_pred8 <- predict(rf_model6, newdata = test_some_years)
mse8 <- mean((my_pred8 - test_some_years$salinity)^2)
r_squared8 <- 1 - sum((test_some_years$salinity - my_pred8)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)

# Predict for 2020 + 2022
# R^2 seems to increase as I predict for a bigger chunk of years
DBSAL_other_subset <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date("2020-01-01") & formatted_date < as.Date("2023-01-01"))


my_pred9 <- predict(rf_model6, newdata = DBSAL_other_subset)
mse9 <- mean((my_pred9 - DBSAL_other_subset$salinity)^2)
r_squared9 <- 1 - sum((DBSAL_other_subset$salinity - my_pred9)^2) / sum((DBSAL_other_subset$salinity - mean(DBSAL_other_subset$salinity))^2)

# Did slightly better than when I trained on 1 year and predicted for the next year


# Subset to 2016-2020
DBSAL_some_years <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date("2016-01-01") & formatted_date < as.Date("2021-01-01"))

sample3 <- sample(1:nrow(DBSAL_some_years), 0.8*nrow(DBSAL_some_years)) 
train_some_years <- DBSAL_some_years[sample3,] 
test_some_years <- DBSAL_some_years[-sample3,]


rf_model_test <- randomForest(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                            B01 + B02 + B03 + B04 + B05 + B06 + B07 +
                            NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI,
                          data = train_some_years,
                          importance = TRUE)

my_pred12 <- predict(rf_model_test, newdata = test_some_years)
mse12 <- mean((my_pred12 - test_some_years$salinity)^2)
r_squared12 <- 1 - sum((test_some_years$salinity - my_pred12)^2) / sum((test_some_years$salinity - mean(test_some_years$salinity))^2)

# Predict for 2021 + 2022
# R^2 seems to increase as I predict for a bigger chunk of years
DBSAL_other_subset <- DBSAL_v2 %>%
  dplyr::filter(formatted_date >= as.Date("2021-01-01") & formatted_date < as.Date("2023-01-01"))


my_pred13 <- predict(rf_model_test, newdata = DBSAL_other_subset)
mse13 <- mean((my_pred13 - DBSAL_other_subset$salinity)^2)
r_squared13 <- 1 - sum((DBSAL_other_subset$salinity - my_pred13)^2) / sum((DBSAL_other_subset$salinity - mean(DBSAL_other_subset$salinity))^2)

# some luck training on 5 years and then predicting for the next 2

## --------------------------------------------- ##
#            Feature Selection -----
## --------------------------------------------- ##

rf_index.sdesign.vsurf <- VSURF(train_all_years[, 8:27], 
                                train_all_years$salinity,
                                ntree = 1000,
                                RFimplem = "randomForest", 
                                clusterType = "PSOCK", 
                                verbose = TRUE,
                                parallel= TRUE)

# VSURF chose distCoast, slope, B03, NDSI, srad, CRSI, SI, B06, B02, tavg
train_all_years[, 8:27][rf_index.sdesign.vsurf$varselect.pred]

predictions <- predict(rf_index.sdesign.vsurf, newdata = test_all_years)

mse10 <- mean((predictions$pred - test_all_years$salinity)^2)
r_squared10 <- 1 - sum((test_all_years$salinity - predictions$pred)^2) / sum((test_all_years$salinity - mean(test_all_years$salinity))^2)

# 0.78

## --------------------------------------------- ##
#              Cross Validation -----
## --------------------------------------------- ##

# Define 5-fold cross-validation
ctrl <- trainControl(method = "cv", number = 5)

# Train the model with CV
model <- train(salinity ~  precip + srad + tavg + elevation + slope + distCoast +
                 B01 + B02 + B03 + B04 + B05 + B06 + B07 +
                 NDVI + SI + NLI + SRSI + S7 + CRSI + NDSI, 
               data = train_all_years, method = "rf", trControl = ctrl)
print(model)

# Best mtry is 12

## --------------------------------------------- ##
#              Gradient Boosting -----
## --------------------------------------------- ##

xgb_model <- xgboost(train_all_years[, 8:27], train_all_years$salinity, 
                     objective = "reg:squarederror",
                     nrounds = 200)

my_pred11 <- predict(xgb_model, test_all_years)
mse11 <- mean((my_pred11 - test_all_years$salinity)^2)
r_squared11 <- 1 - sum((test_all_years$salinity - my_pred11)^2) / sum((test_all_years$salinity - mean(test_all_years$salinity))^2)

# Around same performance as random forest