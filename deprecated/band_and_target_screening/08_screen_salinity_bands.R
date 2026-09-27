# Reproducible exploratory screening; no random row train/test split.
library(data.table)
library(ranger)

root <- "/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model"
spectral_only <- "--spectral-only" %in% commandArgs(trailingOnly = TRUE)
out <- file.path(root, "salinity_band_screening", if (spectral_only) "spectral_all_images" else "context_matched")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
bands <- paste0("B0", 1:7)
context <- c("distCoast", "slope", "precip", "srad", "tavg")
d <- fread(file.path(root, "salinity_targets", "DBSAL_image_targets.csv"))
d <- d[instrument == "Landsat"]
input_rows <- nrow(d)
lookup <- fread("appeears_landsat_lookup/HLSL30-020-Fmask-lookup.csv")
good_masks <- lookup[Cloud == "No" & `Cloud shadow` == "No" &
  `Adjacent to cloud/shadow` == "No" & `Snow/ice` == "No" &
  `Aerosol level` != "High aerosol", Value]
d <- d[Fmask %in% good_masks]
valid_spectra <- Reduce(`&`, lapply(d[, ..bands], function(x) is.finite(x) & x >= 0 & x <= 1))
d <- d[valid_spectra]
clear_rows <- nrow(d)
# Same rows for all predictor sets, so incremental skill is directly comparable.
if (!spectral_only) d <- d[complete.cases(d[, ..context])]
d <- unique(d)
stopifnot(!anyDuplicated(d[, .(station, formatted_date)]))
coords <- unique(fread(file.path(root, "DBHydro_lonlat.csv"))[, .(station, longitude, latitude)])
coords <- coords[station %in% d$station]
stopifnot(!anyDuplicated(coords$station))
set.seed(77)
# Approximate local distances, rather than treating a longitude degree as latitude.
xy <- cbind(coords$longitude * cos(mean(coords$latitude) * pi / 180), coords$latitude)
coords[, spatial_fold := kmeans(xy, centers = 3, nstart = 50)$cluster]
d <- merge(d, coords, by = "station")
d[, year := as.integer(format(formatted_date, "%Y"))]
fwrite(coords, file.path(out, "station_spatial_folds.csv"))
fwrite(data.table(stage = c("Landsat input", "clear plausible spectra", "complete context"),
                 n = c(input_rows, clear_rows, nrow(d))), file.path(out, "sample_flow.csv"))
fwrite(d, file.path(out, "analysis_dataset.csv"))

score <- function(y, p, station) {
  denom <- sum((y - mean(y))^2)
  within <- sum((y - ave(y, station, FUN = mean))^2)
  list(n = length(y), n_stations = uniqueN(station),
       r2 = if (denom > 0) 1 - sum((y-p)^2) / denom else NA_real_,
       within_station_r2 = if (within > 0) 1 - sum((y-p)^2) / within else NA_real_,
       rmse = sqrt(mean((y-p)^2)), mae = mean(abs(y-p)), bias = mean(p-y))
}

models <- c(setNames(lapply(bands, identity), bands),
            list(all_bands = bands, context = context, context_bands = c(context, bands)))
if (spectral_only) models <- models[!names(models) %in% c("context", "context_bands")]
metrics <- list(); predictions <- list(); importance <- list(); skipped <- list()
targets <- as.vector(outer(c("mean", "min", "max", "var"),
                           c("day", "7day", "15day", "30day"), paste, sep = "_"))
for (target in targets) {
  suffix <- sub("^[^_]+_", "", target)
  eligible <- is.finite(d[[target]]) & d[[paste0("coverage_", suffix)]] >= 0.7 &
    d[[paste0("n_high_review_", suffix)]] == 0
  x <- copy(d[eligible])
  x[, response := get(target)]
  for (regime in c("temporal", "spatial")) {
    folds <- if (regime == "temporal") 2023:2025 else 1:3
    for (fold in folds) {
      if (regime == "temporal") {
        test <- x[year == fold]
        if (!nrow(test)) next
        first_test_start <- min(test[[paste0("start_", suffix)]])
        train <- x[year < fold & get(paste0("end_", suffix)) < first_test_start]
      } else {
        train <- x[spatial_fold != fold]
        test <- x[spatial_fold == fold]
      }
      if (nrow(train) && nrow(test)) {
        if (regime == "temporal") stopifnot(max(train[[paste0("end_", suffix)]]) < min(test[[paste0("start_", suffix)]]))
        if (regime == "spatial") stopifnot(!any(train$station %in% test$station))
      }
      if (nrow(train) < 100 || nrow(test) < 30 || uniqueN(train$response) < 2) {
        skipped[[length(skipped)+1L]] <- data.table(target, regime, fold,
          n_train = nrow(train), n_test = nrow(test))
        next
      }
      baseline <- rep(mean(train$response), nrow(test))
      if (regime == "temporal") {
        station_means <- train[, .(baseline = mean(response)), by = station]
        baseline <- station_means$baseline[match(test$station, station_means$station)]
        baseline[is.na(baseline)] <- mean(train$response)
      }
      all_models <- c(list(baseline = character()), models)
      for (model_name in names(all_models)) {
        predictors <- all_models[[model_name]]
        if (model_name == "baseline") {
          pred <- baseline
        } else {
          fit <- ranger(x = as.data.frame(train[, ..predictors]), y = train$response,
                        num.trees = 100, min.node.size = 10, max.depth = 15,
                        mtry = max(1L, floor(sqrt(length(predictors)))),
                        num.threads = 4, seed = 77)
          pred <- predict(fit, data = as.data.frame(test[, ..predictors]))$predictions
          if (model_name == "all_bands") {
            original_mse <- mean((test$response-pred)^2)
            for (band in bands) {
              set.seed(77)
              increments <- replicate(3, {
                permuted <- copy(test[, ..predictors])
                permuted[[band]] <- sample(permuted[[band]])
                perm_pred <- predict(fit, data = as.data.frame(permuted))$predictions
                mean((test$response-perm_pred)^2) - original_mse
              })
              importance[[length(importance)+1L]] <- data.table(target, regime, fold,
                band, delta_mse = mean(increments))
            }
          }
        }
        metrics[[length(metrics)+1L]] <- cbind(data.table(target, regime, fold,
          model = model_name, n_train = nrow(train)), as.data.table(score(test$response, pred, test$station)))
        predictions[[length(predictions)+1L]] <- data.table(target, regime, fold,
          model = model_name, station = test$station, date = test$formatted_date,
          observed = test$response, predicted = pred)
      }
      message(target, " ", regime, " ", fold, ": ", nrow(train), " train, ", nrow(test), " test")
    }
  }
  fwrite(rbindlist(metrics), file.path(out, "fold_metrics.csv"))
}
m <- rbindlist(metrics)
p <- rbindlist(predictions)
imp <- rbindlist(importance)
fwrite(p, file.path(out, "heldout_predictions.csv"))
fwrite(imp, file.path(out, "heldout_band_importance.csv"))
fwrite(rbindlist(skipped), file.path(out, "skipped_folds.csv"))
summary <- m[, .(folds = .N, mean_r2 = mean(r2), worst_r2 = min(r2),
  mean_within_station_r2 = mean(within_station_r2, na.rm = TRUE),
  mean_rmse = mean(rmse), mean_mae = mean(mae)), by = .(target, regime, model)]
setorder(summary, target, regime, -mean_r2)
fwrite(summary, file.path(out, "model_comparison.csv"))
station_scores <- p[, as.data.table(score(observed, predicted, station)),
                    by = .(target, regime, model, fold, station)]
fwrite(station_scores, file.path(out, "station_metrics.csv"))
png(file.path(out, "band_screening.png"), width = 1800, height = 1200, res = 160)
par(mfrow = c(2, 2), mar = c(7, 4, 3, 1))
for (regime_name in c("temporal", "spatial")) {
  for (response_name in c("mean_day", "mean_30day")) {
    z <- summary[target == response_name & regime == regime_name]
    z <- z[match(c("baseline", names(models)), model)]
    barplot(z$mean_r2, names.arg = z$model, las = 2,
            col = c("grey65", rep("#248783", 7), "#b44667", "#647348", "#744b89"),
            main = paste(regime_name, response_name), ylab = "Mean held-out R-squared")
    abline(h = 0, col = "grey40")
  }
}
dev.off()
message("Screening complete: ", out)
