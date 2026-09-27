library(data.table)
root <- "/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/salinity_band_screening"
for (run in c("spectral_all_images", "context_matched")) {
  folder <- file.path(root, run)
  d <- fread(file.path(folder, "analysis_dataset.csv"))
  m <- fread(file.path(folder, "fold_metrics.csv"))
  p <- fread(file.path(folder, "heldout_predictions.csv"))
  stopifnot(uniqueN(m$target) == 16L, all(is.finite(p$predicted)),
            all(m[, .N, by = .(target, regime, model)]$N == 3L),
            !anyDuplicated(p[, .(target, regime, fold, model, station, date)]))
  recomputed <- p[, .(rmse_check = sqrt(mean((observed-predicted)^2)),
                      n_check = .N), by = .(target, regime, fold, model)]
  checked <- merge(m, recomputed, by = c("target", "regime", "fold", "model"))
  stopifnot(nrow(checked) == nrow(m), all(checked$n == checked$n_check),
            max(abs(checked$rmse - checked$rmse_check)) < 1e-8)
  for (target_name in unique(m$target)) {
    suffix <- sub("^[^_]+_", "", target_name)
    x <- d[is.finite(get(target_name)) & get(paste0("coverage_", suffix)) >= 0.7 &
             get(paste0("n_high_review_", suffix)) == 0]
    for (yr in 2023:2025) {
      test <- x[year == yr]
      train <- x[year < yr & get(paste0("end_", suffix)) < min(test[[paste0("start_", suffix)]])]
      stopifnot(max(train[[paste0("end_", suffix)]]) < min(test[[paste0("start_", suffix)]]))
      stopifnot(nrow(train) == m[target == target_name & regime == "temporal" & fold == yr & model == "baseline", n_train])
    }
    for (group in 1:3) {
      stopifnot(!any(x[spatial_fold == group, station] %in% x[spatial_fold != group, station]))
    }
  }
  cat(run, ": all targets, fold counts, prediction metrics and separation checks passed.\n")
}
