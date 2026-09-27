# Build observed salinity targets on actual satellite acquisition dates.
# Usage: Rscript 03b_build_salinity_targets.R [data_folder] [output_folder] [max_salinity]
# No spectral or salinity values are carried forward. Requires data.table.
library(data.table)

clean_salinity <- function(raw, max_salinity = Inf) {
  d <- copy(as.data.table(raw))
  stopifnot(all(c("station", "collectDate", "value", "grab") %in% names(d)),
            length(max_salinity) == 1L, !is.na(max_salinity), max_salinity >= 0)
  d[, source_row := .I]
  original_columns <- setdiff(names(d), "source_row")
  duplicates <- duplicated(d, by = original_columns)
  d[, station := trimws(as.character(station))]
  d[, date := as.IDate(as.character(collectDate), format = "%Y-%m-%d")]
  d[, value := suppressWarnings(as.numeric(as.character(value))) ]
  d[, exclusion_reason := fcase(
    is.na(station) | station == "", "missing_station",
    is.na(date), "invalid_date",
    !is.finite(value), "missing_or_nonfinite_salinity",
    value < 0, "negative_salinity",
    value > max_salinity, "above_configured_maximum",
    duplicates, "exact_duplicate",
    default = "retained")]
  # A review flag is not a scientifically validated rejection threshold.
  d[, high_value_review := is.finite(value) & value > 100]
  d[]
}

summarize_targets <- function(clean, image_dates, widths = c(1L, 7L, 15L, 30L)) {
  observations <- clean[exclusion_reason == "retained"]
  observations[, observation_date := date]
  targets <- unique(copy(as.data.table(image_dates))[, .(station, formatted_date)])
  targets[, formatted_date := as.IDate(formatted_date)]
  setorder(targets, station, formatted_date)
  for (width in widths) {
    suffix <- if (width == 1L) "day" else paste0(width, "day")
    before <- (width - 1L) %/% 2L
    after <- width - 1L - before
    bounds <- targets[, .(station, start = formatted_date - before,
                           end = formatted_date + after)]
    stats <- observations[bounds, on = .(station, date >= start, date <= end),
      .(mean = if (all(is.na(value))) NA_real_ else mean(value, na.rm = TRUE),
        min = if (all(is.na(value))) NA_real_ else min(value, na.rm = TRUE),
        max = if (all(is.na(value))) NA_real_ else max(value, na.rm = TRUE),
        var = if (sum(!is.na(value)) < 2L) NA_real_ else var(value, na.rm = TRUE),
        n_obs = sum(!is.na(value)),
        n_days = uniqueN(observation_date[!is.na(value)]),
        n_grab = sum(grab == 1 & !is.na(value), na.rm = TRUE),
        n_continuous = sum(grab == 0 & !is.na(value), na.rm = TRUE),
        n_high_review = sum(high_value_review, na.rm = TRUE)), by = .EACHI]
    metric_names <- c("mean", "min", "max", "var", "n_obs", "n_days",
                      "n_grab", "n_continuous", "n_high_review")
    for (metric in metric_names) {
      set(targets, j = paste0(metric, "_", suffix), value = stats[[metric]])
    }
    targets[, (paste0("coverage_", suffix)) := stats$n_days / width]
    targets[, (paste0("start_", suffix)) := bounds$start]
    targets[, (paste0("end_", suffix)) := bounds$end]
  }
  targets[]
}

build_salinity_targets <- function(data_folder, output_folder, max_salinity = Inf) {
  raw <- fread(file.path(data_folder, "ENP_DBHydro_salinity_df.csv"))
  clean <- clean_salinity(raw, max_salinity)
  joined <- fread(file.path(data_folder, "DBSAL.csv"))
  # Drop old response and sample-type fields before collapsing repeated matches.
  predictor_columns <- setdiff(names(joined), c("salinity", "grab", "has_salinity"))
  images <- unique(joined[has_surf_ref == 1, ..predictor_columns])
  images[, formatted_date := as.IDate(formatted_date)]
  targets <- summarize_targets(clean, images[, .(station, formatted_date)])
  rebuilt <- merge(images, targets, by = c("station", "formatted_date"), all.x = TRUE)
  daily_dates <- unique(clean[exclusion_reason == "retained",
                              .(station, formatted_date = date)])
  # Reuse the same aggregation implementation for daily summaries.
  daily <- summarize_targets(clean, daily_dates, widths = 1L)
  daily <- daily[, c("station", "formatted_date", grep("_day$", names(daily), value = TRUE)), with = FALSE]
  audit <- clean[, .N, by = exclusion_reason]
  dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)
  fwrite(clean, file.path(output_folder, "salinity_qc_records.csv"))
  fwrite(audit, file.path(output_folder, "salinity_qc_summary.csv"))
  fwrite(daily, file.path(output_folder, "salinity_daily.csv"))
  fwrite(targets, file.path(output_folder, "salinity_image_date_targets.csv"))
  fwrite(rebuilt, file.path(output_folder, "DBSAL_image_targets.csv"))
  print(audit)
  message("Image predictor rows: ", nrow(rebuilt), "; station-date targets: ", nrow(targets))
  message("Retained values flagged above 100: ", clean[exclusion_reason == "retained" & high_value_review, .N])
  invisible(list(targets = targets, rebuilt = rebuilt, audit = audit))
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  data_folder <- if (length(args) >= 1L) args[1L] else "/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model"
  output_folder <- if (length(args) >= 2L) args[2L] else file.path(data_folder, "salinity_targets")
  max_salinity <- if (length(args) >= 3L) as.numeric(args[3L]) else Inf
  build_salinity_targets(data_folder, output_folder, max_salinity)
}
