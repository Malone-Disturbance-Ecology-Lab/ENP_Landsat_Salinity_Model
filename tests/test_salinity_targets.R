source("03b_build_salinity_targets.R")

raw <- data.table(station = c(rep("A", 7), "B", "A", "A", "A"),
  collectDate = c("2025-01-01", "2025-01-04", "2025-01-04", "2025-01-07",
                  "2025-01-08", "2025-01-19", "2024-12-21", "2025-01-04",
                  "2025-01-04", "2025-01-04", "2025-01-04"),
  value = c(1, 3, 5, 7, 9, 11, 13, 99, -1, NA, 1000), grab = 0L)
raw <- rbind(raw, raw[1])
qc <- clean_salinity(raw, max_salinity = 100)
stopifnot(qc[exclusion_reason == "retained", .N] == 8L,
          qc[exclusion_reason == "exact_duplicate", .N] == 1L)
dates <- data.table(station = c("A", "B", "C"), formatted_date = as.IDate("2025-01-04"))
t <- summarize_targets(qc, dates)
stopifnot(t[station == "A", mean_day] == 4,
          t[station == "A", min_day] == 3,
          t[station == "A", max_day] == 5,
          t[station == "A", var_day] == 2,
          t[station == "A", mean_7day] == 4,
          t[station == "A", n_days_7day] == 3,
          t[station == "A", n_obs_15day] == 5,
          t[station == "A", n_obs_30day] == 7,
          t[station == "A", end_30day - start_30day] == 29,
          is.na(t[station == "B", var_day]),
          t[station == "C", n_obs_day] == 0,
          is.na(t[station == "C", mean_day]))
stopifnot(clean_salinity(raw)[value == 1000, exclusion_reason] == "retained")
cat("Salinity target tests passed.\n")
