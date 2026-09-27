# Image-date salinity targets

Run `Rscript 03b_build_salinity_targets.R` from the project directory.
Optional positional arguments are data directory, output directory, and maximum
accepted salinity (default `Inf`). Inputs are `ENP_DBHydro_salinity_df.csv` and
`DBSAL.csv`. Outputs default to the server project's `salinity_targets` directory.
Existing source files and model results are not overwritten.

## Definitions

Targets are matched by station and actual Landsat **or Sentinel** acquisition
date from `DBSAL.csv` (`has_surf_ref == 1`). Filter `instrument == "Landsat"`
for Landsat-only analysis. There is no forward filling. Duplicate predictor
rows caused by repeated salinity matches are collapsed; distinct sensor or
predictor rows are retained. Cloud flags are preserved for subsequent filtering.

For each date, `mean_day`, `min_day`, `max_day`, and `var_day` summarize retained
observations on that day. The equivalent `_7day`, `_15day`, and `_30day` fields
use inclusive windows of -3/+3, -7/+7, and -14/+15 calendar days, respectively.
Variance is sample variance (denominator n-1); fewer than two observations gives
NA. No observations gives NA for all four statistics, not zero.

Window statistics pool all retained source observations, not daily averages.
Days with more observations have more weight. `n_obs`, `n_days`, `n_grab`,
`n_continuous`, `n_high_review`, `coverage`, `start`, and `end` fields accompany
each scale. Coverage is distinct observed days divided by nominal window length.
No minimum coverage is imposed; coverage fields allow later model-specific rules.
Missing dates are not interpolated. Windows near the record edges can be partial.

The inspected continuous DBHydro source contains **daily means**. Therefore these
statistics describe available records, not recovered subdaily extrema or variance.
For a day with one daily mean, min=max=mean and variance is NA. Grab records and
continuous daily summaries are pooled, with counts exposing their contributions.
Original subdaily observations are needed for true within-day statistics.

## Quality control

Exclude missing station, invalid date, missing/nonfinite or negative salinity,
values above an explicitly configured maximum, and exact duplicate input rows.
All source rows and exclusion reasons are saved. Exact duplicates are identified
before normalization; without sample IDs, identical independent samples cannot
be distinguished from duplicate exports.

Values above 100 are flagged for review, **not automatically rejected**. This is
a screening flag, not a validated physical threshold. No outlier correction or
unit conversion is inferred. Source quality flags, units, measurement depths,
and timeseries IDs were discarded by the upstream harmonizer, so agency quality
control cannot be reconstructed. These outputs are provisional pending review
of suspect values and original source metadata.

## Outputs

- `salinity_qc_records.csv`: every input row with QC outcome and review flag.
- `salinity_qc_summary.csv`: counts by mutually exclusive QC outcome.
- `salinity_daily.csv`: daily statistics for retained station-days.
- `salinity_image_date_targets.csv`: one target row per station-image date.
- `DBSAL_image_targets.csv`: image predictors joined to the new target columns.

Existing model scripts still use the original `salinity` response. Select a new
target explicitly and filter missingness for that target and its predictors;
blanket `na.omit()` across these outputs would drop valid rows with undefined
daily variance. Centered windows include future measurements and are intended
for retrospective estimation. Validation must prevent overlapping target windows
from sharing observations across training and test partitions.

Run checks with `Rscript tests/test_salinity_targets.R`.
