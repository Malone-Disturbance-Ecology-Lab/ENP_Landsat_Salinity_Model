# Landsat salinity band screening

Run `Rscript 08_screen_salinity_bands.R --spectral-only` for the primary spectral
screen and `Rscript 08_screen_salinity_bands.R` for context-matched comparisons.
Outputs are in `salinity_band_screening/spectral_all_images` and
`salinity_band_screening/context_matched` under the server project, respectively.

The input is the existing station/date join in `salinity_targets/DBSAL_image_targets.csv`.
This uses extracted harmonized Landsat values, not new raster extraction or filled
imagery. Sentinel rows are excluded. Band definitions: B01 coastal aerosol,
B02 blue, B03 green, B04 red, B05 near infrared, B06 SWIR1, B07 SWIR2.

## Screening rules

- Exclude cloud, shadow, adjacent cloud/shadow, snow/ice, and high aerosol using
  the repository's Fmask lookup. Retain both water and land pixels.
- Require finite nonnegative bands no greater than 1. The context-matched run also
  requires complete terrain/weather predictors for fair incremental comparisons.
  The primary spectral-only run retains images without same-day salinity or
  weather matches, which can still have valid multi-day salinity targets.
- Require at least 70% of calendar days represented within each target window.
- Exclude entire target windows containing values flagged above 100, as a
  conservative sensitivity screen. This is not an agency-validated QC decision;
  the underlying salinity dataset remains unchanged.
- Do not interpret within-day min/max/variance as subdaily statistics when source
  records are daily means. Window statistics pool available observations.

## Models and validation

For every mean/min/max/variance target at 1, 7, 15, and 30 days, compare seven
single-band random forests, all bands, context alone (distance to coast, slope,
precipitation, solar radiation, temperature), and context plus all bands.
Forests have 100 trees, minimum node size 10, maximum depth 15, and fixed seed 77.
This is a fixed exploratory screen without hyperparameter tuning.

Temporal: predict 2023, 2024, and 2025 separately, using only earlier years.
Purge training windows that reach the earliest test-window start date. Baseline
is the training station mean, falling back to the overall training mean for a
new station. Pooled temporal R-squared includes differences between stations;
within-station R-squared uses the same prediction errors divided by variation
around test station means and is a stricter diagnostic, not a trained baseline.

Spatial: hold out each of three geographic station clusters in turn, defined by
k-means on latitude and longitude adjusted for local longitude distance. All
dates for each held-out station stay outside training. Baseline is the overall
training mean. This tests regional transfer, not a simultaneous future-year and
new-region forecast. Nearby sites across cluster boundaries may still be similar.

Skip folds with fewer than 100 training rows, 30 test rows, or a constant training
response. There is no random row holdout. Summary metrics average folds equally;
RMSE is comparable between models on the same target, not across mean and variance.
Min/max/mean comparisons may also have different eligible samples from variance.

Held-out permutation importance shuffles each band three times in the all-band
model. Positive delta MSE indicates that disrupting that band worsens predictions.
Correlated bands share information; importance is not causal or uniquely attributable.
Station-level metrics and held-out predictions are saved for detailed assessment.

Spectral values in the existing extracted dataset are unusually small. The
upstream scripts contain multiple scaling steps that need separate auditing.
Random-forest splits are invariant to a uniform positive rescaling, but any
nonuniform scaling across scenes would affect these results. No physical band
thresholds or salinity spectral mechanisms should be inferred from this screen.

These validation folds are used to compare many candidate targets and models.
Results are exploratory; a final independent holdout is needed after selecting
the target and predictors. Centered targets describe retrospective conditions,
not operational forecasting. No uncertainty intervals are estimated in this run.
