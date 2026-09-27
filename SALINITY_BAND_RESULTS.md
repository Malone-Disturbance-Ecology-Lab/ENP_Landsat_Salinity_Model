# Landsat salinity screening results

## Main findings

The combination of bands carries more predictive information than any band alone.
Green (B03) has the largest held-out permutation importance for daily mean salinity
in both temporal and spatial tests. Blue (B02) is the strongest single band in
temporal tests, but its mean R-squared is only 0.083. Green is the strongest single
band in spatial tests, but its mean R-squared is negative (-0.096).

No individual band provides convincing standalone temporal or regional prediction
in this screen. The all-band model has moderate pooled temporal scores, but does
not beat a station's historical mean for daily or 30-day mean salinity. Negative
within-station diagnostic scores reinforce that pooled performance must not be
read as strong tracking of fluctuations at a station.

## All-band performance

Values below are equally weighted averages of three held-out folds. Temporal
folds are 2023, 2024, and 2025, each trained on earlier nonoverlapping windows.
Spatial folds are three geographic station clusters with all dates held out.

| Target | Temporal R-squared | Spatial R-squared |
| --- | ---: | ---: |
| Daily mean | 0.593 | 0.213 |
| 7-day mean | 0.613 | 0.204 |
| 15-day mean | 0.611 | 0.195 |
| 30-day mean | 0.613 | 0.184 |
| Daily minimum | 0.577 | 0.250 |
| 7-day minimum | 0.587 | 0.272 |
| 15-day minimum | 0.562 | 0.276 |
| 30-day minimum | 0.544 | 0.294 |
| Daily maximum | 0.595 | 0.162 |
| 7-day maximum | 0.625 | 0.079 |
| 15-day maximum | 0.622 | 0.061 |
| 30-day maximum | 0.621 | 0.010 |
| Daily variance | 0.275 | -3.213 |
| 7-day variance | 0.124 | -0.760 |
| 15-day variance | 0.180 | -0.640 |
| 30-day variance | 0.161 | -0.639 |

The historical station-mean baseline reaches temporal R-squared 0.671 for daily
mean and 0.708 for 30-day mean. The all-band model has temporal RMSE 7.58 and
spatial RMSE 10.76 for daily mean, in the input salinity units.

The 30-day minimum transfers best spatially among the all-band targets tested,
but R-squared 0.294 remains modest. Longer averaging windows do not solve the
temporal problem. Variance has no demonstrated regional transfer.

## Band contributions

For daily mean, shuffling B03 in the all-band model increases held-out MSE by
75.83 in temporal tests and 55.03 in spatial tests. The temporal ranking is B03,
B02, B04, B06, B07, B05, B01. The spatial ranking is B03, B07, B02, B06, B01,
B04, B05. These are descriptive rankings, not significance tests. Correlated
bands can substitute for one another; low standalone skill does not imply a band
is useless in a multivariate model.

## Terrain and weather comparison

On the common subset with complete context, daily-mean temporal R-squared rises
from 0.648 for terrain/weather alone to 0.686 with bands added. Spatial R-squared
rises from -0.128 to 0.204. This suggests complementary spectral information,
although uncertainty in that improvement has not been estimated. Do not compare
these directly to the larger spectral-only sample as if rows were identical.

## Scope and limitations

Of 36,784 Landsat rows, 9,194 pass spectral and Fmask screening, spanning 90
stations before target-specific coverage filtering. The context-complete subset
contains 4,529 rows. Require 70% window coverage and exclude windows with salinity
flagged above 100. This latter sensitivity rule is not validated agency QC.
Sentinel data were excluded from this experiment.

Available continuous records are daily means, so daily extrema and variance
describe the available records, not recovered subdaily conditions. Units and
agency quality flags were not preserved upstream. Existing band magnitudes are
unusually small and upstream scaling needs auditing. Uniform rescaling does not
change tree splits, but scene-dependent scaling would compromise inference.

These are exploratory comparisons using fixed random forests, not a final model
selection with a separate untouched test set. Temporal and spatial generalization
are tested separately, not jointly. The sensible next steps are to audit scaling
and source salinity metadata, then test station anomalies and water/vegetation
strata before claiming a transferable spectral salinity model.

## Files and verification

Outputs are under `/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model/salinity_band_screening/`:

- `spectral_all_images/`: primary band-only analysis.
- `context_matched/`: controlled comparisons with terrain/weather predictors.
- Each contains the analysis dataset, fold assignments, fold and station metrics,
  held-out predictions, permutation importance, comparison table, and a PNG figure.

`tests/test_band_screening_outputs.R` verifies all 16 targets, three folds per
comparison, unique predictions, independently recomputed RMSE, nonoverlapping
temporal windows, and station separation. Both runs passed.
See `SALINITY_BAND_SCREENING.md` for reproducible methods and commands.
