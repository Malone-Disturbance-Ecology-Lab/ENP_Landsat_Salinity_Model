> Scripts for the final model are in this folder (`fit_final.py`, `gauge_events.py`, `q_levels.json`) and in `hpc/`. The fitted model (`final_model_onestep.pkl`, `m3_final.so`), CV outputs and sea-level table are on the lab server in `best_two_step/onestep_final/`. Paths in the scripts point to that server layout or to `~/enp_salinity` on Bouchet.

# Final one-step (no-HLS) model, 2026-10-13

This is the model behind the one-step daily 30 m maps (2000-04 to 2025-12). It replaces the released m3 for the one-step record.
- **Predictors (56):** the 53 released one-step predictors plus sea level (msl_1m, msl_3m, msl_12m).
  - Sea level is the NOAA monthly mean sea level averaged over Key West 8724580, Vaca Key 8723970 and Virginia Key 8723214.
    Each station is taken relative to its own 2000–2025 mean (m).
  - msl_1m is the month's value; msl_3m and msl_12m are trailing 3- and 12-month means. Daily values are in msl_by_date.csv.
- **Training data:** 2000-04-01 to 2025-12-31 at all non-anchor stations plus 13 fresh DBHydro field-conductivity sites
  (BCWQA19, EP, NE0, NE1, NP201, P33, P34, P35, P36, P37, S174, TSB, UTS-CN). All rows have equal weight.
  355,156 rows from 97 stations.
- **Model:** HistGradientBoosting with the released hyperparameters (lr 0.05, max_iter 400, max_leaf_nodes 31,
  min_samples_leaf 200, l2 1.0, random_state 0); scikit-learn 1.8.0.

## Skill, 10-block spatial CV (cv_metrics_final.csv; released setup in brackets)
| measure | value [released setup] |
|---|---|
| R² / bias, 2000–2012 | 0.851 / +0.02 [0.81 / +1.04] |
| R² / bias, 2013–2025 | 0.846 / −0.29 [0.83 / +0.07] |
| bias, stations < 2 psu | 2.46 [5.27] |
| bias, stations < 10 psu | 1.77 [3.56] |
| bias, 13 fresh sites (held out) | 1.33 |
| station-mean R² | 0.896 [0.862] |
| long-record gauges (36): trend sign agreement | 64 % [42 %] |
| long-record gauges: Spearman ρ of trends | 0.36 [0.09] |
| median change, 2000–12 to 2013–25 | +0.13 [−0.44]; observed +0.37 |

90 % split-conformal half-widths for predicted 0–5 / 5–15 / 15–30 / >30 psu: 4.34 / 10.30 / 9.97 / 5.85 (q_levels.json).
Coverage: 0.895 (2000–12), 0.905 (2013–25), 0.97 (fresh sites).

## Known limitation
In drought dry seasons (for example May 2020) the unmasked interior does not drop below about 1.5 psu. Most of those cells are dry or
uncertain in the surface-water mask that day. Fresh-site cells, 2019–2021 daily mean: about 1.0–1.2 psu, against about 0.25 observed.

## Files
final_model_onestep.pkl (dict: m3, F, PRM, q_levels, sea_level, train); m3_final.so (tl2cgen, max |compiled − sklearn| = 2e-13);
q_levels.json; msl_by_date.csv; cv_metrics_final.csv; cv_predictions_final.parquet; scripts/fit_final.py;
hpc/daily30_onestep_hpc.py + hpc/run_daily30_onestep.slurm (Bouchet, 309 monthly tasks, output ~/enp_salinity/daily30m_onestep_final).
The worker was tested on a 100-row strip for May 2020: it matches scikit-learn predictions to within the 0.1 psu storage step.
Trial versions v2 and v3 (2026-10-13) were retired and moved to superseded_onestep_trials/.

## Bouchet runs (submitted 2026-10-13)
| job | what | output on Bouchet |
|---|---|---|
| 29464670 | compile m3_final.so (max \|compiled − sklearn\| = 1e-13) | ~/enp_salinity/hpc_inputs/m3_final.so |
| 29464671 [0–308] | daily 30 m maps, 2000-04 to 2025-12 | ~/enp_salinity/daily30m_onestep_final/YYYY/ |
| 29465212 [0–74], after 29464671 | masked 2000–2025 trends (seasonal Kendall / Theil–Sen) + high-salinity events, per tile | ~/enp_salinity/long_masked_final/parts/ |
| 29465213, after 29465212 | reduce step: pixel file with BH-FDR q, region/ecosystem tables, station daily series | ~/enp_salinity/long_masked_final/results/ |
Scripts: hpc/long_masked_hpc.py, long_masked_reduce.py, run_long_*.slurm, submit_long_final.sh. They are the long_masked versions with input switched to
daily30m_onestep_final and output to long_masked_final. The dependent jobs only start if every task of the step before succeeds.
