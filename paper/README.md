# ESSD paper: figure and statistics scripts

The paper presents only the final one-step model (no vegetation spectra; sea level; trained 2000–2025 with 13 freshwater sites)
and its daily 30 m record for April 2000–December 2025. Results focus on spatial and temporal salinity patterns.
Manuscript (working version): [ESSD manuscript, Google Doc](https://docs.google.com/document/d/1XrcMQsavZ3_EjXpiEht6dzYCoFI50tBKd_fPUGjt-tA/edit).

Figures (numbering as in the paper):
Fig02_workflow.png (fig_workflow.py), Fig03_validation.png (fig_validation.py), Fig04_station_series.png (fig_stations.py),
Fig06_spatial_patterns.png (fig_spatial.py), Fig07_temporal_patterns.png (fig_temporal.py), Fig08_trends.png (fig_trends.py),
Fig09_events.png (fig_events.py). Fig 1 (study area) and Fig 5 (surface-water mask) are unchanged from earlier drafts.

Statistics: station_skill_final.csv, skill_by_region_final.csv, interval_coverage_final.csv (validation_final.py);
interval_coverage_crosscal.csv, interval_coverage_crosscal_detail.csv (interval_crosscal.py); area_by_salinity_class.csv, climatology_by_group.csv (fig_spatial.py); annual_drivers.csv,
interannual_driver_correlations.csv. Trend/event/gauge tables are in long_masked_final/results/.
Products made for this version: climatology_final/ (climatology_final.py) and monthly30m_wetdays/ for 2013–2025 rebuilt from
daily30m_onestep_final/ (wetdays_final.py); the earlier 2013–2025 wet-day files from the HLS product are in superseded_hls_wetdays_2013_2025/.

## Drivers (Sect. 5.6, Fig. 10)
- Fig10_drivers.png — made by scripts/fig_drivers.py
- drivers_permutation_importance_groups.csv / _combined.csv — grouped permutation importance with blocks held out (scripts/drivers_perm_importance.py, drivers_perm_combined.py)
- drivers_group_response.csv — each group set to each station's own 5th–95th percentile (scripts/drivers_group_response.py)
- drivers_attribution_by_group.csv — share of temporal variance removed at 500 m nodes when each group is held at the node's mean, 15th of each month 2000–2025 (scripts/drivers_attrib_nodes.py)
- drivers_partial_dependence.csv — single-variable partial dependence (supplementary; not in the figure)

The figures, CSV tables and a text copy of the manuscript are on the lab server in `best_two_step/paper_final/`; the scripts here regenerate them from the released files (`ENP_DATA_ROOT` points to the server folder).
