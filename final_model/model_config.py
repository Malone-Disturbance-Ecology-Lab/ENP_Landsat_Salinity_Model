"""Predictors and hyperparameters of the final one-step model (as fitted on 2026-10-13; scikit-learn 1.8.0).
Single source for harness.py, fit_final.py, onestep_maps.py and the driver scripts."""

# 53 one-step predictors: 27 site variables, 6 climate, 3 water level, 3 flows, 12 anchor-station, 2 season
F0 = ['distkm',
      'p_Open water_1000m',
      'p_Mangrove_1000m',
      'p_Freshwater marsh/prairie_1000m',
      'p_Swamp scrub/shrubland_1000m',
      'p_Salt marsh_1000m',
      'p_Open water_250m',
      'p_Mangrove_250m',
      'dist_to_eden_cell_km',
      'P30',
      'P90',
      'P180',
      'P365',
      'R30',
      'T90',
      'doy_cos',
      'doy_sin',
      'stage_5km_cm30',
      'stage_cell_cm30',
      'stage_5km_cm90',
      'SRS90',
      'TSC11190',
      'chan_dist_soft_km',
      'chan_dist_strict_km',
      'chan_sinuosity',
      'chan_density_500m',
      'chan_density_1km',
      'inland_water_1km',
      'src_dist_SRS_km',
      'src_dist_TaylorSlough_km',
      'src_dist_C111_km',
      'src_dist_BigCypressGulf_km',
      'w_SRS',
      'w_TaylorSlough',
      'w_C111',
      'w_BigCypress',
      'basin_cost_SRS_km',
      'basin_cost_TaylorSlough_km',
      'basin_cost_C111_km',
      'basin_cost_BigCypress_km',
      'src_flow90_w',
      'anc_ENPWP_30',
      'anc_ENPWP_3',
      'anc_ENPSR_30',
      'anc_ENPSR_3',
      'anc_ENPWW_30',
      'anc_ENPWW_3',
      'anc_ENPGB_30',
      'anc_ENPGB_3',
      'anc_TAYLORS3_30',
      'anc_TAYLORS3_3',
      'anc_ENPWB_30',
      'anc_ENPWB_3']

# sea level: NOAA monthly MSL anomaly (m), month value and trailing 3- and 12-month means
MS = ['msl_1m', 'msl_3m', 'msl_12m']

FE = F0 + MS   # the 56 predictors of the final model

# HistGradientBoostingRegressor parameters (get_params() of the fitted model; early_stopping="auto" is active above 10 000 rows)
PRM = {'categorical_features': 'from_dtype', 'early_stopping': 'auto', 'interaction_cst': None, 'l2_regularization': 1.0, 'learning_rate': 0.05, 'loss': 'squared_error', 'max_bins': 255, 'max_depth': None, 'max_features': 1.0, 'max_iter': 400, 'max_leaf_nodes': 31, 'min_samples_leaf': 200, 'monotonic_cst': None, 'n_iter_no_change': 10, 'quantile': None, 'random_state': 0, 'scoring': 'loss', 'tol': 1e-07, 'validation_fraction': 0.1, 'verbose': 0, 'warm_start': False}
