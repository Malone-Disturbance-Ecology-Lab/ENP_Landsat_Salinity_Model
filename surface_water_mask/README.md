> Scripts: `scripts/` builds the daily mask (`mask_static.py`, `wetdry30_hpc.py`, `run_wetdry30.slurm`) and checks it against EDEN gauge dry flags (`validate_gauges.py`). `apply_mask/` writes the released masked daily salinity files. The mask rasters themselves are on the lab server (`best_two_step/wetdry30m/`).

# Surface-water (wet/dry) mask for the 30 m salinity maps, 2000-04 → 2025-12

Made 2026-10-13. Purpose: set salinity to NA where there is no surface water.

## Files (in `best_two_step/`)
| folder | content |
|---|---|
| `wetdry30m/YYYY/ENP_wetdry30m_daily_YYYY_MM.tif` | 309 files; one band per day (band description = date); uint8 codes below; same 30 m grid as `daily30m/` |
| `monthly30m_wetdays/ENP_salinity30m_monthly_wetdays_YYYY_MM.tif` | band 1: monthly mean salinity over days kept (psu × 10; −1 = no days kept); band 2: number of days kept. Source: `daily30m_onestep_final/` (final one-step model) for all months, rebuilt 2026-10-14 by `onestep_final/hpc/derive_final.py` (before 2013) and `wetdays_final.py` (from 2013); the earlier 2013–2025 files from the HLS product are in `superseded_hls_wetdays_2013_2025/` |
| `wetdry/` | scripts, inputs (static classes, 30 m DEM bounds), gauge validation, figures |

## Codes
| code | meaning | salinity |
|---|---|---|
| 0 | dry: EDEN water level below the lowest ground in the cell, even with the lower DEM | dropped (NA) |
| 1 | wet: water level above the upper DEM | kept |
| 2 | uncertain: water level within the DEM lower–upper band (about ±4 cm) | kept, flagged |
| 3 | permanent water (Landsat water frequency ≥ 50%, NPS map water or canal) | kept |
| 4 | tidal/intermittent coastal wetland (mangrove, salt marsh, tidal flat); not assessed daily | kept |
| 5 | not assessed (other wetland outside EDEN coverage, or no EDEN value that day) | kept |
| 9 | never wet: upland vegetation, roads, levees, developed | dropped (NA) |
| 255 | outside the domain | — |

Share of the domain: permanent water 39%, EDEN-assessed 32%, tidal coastal 21%, not assessed 5.5%, never wet 2%.

## Method
- Ground elevation: EDEN DEMs with elevation uncertainty treatment v2 (Simons et al. 2024, ver. 2.0 March 2025, doi:10.5066/P1XTHZ7D), 10 m best/lower/upper, cm NAVD88. Each 30 m cell takes the LOWEST 10 m value inside it, so "wet" means surface water somewhere in the cell (keeps sloughs and creeks wet when ridges dry).
- Water level: EDEN v3 daily water surface (400 m, cm NAVD88), bilinear on native EDEN cells only (no extrapolation).
- Static classes from the NPS vegetation map (30 m) and the Landsat water-frequency layer (`site_variables/`).
- Script: `wetdry/scripts/wetdry30_hpc.py`, run on Bouchet (job 29427867, 309 monthly tasks). Static layers: `mask_static.py`.
- Uncertain days are kept by default; set `DROP_UNCERTAIN=1` to drop them too.

## Validation (EDEN marsh gauges, 2013–2025; `wetdry/wetdry_gauge_validation_*.csv*`)
42 freshwater gauges in EDEN-assessed cells, 199,014 gauge-days, against each gauge's own dry ("D") flag:
- gauge wet days: 92.3% called wet, 3.8% uncertain, 4.0% dry
- gauge dry days: 74.3% called dry, 13.1% uncertain, 12.6% wet
- agreement on days called wet or dry: 94.5%
- The grid water level matches the gauges (median difference 0.3 cm), so the errors come from ground elevation.
- Alternatives tested: the 30 m MEAN ground called the deep slough gauges (P33, P36) dry in the 2015 drought (86% agreement overall, 14% of wet days called dry); the lowest-point rule fixes that at the cost of calling more dry slough-edge days wet.

## Using it
- Daily: `sal[code in (0, 9)] = NA` (also code 2 for a strict mask).
- Monthly: band 1 of `monthly30m_wetdays` keeps a value if at least one day was wet. In dry-season months many marsh cells have only a few wet days (April 2015: median 2 days in the EDEN area), so apply a minimum with band 2 (e.g. ≥ 10 days) for maps of typical conditions.
- Caveats: tidal coastal wetland (code 4) is never masked; the mask cannot say when a mangrove floor is exposed. Outside the EDEN domain other wetlands are not assessed (code 5).
