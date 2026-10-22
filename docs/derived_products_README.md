# Single-model record, derived products (final one-step model), April 2000 – December 2025

Made 2026-10-14 from `daily30m_onestep_final/` (final one-step model: no vegetation spectra, sea level, trained 2000–2025 with 13 freshwater
sites; see `onestep_final/README.md`) by `onestep_final/hpc/derive_final.py` and `reduce_monthly250.py` (Bouchet jobs 29547123, 29547124).
These replace `superseded_released_onestep/hindcast_onestep/` (released one-step model trained on 2013–2025).

| folder | content |
|---|---|
| `monthly30/ENP_salinity30m_onestep_monthly_YYYY_MM.tif` | 309 monthly means of the daily 30 m maps (all days; int16 psu × 10, −1 = no data) |
| `daily250/YYYY/ENP_salinity_onestep_daily_YYYY_MM.tif` | daily maps at the 250 m nodes (30 m cells 8i+4, 8j+4); 480 × 521; int16 psu × 10 |
| `monthly250/ENP_salinity_onestep_monthly_{mean,lower90,upper90}_YYYY.tif` | monthly mean of the daily 250 m maps and of their 90 % bounds (float32 psu) |

The 90 % bounds use split-conformal half-widths by predicted level (0–5, 5–15, 15–30, > 30 psu): 4.34, 10.30, 9.97, 5.85 psu.
They come from 2000–2025 block-CV residuals (held out in space). The monthly maps are exact means of the daily maps (no detail-layer approximation).
The wet-day monthly means before 2013 (`monthly30m_wetdays/`) were rebuilt from the same daily maps.
