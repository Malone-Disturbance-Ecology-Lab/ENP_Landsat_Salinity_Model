# Masked daily 30 m surface-water salinity, final model, April 2000 – December 2025 (release record)

`YYYY/ENP_salinity30m_masked_daily_YYYY_MM.tif`: 309 files, one int16 band per day (band description = date YYYY-MM-DD),
salinity in psu × 10 (PSS-78), nodata −1. Grid: EPSG:4326, 3839 × 4171, upper-left (−81.51875°, 25.88965°), cell 0.00027092271°.
DEFLATE + predictor 2, 256 × 256 tiles.

Made 2026-10-15 by `mask_daily30.py` from `../daily30m_onestep_final/` (final one-step model) and the surface-water mask
`../wetdry30m/`: salinity is set to −1 on days with mask code 0 (dry) or 9 (never wet). Uncertain days (code 2) are kept; use the
companion mask to flag or drop them. This is the same masking used for the wet-day monthly means, climatology, trends and events
in the paper; band 2 of `../monthly30m_wetdays/` equals the count of valid days here (checked on 7 months, exact match).

## Storage and reading
- **Scaled integers.** Salinity is stored as int16 psu × 10 (divide by 10 for psu). This halves the file size relative to
  float32 and compresses better (12.3 GB vs roughly 25 GB), and 0.1 psu resolution is far finer than the model error
  (held-out RMSE 5.3 psu). The scale is recorded in the file tags (`units`, `scale`), not as a GDAL band scale factor, so
  most tools will not apply it automatically.
- **No data = −1.** Integer GeoTIFFs cannot hold NA/NaN; missing values are a declared nodata value in the header. −1 is
  impossible for salinity. GDAL, QGIS, R `terra` and Python `rasterio`/`rioxarray` read the declared nodata as missing
  automatically. −1 marks both cells outside the domain and days without surface water (mask code 0 or 9); use
  `../wetdry30m/` to tell them apart or to drop uncertain days (code 2).
- R: `r <- terra::rast("2015/ENP_salinity30m_masked_daily_2015_04.tif") / 10`   # −1 already NA; layer names are dates
- Python: `rioxarray.open_rasterio(f, masked=True) / 10`

