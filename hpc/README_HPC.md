# Daily 30 m salinity on HPC (Grace)

This produces the full-model **daily** salinity archive on the native HLS/Landsat 30 m grid, covering 2013-01-01 to 2025-12-31.

- **Grid:** 3839 × 4171 cells, EPSG:4326, cell size 0.00027092271°. This is the same grid as `site_variables/ENP_WaterFrequency_L30.tif`.
- **Output:** 156 monthly multi-band GeoTIFFs named `YYYY/ENP_salinity30m_daily_YYYY_MM.tif`.
  - One band per day; band descriptions are the dates.
  - Values are int16 psu×10; −1 = no data.
  - About 2.5 MB per daily band after DEFLATE compression (the 30 m monthly maps average 2.3 MB), so roughly 11–15 GB for all 4,748 days (≈ 70 GB uncompressed).

## Model and inputs

The model and inputs are the same as for the 30 m monthly maps (`best_two_step/salinity30m/monthly/`):

final = clip(0.7·(B + a) + 0.3·v3, 0)

- B is the step-1 12-month spectral mean.
- a is the step-2 daily anomaly.
- v3 is the one-step model.

Each job resamples all inputs to the 30 m grid itself, using NaN-aware bilinear interpolation:

- HLS spectral composites, from the 16-row vegetation sampling grid;
- gridMET (4 km);
- EDEN water level (400 m).

Static layers are native 30 m: distance to coast, water frequency, connectivity, basin weights and the others.

## 1. Prepare `input_dir` (once)

Copy `best_two_step/salinity30m/hpc_inputs/` from the lab server to Grace. It contains:

| file | content |
|---|---|
| `final_models.pkl` | fitted m1 (monthly spectral mean), m2 (daily anomaly), m3 (one-step v3) |
| `m1.so m2.so m3.so` | treelite-compiled models (x86-64 Linux). If they fail to load, rebuild with `compile_models.py` |
| `static30.npy.gz.part*` | 27 static predictors × 7,622,366 cells (float32), gzip, split into parts |
| `spec_grid.npy.gz.part*` | 32 spectral composites × 156 months on the vegetation sampling grid (float16), gzip, split into parts |
| `static30_keys.json`, `spec_grid_keys.json` | band names |
| `gridmet_filled.npz.part*` | gap-filled daily gridMET-derived rasters (split; `cat gridmet_filled.npz.part* > gridmet_filled.npz`) |
| `daily_global_drivers.csv` | daily basin flows, departures and anchor-station series |

Then run:

```bash
cd hpc_inputs
cat static30.npy.gz.part* | gunzip > static30.npy && rm static30.npy.gz.part*
cat spec_grid.npy.gz.part* | gunzip > spec_grid.npy && rm spec_grid.npy.gz.part*
cat gridmet_filled.npz.part* > gridmet_filled.npz && rm gridmet_filled.npz.part*

# grid index and coordinates (px30_index.npy, lonlat30.npy, utm30.npy)
python ../make_grid30.py <SERVER>/site_variables/ENP_WaterFrequency_L30.tif .

# daily EDEN stage fields (~4.6 GB, ~20 min): python/eden_daily.py and python/fillers.py from the repo.
#  - Replace the SERVER/ prefix with the server root.
#  - The script reads data/EDEN_GRIDS_Q/*.nc and writes eden_daily/ into the current directory.
sed "s#SERVER/#<SERVER>/#g" ../../python/eden_daily.py > eden_daily_local.py
cp ../../python/fillers.py .
python eden_daily_local.py
```

## 2. Run

```bash
module load miniconda; conda create -n sal30 -c conda-forge python=3.11 numpy pandas scikit-learn rasterio pyproj scipy
conda activate sal30; pip install treelite==4.1.2 tl2cgen   # optional, ~12x faster; identical predictions
# edit IN / OUT (and the conda env name) in run_daily30.slurm, then:
sbatch run_daily30.slurm           # array 0-155, one month per task
```

**Per task:** 4 CPUs, about 20 GB RAM. With the compiled models and one thread, the test ran at about 1 s per row for a 31-day month (≈ 1–1.5 h per month). Extra threads speed up the model predictions but not the NumPy resampling, so expect somewhat less with 4 CPUs; `--time=06:00:00` leaves ample margin.

- The scikit-learn fallback is about 6× slower, so raise `--time` if you use it.
- `TILE_ROWS` (default 64) trades memory for speed.
- Months are independent. Re-submit failed months with `sbatch --array=<ids> run_daily30.slurm`.

## 3. Check

- Each file has one band per day. For example, `rasterio.open('OUT/2020/ENP_salinity30m_daily_2020_05.tif').count` should be 31.
- The mean of a month's bands should match the 30 m monthly map for that month. A test on May 2020, rows 1500–1531 (76,055 cells), gave mean absolute difference 0.04 psu (95th percentile 0.1 psu, r = 1.000).
- `ROWS=r0,r1` runs a row subset for quick tests, e.g. `ROWS=1500,1532 TILE_ROWS=32 python daily30_hpc.py 88 IN OUT`.
