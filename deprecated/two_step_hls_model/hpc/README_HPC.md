# Daily 30 m salinity on HPC (Bouchet)

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

**On Bouchet:** after copying the folders, run `cd ~/enp_salinity/hpc && sbatch setup_bouchet.slurm`. It builds the `sal30` conda env, rebuilds the split files, makes the grid index, checks the compiled models and builds the daily EDEN fields (the steps below). `eden_daily.py` and `fillers.py` are copied into `hpc/` for this.

Copy `best_two_step/salinity30m/hpc/`, `hpc_inputs/`, `data/EDEN_GRIDS_Q/` and `site_variables/ENP_WaterFrequency_L30.tif` from the lab server to Bouchet (`transfer-bouchet.ycrc.yale.edu`, on campus or VPN), e.g. into `~/enp_salinity/`. Home is 125 GiB, enough for the inputs (~11 GB with daily EDEN) and outputs (~15 GB); `~/scratch_pi_<netid>` is purged after 30 days. It contains:

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

## 4. Daily 30 m hindcast, 2000-04 → 2012-12 (one-step model, no spectra)

`daily30_onestep_hpc.py` + `run_daily30_onestep.slurm` run the one-step model m3 (no vegetation spectra) at every 30 m cell for every day. Array index 0 = 2000-04 … 152 = 2012-12. Indices 153–308 (2013-01 … 2025-12) are also allowed, which gives a one-step daily record for 2000–2025 with no model switch at 2013.

- Output: `daily30m_onestep/YYYY/ENP_salinity30m_onestep_daily_YYYY_MM.tif`. Each file has one band per day; values are int16 psu×10; −1 = no data.
- Inputs: the same `hpc_inputs/` as above, plus three files from `dbhydro_hindcast/`:
  - `features/gridmet_filled_1998_2013.npz`
  - `features/daily_drivers_1998_2025.csv`
  - `eden_anyv.npy`

  The EDEN quarterly files for 2000–2012 also go in `data/EDEN_GRIDS_Q/`. Months from 2013 on use `gridmet_filled.npz`.
- Cost: a test on May 2005 rows 2000–2040 (89,123 cells × 31 days) took 28 s of prediction on 2 threads. A full month is therefore about 40 min on 2 threads, and less with 4. Each task has 4 CPUs, 24 GB and a 4 h limit.
- Check (May 2005, rows 1500–1563): at the 250 m nodes the daily 30 m values equal the 250 m hindcast maps exactly (difference 0.0 psu). On the 15th they match the stored 30 m full-model run within rounding (mean 0.02 psu).
