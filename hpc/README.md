# HPC setup (Bouchet, Yale)

One-time preparation for the daily 30 m runs. The model jobs themselves are in [`../final_model/hpc/`](../final_model/hpc/).

- Copy to `~/enp_salinity/`: `hpc/` (these scripts plus `final_model/hpc/*`, `python/eden_daily.py` and `python/fillers.py`), `hpc_inputs/` (from `best_two_step/salinity30m/hpc_inputs/` on the lab server, plus `final_model_onestep.pkl`, `msl_by_date.csv` and the pre-2013 inputs from `dbhydro_hindcast/features/`), `data/EDEN_GRIDS_Q/` and `site_variables/ENP_WaterFrequency_L30.tif`.
- `make_hpc_inputs.py`: filled gridMET fields and daily global drivers (run in the pipeline working directory).
- `make_grid30.py`: the 30 m grid index from the water-frequency raster.
- `setup_bouchet.slurm`: builds the `sal30` conda environment (Python 3.11, scikit-learn, rasterio, treelite 4.1.2, tl2cgen), unpacks the inputs and builds the daily EDEN fields (about 20 min). Then run `final_model/hpc/submit_final.sh`.
