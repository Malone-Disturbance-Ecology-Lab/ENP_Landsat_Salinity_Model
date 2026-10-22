#!/bin/bash
# One-time preparation on Bouchet. Expects ~/enp_salinity/ to hold:
#   hpc/  hpc_inputs/  data/EDEN_GRIDS_Q/  site_variables/ENP_WaterFrequency_L30.tif
# Run from ~/enp_salinity/hpc:  sbatch setup_bouchet.slurm
#SBATCH --job-name=enp_setup
#SBATCH --partition=day
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --output=setup_%j.out
set -eo pipefail
R=$HOME/enp_salinity
module load miniconda
if ! conda env list | grep -q '^sal30 '; then
  conda create -y -n sal30 -c conda-forge python=3.11 numpy pandas scikit-learn=1.8.0 rasterio pyproj scipy pyarrow
fi
conda activate sal30
pip install -q treelite==4.1.2 tl2cgen
cd $R/hpc_inputs
[ -f static30.npy ]        || { cat static30.npy.gz.part* | gunzip > static30.npy && rm static30.npy.gz.part*; }
[ -f gridmet_filled.npz ]  || { cat gridmet_filled.npz.part* > gridmet_filled.npz && rm gridmet_filled.npz.part*; }
python $R/hpc/make_grid30.py $R/site_variables/ENP_WaterFrequency_L30.tif .
# the final model is compiled on Bouchet by final_model/hpc/compile_final.slurm (submit_final.sh)
# daily EDEN fields (~4.6 GB, ~20 min) -> hpc_inputs/eden_daily/
sed "s#SERVER/#$R/#g" $R/hpc/eden_daily.py > eden_daily_local.py
cp $R/hpc/fillers.py .
python eden_daily_local.py
ls eden_daily; echo "setup done"
