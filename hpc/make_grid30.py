"""Regenerate px30_index.npy, lonlat30.npy, utm30.npy (30 m HLS grid cells inside the domain) from ENP_WaterFrequency_L30.tif.
Usage: python make_grid30.py <path/to/ENP_WaterFrequency_L30.tif> <input_dir>"""
import sys, numpy as np, rasterio
from pyproj import Transformer
wfr=rasterio.open(sys.argv[1]); O=sys.argv[2]; H,W=wfr.height,wfr.width; T=wfr.transform
idx=np.flatnonzero((wfr.read(1)>-1).ravel()); np.save(f'{O}/px30_index.npy',idx.astype(np.int64))
rr,cc=np.unravel_index(idx,(H,W)); lon,lat=rasterio.transform.xy(T,rr,cc); lon=np.asarray(lon,'float64'); lat=np.asarray(lat,'float64')
np.save(f'{O}/lonlat30.npy',np.stack([lon,lat]))
X,Y=Transformer.from_crs(4326,26917,always_xy=True).transform(lon,lat); np.save(f'{O}/utm30.npy',np.stack([X,Y]).astype('float64'))
print('inside cells',len(idx))
