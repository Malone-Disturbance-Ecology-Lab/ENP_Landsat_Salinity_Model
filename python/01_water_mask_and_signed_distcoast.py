# Step 1a (ran on the Mac-side VM because the Fmask stack is 1.5 GB): water frequency from the L30 Fmask stack.
# The stack is pixel-interleaved, so read row blocks across all scenes.
import rasterio, numpy as np, os, shapefile
from rasterio import features
from rasterio.windows import Window
from pyproj import Transformer
from scipy import ndimage as ndi

PROJ = '/Volumes/malonelab/Research/ENP_Landsat_Salinity_Model'
FMASK = f'{PROJ}/harmonized_appeears_landsat_data/L30/L30_ENP_Fmask.tif'
FLORIDA = '/Volumes/malonelab/Research/ENP/shapefiles/Florida'   # NAD83 UTM17N land polygon
OUT = f'{PROJ}/site_variables'

r = rasterio.open(FMASK); H, W = r.height, r.width
clear = np.zeros((H, W), np.uint16); water = np.zeros_like(clear)
for row in range(0, H, 100):
    n = min(100, H - row); a = r.read(window=Window(0, row, W, n))
    valid = (a != 255) & ((a & 31) == 0)            # no cirrus/cloud/adjacent/shadow/snow
    clear[row:row+n] = valid.sum(0); water[row:row+n] = (valid & ((a & 32) > 0)).sum(0)
wf = np.where(clear >= 20, water / np.maximum(clear, 1), np.nan)

# Florida land polygon for pixels outside the Landsat footprint (parity rasterization handles holes)
tr = Transformer.from_crs(26917, 4326, always_xy=True)
geoms = []
for shp in shapefile.Reader(FLORIDA).shapes():
    pts = np.array(shp.points); parts = list(shp.parts) + [len(pts)]
    for i in range(len(parts) - 1):
        p = pts[parts[i]:parts[i+1]]
        if len(p) < 4: continue
        x, y = tr.transform(p[:, 0], p[:, 1])
        if x.max() < r.bounds.left or x.min() > r.bounds.right or y.max() < r.bounds.bottom or y.min() > r.bounds.top: continue
        geoms.append(({'type': 'Polygon', 'coordinates': [list(zip(x, y))]}, 1))
fl_land = (features.rasterize(geoms, out_shape=(H, W), transform=r.transform, fill=0, dtype='uint16',
                              merge_alg=features.MergeAlg.add) % 2) == 1

nod = np.isnan(wf)
wet = np.where(nod, ~fl_land, wf >= 0.5)
wet_open = ndi.binary_opening(wet, structure=ndi.generate_binary_structure(2, 1), iterations=5)
lab, n = ndi.label(wet_open); ocean = lab == (np.argmax(ndi.sum(wet_open, lab, range(1, n+1))) + 1)
lab2, n2 = ndi.label(~ocean, structure=np.ones((3, 3)))
mainland = lab2 == (np.argmax(ndi.sum(~ocean, lab2, range(1, n2+1))) + 1)

lat = (r.bounds.top + r.bounds.bottom) / 2
dy = r.res[1] * 110574; dx = r.res[0] * 111320 * np.cos(np.radians(lat))
sd = np.where(mainland, ndi.distance_transform_edt(mainland, sampling=(dy, dx)),
              -ndi.distance_transform_edt(~mainland, sampling=(dy, dx))).astype('float32')

p = r.profile.copy(); p.update(count=1, interleave='band', tiled=True, blockxsize=256, blockysize=256)
with rasterio.open(f'{OUT}/ENP_DistCoast_signed.tif', 'w', **{**p, 'dtype': 'float32', 'nodata': -99999, 'compress': 'deflate', 'predictor': 3}) as o: o.write(sd, 1)
with rasterio.open(f'{OUT}/ENP_WaterFrequency_L30.tif', 'w', **{**p, 'dtype': 'float32', 'nodata': -9999, 'compress': 'lzw'}) as o: o.write(np.where(nod, -9999, wf).astype('float32'), 1)
cls = np.where(mainland, 1, np.where(ocean, 0, 2)).astype('uint8')
with rasterio.open(f'{OUT}/ENP_MainlandMask.tif', 'w', **{**p, 'dtype': 'uint8', 'nodata': 255, 'compress': 'lzw'}) as o: o.write(cls, 1)
