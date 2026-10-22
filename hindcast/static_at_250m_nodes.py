"""Static predictors (static30) at the 250 m map nodes = 30 m cells (8i+4, 8j+4) inside the domain (119,019 cells).
Output: nodes250_static.npz (row250, col250, lon, lat, X, Y (UTM17N), and the 27 static layers)."""
import numpy as np, rasterio, json, glob, gzip, sys
from pyproj import Transformer
R = sys.argv[1]; OUT = sys.argv[2]; HPC = f'{R}/best_two_step/salinity30m/hpc_inputs'
W30 = 4171; D30 = 0.00027092271166751034; LEFT, TOP = -81.51874748177192, 25.889645249658873
wf = rasterio.open(f'{R}/site_variables/ENP_WaterFrequency_L30.tif').read(1); dom = wf > -1
idx = np.flatnonzero(dom.ravel())
ii, jj = np.nonzero(dom[4::8, 4::8][:480, :521]); rr, cc = ii * 8 + 4, jj * 8 + 4
pos = np.searchsorted(idx, rr * W30 + cc); assert (idx[pos] == rr * W30 + cc).all()
keys = json.load(open(f'{HPC}/static30_keys.json')); NP = len(idx)
class Cat:
    def __init__(s, fs): s.fs = iter(fs); s.f = open(next(s.fs), 'rb')
    def read(s, n=-1):
        b = s.f.read(n)
        while n < 0 or len(b) < n:
            try: s.f = open(next(s.fs), 'rb')
            except StopIteration: break
            b += s.f.read(-1 if n < 0 else n - len(b))
        return b
g = gzip.GzipFile(fileobj=Cat(sorted(glob.glob(f'{HPC}/static30.npy.gz.part*')))); ver = np.lib.format.read_magic(g); shp, fo, dt = np.lib.format._read_array_header(g, ver)
out = {}
for k in keys:
    a = np.frombuffer(g.read(NP * 4), dtype=dt, count=NP); out[k] = a[pos].astype('float32')
lon = LEFT + (cc + 0.5) * D30; lat = TOP - (rr + 0.5) * D30
X, Y = Transformer.from_crs(4326, 26917, always_xy=True).transform(lon, lat)
np.savez_compressed(OUT, row250=ii, col250=jj, lon=lon, lat=lat, X=X, Y=Y, keys=np.array(keys), **out)
print('nodes', len(ii), {k: round(float(np.nanmean(v)), 2) for k, v in list(out.items())[:4]})
