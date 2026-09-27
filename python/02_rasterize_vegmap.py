import rasterio, numpy as np, pandas as pd, shapefile, glob
from rasterio import features
from pyproj import Transformer
D='/Volumes/malonelab/Research/ENP/VegetationMap/Data/'; U='/Volumes/malonelab/Research/ENP_Landsat_Salinity_Model/'
m=rasterio.open(U+'site_variables/ENP_MainlandMask.tif'); T=m.transform; shape=m.shape
tr=Transformer.from_crs(26917,4326,always_xy=True)
def cls_of(layer,code,name):
    if layer=='Other':
        if code=='ONW': return 'Open water'
        if code in ('ONT','ONS','ONB'): return 'Tidal/salt flat, beach'
        return 'Anthropogenic/other'
    if layer in ('MangroveForest','MangroveScrub','MangroveShrubland'): return 'Mangrove'
    if layer=='Woodlands':
        if code.startswith('WM'): return 'Mangrove'
        if code.startswith('WU'): return 'Upland'
        return 'Swamp/cypress woodland'
    if layer=='SaltMarsh': return 'Salt marsh'
    if layer in ('FreshwaterMarsh','Prairie'): return 'Freshwater marsh/prairie'
    if layer in ('SwampScrub','SwampShrubland'): return 'Swamp scrub/shrubland'
    if layer in ('UplandForest','UplandScrub','UplandShrubland','Dune'): return 'Upland'
    if layer=='InvasiveSpecies': return 'Exotic'
codes=[]; ras=np.zeros(shape,np.uint16)
for f in sorted(glob.glob(D+'*.shp')):
    layer=f.split('/')[-1][:-4]; s=shapefile.Reader(f[:-4])
    for shp,rec in zip(s.shapes(),s.records()):
        codes.append(dict(id=len(codes)+1,layer=layer,VegCode=rec[0],name=rec[3],cls=cls_of(layer,rec[0],rec[3])))
        pts=np.array(shp.points); parts=list(shp.parts)+[len(pts)]
        rings=[]
        for i in range(len(parts)-1):
            p=pts[parts[i]:parts[i+1]]
            if len(p)<4: continue
            x,y=tr.transform(p[:,0],p[:,1]); rings.append(list(zip(x,y)))
        # parity rasterization per record (handles holes/multipart)
        if not rings: continue
        cnt=features.rasterize([({'type':'Polygon','coordinates':[rg]},1) for rg in rings],out_shape=shape,transform=T,fill=0,dtype='uint16',merge_alg=features.MergeAlg.add)
        ras[(cnt%2)==1]=len(codes)
    print(layer,len(s))
cd=pd.DataFrame(codes); cd.to_csv(U+'site_variables/vegmap_codes.csv',index=False)
p=m.profile.copy(); p.update(dtype='uint16',nodata=0,compress='deflate',tiled=True,blockxsize=256,blockysize=256)
with rasterio.open(U+'site_variables/ENP_VegMap_30m.tif','w',**p) as o: o.write(ras,1)
print('mapped frac',(ras>0).mean().round(3))
