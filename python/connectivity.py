import numpy as np, rasterio, pandas as pd, time
from scipy import ndimage as ndi
from skimage.graph import MCP_Geometric
U='SERVER/site_variables/'
r=rasterio.open(U+'ENP_WaterFrequency_L30.tif'); wf=r.read(1); T=r.transform; prof=r.profile
cls=rasterio.open(U+'ENP_MainlandMask.tif').read(1)       # 0 ocean, 1 mainland, 2 islands/sea-side
nod=wf<-1
water=np.where(nod, cls==0, wf>=0.4)
ocean=(cls==0)&water
dy=r.res[1]*110574; dx=r.res[0]*111320*np.cos(np.radians(25.37)); px=(dx+dy)/2
print('grid',wf.shape,'px m',round(px,1),'water frac',water.mean().round(3))
out={}
t0=time.time()
for nm,landcost in [('strict',1000.0),('soft',10.0)]:
    cost=np.where(water,1.0,landcost).astype('float32')
    m=MCP_Geometric(cost,fully_connected=True)
    starts=np.argwhere(ocean)
    cum,_=m.find_costs(starts[::1])
    out[f'chan_dist_{nm}_km']=(cum*px/1000).astype('float32'); print(nm,'done',round(time.time()-t0))
# narrow channels: water that disappears after opening with ~200 m structuring element, inside mainland/islands (not ocean)
wo=ndi.binary_opening(water,structure=np.ones((3,3)),iterations=3)
narrow=water&~wo&~ocean
inland_water=water&~ocean
def density(mask,rad_m):
    r_=int(round(rad_m/px)); yy,xx=np.mgrid[-r_:r_+1,-r_:r_+1]; k=((yy**2+xx**2)<=r_**2).astype('float32'); k/=k.sum()
    return ndi.convolve(mask.astype('float32'),k,mode='constant').astype('float32')
out['chan_density_500m']=density(narrow,500); out['chan_density_1km']=density(narrow,1000); out['inland_water_1km']=density(inland_water,1000)
print('density done',round(time.time()-t0))
# source distances through mainland (cost 1 on mainland, very high elsewhere)
seeds={'SRS':[(25.7616,-80.8213),(25.7617,-80.7694),(25.7617,-80.7269),(25.7617,-80.6819),(25.7619,-80.6742)],
       'TaylorSlough':[(25.4831,-80.5642),(25.4181,-80.5739),(25.4030,-80.6073)],
       'C111':[(25.3303,-80.5253)],
       'BigCypressGulf':[(25.8873,-81.2618),(25.8514,-80.9806)]}
cost=np.where(cls==1,1.0,50.0).astype('float32')
for k,pts in seeds.items():
    rc=[rasterio.transform.rowcol(T,lo,la) for la,lo in pts]
    rc=[(min(max(a,0),wf.shape[0]-1),min(max(b,0),wf.shape[1]-1)) for a,b in rc]
    m=MCP_Geometric(cost,fully_connected=True); cum,_=m.find_costs(rc)
    out[f'src_dist_{k}_km']=(cum*px/1000).astype('float32'); print(k,'done',round(time.time()-t0))
np.savez_compressed('connectivity_rasters.npz',**out)
# station extraction (min over 3x3 for channel distance: stations sit in creeks)
st=pd.read_csv('sitevars/stations_coast.csv')
rows=[]
for _,s in st.iterrows():
    rr,cc=rasterio.transform.rowcol(T,s.longitude,s.latitude); d=dict(station=s.station)
    for k,v in out.items():
        w=v[max(rr-3,0):rr+4,max(cc-3,0):cc+4]
        d[k]=float(np.nanmin(w)) if k.startswith('chan_dist') else float(v[rr,cc])
    rows.append(d)
S=pd.DataFrame(rows).set_index('station')
S['src_nearest']=S[[c for c in S if c.startswith('src_dist')]].idxmin(axis=1).str.replace('src_dist_','').str.replace('_km','')
S['euclid_coast_km']=st.set_index('station').distCoast_signed_m.clip(lower=0)/1000
S['chan_sinuosity']=S.chan_dist_soft_km/np.maximum(S.euclid_coast_km,0.1)
S.round(3).to_csv('station_connectivity.csv'); print(S.round(2).to_string())
