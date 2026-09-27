import numpy as np, rasterio, pandas as pd
from scipy import ndimage as ndi
sd=rasterio.open('sitevars/ENP_DistCoast_signed.tif'); dist=sd.read(1)/1000; T=sd.transform
C=np.load('connectivity_rasters.npz'); B=np.load('basins_30m.npz')
names=['SRS','TaylorSlough','C111','BigCypress']
fields={f'src_dist_{k}_km':C[f'src_dist_{k}_km'] for k in ['SRS','TaylorSlough','C111','BigCypressGulf']}
fields.update({f'cost_{n}':B[f'cost_{n}'] for n in names}); fields.update({f'w_{n}':B[f'w_{n}'] for n in names})
off=dist<0; alpha=np.clip(-dist/3.0,0,1).astype('float32')     # 0 at shoreline -> 1 at 1 km offshore
sig=150  # px (~4.3 km)
m=off.astype('float32'); den=ndi.gaussian_filter(m,sig)
out={}
for k,v in fields.items():
    v=v.astype('float32'); v=np.where(np.isfinite(v),v,0)
    sm=ndi.gaussian_filter(v*m,sig)/np.maximum(den,1e-6)
    out[k]=np.where(off,(1-alpha)*v+alpha*sm,v).astype('float32'); print(k,flush=True)
S=sum(out[f'w_{n}'] for n in names)
for n in names: out[f'w_{n}']=out[f'w_{n}']/np.maximum(S,1e-6)
out['basin']=(np.argmax(np.stack([out[f'w_{n}'] for n in names]),0)+1).astype('int8')
np.savez_compressed('offshore_smoothed_v3.npz',**out)
st=pd.read_csv('sitevars/stations_coast.csv'); rows=[]
for _,s in st.iterrows():
    r_,c_=rasterio.transform.rowcol(T,s.longitude,s.latitude)
    d=dict(station=s.station,basin=names[int(out['basin'][r_,c_])-1])
    for n in names: d[f'w_{n}']=float(out[f'w_{n}'][r_,c_]); d[f'basin_cost_{n}_km']=float(out[f'cost_{n}'][r_,c_])
    for k in ['SRS','TaylorSlough','C111','BigCypressGulf']: d[f'src_dist_{k}_km']=float(out[f'src_dist_{k}_km'][r_,c_])
    rows.append(d)
pd.DataFrame(rows).set_index('station').round(3).to_csv('station_basins_v3.csv'); print('done')
