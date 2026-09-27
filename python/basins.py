import numpy as np, rasterio, pandas as pd, time
from pyproj import Transformer
from skimage.graph import MCP_Geometric
E=np.load('downloads/eden_mean_grids_2025.npz')
Z=np.nanmean(np.stack([E['2025-05-15_m90'],E['2025-10-15_m90']]),0)   # (405,287), row 0 = south
valid=np.isfinite(Z); ny,nx=Z.shape
xs=463400+np.arange(nx)*400; ys=2790200+np.arange(ny)*400
tLL=Transformer.from_crs(26917,4326,always_xy=True)
XX,YY=np.meshgrid(xs,ys); LON,LAT=tLL.transform(XX,YY)
# steepest-ascent trace; basin = first inflow-structure zone the trace enters, else terminus rule
zones=[(2,25.483,-80.564,3),(2,25.418,-80.574,3),(2,25.403,-80.607,3),(2,25.437,-80.590,3),
       (3,25.330,-80.525,4),
       (1,25.7616,-80.8213,4),(1,25.7617,-80.7694,4),(1,25.7617,-80.7269,4),(1,25.7617,-80.6819,4),(1,25.7619,-80.6742,4)]
Zlab=np.zeros(Z.shape,np.int8)
for k,la,lo,rk in zones:
    d=np.hypot((LAT-la)*111,(LON-lo)*100.5); Zlab[(d<=rk)&(Zlab==0)]=k
nb=[(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]
def src(la,lo):
    if la>=25.70 and lo>-80.92: return 1
    if lo<=-80.92: return 4 if la>=25.60 else 1
    if la<25.40 and lo>-80.60: return 3
    return 2
B=np.zeros(Z.shape,np.int8); memo={}
for r0 in range(ny):
    for c0 in range(nx):
        if not valid[r0,c0]: continue
        path=[]; r,c=r0,c0; lab=0
        for _ in range(2000):
            if (r,c) in memo: lab=memo[(r,c)]; break
            path.append((r,c))
            if Zlab[r,c]>0: lab=int(Zlab[r,c]); break
            bz=Z[r,c]; br,bc=r,c
            for dr,dc in nb:
                rr,cc=r+dr,c+dc
                if 0<=rr<ny and 0<=cc<nx and valid[rr,cc] and Z[rr,cc]>bz+1e-6: bz=Z[rr,cc]; br,bc=rr,cc
            if (br,bc)==(r,c): lab=src(LAT[r,c],LON[r,c]); break
            r,c=br,bc
        for p in path: memo[p]=lab
        B[r0,c0]=lab
names={1:'SRS',2:'TaylorSlough',3:'C111',4:'BigCypress'}
print({names[k]:int((B==k).sum()) for k in names})
# map to 30 m grid
U='SERVER/site_variables/'
wr=rasterio.open(U+'ENP_WaterFrequency_L30.tif'); wf=wr.read(1); T=wr.transform
cls=rasterio.open(U+'ENP_MainlandMask.tif').read(1)
H,W=wf.shape; rr,cc=np.mgrid[0:H,0:W]
lon,lat=rasterio.transform.xy(T,rr.ravel()[::1],cc.ravel()[::1]); lon=np.array(lon).reshape(H,W); lat=np.array(lat).reshape(H,W)
tU=Transformer.from_crs(4326,26917,always_xy=True); X,Y=tU.transform(lon,lat)
ex=np.round((X-463400)/400).astype(int); ey=np.round((Y-2790200)/400).astype(int)
ok=(ex>=0)&(ex<nx)&(ey>=0)&(ey<ny)
B30=np.zeros((H,W),np.int8); B30[ok]=B[ey[ok],ex[ok]]
water=np.where(wf<-1,cls==0,wf>=0.4)
cost=np.where(water,1.0,np.where(cls==1,2.0,5.0)).astype('float32')
px=28.6; costs={}
t0=time.time()
for k in names:
    seeds=np.argwhere(B30==k)
    m=MCP_Geometric(cost,fully_connected=True); cum,_=m.find_costs(seeds[::3])
    costs[k]=(cum*px/1000).astype('float32'); print(names[k],round(time.time()-t0),flush=True)
Cst=np.stack([costs[k] for k in names]); basin=(np.argmin(Cst,0)+1).astype(np.int8)
L=10.0; Wt=np.exp(-(Cst-Cst.min(0))/L); Wt/=Wt.sum(0)
np.savez_compressed('basins_30m.npz',basin=basin,**{f'w_{names[k]}':Wt[i].astype('float32') for i,k in enumerate(names)},**{f'cost_{names[k]}':Cst[i] for i,k in enumerate(names)})
np.save('eden_basin_400m.npy',B)
# station values
st=pd.read_csv('sitevars/stations_coast.csv')
rows=[]
for _,s in st.iterrows():
    r_,c_=rasterio.transform.rowcol(T,s.longitude,s.latitude)
    rows.append(dict(station=s.station,basin=names[int(basin[r_,c_])],**{f'w_{names[k]}':float(Wt[i,r_,c_]) for i,k in enumerate(names)},**{f'basin_cost_{names[k]}_km':float(Cst[i,r_,c_]) for i,k in enumerate(names)}))
S=pd.DataFrame(rows).set_index('station'); S.round(3).to_csv('station_basins.csv'); print(S.basin.value_counts().to_dict())
old=pd.read_csv('station_connectivity_with_region.csv',index_col=0).region
print(pd.crosstab(old,S.basin))
