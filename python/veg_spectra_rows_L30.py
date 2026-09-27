import rasterio, numpy as np, sys, os, datetime, time, warnings; warnings.filterwarnings('ignore')
from rasterio.windows import Window
k,nw=int(sys.argv[1]),int(sys.argv[2])
D=os.path.expanduser('SERVER/harmonized_appeears_landsat_data/L30/')
WF=os.path.expanduser('SERVER/site_variables/ENP_WaterFrequency_L30.tif')
OUT=os.path.expanduser('veg_rows'); os.makedirs(OUT,exist_ok=True)
R={b:rasterio.open(D+f'L30_ENP_{b}.tif') for b in ['B04','B05','B06','Fmask']}
desc=R['B04'].descriptions
dts=[datetime.date(int(s.split('doy')[1][:4]),1,1)+datetime.timedelta(int(s.split('doy')[1][4:7])-1) for s in desc]
mi=np.array([(d.year-2013)*12+d.month-1 for d in dts]); NM=156
wf=rasterio.open(WF); W=4171; K=8; nb=W//K
rows=[r for r in range(K//2,3839,2*K)]; rows=rows[k::nw]
for r in rows:
    fn=f'{OUT}/row_{r:04d}.npz'
    if os.path.exists(fn): continue
    t=time.time()
    a={b:R[b].read(window=Window(0,r,W,1))[:,0,:] for b in R}
    wfr=wf.read(1,window=Window(0,r,W,1))[0]
    veg=(wfr>=0)&(wfr<0.2)
    f=a['Fmask'].astype(np.int32); clear=((f&31)==0)&(f!=255)
    Rd,N,S1=[a[b].astype('float32') for b in ['B04','B05','B06']]; G=np.full_like(N,np.nan); S2=G
    ok=clear&np.isfinite(N)&np.isfinite(Rd)&(N>0)&(N<1)&veg[None,:]
    with np.errstate(all='ignore'):
        I=np.stack([(N-Rd)/(N+Rd),2.5*(N-Rd)/(N+2.4*Rd+1),(N-S1)/(N+S1),(S1-S2)/(S1+S2),(G-N)/(G+N),S1])
    I[:,~ok]=np.nan
    # block mean per obs (require >=25% of veg pixels clear)
    Ib=I[:,:,:nb*K].reshape(6,len(mi),nb,K); okb=ok[:,:nb*K].reshape(len(mi),nb,K)
    nveg=veg[:nb*K].reshape(nb,K).sum(1)
    cnt=okb.sum(2); bm=np.nanmean(Ib,3)
    bm[:,cnt<np.maximum(1,0.25*nveg)]=np.nan
    M=np.full((NM,6,nb),np.nan,'float32'); C=np.zeros((NM,nb),'int16')
    for m in range(NM):
        s=mi==m
        if s.any():
            M[m]=np.nanmedian(bm[:,s,:],1); C[m]=np.isfinite(bm[0,s,:]).sum(0)
    np.savez_compressed(fn+'.tmp.npz',M=M.astype('float16'),C=C,nveg=nveg.astype('int16'))
    os.replace(fn+'.tmp.npz',fn); print(r,round(time.time()-t,1),flush=True)
print('worker done',k)
