"""Daily EDEN stage fields for mapping, 2013-01-01..2025-12-31, matching the station features:
cell (daily), cell30, R365 of cell (centered), 5km30, 5km90, R365 of 5km.
Outside the EDEN domain the fields are extended by harmonic (Laplace) extrapolation (as in v3),
then cropped to the map domain. Stored as float16 memmaps in eden_daily/."""
import numpy as np, scipy.io as sio, os, glob, time, rasterio
from scipy import ndimage as ndi
from scipy.signal import fftconvolve
from pyproj import Transformer
from fillers import laplace_solver, nearest_inside
OUT='eden_daily'; os.makedirs(OUT,exist_ok=True)
files=sorted(glob.glob('SERVER/data/EDEN_GRIDS_Q/20[12][0-9]_q[1-4].nc')); files=[f for f in files if 2013<=int(os.path.basename(f)[:4])<=2025]
assert len(files)==52, len(files)
NT=4748; H,W=405,287
# ---- pass 0: daily cell and 5-km fields on the native grid
t0=time.time(); anyv=np.zeros((H,W),bool); stacks=[]
for f in files[:8]:
    s=sio.netcdf_file(f,'r',mmap=False).variables['stage'][:].astype('float32'); s[(s<-900)|(s>9000)]=np.nan; anyv|=np.isfinite(s).any(0)
r_=12; yy,xx=np.mgrid[-r_:r_+1,-r_:r_+1]; kk=((yy**2+xx**2)<=12.5**2).astype('float32')
reg5=ndi.convolve(anyv.astype('float32'),kk,mode='constant')>0
den5=ndi.convolve(anyv.astype('float32'),kk,mode='constant')
CELL=np.lib.format.open_memmap(f'{OUT}/_cell_native.npy','w+','float16',(NT,H,W)); FIVE=np.lib.format.open_memmap(f'{OUT}/_5km_native.npy','w+','float16',(NT,H,W))
t=0
for f in files:
    s=sio.netcdf_file(f,'r',mmap=False).variables['stage'][:].astype('float32'); s[(s<-900)|(s>9000)]=np.nan
    for d in range(s.shape[0]): s[d]=nearest_inside(s[d],anyv)
    num=fftconvolve(np.nan_to_num(s),kk[None],mode='same',axes=(1,2))
    five=np.where(reg5[None],num/np.maximum(den5,1)[None],np.nan)
    n=s.shape[0]; CELL[t:t+n]=s; FIVE[t:t+n]=five; t+=n
print('pass0 days',t,round(time.time()-t0),'s',flush=True); assert t==NT
# ---- pass A: temporal windows per row block (trailing 30/90-day means; centered 365-day mean)
def trailing(x,w,minp):   # x (T,...) float32 with NaN
    v=np.isfinite(x); c=np.cumsum(np.where(v,x,0),0,dtype='float64'); n=np.cumsum(v,0)
    cs=np.concatenate([np.zeros((1,)+x.shape[1:]),c]); ns=np.concatenate([np.zeros((1,)+x.shape[1:]),n])
    i=np.arange(1,len(x)+1); j=np.maximum(i-w,0)
    S=cs[i]-cs[j]; N=ns[i]-ns[j]; return np.where(N>=minp,S/np.maximum(N,1),np.nan).astype('float32')
def centered(x,w,minp):
    v=np.isfinite(x); c=np.cumsum(np.where(v,x,0),0,dtype='float64'); n=np.cumsum(v,0)
    cs=np.concatenate([np.zeros((1,)+x.shape[1:]),c]); ns=np.concatenate([np.zeros((1,)+x.shape[1:]),n])
    T=len(x); i=np.arange(T); hi=np.minimum(i+w//2+1,T); lo=np.maximum(i-w//2,0)
    S=cs[hi]-cs[lo]; N=ns[hi]-ns[lo]; return np.where(N>=minp,S/np.maximum(N,1),np.nan).astype('float32')
NAMES=['cell','cell30','Rcell','5km30','5km90','R5']
NAT={k:np.lib.format.open_memmap(f'{OUT}/_{k}_nat.npy','w+','float16',(NT,H,W)) for k in NAMES}
for r0 in range(0,H,45):
    c=np.array(CELL[:,r0:r0+45],dtype='float32'); f5=np.array(FIVE[:,r0:r0+45],dtype='float32')
    NAT['cell'][:,r0:r0+45]=c; NAT['cell30'][:,r0:r0+45]=trailing(c,30,10); NAT['Rcell'][:,r0:r0+45]=centered(c,365,60)
    NAT['5km30'][:,r0:r0+45]=trailing(f5,30,10); NAT['5km90'][:,r0:r0+45]=trailing(f5,90,30); NAT['R5'][:,r0:r0+45]=centered(f5,365,60)
print('passA',round(time.time()-t0),'s',flush=True)
del CELL,FIVE; os.remove(f'{OUT}/_cell_native.npy'); os.remove(f'{OUT}/_5km_native.npy')
# ---- pass B: harmonic extrapolation on the padded grid, crop to map domain
PS,PW,PE,PN=140,70,10,10; y0p=2790200-PS*400; x0p=463400-PW*400
def pad(g): return np.pad(g,((0,0),(PS,PN),(PW,PE)),constant_values=np.nan)
fill_c=laplace_solver(np.pad(anyv,((PS,PN),(PW,PE)))); fill_5=laplace_solver(np.pad(reg5,((PS,PN),(PW,PE))))
wf=rasterio.open('SERVER/site_variables/ENP_WaterFrequency_L30.tif'); b=wf.bounds
tU=Transformer.from_crs(4326,26917,always_xy=True)
xs,ys=tU.transform([b.left,b.left,b.right,b.right],[b.bottom,b.top,b.bottom,b.top])
c0=int(np.floor((min(xs)-x0p)/400))-2; c1=int(np.ceil((max(xs)-x0p)/400))+3; r0=int(np.floor((min(ys)-y0p)/400))-2; r1=int(np.ceil((max(ys)-y0p)/400))+3
c0,r0=max(c0,0),max(r0,0); c1=min(c1,W+PW+PE); r1=min(r1,H+PS+PN)
np.save(f'{OUT}/crop_meta.npy',np.array([y0p+r0*400,x0p+c0*400,r1-r0,c1-c0],dtype='float64'))
print('crop rows',r0,r1,'cols',c0,c1,flush=True)
for k in NAMES:
    O=np.lib.format.open_memmap(f'{OUT}/{k}.npy','w+','float16',(NT,r1-r0,c1-c0)); fill=fill_5 if k in ('5km30','5km90','R5') else fill_c
    for d0 in range(0,NT,200):
        g=pad(np.array(NAT[k][d0:d0+200],dtype='float32')); ok=np.isfinite(g).any(axis=(1,2))
        known=fill.__closure__ if False else None
        gg=np.where(np.isfinite(g),g,0)
        # days where the field is undefined (start of record for windows) stay NaN
        res=fill(gg).astype('float32'); res[~ok]=np.nan
        # cells inside the known mask that are NaN for that day (window warm-up) -> NaN
        nanin=~np.isfinite(g)&np.isfinite(res)&(np.pad(anyv if fill is fill_c else reg5,((PS,PN),(PW,PE)))[None])
        res[nanin]=np.nan
        O[d0:d0+200]=res[:,r0:r1,c0:c1]
    O.flush(); print('passB',k,round(time.time()-t0),'s',flush=True)
for k in NAMES: os.remove(f'{OUT}/_{k}_nat.npy')
print('done')
