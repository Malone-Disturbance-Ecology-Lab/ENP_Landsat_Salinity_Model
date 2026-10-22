"""High-salinity events at 30 m from daily 30 m salinity reconstructed on the fly:
daily30(t) = bilinear(250 m daily map, t) + C(t)  (C = 30 m full-model detail, linear in time between mid-month runs).
Event = >= 7 consecutive days above the pixel's own 95th percentile. Also days > 40 psu, and daily series at stations."""
import numpy as np, pandas as pd, rasterio, json, time, os
t0=time.time(); H,W=3839,4171; NM=156; D='change30'; os.makedirs(D,exist_ok=True)
MED=np.load('daily_maps/median_int16.npy',mmap_mode='r'); C=np.load('map30/C_int16.npy',mmap_mode='r'); idx=np.load('px30_index.npy')
eco=np.load(f'{D}/eco30.npy'); reg=np.load(f'{D}/region30.npy'); nE,nR=7,6
edates=pd.date_range('2013-01-01','2025-12-31'); T=len(edates); yrs=edates.year.values-2013; mon=edates.month.values
months=pd.date_range('2013-01-01',periods=NM,freq='MS'); midn=(months+pd.Timedelta(days=14)).values.astype('datetime64[D]').astype(float)
dn=edates.values.astype('datetime64[D]').astype(float); jj=np.clip(np.searchsorted(midn,dn)-1,0,NM-2); ff=np.clip((dn-midn[jj])/(midn[jj+1]-midn[jj]),0,1).astype('float32')[:,None]
# station pixels
wfr=rasterio.open('/mnt/user-data/uploads/ENP_Landsat_Salinity_Model/site_variables/ENP_WaterFrequency_L30.tif'); Tr=wfr.transform
st=pd.read_csv('sitevars/stations_coast.csv'); spix={}
for _,s in st.iterrows():
    r,c=rasterio.transform.rowcol(Tr,s.longitude,s.latitude); spix[r*W+c]=s.station
SSER={}
# accumulators
PB=np.arange(0,80.01,0.1); XB=np.arange(-5,60.01,0.1); DB=np.arange(7,400)
hist={g:{k:np.zeros((n,len(b)-1)) for k,b in (('peak',PB),('excess',XB),('dur',DB))} for g,n in (('reg',nR),('eco',nE))}
ev_year={g:np.zeros((n,13)) for g,n in (('reg',nR),('eco',nE))}; evd_year={g:np.zeros((n,13)) for g,n in (('reg',nR),('eco',nE))}
px_event_year={g:np.zeros((n,13)) for g,n in (('reg',nR),('eco',nE))}; hyp_year={g:np.zeros((n,13)) for g,n in (('reg',nR),('eco',nE))}
hyp30={g:np.zeros((n,13)) for g,n in (('reg',nR),('eco',nE))}; startmon={g:np.zeros((n,12)) for g,n in (('reg',nR),('eco',nE))}
npx={'reg':np.bincount(reg[reg>=0],minlength=nR).astype(float),'eco':np.bincount(eco,minlength=nE).astype(float)}
NP=len(idx); P95=np.full(NP,np.nan,'float32'); P50=P95.copy(); EXC=P95.copy(); EVY=P95.copy(); HYP=P95.copy()
TR=6
for r0 in range(0,H,TR):
    r1=min(r0+TR,H); a,b=np.searchsorted(idx,[r0*W,r1*W]); n=b-a
    if n==0: continue
    p=idx[a:b]; rr,cc=p//W,p%W
    fi=(rr-4)/8.0; iy=np.clip(np.floor(fi).astype(int),0,478); fy=np.clip(fi-iy,0,1).astype('float32'); fj=(cc-4)/8.0; ix=np.clip(np.floor(fj).astype(int),0,519); fx=np.clip(fj-ix,0,1).astype('float32')
    ia=iy.min(); ib=iy.max()+2; G=np.asarray(MED[:,ia:ib,:]).astype('float32')/10; G[G<0]=np.nan; iy2=iy-ia
    num=np.zeros((T,n),'float32'); den=np.zeros((T,n),'float32')
    for (dy,dx,w) in ((0,0,(1-fy)*(1-fx)),(0,1,(1-fy)*fx),(1,0,fy*(1-fx)),(1,1,fy*fx)):
        v=G[:,iy2+dy,ix+dx]; ok=np.isfinite(v); num+=np.where(ok,v*w,0); den+=np.where(ok,w,0)
    X=np.where(den>0,num/np.maximum(den,1e-9),np.nan); del num,den,G
    Ct=np.asarray(C[:,rr,cc]).astype('float32'); Ct=np.where(Ct==-32768,0,Ct/10)
    X+=Ct[jj]*(1-ff)+Ct[jj+1]*ff; np.clip(X,0,None,out=X); del Ct
    for q,s in enumerate(p):
        if s in spix: SSER[spix[s]]=X[:,q].copy()
    q95=np.nanpercentile(X,95,axis=0); q50=np.nanmedian(X,axis=0); P95[a:b]=q95; P50[a:b]=q50
    hyp=np.stack([(X[yrs==y]>40).sum(0) for y in range(13)]); HYP[a:b]=hyp.mean(0)
    ex=(X>q95[None]).astype(np.int8); ex[np.isnan(X)]=0
    d=np.diff(np.vstack([np.zeros((1,n),np.int8),ex,np.zeros((1,n),np.int8)]),axis=0)
    sp,stt=np.nonzero(d.T==1); _,ent=np.nonzero(d.T==-1); dur=ent-stt; k=dur>=7; sp,stt,ent,dur=sp[k],stt[k],ent[k],dur[k]
    ec=eco[a:b]; rc=reg[a:b]
    if len(sp):
        AT=np.append(np.nan_to_num(X.T,nan=-1).ravel(),-1); off=sp*T; ii=np.empty(2*len(sp),np.int64); ii[0::2]=off+stt; ii[1::2]=off+ent
        pk=np.maximum.reduceat(AT,ii)[0::2]; exs=pk-q50[sp]; y0=yrs[stt]; m0=mon[stt]-1
        EXC[a:b]=pd.Series(exs).groupby(sp).mean().reindex(range(n)).values; EVY[a:b]=np.bincount(sp,minlength=n)/13.0
        for g,cat,nn in (('reg',rc,nR),('eco',ec,nE)):
            cg=cat[sp]; okg=cg>=0
            for kk,vals,bins in (('peak',pk,PB),('excess',exs,XB),('dur',dur,DB)):
                h,_,_=np.histogram2d(cg[okg],vals[okg],bins=[np.arange(nn+1),bins]); hist[g][kk]+=h
            np.add.at(ev_year[g],(cg[okg],y0[okg]),1); np.add.at(evd_year[g],(cg[okg],y0[okg]),dur[okg]); np.add.at(startmon[g],(cg[okg],m0[okg]),1)
            pe=np.zeros((n,13),bool); pe[sp,y0]=True
            for c_ in range(nn):
                mm=cat==c_
                if mm.any(): px_event_year[g][c_]+=pe[mm].sum(0)
    else:
        EVY[a:b]=0
    for g,cat,nn in (('reg',rc,nR),('eco',ec,nE)):
        for c_ in range(nn):
            mm=cat==c_
            if mm.any(): hyp_year[g][c_]+=hyp[:,mm].sum(1); hyp30[g][c_]+=(hyp[:,mm]>=30).sum(1)
    if r0%480==0: print('rows',r0,round(time.time()-t0),'s',flush=True)
np.savez_compressed(f'{D}/pixel_events30.npz',p95=P95,p50=P50,mean_excess=EXC,events_per_yr=EVY,hyper_days_per_yr=HYP)
REG=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']; ECO=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay','Other']
def qh(h,bins,q):
    c=np.cumsum(h); return bins[np.searchsorted(c,q*c[-1])+1] if c[-1]>0 else np.nan
for g,labels in (('reg',REG),('eco',ECO)):
    rows=[]
    for i,l in enumerate(labels):
        H_=hist[g]; rows.append(dict(group=l,events=H_['dur'][i].sum(),events_per_pixel_per_yr=H_['dur'][i].sum()/npx[g][i]/13,median_duration_d=qh(H_['dur'][i],DB,.5),p90_duration_d=qh(H_['dur'][i],DB,.9),
            median_peak_psu=qh(H_['peak'][i],PB,.5),p90_peak_psu=qh(H_['peak'][i],PB,.9),max_peak_psu=PB[np.flatnonzero(H_['peak'][i])[-1]+1] if H_['peak'][i].sum() else np.nan,
            median_excess_psu=qh(H_['excess'][i],XB,.5),p90_excess_psu=qh(H_['excess'][i],XB,.9),hypersaline_days_per_yr_mean=hyp_year[g][i].sum()/npx[g][i]/13))
    pd.DataFrame(rows).round(2).to_csv(f'{D}/high_salinity_events_by_{"region" if g=="reg" else "ecosystem"}.csv',index=False)
    Y=[]
    for i,l in enumerate(labels):
        for y in range(13): Y.append(dict(group=l,year=2013+y,events_per_pixel=ev_year[g][i,y]/npx[g][i],event_days_per_pixel=evd_year[g][i,y]/npx[g][i],pct_area_with_event=100*px_event_year[g][i,y]/npx[g][i],
                                         hyper_days_mean=hyp_year[g][i,y]/npx[g][i],pct_area_30plus_hyper_days=100*hyp30[g][i,y]/npx[g][i]))
    pd.DataFrame(Y).round(3).to_csv(f'{D}/high_salinity_events_by_{"region" if g=="reg" else "ecosystem"}_year.csv',index=False)
    sm=pd.DataFrame(startmon[g],index=labels,columns=range(1,13)); (sm.div(sm.sum(1),axis=0)*100).round(1).to_csv(f'{D}/event_start_month_pct_{"region" if g=="reg" else "ecosystem"}.csv')
S=pd.DataFrame(SSER,index=edates); S.index.name='date'; S=S.stack().rename('model_salinity_30m').reset_index().rename(columns={'level_1':'station'})
S.to_csv(f'{D}/station_daily_model_series_30m.csv.gz',index=False)
print(pd.read_csv(f'{D}/high_salinity_events_by_region.csv').to_string()); print('done',round(time.time()-t0))
