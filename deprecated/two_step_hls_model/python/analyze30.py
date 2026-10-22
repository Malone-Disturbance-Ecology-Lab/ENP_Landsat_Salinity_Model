"""30 m analyses from monthly 30 m maps: region x ecosystem means, seasonal cycle, annual/seasonal series,
Theil-Sen + Mann-Kendall trends per pixel with BH-FDR, 2020-25 vs 2013-18 change."""
import numpy as np, pandas as pd, rasterio, json, time
from scipy.stats import norm
t0=time.time(); H,W=3839,4171; NM=156; D='change30'; import os; os.makedirs(D,exist_ok=True)
MO=np.load('map30/monthly30_int16.npy',mmap_mode='r'); idx=np.load('px30_index.npy')
STK=json.load(open('static30_keys.json')); STA=np.load('static30.npy',mmap_mode='r'); distkm=np.asarray(STA[STK.index('distkm')])
LL=np.load('lonlat30.npy',mmap_mode='r'); lon=np.asarray(LL[0]); lat=np.asarray(LL[1])
V=np.load('vegclass30.npy'); cd=pd.read_csv('sitevars/vegmap_codes.csv').set_index('id')
cls_of=np.array(['']+list(cd.cls.reindex(range(1,cd.index.max()+1)).fillna('')),dtype=object)
GROUP={'Freshwater marsh/prairie':'Freshwater marsh & prairie','Swamp scrub/shrubland':'Swamp forest & scrub','Swamp/cypress woodland':'Swamp forest & scrub',
       'Mangrove':'Mangrove','Salt marsh':'Salt marsh','Open water':'Inland & estuarine water'}
ECO=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
eco=pd.Series(cls_of[V]).map(GROUP).fillna('Other').values.astype(object)
unm=(V==0); eco[unm&(distkm<0)]='Open bay'; eco[unm&(distkm>=0)]='Inland & estuarine water'
names=np.array(['Shark River Slough','Taylor Slough','C-111','Big Cypress']); region=names[np.load('basin30.npy')-1].astype(object)
off=distkm<-0.5; region[off]=np.where((lon[off]<-80.95)&(lat[off]>25.2),'Gulf coastal waters','Florida Bay')
REG=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
eco_c=pd.Categorical(eco,categories=ECO+['Other']).codes.astype('int8'); reg_c=pd.Categorical(region,categories=REG).codes.astype('int8')
np.save(f'{D}/eco30.npy',eco_c); np.save(f'{D}/region30.npy',reg_c)
AREA=(abs(0.00027092271)*111.0)*(0.00027092271*111.0*np.cos(np.deg2rad(25.3)))   # km2 per 30 m pixel
months=pd.period_range('2013-01','2025-12',freq='M'); yr=np.array([p.year for p in months]); mon=np.array([p.month for p in months]); yrs=np.arange(2013,2026)
sel={'annual':[yr==y for y in yrs],'dry':[((yr==y)&(mon<=5))|((yr==y-1)&(mon==12)) for y in yrs],'wet':[(yr==y)&(mon>=6)&(mon<=11) for y in yrs]}
i_,j_=np.triu_indices(13,1); dt=(yrs[j_]-yrs[i_]).astype('float32'); var=13*12*31/18.0
NP=len(idx); OUTK=['mean','mean_dry','mean_wet','sd_annual']+[f'{p}_{k}' for k in sel for p in ('slope','p','delta')]
R={k:np.full(NP,np.nan,'float32') for k in OUTK}
ne,nr=len(ECO)+1,len(REG)
acc_clim=np.zeros((ne,12)); acc_ms=np.zeros((ne,NM)); cnt=np.zeros(ne)
acc_sea={k:np.zeros((ne,13)) for k in sel}; acc_sea_r={k:np.zeros((nr,13)) for k in sel}; cnt_r=np.zeros(nr)
CH=400000
for a in range(0,NP,CH):
    b=min(a+CH,NP); p=idx[a:b]; rr,cc=p//W,p%W
    X=np.asarray(MO[:,rr,cc],dtype='float32')/10; X[X<0]=np.nan          # (156,n)
    S={k:np.stack([np.nanmean(X[s],0) for s in v]) for k,v in sel.items()}
    R['mean'][a:b]=np.nanmean(S['annual'],0); R['mean_dry'][a:b]=np.nanmean(S['dry'],0); R['mean_wet'][a:b]=np.nanmean(S['wet'],0); R['sd_annual'][a:b]=np.nanstd(S['annual'],0)
    for k,Y in S.items():
        d=Y[j_]-Y[i_]; R[f'slope_{k}'][a:b]=np.nanmedian(d/dt[:,None],0); Ssum=np.nansum(np.sign(d),0)
        z=np.where(Ssum>0,(Ssum-1)/np.sqrt(var),np.where(Ssum<0,(Ssum+1)/np.sqrt(var),0)); R[f'p_{k}'][a:b]=2*norm.sf(np.abs(z))
        R[f'delta_{k}'][a:b]=np.nanmean(Y[yrs>=2020],0)-np.nanmean(Y[yrs<=2018],0)
    ec=eco_c[a:b]; rc=reg_c[a:b]
    for e in range(ne):
        m=ec==e
        if not m.any(): continue
        cnt[e]+=m.sum(); acc_ms[e]+=np.nansum(X[:,m],1)
        for mm in range(12): acc_clim[e,mm]+=np.nansum(np.nanmean(X[mon==mm+1][:,m],0))
        for k in sel: acc_sea[k][e]+=np.nansum(S[k][:,m],1)
    for r in range(nr):
        m=rc==r
        if not m.any(): continue
        cnt_r[r]+=m.sum()
        for k in sel: acc_sea_r[k][r]+=np.nansum(S[k][:,m],1)
    print('chunk',a,round(time.time()-t0),flush=True)
def bh(p):
    ok=np.isfinite(p); q=np.full(p.shape,np.nan,'float32'); ps=p[ok]; o=np.argsort(ps); n=len(ps)
    qq=ps[o]*n/np.arange(1,n+1); qq=np.minimum.accumulate(qq[::-1])[::-1]; t=np.empty(n,'float32'); t[o]=np.minimum(qq,1); q[ok]=t; return q
for k in sel: R[f'q_{k}']=bh(R[f'p_{k}'])
np.savez_compressed(f'{D}/pixel_results30.npz',**{k:v.astype('float32') for k,v in R.items()})
pix=pd.DataFrame({'eco':eco,'region':region,**{k:R[k] for k in ['mean','mean_dry','mean_wet','sd_annual']+[f'{p}_{k}' for k in sel for p in ('slope','q','p','delta')]}})
g=pix[pix.eco!='Other']
g.groupby(['eco','region'])['mean'].mean().unstack().reindex(index=ECO,columns=REG).round(1).to_csv(f'{D}/mean_salinity_region_x_ecosystem.csv')
(g.groupby(['eco','region']).size().unstack()*AREA).reindex(index=ECO,columns=REG).round(1).to_csv(f'{D}/area_km2_region_x_ecosystem.csv')
rows=[]
for e in ECO:
    s=g[g.eco==e]; d=dict(ecosystem=e,pixels=len(s),area_km2=round(len(s)*AREA),mean=s['mean'].mean(),dry_mean=s.mean_dry.mean(),wet_mean=s.mean_wet.mean(),sd_between_pixels=s['mean'].std(),mean_interannual_sd=s.sd_annual.mean())
    for k in sel:
        d[f'{k}_median_slope']=s[f'slope_{k}'].median(); d[f'{k}_pct_sig_up']=100*((s[f'q_{k}']<0.05)&(s[f'slope_{k}']>0)).mean(); d[f'{k}_pct_sig_down']=100*((s[f'q_{k}']<0.05)&(s[f'slope_{k}']<0)).mean()
        d[f'{k}_pct_p05']=100*(s[f'p_{k}']<0.05).mean(); d[f'{k}_delta_2020_25_vs_2013_18']=s[f'delta_{k}'].mean()
    rows.append(d)
pd.DataFrame(rows).round(3).to_csv(f'{D}/ecosystem_summary_and_trends.csv',index=False)
rows=[]
for r in REG:
    s=g[g.region==r]; rows.append(dict(region=r,pixels=len(s),area_km2=round(len(s)*AREA),mean=s['mean'].mean(),sd_annual_mean=s.sd_annual.mean(),sd_annual_p95=s.sd_annual.quantile(.95),
        **{f'{k}_median_slope':s[f'slope_{k}'].median() for k in sel},**{f'{k}_pct_sig':100*(s[f'q_{k}']<0.05).mean() for k in sel}))
pd.DataFrame(rows).round(3).to_csv(f'{D}/region_summary_and_trends.csv',index=False)
cl=[dict(ecosystem=e,month=mm+1,mean=acc_clim[i,mm]/cnt[i]) for i,e in enumerate(ECO) if cnt[i]>0 for mm in range(12)]; pd.DataFrame(cl).to_csv(f'{D}/monthly_climatology_by_ecosystem.csv',index=False)
pd.DataFrame({e:acc_ms[i]/cnt[i] for i,e in enumerate(ECO) if cnt[i]>0},index=months.to_timestamp()).to_csv(f'{D}/monthly_series_by_ecosystem.csv')
ts=[dict(group='ecosystem',name=e,season=k,year=y,mean=acc_sea[k][i,jj]/cnt[i]) for k in sel for i,e in enumerate(ECO) if cnt[i]>0 for jj,y in enumerate(yrs)]
ts+=[dict(group='region',name=r,season=k,year=y,mean=acc_sea_r[k][i,jj]/cnt_r[i]) for k in sel for i,r in enumerate(REG) for jj,y in enumerate(yrs)]
pd.DataFrame(ts).to_csv(f'{D}/annual_seasonal_means_by_group.csv',index=False)
print(pd.read_csv(f'{D}/ecosystem_summary_and_trends.csv').round(2).T.to_string()); print(pd.read_csv(f'{D}/mean_salinity_region_x_ecosystem.csv').to_string()); print('done',round(time.time()-t0))
