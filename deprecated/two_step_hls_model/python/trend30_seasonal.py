"""Trends in monthly salinity anomalies at every 30 m pixel (deseasonalized trend).
Seasonal Kendall test (Hirsch et al. 1982) with the Hirsch & Slack (1984) covariance correction for
dependence between months, and the seasonal Theil-Sen slope (median of all within-calendar-month pairwise slopes).
Because only same-month pairs are compared, the test is identical on raw monthly values and on monthly anomalies
(value minus that pixel's 2013-2025 mean for the calendar month).
Sets: all 12 months, dry-season months (Dec-May), wet-season months (Jun-Nov). BH-FDR across pixels."""
import numpy as np, pandas as pd, time, json
from scipy.stats import norm, rankdata
t0=time.time(); H,W=3839,4171; D='change30'; NY=13
MO=np.load('map30/monthly30_int16.npy',mmap_mode='r'); idx=np.load('px30_index.npy'); NP=len(idx)
iu,ju=np.triu_indices(NY,1); DT=(ju-iu).astype('float32')             # 78 year pairs
SETS={'monthly':list(range(12)),'drymonths':[11,0,1,2,3,4],'wetmonths':[5,6,7,8,9,10]}

def seasonal_kendall(X):
    """X: (NY, 12, n) float32 with NaN for missing. Returns dict of slope, S, var, z, p per set."""
    n=X.shape[2]; d=X[ju]-X[iu]                                          # (78,12,n)
    sg=np.sign(d); sg[np.isnan(sg)]=0; sg=sg.astype(np.int8)
    Sg=sg.sum(0,dtype=np.int32)                                          # (12,n)
    ng=np.isfinite(X).sum(0)                                             # (12,n)
    # Mann-Kendall variance per month with the tie correction
    Xf=np.where(np.isfinite(X),X,np.inf); srt=np.sort(Xf,0)
    eq=(srt[1:]==srt[:-1])&np.isfinite(srt[1:])
    run=np.zeros(X.shape[1:],'float32'); tie=np.zeros(X.shape[1:],'float32')
    for k in range(eq.shape[0]):
        e=eq[k]; t=run+1; tie+=np.where((~e)&(run>0),t*(t-1)*(2*t+5),0); run=np.where(e,run+1,0)
    t=run+1; tie+=np.where(run>0,t*(t-1)*(2*t+5),0)
    varg=(ng*(ng-1)*(2*ng+5)-tie)/18.0
    # ranks within each month's series (ties averaged); missing -> (ng+1)/2
    r=rankdata(Xf,axis=0).astype('float32')
    R=np.where(np.isfinite(X),r,((ng+1)/2.0)[None])
    # K_gh = sum over pairs of sgn(d_g) sgn(d_h)
    K=np.empty((12,12,n),'float32')
    for g in range(12):
        K[g]=(sg[:,g][:,None,:]*sg).sum(0,dtype=np.int32)
    RR=np.einsum('ign,ihn->ghn',R,R)
    ng1=(ng+1).astype('float32'); cov=(K+4*RR-NY*ng1[:,None,:]*ng1[None,:,:])/3.0   # Hirsch & Slack (1984) with missing values
    for g in range(12): cov[g,g]=varg[g]
    slopes=d/DT[:,None,None]                                             # (78,12,n)
    out={}
    for k,ms in SETS.items():
        S=Sg[ms].sum(0).astype('float32'); V=cov[np.ix_(ms,ms)].sum((0,1))
        V=np.maximum(V,1e-6); z=np.where(S>0,(S-1)/np.sqrt(V),np.where(S<0,(S+1)/np.sqrt(V),0))
        sl=slopes[:,ms].reshape(-1,n)
        med=np.median(sl,0) if not np.isnan(sl).any() else np.nanmedian(sl,0)
        out[k]=dict(slope=med.astype('float32'),z=z.astype('float32'),p=(2*norm.sf(np.abs(z))).astype('float32'))
    return out

if __name__=='__main__':
    R={f'{q}_{k}':np.full(NP,np.nan,'float32') for k in SETS for q in ('slope','z','p')}
    CH=60000
    for a in range(0,NP,CH):
        b=min(a+CH,NP); p=idx[a:b]; rr,cc=p//W,p%W
        X=np.asarray(MO[:,rr,cc],dtype='float32')/10; X[X<0]=np.nan
        o=seasonal_kendall(X.reshape(NY,12,-1))
        for k,v in o.items():
            for q in ('slope','z','p'): R[f'{q}_{k}'][a:b]=v[q]
        if (a//CH)%10==0: print('px',a,round(time.time()-t0),'s',flush=True)
    def bh(p):
        ok=np.isfinite(p); q=np.full(p.shape,np.nan,'float32'); ps=p[ok]; o=np.argsort(ps); n=len(ps)
        qq=ps[o]*n/np.arange(1,n+1); qq=np.minimum.accumulate(qq[::-1])[::-1]; t=np.empty(n,'float32'); t[o]=np.minimum(qq,1); q[ok]=t; return q
    for k in SETS: R[f'q_{k}']=bh(R[f'p_{k}'])
    np.savez_compressed(f'{D}/pixel_trends_monthly_anomaly30.npz',**R)
    eco=np.load(f'{D}/eco30.npy'); reg=np.load(f'{D}/region30.npy')
    ECO=['Freshwater marsh & prairie','Swamp forest & scrub','Salt marsh','Mangrove','Inland & estuarine water','Open bay']
    REG=['Big Cypress','Shark River Slough','Taylor Slough','C-111','Gulf coastal waters','Florida Bay']
    A=0.000818; rows=[]
    for gname,cat,labels in (('ecosystem',eco,ECO),('region',reg,REG)):
        for i,l in enumerate(labels):
            m=cat==i; d=dict(group=gname,name=l,pixels=int(m.sum()),area_km2=round(m.sum()*A))
            for k in SETS:
                s=R[f'slope_{k}'][m]; q=R[f'q_{k}'][m]; pp=R[f'p_{k}'][m]
                d[f'{k}_median_slope']=float(np.nanmedian(s)); d[f'{k}_p10_slope']=float(np.nanpercentile(s,10)); d[f'{k}_p90_slope']=float(np.nanpercentile(s,90))
                d[f'{k}_pct_sig_up']=100*float(((q<0.05)&(s>0)).mean()); d[f'{k}_pct_sig_down']=100*float(((q<0.05)&(s<0)).mean())
                d[f'{k}_km2_sig_up']=round(float(((q<0.05)&(s>0)).sum()*A),1); d[f'{k}_km2_sig_down']=round(float(((q<0.05)&(s<0)).sum()*A),1)
                d[f'{k}_pct_p05']=100*float((pp<0.05).mean())
            rows.append(d)
    T=pd.DataFrame(rows); T.round(3).to_csv(f'{D}/monthly_anomaly_trends_by_group.csv',index=False)
    print(T.round(2).T.to_string()); print('done',round(time.time()-t0))
