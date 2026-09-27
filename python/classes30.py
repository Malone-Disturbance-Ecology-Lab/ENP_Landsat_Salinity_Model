"""30 m ecosystem/region classes (from from monthly 30 m maps: region x ecosystem means, seasonal cycle, annual/seasonal series,
Theil-Sen + Mann-Kendall trends per pixel with BH-FDR, 2020-25 vs 2013-18 change."""
import numpy as np, pandas as pd, rasterio, json, time
from scipy.stats import norm
import os
t0=time.time(); H,W=3839,4171; NM=156; D=os.environ.get('D30','change30'); import os; os.makedirs(D,exist_ok=True)
idx=np.load('px30_index.npy')
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
