exec(open('ts_setup.py').read())
from sklearn.ensemble import HistGradientBoostingRegressor as H
f1=STATIC+SPEC+SLOW
cfg0={'small(15 leaves,min40)':dict(max_iter=300,learning_rate=0.05,max_leaf_nodes=15,min_samples_leaf=40,l2_regularization=1.0),
     'tiny(7 leaves,min80)':dict(max_iter=400,learning_rate=0.03,max_leaf_nodes=7,min_samples_leaf=80,l2_regularization=2.0),
     'mid(31 leaves,min20)':dict(max_iter=300,learning_rate=0.05,max_leaf_nodes=31,min_samples_leaf=20,l2_regularization=1.0),
     'small+monotone-free lr.02 800':dict(max_iter=800,learning_rate=0.02,max_leaf_nodes=15,min_samples_leaf=40,l2_regularization=1.0)}
cfg={'4 leaves min250':dict(max_iter=600,learning_rate=0.03,max_leaf_nodes=4,min_samples_leaf=250,l2_regularization=3.0),
 '3 leaves min300':dict(max_iter=800,learning_rate=0.03,max_leaf_nodes=3,min_samples_leaf=300,l2_regularization=3.0),
 '4 leaves min400':dict(max_iter=700,learning_rate=0.03,max_leaf_nodes=4,min_samples_leaf=400,l2_regularization=3.0),
 '2 leaves (stumps) min300':dict(max_iter=1500,learning_rate=0.03,max_leaf_nodes=2,min_samples_leaf=300,l2_regularization=3.0)}
for fs_name,fs in {'static+spectra+slow':f1,'static+slow (NO spectra)':STATIC+SLOW}.items():
  for cn,c in cfg.items():
    pass
    sA=np.full(len(TA),np.nan); sB=np.full(len(TB),np.nan)
    for g in sorted(TA.block.unique()):
        tr=TA[(TA.block!=g)&TA.B.notna()]; m=H(random_state=0,**c).fit(tr[fs],tr.B,sample_weight=np.sqrt(tr.n))
        sA[(TA.block==g).values]=m.predict(TA.loc[TA.block==g,fs]); sB[(TB.block==g).values]=m.predict(TB.loc[TB.block==g,fs])
    TA['bh']=sA; TB['bh']=sB; both=pd.concat([TA[['station','month','bh']],TB[['station','month','bh']]]).reset_index(drop=True); both['s']=smooth(both,'bh')
    a=both.s.values[:len(TA)]; b=both.s.values[len(TA):]; o=TA.B.notna().values; ob=TB.B.notna().values
    print(f'{fs_name:22s} {cn:30s} CV R2 {s1r(TA.B.values[o],a[o]):.3f}  2025 R2 {s1r(TB.B.values[ob],b[ob]):.3f}',flush=True)
