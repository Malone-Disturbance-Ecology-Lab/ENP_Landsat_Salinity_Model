"""Select the best two-step combination (spatial-block CV + 2025 hold-out) and save daily OOF predictions."""
exec(open('ts_setup.py').read())
from sklearn.ensemble import HistGradientBoostingRegressor as H
S1CFG=dict(max_iter=700,learning_rate=0.03,max_leaf_nodes=4,min_samples_leaf=400,l2_regularization=3.0,random_state=0)
F1=STATIC+SPEC+SLOW
def fit_s1(tr,f=F1): return H(**S1CFG).fit(tr[f],tr.B,sample_weight=np.sqrt(tr.n))
blocks=sorted(A.block.unique())
def smooth_pair(a,b):
    both=pd.concat([TA[['station','month']].assign(bh=a),TB[['station','month']].assign(bh=b)]).reset_index(drop=True); both['s']=smooth(both,'bh')
    return both.s.values[:len(TA)],both.s.values[len(TA):]
P={k:(np.full(len(A),np.nan),np.full(len(B_),np.nan)) for k in ['v3 one-step','K1 two-step','K2 two-step + Bhat in step 2','K3 stacked v3 + Bhat']}
BH=(np.full(len(A),np.nan),np.full(len(B_),np.nan))
for g in blocks:
    ha=(A.block==g).values; hb=(B_.block==g).values; hA=(TA.block==g).values; hB=(TB.block==g).values
    # outer step-1 prediction for the held-out block
    m1=fit_s1(TA[(TA.block!=g)&TA.B.notna()]); a=np.full(len(TA),np.nan); b=np.full(len(TB),np.nan); a[hA]=m1.predict(TA.loc[hA,F1]); b[hB]=m1.predict(TB.loc[hB,F1])
    # nested step-1 predictions for the training blocks
    for g2 in [x for x in blocks if x!=g]:
        m=fit_s1(TA[(~TA.block.isin([g,g2]))&TA.B.notna()]); h=(TA.block==g2).values; a[h]=m.predict(TA.loc[h,F1])
    sa,sb=smooth_pair(a,b); A['Bhat']=attach(A,TA,sa); B_['Bhat']=attach(B_,TB,sb)
    BH[0][ha]=A.Bhat[ha]; BH[1][hb]=B_.Bhat[hb]
    tr=A[A.block!=g]; trb=tr[tr.B.notna()]
    fits={'v3 one-step':(FEATS_V3,tr,tr.salinity,False),'K1 two-step':(S2F,trb,trb.salinity-trb.B,True),
          'K2 two-step + Bhat in step 2':(S2F+['Bhat'],trb,trb.salinity-trb.Bhat,True),'K3 stacked v3 + Bhat':(FEATS_V3+['Bhat'],tr,tr.salinity,False)}
    for k,(F,d,y,add) in fits.items():
        m=gbm().fit(d[F],y); pa=m.predict(A.loc[ha,F]); P[k][0][ha]=pa+(A.Bhat[ha] if add else 0)
        if hb.any(): pb=m.predict(B_.loc[hb,F]); P[k][1][hb]=pb+(B_.Bhat[hb] if add else 0)
    print('fold',g,flush=True)
P['K4 ensemble K2+K3']=((P['K2 two-step + Bhat in step 2'][0]+P['K3 stacked v3 + Bhat'][0])/2,(P['K2 two-step + Bhat in step 2'][1]+P['K3 stacked v3 + Bhat'][1])/2)
P['K5 ensemble K1+K2+K3']=tuple((P['K1 two-step'][i]+P['K2 two-step + Bhat in step 2'][i]+P['K3 stacked v3 + Bhat'][i])/3 for i in (0,1))
res=[]
for k,(pa,pb) in P.items():
    m=metrics(A,pa); mb=metrics(B_,pb)
    res.append(dict(model=k,R2=m['R2'],RMSE=m['RMSE'],R2_2025=mb['R2'],RMSE_2025=mb['RMSE'],station_r=m['median_station_r'],station_r_2025=mb['median_station_r'],station_mean_R2=m['station_mean_R2'],station_mean_R2_2025=mb['station_mean_R2'],med_abs_bias=m['median_abs_bias']))
R=pd.DataFrame(res); R.to_csv('best_two_step_cv.csv',index=False); print(R.round(3).to_string())
pickle.dump(dict(P=P,BH=BH,keyA=A[['station','date','salinity','block','basin','B']],keyB=B_[['station','date','salinity','block','basin','B']]),open('best_two_step_oof.pkl','wb'))
