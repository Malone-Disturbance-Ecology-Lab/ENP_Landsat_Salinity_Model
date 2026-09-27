import rasterio, numpy as np, pandas as pd
U='/Volumes/malonelab/Research/ENP_Landsat_Salinity_Model/'
v=rasterio.open(U+'site_variables/ENP_VegMap_30m.tif'); V=v.read(1); T=v.transform
cd=pd.read_csv(U+'site_variables/vegmap_codes.csv').set_index('id')
classes=sorted(cd.cls.unique())
st=pd.read_csv(U+'DBHydro_lonlat.csv').drop_duplicates('station')
sdr=rasterio.open(U+'site_variables/ENP_DistCoast_signed.tif'); SD=sdr.read(1)
mk=rasterio.open(U+'site_variables/ENP_MainlandMask.tif'); MK=mk.read(1)
rr,cc=rasterio.transform.rowcol(sdr.transform,st.longitude.values,st.latitude.values)
rr=np.array(rr); cc=np.array(cc); inb=(rr>=0)&(rr<SD.shape[0])&(cc>=0)&(cc<SD.shape[1])
st['in_scene']=inb
st['distCoast_signed_m']=np.where(inb,SD[np.clip(rr,0,SD.shape[0]-1),np.clip(cc,0,SD.shape[1]-1)],np.nan)
st['mask_class']=np.where(inb,MK[np.clip(rr,0,SD.shape[0]-1),np.clip(cc,0,SD.shape[1]-1)],255)
el=rasterio.open(U+'ENP_Elevation.tif'); E=el.read(1)
oc=rasterio.open(U+'ENP_DistCoast.tif'); O=oc.read(1)
dy=v.res[1]*110574; dx=v.res[0]*111320*np.cos(np.radians(25.37))
out=[]
for _,s in st.iterrows():
    r,c=rasterio.transform.rowcol(T,s.longitude,s.latitude)
    rec=dict(station=s.station)
    ok=0<=r<V.shape[0] and 0<=c<V.shape[1]
    rec['elev_elevatr']=float(E[r,c]) if ok else np.nan
    rec['distCoast_old']=float(O[r,c]) if ok else np.nan
    code=V[r,c] if ok else 0
    rec['veg_point_code']=cd.VegCode.get(code,np.nan) if code else np.nan
    rec['veg_point_name']=cd.name.get(code,np.nan) if code else np.nan
    rec['veg_point_class']=cd.cls.get(code,np.nan) if code else np.nan
    for rad in [250,1000]:
        ny,nx=int(np.ceil(rad/dy)),int(np.ceil(rad/dx))
        yy,xx=np.mgrid[-ny:ny+1,-nx:nx+1]; disk=(yy*dy)**2+(xx*dx)**2<=rad**2
        r0,c0=r-ny,c-nx
        if not ok: continue
        win=np.zeros(disk.shape,np.uint16)
        rs=slice(max(r0,0),min(r+ny+1,V.shape[0])); cs=slice(max(c0,0),min(c+nx+1,V.shape[1]))
        win[rs.start-r0:rs.stop-r0, cs.start-c0:cs.stop-c0]=V[rs,cs]
        vals=win[disk]; mapped=vals[vals>0]
        rec[f'mapped_frac_{rad}m']=len(mapped)/disk.sum()
        if len(mapped)==0: continue
        cl=cd.cls.reindex(mapped).values
        for k in classes: rec[f'p_{k}_{rad}m']=(cl==k).mean()
        nonw=cl[cl!='Open water']
        rec[f'veg_dominant_{rad}m']=pd.Series(nonw).mode()[0] if len(nonw) else 'Open water'
    out.append(rec)
sv=st.merge(pd.DataFrame(out),on='station')
sv.to_csv(U+'site_variables/station_site_variables.csv',index=False)
pd.set_option('display.width',250)
print(sv[['station','distCoast_signed_m','distCoast_old','elev_elevatr','veg_point_class','veg_dominant_1000m','mapped_frac_1000m']].round(0).to_string())
sv[['station','longitude','latitude','grab','in_scene','distCoast_signed_m','mask_class']].to_csv(U+'site_variables/stations_coast.csv',index=False)   # compact station table read by the model scripts
