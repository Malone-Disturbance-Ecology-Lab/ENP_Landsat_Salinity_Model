import numpy as np, scipy.sparse as sp
from scipy.sparse.linalg import splu
from scipy import ndimage as ndi
def laplace_solver(known):
    """Harmonic (Laplace) fill of unknown cells, Dirichlet at known cells, Neumann at grid edges."""
    H,W=known.shape; unk=~known; idx=-np.ones((H,W),int); uy,ux=np.nonzero(unk); n=len(uy); idx[uy,ux]=np.arange(n)
    rows=[];cols=[];vals=[]; krow=[];kcol=[]; diag=np.zeros(n)
    ky,kx=np.nonzero(known); kidx=-np.ones((H,W),int); kidx[ky,kx]=np.arange(len(ky))
    for dy,dx in ((1,0),(-1,0),(0,1),(0,-1)):
        ny,nx=uy+dy,ux+dx; ok=(ny>=0)&(ny<H)&(nx>=0)&(nx<W); diag[ok]+=1
        i=np.nonzero(ok)[0]; nyy,nxx=ny[i],nx[i]; u=unk[nyy,nxx]
        rows+=list(i[u]); cols+=list(idx[nyy[u],nxx[u]]); vals+=[-1.0]*int(u.sum())
        krow+=list(i[~u]); kcol+=list(kidx[nyy[~u],nxx[~u]])
    A=sp.csc_matrix((np.r_[vals,diag],(np.r_[rows,np.arange(n)],np.r_[cols,np.arange(n)])),shape=(n,n))
    Bk=sp.csr_matrix((np.ones(len(krow)),(krow,kcol)),shape=(n,len(ky))); lu=splu(A)
    def fill(g):  # g: (H,W) or (T,H,W)
        g=np.asarray(g,float); single=g.ndim==2; G=g[None] if single else g
        kv=G[:,ky,kx].T; out=G.copy(); out[:,uy,ux]=lu.solve(Bk@kv).T
        return out[0] if single else out
    return fill
def nearest_inside(g,region):
    v=np.isfinite(g)&region
    if (v==region).all(): return np.where(region,g,np.nan)
    _,(iy,ix)=ndi.distance_transform_edt(~v,return_indices=True); return np.where(region,g[iy,ix],np.nan)
