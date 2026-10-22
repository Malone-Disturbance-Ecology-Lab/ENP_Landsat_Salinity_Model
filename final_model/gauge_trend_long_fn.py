import numpy as np, pandas as pd
from scipy.stats import norm, rankdata
NY = 26; iu, ju = np.triu_indices(NY, 1); DT = (ju - iu).astype(float)
def sk(series):
    """series: pd.Series indexed by Period[M]; returns slope (per yr), z, p (Hirsch-Slack corrected)"""
    X = np.full((NY, 12), np.nan)
    for per, v in series.items(): X[per.year - 2000, per.month - 1] = v
    d = X[ju] - X[iu]; sg = np.nan_to_num(np.sign(d)); S = sg.sum(); ng = np.isfinite(X).sum(0)
    var = 0.0; R = np.zeros_like(X)
    for g in range(12):
        x = X[:, g]; ok = np.isfinite(x); r = np.full(NY, (ng[g] + 1) / 2.0); r[ok] = rankdata(x[ok]); R[:, g] = r
    K = sg.T @ sg; RR = R.T @ R; cov = (K + 4 * RR - NY * np.outer(ng + 1, ng + 1)) / 3.0
    for g in range(12):
        x = X[:, g][np.isfinite(X[:, g])]; n = len(x); _, c = np.unique(x, return_counts=True); cov[g, g] = (n * (n - 1) * (2 * n + 5) - (c * (c - 1) * (2 * c + 5)).sum()) / 18.0
    V = max(cov.sum(), 1e-9); z = (S - np.sign(S)) / np.sqrt(V) if S != 0 else 0.0
    sl = (d / DT[:, None]).ravel(); return np.nanmedian(sl), z, 2 * norm.sf(abs(z))
