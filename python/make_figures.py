"""Make every figure for the ENP 30 m salinity work from the released outputs on the lab server.

    python make_figures.py --set ms            # the 7 manuscript figures (Scientific Data Data Descriptor)
    python make_figures.py --set ms 3 6        # only manuscript Figures 3 and 6
    python make_figures.py --set doc           # methods/results document figures (legacy, see make_figures_doc.py)
    python make_figures.py                     # both sets

Options
  --root PATH   project data root (default: $ENP_DATA_ROOT, else /Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model)
  --out PATH    output folder for the manuscript set (default: <root>/best_two_step/salinity30m/figures/manuscript)
  --refresh     recompute the cached map summaries (30 m means and the monthly-anomaly trend map)

Manuscript figures (PNG and PDF, 300 dpi) and the files each one reads (relative to --root):
  Fig 1  study area          best_two_step/salinity30m/ENP_region_ecosystem_classes_30m.tif, site_variables/station_site_variables.csv
  Fig 2  workflow schematic  (no inputs)
  Fig 3  example products    best_two_step/salinity30m/monthly/ENP_salinity30m_monthly_YYYY_MM.tif (156)
  Fig 4  validation          best_two_step/timeseries/station_daily_model_vs_gauge.csv.gz (held-out predictions),
                             best_two_step/timeseries/station_timeseries_metrics.csv
  Fig 5  station series      same as Fig 4
  Fig 6  250 m vs 30 m       best_two_step/daily_maps/monthly/ENP_salinity_monthly_mean_2020.tif, 30 m monthly 2020_05
  Fig 7  usage examples      best_two_step/salinity30m/change30/{monthly_climatology_by_ecosystem, high_salinity_events_by_region_year,
                             event_days_gauge_vs_model_30m}.csv and the 30 m monthly maps (trend map)

The 90% intervals are rebuilt from the held-out predictions exactly as in the paper: split-conformal half-widths = the
conformal 90th percentile of |held-out residual| in 2013-2024 within bins of predicted salinity (0-5, 5-15, 15-30, >30 psu).

Colour convention: wherever colour encodes salinity, red = more saline and blue = fresher.
"""
import os as _os; ROOT = _os.environ.get('ENP_DATA_ROOT', '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'); WORK = _os.environ.get('ENP_WORK', 'work')   # data root and local work folder (README)
import os, sys, argparse, glob, numpy as np, pandas as pd
import matplotlib, matplotlib.ticker
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm, TwoSlopeNorm, to_rgb
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
from matplotlib.lines import Line2D

# ------------------------------------------------------------------ style
INK = '#3d3d3a'; MUTED = '#8a897f'; GRIDC = '#e6e5df'; OBS = '#e0592a'; MOD = '#1f5fae'; BAND = '#c9d7ea'
SAL = 'RdYlBu_r'; DIV = 'RdBu_r'
RAMP = ['#1f5fae', '#56b1e3', '#2f9e6b', '#e7a33a', '#e0592a', '#a3172c']      # categorical, ordered fresh -> saline
ECO = ['Freshwater marsh & prairie', 'Swamp forest & scrub', 'Salt marsh', 'Mangrove', 'Inland & estuarine water', 'Open bay']
REG = ['Big Cypress', 'Shark River Slough', 'Taylor Slough', 'C-111', 'Gulf coastal waters', 'Florida Bay']   # band-1 codes 0-5
REG_BY_SAL = ['Taylor Slough', 'Shark River Slough', 'Big Cypress', 'C-111', 'Gulf coastal waters', 'Florida Bay']
COL_E = dict(zip(ECO, RAMP)); COL_R = dict(zip(REG_BY_SAL, RAMP))
BASIN = {'SRS': 'Shark River Slough', 'TaylorSlough': 'Taylor Slough', 'C111': 'C-111', 'BigCypress': 'Big Cypress'}
ANCHORS = ['ENPWP', 'ENPSR', 'ENPWW', 'ENPGB', 'TAYLORS3', 'ENPWB']
TS_STATIONS = ['ENPCN', 'ENPNR', 'ENPHR', 'ENPCA', 'ENPLB', 'ENPJK']             # fresh -> saline
LEVELS = [-1, 5, 15, 30, 1e9]; LEVEL_LABELS = ['0–5 psu', '5–15 psu', '15–30 psu', '> 30 psu']
ZOOM = (-81.20, -80.90, 25.28, 25.47)                                              # Shark River estuary / Whitewater Bay
FULL_W = 7.09                                                                      # Scientific Data full page width (180 mm)
DPI = 300
plt.rcParams.update({'font.size': 7.5, 'font.family': 'DejaVu Sans', 'axes.edgecolor': GRIDC, 'axes.labelcolor': INK,
                     'xtick.color': INK, 'ytick.color': INK, 'axes.titlecolor': INK, 'text.color': INK,
                     'axes.titlesize': 8, 'axes.titlelocation': 'left', 'legend.frameon': False,
                     'pdf.fonttype': 42, 'savefig.dpi': DPI})

# ------------------------------------------------------------------ paths
A = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
A.add_argument('figs', nargs='*', type=int)
A.add_argument('--set', choices=['ms', 'doc'], default=None)
A.add_argument('--root', default=ROOT)
A.add_argument('--out', default=None)
A.add_argument('--refresh', action='store_true')
ARGS = A.parse_args() if __name__ == '__main__' else A.parse_args([])
ROOT = ARGS.root
S30 = f'{ROOT}/best_two_step/salinity30m'; TS = f'{ROOT}/best_two_step/timeseries'; DM = f'{ROOT}/best_two_step/daily_maps'
CH = f'{S30}/change30'; SITE = f'{ROOT}/site_variables'
OUT = ARGS.out or f'{S30}/figures/manuscript'; CACHE = f'{S30}/ms_cache'


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    for ext in ('png', 'pdf'):
        fig.savefig(f'{OUT}/{name}.{ext}', dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close(fig); print('  wrote', name)


def clean(a):
    for s in ('top', 'right'): a.spines[s].set_visible(False)


def ink_on(cmap, norm, v):
    r, g, b = to_rgb(plt.get_cmap(cmap)(norm(v))); return 'white' if 0.299 * r + 0.587 * g + 0.114 * b < 0.5 else INK


# ------------------------------------------------------------------ grid, classes and map helpers
_G = {}
def grid():
    if not _G:
        import rasterio
        with rasterio.open(f'{S30}/ENP_region_ecosystem_classes_30m.tif') as r:
            reg = r.read(1); eco = r.read(2); b = r.bounds; _G['tr'] = r.transform
        _G['ext'] = (b.left, b.right, b.bottom, b.top); _G['reg'] = reg; _G['eco'] = eco; _G['dom'] = reg != 255
        _G['H'], _G['W'] = reg.shape; _G['d'] = r.transform.a
    return _G


def outlines(a, step=3, lw=0.4, domain=True):
    g = grid(); x0, x1, y0, y1 = g['ext']
    reg = g['reg'][::step, ::step].astype('float32'); reg[reg == 255] = np.nan
    H, W = reg.shape; xs = x0 + (np.arange(W) + 0.5) * (x1 - x0) / W; ys = y1 - (np.arange(H) + 0.5) * (y1 - y0) / H
    for k in range(6):
        a.contour(xs, ys, (reg == k).astype('float32'), levels=[0.5], colors=INK, linewidths=lw)
    if domain:
        a.contour(xs, ys, np.isfinite(reg).astype('float32'), levels=[0.5], colors=INK, linewidths=lw * 1.2)


def mapax(a, ticks=False):
    a.set_aspect(1 / np.cos(np.deg2rad(25.37)))
    if not ticks:
        a.set_xticks([]); a.set_yticks([]); [a.spines[s].set_visible(False) for s in a.spines]


def monthly_files():
    f = sorted(glob.glob(f'{S30}/monthly/ENP_salinity30m_monthly_*_*.tif'))
    f = [x for x in f if os.path.basename(x).count('_') == 4]
    if len(f) != 156: print(f'  warning: found {len(f)} monthly 30 m files (expected 156)')
    return f


def read_month(path, window=None, step=1):
    import rasterio
    with rasterio.open(path) as r:
        v = r.read(1, window=window)
    v = v[::step, ::step].astype('float32'); v[v < 0] = np.nan; return v / 10


def map_summaries():
    """30 m 2013-2025 mean and the May and October climatologies (cached as int16 x 10)."""
    fn = f'{CACHE}/mean_may_oct_30m.npz'
    if os.path.exists(fn) and not ARGS.refresh:
        z = np.load(fn); return {k: np.where(z[k] < 0, np.nan, z[k] / 10.0).astype('float32') for k in z.files}
    print('  computing 30 m mean and May/October climatology from the monthly maps ...')
    g = grid(); H, W = g['H'], g['W']
    acc = {k: (np.zeros((H, W), 'float32'), np.zeros((H, W), 'uint8')) for k in ('mean', 'may', 'oct')}
    for f in monthly_files():
        m = int(os.path.basename(f)[-6:-4]); v = read_month(f); ok = np.isfinite(v); v0 = np.where(ok, v, 0)
        for k in ['mean'] + (['may'] if m == 5 else []) + (['oct'] if m == 10 else []):
            acc[k][0][:] += v0; acc[k][1][:] += ok
    out = {}
    for k, (s, n) in acc.items():
        with np.errstate(invalid='ignore', divide='ignore'): out[k] = np.where(n > 0, s / n, np.nan).astype('float32')
    os.makedirs(CACHE, exist_ok=True)
    np.savez_compressed(fn, **{k: np.where(np.isfinite(v), np.round(v * 10), -1).astype('int16') for k, v in out.items()})
    return out


def trend_map(step=4):
    """Seasonal Theil-Sen slope of monthly salinity (psu/yr) on every `step`-th 30 m cell (cached).
    Same estimator as trend30_seasonal.py: median of all same-calendar-month pairwise slopes, 2013-2025."""
    fn = f'{CACHE}/trend_monthly_anomaly_step{step}.npy'
    if os.path.exists(fn) and not ARGS.refresh: return np.load(fn)
    print(f'  computing the monthly-anomaly trend map (every {step}th cell) ...')
    files = monthly_files(); X = None
    for i, f in enumerate(files):
        v = read_month(f, step=step)
        if X is None: X = np.full((len(files),) + v.shape, -1, 'int16')
        X[i] = np.where(np.isfinite(v), np.round(v * 10), -1)
    h, w = X.shape[1:]; NY = len(files) // 12
    Xf = X.reshape(NY, 12, h * w); valid = (Xf >= 0).sum((0, 1)) >= NY * 12 * 0.8; px = np.flatnonzero(valid)
    iu, ju = np.triu_indices(NY, 1); dt = (ju - iu).astype('float32')
    slope = np.full(h * w, np.nan, 'float32')
    for a in range(0, len(px), 40000):
        p = px[a:a + 40000]; Y = Xf[:, :, p].astype('float32') / 10; Y[Y < 0] = np.nan
        d = (Y[ju] - Y[iu]) / dt[:, None, None]
        slope[p] = np.nanmedian(d.reshape(-1, len(p)), 0)
    slope = slope.reshape(h, w); os.makedirs(CACHE, exist_ok=True); np.save(fn, slope); return slope


# ------------------------------------------------------------------ held-out predictions and intervals
_V = {}
def validation():
    if not _V:
        d = pd.read_csv(f'{TS}/station_daily_model_vs_gauge.csv.gz', parse_dates=['date']).dropna(subset=['salinity', 'model'])
        d['p'] = d.model.clip(lower=0); d['lev'] = pd.cut(d.p, LEVELS, labels=False); d['y25'] = d.date.dt.year == 2025
        cal = d[~d.y25]; q = []
        for k in range(4):                                          # split-conformal 90% half-widths
            r = np.sort((cal.salinity - cal.p).abs()[cal.lev == k].values); n = len(r)
            q.append(r[min(n - 1, int(np.ceil((n + 1) * 0.9)) - 1)])
        q = np.array(q); d['lo'] = (d.p - q[d.lev]).clip(lower=0); d['hi'] = d.p + q[d.lev]
        d['inside'] = (d.salinity >= d.lo) & (d.salinity <= d.hi)
        _V.update(d=d, q=q, st=pd.read_csv(f'{TS}/station_timeseries_metrics.csv'))
    return _V


# ================================================================== Figure 1: study area
def fig1():
    g = grid(); step = 3
    eco = g['eco'][::step, ::step].astype('float32'); eco[eco == 255] = np.nan; eco[eco == 6] = 6
    cmap = ListedColormap(RAMP + ['#d9d8d2']); norm = BoundaryNorm(np.arange(-0.5, 7.5), 7)
    s = pd.read_csv(f'{SITE}/station_site_variables.csv')
    s['role'] = np.where(s.station.isin(ANCHORS), 'anchor', np.where(s.grab == 1, 'grab', 'sonde'))
    fig, a = plt.subplots(figsize=(FULL_W, FULL_W * 0.66))
    a.imshow(eco, extent=g['ext'], cmap=cmap, norm=norm, interpolation='nearest')
    outlines(a, step=step)
    so = s[s.role == 'sonde']; gr = s[s.role == 'grab']; an = s[s.role == 'anchor']
    a.scatter(so.longitude, so.latitude, s=16, facecolor='white', edgecolor=INK, lw=0.8, zorder=4)
    a.scatter(gr.longitude, gr.latitude, s=12, facecolor=INK, edgecolor='white', lw=0.5, zorder=4)
    a.scatter(an.longitude, an.latitude, s=80, marker='*', facecolor='#f2c12e', edgecolor=INK, lw=0.6, zorder=5)
    for _, r in an.iterrows():
        a.annotate(r.station, (r.longitude, r.latitude), xytext=(5, 4), textcoords='offset points', fontsize=6,
                   bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.75), zorder=6)
    lab = {'Big Cypress': (-81.40, 25.64), 'Shark River\nSlough': (-80.82, 25.70), 'Taylor\nSlough': (-80.62, 25.43),
           'C-111': (-80.48, 25.37), 'Gulf\nwaters': (-81.44, 25.48), 'Florida Bay': (-80.56, 24.96)}
    for t, (x, y) in lab.items():
        a.text(x, y, t, ha='center', va='center', fontsize=7.5, bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7))
    h = [Patch(fc=c, label=e) for e, c in zip(ECO, RAMP)]
    h += [Line2D([], [], ls='', marker='o', mfc='white', mec=INK, ms=5, label=f'Model station, sonde (n = {len(so)})'),
          Line2D([], [], ls='', marker='o', mfc=INK, mec='white', ms=4.5, label=f'Model station, grab (n = {len(gr)})'),
          Line2D([], [], ls='', marker='*', mfc='#f2c12e', mec=INK, ms=10, label=f'Anchor station (n = {len(an)})')]
    a.legend(handles=h, loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=6.5, title='Ecosystem (30 m) and stations', title_fontsize=6.5,
             frameon=True, facecolor='white', edgecolor=GRIDC, framealpha=0.95)
    mapax(a, ticks=True); clean(a); a.set_xlabel('Longitude (°)'); a.set_ylabel('Latitude (°)')
    a.set_xticks([-81.5, -81.25, -81.0, -80.75, -80.5]); a.set_yticks([25.0, 25.25, 25.5, 25.75])
    save(fig, 'Fig1_study_area')


# ================================================================== Figure 2: workflow
def fig2():
    fig = plt.figure(figsize=(FULL_W, FULL_W * 0.74)); a = fig.add_axes([0, 0, 1, 1]); a.set_xlim(0, 100); a.set_ylim(0, 74); a.axis('off')
    def box(x, y, w, h, title, body, fc='white'):
        a.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=1.2', fc=fc, ec=INK, lw=1))
        a.text(x + 1.5, y + h - 1.6, title, fontsize=8.5, weight='bold', va='top')
        a.text(x + 1.5, y + h - 5.2, body, fontsize=7, va='top', linespacing=1.45)
    def arr(x0, y0, x1, y1): a.annotate('', (x1, y1), (x0, y0), arrowprops=dict(arrowstyle='-|>', color=MUTED, lw=1.4, mutation_scale=9))
    B1 = '#eef4fb'; B2 = '#fbeee8'
    box(2, 55, 31, 17, 'Salinity observations', 'DBHydro, 90 stations, 2013–2025\n84 model stations (train, validate)\n6 anchor stations (daily predictors)\nQC: 0–80 psu; daily means')
    box(2, 36.5, 31, 17, 'Vegetation spectra (HLS L30)', '737 clear scenes, 30 m pixels\nNDVI, EVI2, NDMI, SWIR1\nmeans within 1 km and 3 km\n12-mo, 3-mo, median, departure')
    box(2, 18, 31, 17, 'Site variables (30 m)', 'signed distance to coast,\nwater frequency, vegetation\nfractions, connectivity, basins')
    box(2, 1, 31, 15.5, 'Daily drivers', 'EDEN stage (400 m), gridMET\n(4 km), structure flows, anchor\nstations, 365-day departures')
    box(38, 46, 26, 26, 'Step 1: site mean B', '12-month centred mean\nshallow gradient boosting\n(4 leaves, ≥ 400 per leaf)\nspectra + site + slow drivers\nsmoothed over 12 months', B1)
    box(38, 20, 26, 23, 'Step 2: daily anomaly A', 'A = S − B\ngradient boosting (31 leaves)\ndaily drivers + departures,\nday of year, site variables', B1)
    box(38, 1, 26, 16, 'Blend', '$\\hat{S} = \\max(0,\\ 0.7(\\hat{B}+\\hat{A})$\n$\\quad +\\ 0.3\\,\\hat{S}_{\\mathrm{one\\text{-}step}})$', B1)
    box(69, 47.5, 29, 24.5, 'Validation', 'spatial-block CV, 10 blocks\n2013–2024 held out in space\n2025 held out in space and time\nsplit-conformal 90% intervals')
    box(69, 1, 29, 43.5, 'Data records', 'Daily 250 m maps, 2013–2025,\n  with monthly 90% bounds\nMonthly 30 m maps (Landsat\n  grid), 156 GeoTIFFs\nDaily 30 m maps (full model),\n  156 GeoTIFFs\nRegion and ecosystem classes\nStation daily series\nHeld-out predictions\nCode and figure script', B2)
    arr(33, 63, 37.8, 61); arr(33, 45, 37.8, 51); arr(33, 26.5, 37.8, 31); arr(33, 9, 37.8, 28)
    arr(51, 46, 51, 43.3); arr(51, 20, 51, 17.3); arr(64, 59, 68.8, 59); arr(83.5, 47.5, 83.5, 44.8); arr(64, 9, 68.8, 9)
    save(fig, 'Fig2_workflow')


# ================================================================== Figure 3: example products
def fig3():
    g = grid(); M = map_summaries(); st = 3
    fig = plt.figure(figsize=(FULL_W, FULL_W * 0.93))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 0.05], wspace=0.04, hspace=0.12)
    kw = dict(extent=g['ext'], cmap=SAL, vmin=0, vmax=40, interpolation='nearest')
    for (i, j), k, t in [((0, 0), 'mean', 'a  Mean, 2013–2025'), ((0, 1), 'may', 'b  May (end of dry season)'),
                         ((1, 0), 'oct', 'c  October (end of wet season)')]:
        a = fig.add_subplot(gs[i, j]); im = a.imshow(M[k][::st, ::st], **kw); outlines(a, lw=0.35); mapax(a); a.set_title(t)
        if k == 'mean':
            a.add_patch(Rectangle((ZOOM[0], ZOOM[2]), ZOOM[1] - ZOOM[0], ZOOM[3] - ZOOM[2], fill=False, ec=INK, lw=0.9, ls='--'))
    a = fig.add_subplot(gs[1, 1]); x0, d = g['ext'][0], g['d']; y1 = g['ext'][3]
    c0, c1 = int((ZOOM[0] - x0) / d), int((ZOOM[1] - x0) / d); r0, r1 = int((y1 - ZOOM[3]) / d), int((y1 - ZOOM[2]) / d)
    a.imshow(M['mean'][r0:r1, c0:c1], extent=(x0 + c0 * d, x0 + c1 * d, y1 - r1 * d, y1 - r0 * d), **{k: v for k, v in kw.items() if k != 'extent'})
    mapax(a); [a.spines[s].set_visible(True) for s in a.spines]; a.set_title('d  30 m detail (box in a)')
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[:, 2])); cb.set_label('Surface-water salinity (psu)'); cb.outline.set_edgecolor(GRIDC)
    fig.get_axes()[-1].set_box_aspect(14)
    save(fig, 'Fig3_example_products')


# ================================================================== Figure 4: validation
def fig4():
    V = validation(); d = V['d']; ST = V['st']; ok = ST.n >= 50
    y, p = d.salinity, d.p; r2 = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum(); rmse = np.sqrt(((y - p) ** 2).mean())
    fig, ax = plt.subplots(2, 2, figsize=(FULL_W, FULL_W * 0.82)); fig.subplots_adjust(wspace=0.3, hspace=0.38)
    a = ax[0, 0]; a.hexbin(y, p, gridsize=70, bins='log', cmap='Greys', mincnt=1, extent=(0, 70, 0, 70), lw=0)
    a.plot([0, 70], [0, 70], color=OBS, lw=0.9); a.set_xlim(-2, 72); a.set_ylim(-2, 72); a.set_xlabel('Gauge salinity (psu)'); a.set_ylabel('Held-out model salinity (psu)')
    a.text(3, 66, f'2013–2025 pooled\nR² = {r2:.2f}\nRMSE = {rmse:.1f} psu\nn = {len(d):,} station-days', va='top', fontsize=7)
    a.set_title('a  Held-out station-days'); clean(a)
    a = ax[0, 1]; a.hist(ST.r[ok], bins=np.arange(0, 1.01, 0.05), color=MOD, edgecolor='white'); a.axvline(ST.r[ok].median(), color=INK, lw=1)
    a.set_xlabel('Within-station correlation r'); a.set_ylabel('Stations'); a.set_title(f'b  Timing: median r = {ST.r[ok].median():.2f}'); clean(a)
    a = ax[1, 0]; k = ST.kge[ok].clip(lower=-1)
    a.hist(k, bins=np.arange(-1.05, 1.01, 0.1), color=MOD, edgecolor='white'); a.axvline(ST.kge[ok].median(), color=INK, lw=1)
    a.set_xlabel('Kling–Gupta efficiency (clipped at −1)'); a.set_ylabel('Stations'); a.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True)); a.set_title(f'c  Skill: median KGE = {ST.kge[ok].median():.2f}'); clean(a)
    a = ax[1, 1]; cov = d.groupby(['y25', 'lev']).inside.mean().unstack(0) * 100; x = np.arange(4)
    a.axhline(90, color=OBS, ls='--', lw=1); a.text(-0.45, 90.4, 'target 90%', color=OBS, ha='left', fontsize=6.5)
    a.scatter(x - 0.1, cov[False], s=30, color=MUTED, label='2013–2024 (calibration)', zorder=3)
    a.scatter(x + 0.1, cov[True], s=30, marker='s', color=MOD, label='2025 (independent)', zorder=3)
    a.set_xticks(x); a.set_xticklabels(LEVEL_LABELS); a.set_ylim(80, 100); a.set_xlim(-0.5, 3.5)
    a.set_ylabel('Observations inside the 90% interval (%)'); a.set_xlabel('Predicted salinity'); a.grid(axis='y', color=GRIDC, lw=0.6)
    a.legend(loc='lower left', fontsize=6.5); a.set_title('d  Interval coverage by predicted level'); clean(a)
    save(fig, 'Fig4_validation')
    print('   half-widths', np.round(V['q'], 2), '| 2025 coverage', np.round(cov[True].values, 1), '| overall 2025', round(100 * d[d.y25].inside.mean(), 1))


# ================================================================== Figure 5: station time series
def fig5():
    V = validation(); d = V['d']; ST = V['st'].set_index('station')
    fig, ax = plt.subplots(len(TS_STATIONS), 1, figsize=(FULL_W, FULL_W * 1.2), sharex=True); fig.subplots_adjust(hspace=0.42)
    for i, (a, s) in enumerate(zip(ax, TS_STATIONS)):
        g = d[d.station == s].set_index('date').sort_index(); full = pd.date_range(g.index.min(), g.index.max())
        m = g.reindex(full)
        a.axvspan(pd.Timestamp('2025-01-01'), pd.Timestamp('2025-12-31'), color='#eeede8', lw=0, zorder=0)
        a.fill_between(full, m.lo, m.hi, color=BAND, lw=0, label='90% interval', zorder=1)
        a.plot(full, m.p, color=MOD, lw=0.7, label='Model (held out)', zorder=2)
        a.scatter(g.index, g.salinity, s=0.25, lw=0, color=OBS, label='Gauge (daily mean)', zorder=3, rasterized=True)
        r = ST.loc[s]
        a.set_title(f'{s} ({BASIN.get(r.basin, r.basin)}, mean {r.obs_mean:.1f} psu)   r {r.r:.2f} · KGE {r.kge:.2f} · RMSE {r.rmse:.1f} · bias {r.bias:+.1f}', fontsize=7)
        a.set_ylabel('psu'); a.set_ylim(bottom=0 if g.salinity.min() < 15 else None); a.grid(axis='y', color=GRIDC, lw=0.5); clean(a)
        if i == 0: a.legend(loc='upper left', ncol=3, fontsize=6.5, markerscale=6)
        for ys in (a.get_yticklabels()): ys.set_fontsize(6.5)
    save(fig, 'Fig5_station_timeseries')


# ================================================================== Figure 6: 250 m vs 30 m
def fig6():
    import rasterio
    from scipy.ndimage import map_coordinates
    g = grid(); x0, y1, d = g['ext'][0], g['ext'][3], g['d']
    c0, c1 = int((ZOOM[0] - x0) / d), int((ZOOM[1] - x0) / d); r0, r1 = int((y1 - ZOOM[3]) / d), int((y1 - ZOOM[2]) / d)
    with rasterio.open(f'{DM}/monthly/ENP_salinity_monthly_mean_2020.tif') as r:
        g250 = r.read(5).astype('float32')                     # band 5 = 2020-05
    g250[~np.isfinite(g250) | (g250 < 0)] = np.nan
    m30 = read_month(f'{S30}/monthly/ENP_salinity30m_monthly_2020_05.tif')[r0:r1, c0:c1]
    # each 250 m cell covers an 8 x 8 block of 30 m cells and its value sits at the block centre (30 m row/col 8k + 4);
    # interpolate bilinearly to the 30 m cells (NaN-aware, weights renormalised over valid nodes)
    rr, cc = np.meshgrid((np.arange(r0, r1) - 4) / 8.0, (np.arange(c0, c1) - 4) / 8.0, indexing='ij')
    ok = np.isfinite(g250); num = map_coordinates(np.where(ok, g250, 0), [rr, cc], order=1, mode='nearest')
    den = map_coordinates(ok.astype('float32'), [rr, cc], order=1, mode='nearest')
    with np.errstate(invalid='ignore', divide='ignore'): i250 = np.where(den > 0.01, num / den, np.nan)
    det = m30 - i250
    e250 = (x0, x0 + g250.shape[1] * 8 * d, y1 - g250.shape[0] * 8 * d, y1)
    ez = (x0 + c0 * d, x0 + c1 * d, y1 - r1 * d, y1 - r0 * d)
    fig = plt.figure(figsize=(FULL_W, FULL_W * 0.30)); gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.045], hspace=0.08, wspace=0.04)
    ax = [fig.add_subplot(gs[0, j]) for j in range(3)]
    im = ax[0].imshow(g250, extent=e250, cmap=SAL, vmin=0, vmax=40, interpolation='nearest'); ax[0].set_title('a  250 m, May 2020')
    ax[1].imshow(m30, extent=ez, cmap=SAL, vmin=0, vmax=40, interpolation='nearest'); ax[1].set_title('b  30 m, May 2020')
    im2 = ax[2].imshow(det, extent=ez, cmap=DIV, norm=TwoSlopeNorm(0, -3, 3), interpolation='nearest'); ax[2].set_title('c  30 m detail, May 2020')
    for a in ax:
        a.set_xlim(ZOOM[0], ZOOM[1]); a.set_ylim(ZOOM[2], ZOOM[3]); mapax(a); [a.spines[s].set_visible(True) for s in a.spines]
    cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, :2]), orientation='horizontal'); cb.set_label('Salinity (psu)'); cb.outline.set_edgecolor(GRIDC)
    cb2 = fig.colorbar(im2, cax=fig.add_subplot(gs[1, 2]), orientation='horizontal'); cb2.set_label('30 m − 250 m (psu)'); cb2.outline.set_edgecolor(GRIDC)
    save(fig, 'Fig6_250m_vs_30m')
    v = det[np.isfinite(det)]; print(f'   detail layer May 2020 (zoom): mean {v.mean():+.2f}, sd {v.std():.2f}, |99th| {np.percentile(np.abs(v), 99):.2f} psu')


# ================================================================== Figure 7: usage examples
def fig7():
    g = grid(); C = pd.read_csv(f'{CH}/monthly_climatology_by_ecosystem.csv')
    E = pd.read_csv(f'{CH}/high_salinity_events_by_region_year.csv'); GV = pd.read_csv(f'{CH}/event_days_gauge_vs_model_30m.csv')
    T = trend_map(4)
    fig = plt.figure(figsize=(FULL_W, FULL_W * 0.95)); gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.32, wspace=0.18)
    a = fig.add_subplot(gs[0, 0])
    for e in ECO:
        c = C[C.ecosystem == e].sort_values('month'); a.plot(c.month, c['mean'], color=COL_E[e], lw=1.6, marker='o', ms=2.8, label=e)
    a.set_xticks(range(1, 13)); a.set_xticklabels(list('JFMAMJJASOND')); a.set_ylabel('Salinity (psu)'); a.set_ylim(0, 52)
    a.grid(axis='y', color=GRIDC, lw=0.5); a.legend(ncol=2, fontsize=6, loc='upper left', handlelength=1.5, columnspacing=0.8)
    a.set_title('a  Seasonal cycle by ecosystem'); clean(a)
    a = fig.add_subplot(gs[0, 1])
    im = a.imshow(T, extent=g['ext'], cmap=DIV, norm=TwoSlopeNorm(0, -0.5, 0.5), interpolation='nearest')
    outlines(a, lw=0.35); mapax(a); a.set_title('b  Trend in monthly anomalies')
    cb = fig.colorbar(im, ax=a, shrink=0.8, pad=0.02); cb.set_label('psu yr$^{-1}$'); cb.outline.set_edgecolor(GRIDC)
    a = fig.add_subplot(gs[1, :])
    for rgn in REG_BY_SAL:
        e = E[E.group == rgn].sort_values('year'); a.plot(e.year, e.event_days_per_pixel, color=COL_R[rgn], lw=1.6, marker='o', ms=3, label=rgn)
    a.set_xticks(sorted(E.year.unique())); a.set_ylabel('Event days per pixel'); a.grid(axis='y', color=GRIDC, lw=0.5)
    a.legend(ncol=3, fontsize=6.5, loc='upper center', bbox_to_anchor=(0.63, 1.0)); clean(a)
    a.set_title('c  Days in high-salinity events (≥ 7 days above the local 95th percentile)')
    gv = GV.groupby('year')[['obs', 'mod']].mean()
    ins = a.inset_axes([0.69, 0.32, 0.29, 0.40]); ins.set_facecolor('white'); ins.set_xticks([2015, 2020, 2025])
    ins.plot(gv.index, gv.obs, color=INK, lw=1.2, marker='o', ms=2.5, label='gauges'); ins.plot(gv.index, gv['mod'], color=OBS, lw=1.2, marker='s', ms=2.5, label='30 m model')
    ins.set_title(f'{GV.station.nunique()} gauges, event days per station', fontsize=6.5); ins.tick_params(labelsize=6); ins.legend(fontsize=6); clean(ins)
    save(fig, 'Fig7_usage_examples')
    r = np.corrcoef(GV.obs, GV['mod'])[0, 1]; print(f'   gauge check: {GV.station.nunique()} stations, {len(GV)} station-years, r = {r:.2f}')


MS = {1: fig1, 2: fig2, 3: fig3, 4: fig4, 5: fig5, 6: fig6, 7: fig7}

if __name__ == '__main__':
    sets = [ARGS.set] if ARGS.set else ['ms', 'doc']
    if 'ms' in sets:
        want = sorted(set(ARGS.figs) or MS)
        print(f'manuscript figures {want} -> {OUT}')
        for k in want: MS[k]()
    if 'doc' in sets:
        here = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, here)
        try:
            import runpy; sys.argv = [f'{here}/../deprecated/two_step_hls_model/python/make_figures_doc.py'] + [str(f) for f in ARGS.figs]
            runpy.run_path(f'{here}/../deprecated/two_step_hls_model/python/make_figures_doc.py', run_name='__main__')
        except (FileNotFoundError, OSError) as e:
            print('document figures skipped:', e, '\n  (they need the intermediate model files listed in make_figures_doc.py; run from the folder that holds them)')
