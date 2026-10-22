"""ESSD Fig. 1: study area with the 84 model stations, 6 anchor stations and the 13 DBHydro freshwater sites.
Reuses the drawing helpers of salinity30m/scripts/make_figures.py. Usage: python fig01_study_area.py [root]"""
import os, sys, numpy as np, pandas as pd
ROOT = sys.argv[1] if len(sys.argv) > 1 else '/Volumes/MaloneLab/Research/ENP_Landsat_Salinity_Model'
OUT = f'{ROOT}/best_two_step/paper_final'
os.environ['ENP_DATA_ROOT'] = ROOT; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'python'))
import make_figures as M; M.OUT = OUT
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
FRESH = '#2a78d6'
g = M.grid(); step = 3
eco = g['eco'][::step, ::step].astype('float32'); eco[eco == 255] = np.nan
cmap = M.ListedColormap(M.RAMP + ['#d9d8d2']); norm = M.BoundaryNorm(np.arange(-0.5, 7.5), 7)
s = pd.read_csv(f'{M.SITE}/station_site_variables.csv')
s['role'] = np.where(s.station.isin(M.ANCHORS), 'anchor', np.where(s.grab == 1, 'grab', 'sonde'))
f = pd.read_csv(f'{OUT}/fresh13_coords.csv')
fig, a = plt.subplots(figsize=(M.FULL_W, M.FULL_W * 0.66))
a.imshow(eco, extent=g['ext'], cmap=cmap, norm=norm, interpolation='nearest')
M.outlines(a, step=step)
so = s[s.role == 'sonde']; gr = s[s.role == 'grab']; an = s[s.role == 'anchor']
a.scatter(so.longitude, so.latitude, s=16, facecolor='white', edgecolor=M.INK, lw=0.8, zorder=4)
a.scatter(gr.longitude, gr.latitude, s=12, facecolor=M.INK, edgecolor='white', lw=0.5, zorder=4)
a.scatter(f.lon, f.lat, s=30, marker='D', facecolor='white', edgecolor=M.INK, lw=1.1, zorder=5)
a.scatter(an.longitude, an.latitude, s=80, marker='*', facecolor='#f2c12e', edgecolor=M.INK, lw=0.6, zorder=6)
for _, r in an.iterrows():
    a.annotate(r.station, (r.longitude, r.latitude), xytext=(5, 4), textcoords='offset points', fontsize=6,
               bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.75), zorder=7)
lab = {'Big Cypress': (-81.40, 25.64), 'Shark River\nSlough': (-80.97, 25.68), 'Taylor\nSlough': (-80.625, 25.335),
       'C-111': (-80.48, 25.37), 'Gulf\nwaters': (-81.44, 25.48), 'Florida Bay': (-80.56, 24.96)}
for t, (x, y) in lab.items():
    a.text(x, y, t, ha='center', va='center', fontsize=7.5, bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.7))
h = [Patch(fc=c, label=e) for e, c in zip(M.ECO, M.RAMP)]
h += [Line2D([], [], ls='', marker='o', mfc='white', mec=M.INK, ms=5, label=f'Station, sonde (n = {len(so)})'),
      Line2D([], [], ls='', marker='o', mfc=M.INK, mec='white', ms=4.5, label=f'Station, grab sample (n = {len(gr)})'),
      Line2D([], [], ls='', marker='D', mfc='white', mec=M.INK, mew=1.1, ms=5, label=f'Freshwater site, field conductivity (n = {len(f)})'),
      Line2D([], [], ls='', marker='*', mfc='#f2c12e', mec=M.INK, ms=10, label=f'Anchor station (n = {len(an)})')]
a.legend(handles=h, loc='center left', bbox_to_anchor=(1.01, 0.5), fontsize=6.5, title='Ecosystem (30 m) and stations', title_fontsize=6.5,
         frameon=True, facecolor='white', edgecolor=M.GRIDC, framealpha=0.95)
M.mapax(a, ticks=True); M.clean(a); a.set_xlabel('Longitude (°)'); a.set_ylabel('Latitude (°)')
a.set_xticks([-81.5, -81.25, -81.0, -80.75, -80.5]); a.set_yticks([25.0, 25.25, 25.5, 25.75])
M.save(fig, 'Fig01_study_area')
print('sites plotted', len(so), len(gr), len(an), len(f))
