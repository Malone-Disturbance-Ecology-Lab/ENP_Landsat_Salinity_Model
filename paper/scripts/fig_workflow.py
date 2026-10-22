"""Workflow of the final model and products."""
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
INK='#3d3d3a'; MUTED='#77756e'; BLUE='#2a78d6'; LB='#e3eefb'; OR='#eb6834'; LO='#fde8df'; GR='#f1f0ec'
fig,ax=plt.subplots(figsize=(7.09,4.2)); ax.set_xlim(0,100); ax.set_ylim(0,60); ax.axis('off')
def box(x,y,w,h,title,body,fc,ec):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.4,rounding_size=1.2',fc=fc,ec=ec,lw=1))
    ax.text(x+1,y+h-1.6,title,fontsize=7.2,weight='bold',color=INK,va='top'); ax.text(x+1,y+h-5,body,fontsize=6.1,color=INK,va='top',linespacing=1.35)
def arrow(x0,y0,x1,y1): ax.annotate('',(x1,y1),(x0,y0),arrowprops=dict(arrowstyle='-|>',color=MUTED,lw=1.1))
box(1,41,31,17,'Observations','DBHydro salinity, 2000–2025:\n84 model stations (352 714 days)\n+ 13 freshwater sites (2442 days)\n6 anchor stations (daily predictors)',GR,MUTED)
box(1,21,31,17,'Static 30 m site variables (27)','coast distance, land cover, channels,\nflow-source distances and weights,\nbasin connectivity',GR,MUTED)
box(1,1,31,17,'Daily drivers','EDEN water levels (30- and 90-day),\ngridMET rain and temperature,\nstructure flows, anchor stations,\nsea level (NOAA tide gauges)',GR,MUTED)
box(37,17,25,27,'Boosted-tree model','56 predictors\none model, 2000–2025\n\nvalidation: 10-block\nspatial cross-validation\n\nsplit-conformal\n90 % intervals',LB,BLUE)
for y in (49,29,9): arrow(32.6,y,36.6,30.5)
box(67,40,32,18,'Daily maps, 2000–2025','every 30 m cell, every day\n(April 2000–December 2025)\n250 m nodes with 90 % bounds\nmonthly means',LB,BLUE)
box(67,20,32,17,'Daily surface-water mask','EDEN water level vs DEM bounds\n(lowest 10 m point per 30 m cell):\ndry, wet or uncertain;\nstatic classes elsewhere',GR,MUTED)
box(67,1,32,16,'Masked products and analyses','wet-day monthly means, climatology,\nseasonal and interannual patterns,\ntrends, high-salinity events',LO,OR)
arrow(62.6,33,66.6,48); arrow(83,39.5,83,37.5); arrow(83,19.5,83,17.5)
fig.savefig('Fig02_workflow.png',dpi=300,bbox_inches='tight',facecolor='white'); print('ok')
