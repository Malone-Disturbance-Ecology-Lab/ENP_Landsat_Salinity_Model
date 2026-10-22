# ESSD manuscript: final-model version (2026-10-14)

Decision (Sparkle, 2026-10-14): the paper presents only the final one-step model (no vegetation spectra; sea level; trained
2000–2025 with 13 freshwater sites) and its daily 30 m record, April 2000–December 2025. No comparisons with past or
alternative models. Results focus on spatial and temporal salinity patterns. The HLS two-step product is not in the paper.

Doc: https://claude.ai/code/artifact/915b249d-e85c-4128-abab-7b59628dc098 (rev 89). Text copy and all figure/statistics scripts:
best_two_step/paper_final/ on the lab server.

Structure: 1 Introduction · 2 Study area · 3 Data and methods (observations, 56 predictors, model, spatial-block CV and
split-conformal intervals, mapping, surface-water mask, pattern analyses incl. drivers) · 4 Technical validation (held-out skill,
intervals, mask vs gauges, gauge trends and events, limitations) · 5 Spatial and temporal patterns (spatial structure, seasonal
cycle, interannual variation, trends, events, 5.6 drivers of modelled salinity, 5.7 further uses) · 6–8 data/code availability,
conclusions. Tables 1–8, Figs 1–10.

Key numbers: held-out R² 0.85, RMSE 5.3 psu, bias −0.1 (352 714 station-days, 84 stations); R² 0.82–0.86 by period;
station median r 0.87, KGE 0.71; cross-calibrated 90 % interval coverage 89.5 %. Median wet-day salinity 0.3 (freshwater marsh),
12.2 (mangrove), 35.1 psu (open bay). Peak in May everywhere; seasonal range up to 13.9 psu (mangrove). Saltiest years 2015, 2011;
May salinity vs Dec–May rain r −0.65 to −0.69. Trends: 12.7 % of cells significantly increasing, mostly freshwater marsh at
< 0.01 psu/yr. Events: 2015 most extensive (Florida Bay 103 event days per pixel); gauge event days r 0.88.

Drivers (Sect. 5.6, Fig. 10): grouped permutation importance (held out, base RMSE 5.24): site variables +4.25 psu, anchor
stations +4.05, climate +0.33, water level +0.21, sea level +0.12, flows +0.03, season +0.02; hydroclimate together +0.71; all
daily drivers together +5.83 (groups overlap). Group response 5th→95th percentile at 5–25 psu stations: anchors +22.2, climate
−4.5, water level −3.5, sea level +2.9, flows −0.4 psu. Spatial attribution (500 m nodes, 15th of each month): anchors dominate
86 % of nodes with SD ≥ 0.5 psu, water level 13 % (freshwater marsh / Taylor Slough), climate 1 %. Files: paper_final/drivers_*.

Released daily record (decision 2026-10-15: publish the MASKED daily 30 m data): best_two_step/daily30m_masked_final/YYYY/
ENP_salinity30m_masked_daily_YYYY_MM.tif, 309 files, 12.3 GB, int16 psu×10, −1 = no data or mask code 0/9; uncertain (2) kept.
Made on the lab server by daily30m_masked_final/scripts/mask_daily30.py; checked exactly against monthly30m_wetdays (7 months).
Table 1 now lists this file set. Unmasked daily30m_onestep_final/ stays on the server (not in Table 1). _to_delete/ in that
folder holds 1.4 GB of interrupted .part files.

Other products: climatology_final/, monthly30m_wetdays/ (rebuilt from final daily maps). Open CONFIRM items: co-authors, DOIs/
review link, funding, acknowledgements, NPS and DBHydro citations, file renaming ("onestep" in the monthly/250 m names).
