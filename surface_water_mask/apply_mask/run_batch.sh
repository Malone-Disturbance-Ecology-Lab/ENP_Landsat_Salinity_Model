#!/bin/bash
# runs the remaining months, 4 at a time, for up to $1 seconds
D="$(cd "$(dirname "$0")/../.." && pwd)"; T=${1:-160}
todo=$(python3 -c "
import pandas as pd, os
M=pd.date_range('2000-04-01','2025-12-01',freq='MS')
print(' '.join(str(i) for i,m in enumerate(M) if not os.path.exists(f'$D/daily30m_masked_final/{m.year}/ENP_salinity30m_masked_daily_{m.year}_{m.month:02d}.tif')))")
echo $todo | tr ' ' '\n' | grep . | timeout $T xargs -P 6 -n 1 python3 "$D/daily30m_masked_final/scripts/mask_daily30.py" "$D" | grep -c done
echo left $(( $(echo $todo | wc -w) )) before this round
