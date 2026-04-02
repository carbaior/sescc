#!/bin/bash
# Date the Tycho Brahe catalogue by ecliptic longitudes (SESCC-pairs)
# Source: Verbunt & van Gent (2010), A&A 516, A28 (KeplerE edition)
# True epoch: ~1580 CE — Expected result: ~1547 CE
# Filters: positional error < 60', excluding Keid (HIP 19849) and alpha Cen (HIP 71681)

awk -F';' 'substr($1,1,1)=="#" || $8+0 < 60' catalogs/brahe_verbunt.csv | \
  grep -v "^19849;" | grep -v "^71681;" | \
  python3 sescc.py --lon --dlon 30
