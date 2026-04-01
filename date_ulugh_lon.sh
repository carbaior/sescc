#!/bin/bash
# Date the Ulugh Beg catalogue by ecliptic longitudes (SESCC-pairs)
# Source: Verbunt & van Gent (2012), A&A 544, A31
# True epoch: 1437 CE — Expected result: ~1452 CE
# Filters: positional error < 60', excluding Keid (HIP 19849) and alpha Cen (HIP 71681)

awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/ulughbeg_verbunt.csv | \
  grep -v "^19849;" | grep -v "^71681;" | \
  python3 sescc_pairs.py --dlon 30
