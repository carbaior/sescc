#!/bin/bash
# Date the Ulugh Beg catalogue by ecliptic latitudes (SESCC)
# Source: Verbunt & van Gent (2012), A&A 544, A31
# True epoch: 1437 CE — Expected result: ~1177 CE

awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/ulughbeg_verbunt.csv | \
  grep -v "^19849;" | grep -v "^71681;" | \
  python3 sescc.py
