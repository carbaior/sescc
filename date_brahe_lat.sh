#!/bin/bash
# Date the Tycho Brahe catalogue by ecliptic latitudes (SESCC)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Expected result (paper v2): T = 1576 (observations 1576-1597)
awk -F';' 'substr($1,1,1)=="#" || $8+0 < 60' catalogs/brahe_verbunt.csv | \
  python3 sescc.py --lat "$@"
