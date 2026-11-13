#!/bin/bash
# Date the Ulugh Beg catalogue by ecliptic longitudes (SESCC-pairs)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Expected result (paper v2): T = 1462 (true epoch 1437)
awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/ulughbeg_verbunt.csv | \
  python3 sescc.py --lon --dlon 30 --pairnorm 0.5 "$@"
