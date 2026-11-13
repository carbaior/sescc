#!/bin/bash
# Date the Ulugh Beg catalogue by ecliptic longitudes (SESCC-pairs)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Pairs: angular separation < 15 deg, residuals as arcs (--cosb), signed speeds (--signed).
# Expected result (paper v2): T = 1457 (bootstrap median 1451)
awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/ulughbeg_verbunt.csv | \
  python3 sescc.py --lon --sep 15 --cosb --signed --pairnorm 0.5 "$@"
