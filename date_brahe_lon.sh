#!/bin/bash
# Date the Tycho Brahe catalogue by ecliptic longitudes (SESCC-pairs)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Pairs: angular separation < 15 deg, residuals as arcs (--cosb), signed speeds (--signed).
# Expected result (paper v2): T = 1572 (bootstrap median 1568)
awk -F';' 'substr($1,1,1)=="#" || $8+0 < 60' catalogs/brahe_verbunt.csv | \
  python3 sescc.py --lon --sep 15 --cosb --signed --pairnorm 0.5 "$@"
