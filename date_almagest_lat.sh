#!/bin/bash
# Date the Almagest by ecliptic latitudes (SESCC)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Expected result (paper v2): T = -49 (bootstrap median -49, 76% pre-CE)
awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/almagest_verbunt.csv | \
  python3 sescc.py --lat "$@"
