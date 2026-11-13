#!/bin/bash
# Date the Almagest by ecliptic longitudes (SESCC-pairs, normalised pair weights)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Expected result (paper v2): T = -244 (bootstrap median -174, 90% pre-CE)
awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/almagest_verbunt.csv | \
  python3 sescc.py --lon --dlon 30 --pairnorm 0.5 "$@"
