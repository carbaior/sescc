#!/bin/bash
# Date the Almagest by ecliptic longitudes (SESCC-pairs, normalised pair weights)
# Selection rule (same for all catalogues, paper v2): positional error < 60', no individual exclusions.
# Pairs: angular separation < 15 deg, residuals as arcs (--cosb), signed speeds (--signed).
# Expected result (paper v2): T = -240 (bootstrap median -227, 86% pre-CE)
awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/almagest_verbunt.csv | \
  python3 sescc.py --lon --sep 15 --cosb --signed --pairnorm 0.5 "$@"
