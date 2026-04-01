#!/bin/bash
# Date the Almagest by ecliptic latitudes (SESCC)
# Source: Verbunt & van Gent (2012), A&A 544, A31
# Expected result: ~-49 BCE

cat catalogs/almagest_verbunt.csv | python3 sescc.py
