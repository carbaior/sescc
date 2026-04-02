#!/bin/bash
# Date the Tycho Brahe catalogue by ecliptic latitudes (SESCC)
# Source: Verbunt & van Gent (2010), A&A 516, A28 (KeplerE edition)
# True epoch: ~1580 CE — Expected result: ~1570 CE

cat catalogs/brahe_verbunt.csv | python3 sescc.py --lat
