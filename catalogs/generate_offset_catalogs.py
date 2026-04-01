#!/usr/bin/env python3
# generate_offset_catalogs.py
# Generates 13 versions of a star catalog with longitude offsets
# from -6 to +6 degrees in steps of 1 degree.
# Offsets are applied only to full degrees to preserve fractional parts.
#
# Usage: python3 generate_offset_catalogs.py catalogs/almagest.csv
# Output: catalogs/almagest_offset_-6.csv ... almagest_offset_+6.csv

import sys
import os

if len(sys.argv) != 2:
    print("Usage: generate_offset_catalogs.py INPUT_CATALOG")
    sys.exit(1)

infile = sys.argv[1]
basedir = os.path.dirname(infile)
basename = os.path.splitext(os.path.basename(infile))[0]

# Read all lines
with open(infile, 'r', encoding='utf-8') as f:
    lines = f.readlines()

offsets = list(range(-6, 7))  # -6 to +6

for offset in offsets:
    sign = '+' if offset >= 0 else ''
    outname = f"{basename}_offset_{sign}{offset}.csv"
    outpath = os.path.join(basedir, outname)

    with open(outpath, 'w', encoding='utf-8') as out:
        for line in lines:
            stripped = line.strip()
            if not stripped:
                continue

            parts = stripped.split(';')

            # Try to parse HIP — skip header
            try:
                hip = int(parts[0])
            except ValueError:
                out.write(line)
                continue

            # Skip unidentified entries (hip=0)
            if hip == 0:
                out.write(line)
                continue

            # Parse longitude (column 2)
            try:
                lon = float(parts[2])
            except (ValueError, IndexError):
                out.write(line)
                continue

            # Apply offset to longitude (full degrees only, preserves fractional part)
            new_lon = lon + offset

            # Wrap to [0, 360)
            new_lon = new_lon % 360.0

            # Reconstruct line with modified longitude
            parts[2] = str(new_lon)
            out.write(';'.join(parts) + '\n')

    print(f"Generated: {outpath}")

print(f"\nDone. Generated {len(offsets)} catalogs.")
print(f"Offsets: {offsets} degrees")
print(f"Granularity: ~72 years per degree (at 1.4°/century precession rate)")
