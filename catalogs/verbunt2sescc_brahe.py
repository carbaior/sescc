#!/usr/bin/env python3
"""
verbunt2sescc_brahe.py - Convert KeplerE.dat (Verbunt & van Gent 2010, A&A 516, A28)
to SESCC input format (semicolon-separated CSV).

Source: Tycho Brahe star catalogue, Kepler (1627) emended edition.
Coordinates are ecliptic, equinox AD 1601.0.

Format of KeplerE.dat (fixed-width, byte positions from ReadMe):
  1-  4  M        Sequence number in Manuscript
  6-  8  B        Sequence number in Brahe 1602 edition
 10- 13  K        Sequence number in Kepler 1627 edition
 15- 16  C        Constellation number
     18  [=]
 19- 21  cst      Constellation abbreviation
 23- 24  N        Star number in constellation
 26- 27  Elon.Z   Zodiacal sign (1-12), lambda = (Z-1)*30 + deg + min/60
 29- 30  Elon.d   Degrees of ecliptic longitude
 32- 35  Elon.m   Arcminutes of ecliptic longitude
 37- 38  Elat.d   Degrees of ecliptic latitude
 40- 43  Elat.m   Arcminutes of ecliptic latitude
     45  Elat.-   B=Borealis (+), A=Australis (-)
     47  VB       Magnitude
 50- 55  HIP      Hipparcos number
     58  I        Identification quality (1=secure, 2=secure alt, 3=probable,
                  4=possible, 5=not identified, 6=double entry)
 65- 68  Vmag     Hipparcos V magnitude
 70- 75  Dlon     Longitude difference Hipparcos-Brahe (arcmin)
 77- 82  Dlat     Latitude difference Hipparcos-Brahe (arcmin)
 84- 89  Delta    Total offset (arcmin)

Output CSV columns:
  hip ; lat_deg ; lon_deg ; mag_brahe ; vmag ; dlon_arcmin ; dlat_arcmin ; dist_arcmin ; id_flag

Usage:
  python3 verbunt2sescc_brahe.py keplere.dat > catalogs/brahe_verbunt.csv

copyleft (GPLv3) 2024 Carlos Baiget Orts (asinfreedom@gmail.com)
"""

import sys

def parse_keplere(filename):
    results = []
    with open(filename, 'r', encoding='latin-1') as f:
        for lineno, line in enumerate(f, 1):
            if len(line) < 55:
                continue

            # Skip lines with no HIP (unidentified)
            hip_str = line[49:55].strip()
            if not hip_str or hip_str == '0':
                continue
            try:
                hip = int(hip_str)
            except ValueError:
                continue
            if hip == 0:
                continue

            # Identification quality flag (byte 58, 0-indexed 57)
            try:
                id_flag = int(line[57].strip())
            except (ValueError, IndexError):
                id_flag = 0

            # Skip not identified (5) or double entry (6)
            if id_flag in (5, 6):
                continue

            # Ecliptic longitude
            try:
                z    = int(line[25:27].strip())   # zodiacal sign 1-12
                elon_d = int(line[28:30].strip())  # degrees
                elon_m = float(line[31:35].strip()) # arcminutes
            except (ValueError, IndexError):
                continue
            lon_deg = (z - 1) * 30.0 + elon_d + elon_m / 60.0

            # Ecliptic latitude
            try:
                elat_d = int(line[36:38].strip())
                elat_m_str = line[39:43].strip()
                elat_m = float(elat_m_str) if elat_m_str else 0.0
                sign_char = line[44].strip()
            except (ValueError, IndexError):
                continue
            sign = -1.0 if sign_char == 'A' else 1.0
            lat_deg = sign * (elat_d + elat_m / 60.0)

            # Brahe magnitude
            try:
                mag_brahe = int(line[46].strip())
            except (ValueError, IndexError):
                mag_brahe = 0

            # Hipparcos magnitude
            try:
                vmag_str = line[64:68].strip()
                vmag = float(vmag_str) if vmag_str and vmag_str != '0' else 0.0
            except (ValueError, IndexError):
                vmag = 0.0

            # Positional differences
            try:
                dlon_str = line[69:75].strip()
                dlon = float(dlon_str) if dlon_str else 0.0
            except (ValueError, IndexError):
                dlon = 0.0

            try:
                dlat_str = line[76:82].strip()
                dlat = float(dlat_str) if dlat_str else 0.0
            except (ValueError, IndexError):
                dlat = 0.0

            try:
                delta_str = line[83:89].strip()
                delta = float(delta_str) if delta_str else 0.0
            except (ValueError, IndexError):
                delta = 0.0

            results.append((hip, lat_deg, lon_deg, mag_brahe, vmag,
                            dlon, dlat, delta, id_flag))
    return results

def main():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} keplere.dat", file=sys.stderr)
        sys.exit(1)

    filename = sys.argv[1]
    stars = parse_keplere(filename)

    print('# Tycho Brahe star catalogue (Kepler 1627 emended edition)')
    print('# Source: Verbunt & van Gent (2010), A&A 516, A28')
    print('# Coordinates: ecliptic, equinox AD 1601.0')
    print('# Columns: hip;lat_deg;lon_deg;mag_brahe;vmag;dlon_arcmin;dlat_arcmin;dist_arcmin;id_flag')
    print(f'# Stars written: {len(stars)}')
    print('#')

    for (hip, lat, lon, mag_b, vmag, dlon, dlat, delta, idflag) in stars:
        print(f'{hip};{lat:.4f};{lon:.4f};{mag_b};{vmag:.1f};'
              f'{dlon:.1f};{dlat:.1f};{delta:.1f};{idflag}')

if __name__ == '__main__':
    main()
