#!/usr/bin/env python3
# verbunt2sescc.py
#
# Converts Verbunt & van Gent (2012) catalogue files to SESCC input format.
#
# Usage:
#   python3 verbunt2sescc.py ptolema.dat   > almagest_verbunt.csv
#   python3 verbunt2sescc.py ulughbeg.dat  > ulughbeg_verbunt.csv
#
# Output format (semicolon-separated):
#   hip ; lat_deg ; lon_deg ; seq ; cst ; name ; mag ; id_flag ; dist_arcmin
#
# Notes:
#   - Stars with HIP=0 (unidentified) are skipped
#   - Stars with id_flag=5 (not identified) are skipped
#   - Longitude in ptolema.dat uses equinox -128 (Hipparcos epoch)
#     per Verbunt note (1): Elon = (Elon.Z - 1)*30 + Elon.d + Elon.m/60
#   - Longitude in ulughbeg.dat uses equinox 1437
#     per Verbunt note (1): Elon = Elon.Z*30 + Elon.d + Elon.m/60
#   - Latitude: positive for boreal (B), negative for austral (A)

import sys, os

def parse_ptolema(line):
    """Parse one line of ptolema.dat. Returns dict or None to skip."""
    if len(line) < 48:
        return None

    try:
        seq    = int(line[0:4])
        c      = int(line[5:7])
        cst    = line[8:12].strip()
        n_star = line[13:16].strip()

        # Longitude: (Elon.Z - 1)*30 + Elon.d + Elon.m/60
        elon_z = int(line[18:20])
        elon_d = int(line[21:23])
        elon_m = int(line[24:26])
        lon    = (elon_z - 1) * 30.0 + elon_d + elon_m / 60.0

        # Latitude: sign from col 35 (B=boreal=+, A=austral=-)
        elat_d = int(line[28:30])
        elat_m = int(line[31:33])
        elat_s = line[34].strip()
        lat    = elat_d + elat_m / 60.0
        if elat_s == 'A':
            lat = -lat

        # Magnitude
        try:
            mag_str = line[37].strip()
            mag_q   = line[38].strip() if len(line) > 38 else ''
            mag     = int(mag_str) if mag_str.isdigit() else 0
        except:
            mag = 0

        # HIP
        hip_str = line[40:46].strip()
        hip     = int(hip_str) if hip_str and hip_str != '0' else 0

        # ID flag
        id_flag_str = line[47].strip() if len(line) > 47 else '0'
        id_flag = int(id_flag_str) if id_flag_str.isdigit() else 0

        # Vmag
        try:
            vmag = float(line[51:55])
        except:
            vmag = 0.0

        # Positional errors
        try:
            dl   = float(line[56:62])
            db   = float(line[63:69])
            dist = float(line[70:76])
        except:
            dl = db = dist = 0.0

        name = f"{n_star}{cst}"

        return {
            'seq': seq, 'cst': cst, 'name': name,
            'lon': lon, 'lat': lat,
            'hip': hip, 'id_flag': id_flag,
            'mag': mag, 'vmag': vmag,
            'dl': dl, 'db': db, 'dist': dist,
        }
    except Exception as e:
        return None


def parse_ulughbeg(line):
    """Parse one line of ulughbeg.dat. Returns dict or None to skip."""
    if len(line) < 56:
        return None

    try:
        seq    = int(line[0:4])
        n_u    = line[4].strip()    # 'c' if copied from Ptolemaios via al-Sufi
        p_seq_str = line[6:10].strip()
        p_seq  = int(p_seq_str) if p_seq_str else 0
        c      = int(line[12:14])
        cst    = line[15:19].strip()
        n_star = line[20:23].strip()

        # Longitude: Elon.Z*30 + Elon.d + Elon.m/60
        elon_z = int(line[25:27])
        elon_d = int(line[28:30])
        elon_m = int(line[31:33])
        lon    = elon_z * 30.0 + elon_d + elon_m / 60.0

        # Latitude
        elat_d = int(line[36:38])
        elat_m = int(line[39:41])
        elat_s = line[42].strip()
        lat    = elat_d + elat_m / 60.0
        if elat_s == 'A':
            lat = -lat

        # Magnitude
        try:
            mag_str = line[45].strip()
            mag     = int(mag_str) if mag_str.isdigit() else 0
        except:
            mag = 0

        # HIP
        hip_str = line[48:54].strip()
        hip     = int(hip_str) if hip_str and hip_str != '0' else 0

        # ID flag
        id_flag_str = line[55].strip() if len(line) > 55 else '0'
        id_flag = int(id_flag_str) if id_flag_str.isdigit() else 0

        # Vmag
        try:
            vmag = float(line[61:65])
        except:
            vmag = 0.0

        # Positional errors
        try:
            dl   = float(line[66:72])
            db   = float(line[73:79])
            dist = float(line[80:86])
        except:
            dl = db = dist = 0.0

        copied = (n_u == 'c')
        name   = f"{n_star}{cst}"

        return {
            'seq': seq, 'p_seq': p_seq, 'cst': cst, 'name': name,
            'lon': lon, 'lat': lat,
            'hip': hip, 'id_flag': id_flag,
            'mag': mag, 'vmag': vmag,
            'dl': dl, 'db': db, 'dist': dist,
            'copied': copied,
        }
    except Exception as e:
        return None


def convert(filename):
    basename = os.path.basename(filename).lower()
    is_ptolema   = 'ptolema'  in basename
    is_ulughbeg  = 'ulugh'    in basename or 'ulug' in basename

    if not (is_ptolema or is_ulughbeg):
        # Try to detect from content
        with open(filename) as f:
            first = f.readline()
        is_ptolema  = len(first.rstrip()) <= 78
        is_ulughbeg = not is_ptolema

    parser = parse_ptolema if is_ptolema else parse_ulughbeg
    cat    = 'Ptolemaios (equinox -128)' if is_ptolema else 'Ulugh Beg (equinox 1437)'

    # Header comment
    print(f"# Converted from {filename} — {cat}")
    print(f"# Source: Verbunt & van Gent (2012) A&A 544 A31")
    print(f"# id_flag: 1=secure(nearest) 2=secure(not nearest) 3=probable 4=possible 5=unident 6=repeat")
    print(f"# Skipped: HIP=0 (unidentified) and id_flag=5 (not identified)")
    print(f"# hip;lat_deg;lon_deg;seq;cst;name;mag_ptolema;vmag_hip;id_flag;dist_arcmin")

    n_ok = n_skip_hip = n_skip_flag = n_parse_err = 0

    with open(filename, 'r', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line.strip():
                continue

            record = parser(line)

            if record is None:
                n_parse_err += 1
                continue

            # Skip unidentified
            if record['hip'] == 0:
                n_skip_hip += 1
                continue
            if record['id_flag'] == 5:
                n_skip_flag += 1
                continue

            # Format output
            # Required: hip ; lat ; lon  (first three fields)
            # Then: seq ; cst ; name ; mag ; vmag ; id_flag ; dist
            print(
                f"{record['hip']};"
                f"{record['lat']:.4f};"
                f"{record['lon']:.4f};"
                f"{record['seq']};"
                f"{record['cst']};"
                f"\"{record['name']}\";"
                f"{record['mag']};"
                f"{record['vmag']:.1f};"
                f"{record['id_flag']};"
                f"{record['dist']:.1f}"
            )
            n_ok += 1

    print(f"# Written: {n_ok} stars  |  Skipped: {n_skip_hip} (HIP=0)  "
          f"{n_skip_flag} (unidentified)  {n_parse_err} (parse errors)",
          file=sys.stderr)


if __name__ == '__main__':
    if len(sys.argv) < 2 or sys.argv[1] in ('-h', '--help'):
        print("Usage: python3 verbunt2sescc.py ptolema.dat   > almagest_verbunt.csv")
        print("       python3 verbunt2sescc.py ulughbeg.dat  > ulughbeg_verbunt.csv")
        sys.exit(0)
    convert(sys.argv[1])
