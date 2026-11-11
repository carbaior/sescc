#!/usr/bin/env python3
"""
¿El error de longitud del catálogo es constante en grados de λ o en ángulo
sobre la esfera (es decir, crece como 1/cos β)?

Usa las diferencias catálogo − Hipparcos de Verbunt & van Gent:
  Δλ (en λ, sin escalar), Δβ y la distancia total, que cumple
  dist ≈ sqrt((Δλ cos β)^2 + Δβ^2).
No interviene ninguna datación.

Uso (desde sescc-main):
  python3 check_noise.py catalogs/almagest.dat catalogs/almagest_verbunt.csv
  python3 check_noise.py catalogs/ulughbeg.dat catalogs/ulughbeg_verbunt.csv
  python3 check_noise.py --brahe catalogs/brahe_verbunt.csv
"""
import sys, math, statistics

def rsig(v):
    """Desviación robusta: 1.4826 × mediana de |x − mediana|."""
    if len(v) < 3: return float('nan')
    m = statistics.median(v)
    return 1.4826 * statistics.median([abs(x - m) for x in v])

def load_dat(path):
    """Últimos tres números de cada fila del .dat: Δλ, Δβ, dist."""
    d = {}
    for line in open(path, encoding='utf-8', errors='replace'):
        t = line.split()
        if not t or not t[0].isdigit(): continue
        nums = []
        for x in reversed(t):
            x = x.rstrip('*')
            if not x: continue
            try: nums.append(float(x))
            except ValueError: break
            if len(nums) == 3: break
        if len(nums) == 3:
            d[int(t[0])] = (nums[2], nums[1], nums[0])
    return d

rows = []                                   # (beta, dlon, dlat, dist)
if sys.argv[1] == '--brahe':
    for line in open(sys.argv[2], encoding='utf-8'):
        if line.startswith('#') or not line.strip(): continue
        f = line.strip().split(';')
        rows.append((float(f[1]), float(f[5]), float(f[6]), float(f[7])))
else:
    dat = load_dat(sys.argv[1])
    for line in open(sys.argv[2], encoding='utf-8'):
        if line.startswith('#') or not line.strip(): continue
        f = line.strip().split(';')
        seq = int(f[3])
        if seq in dat:
            dl, db, dist = dat[seq]
            rows.append((float(f[1]), dl, db, dist))

# Comprobación de que las columnas se han leído bien
ok = sum(1 for b, dl, db, d in rows
         if abs(d - math.hypot(dl*math.cos(math.radians(b)), db)) < 1.0)
print(f"{len(rows)} filas; columnas coherentes (dist = hypot(Δλ cos β, Δβ)) en {ok}")

rows = [r for r in rows if r[3] < 60]
print(f"{len(rows)} filas con error < 60'\n")
print(f"{'|β|':>9} {'N':>5} {'σ(Δλ)':>8} {'σ(Δλ·cosβ)':>11} {'σ(Δβ)':>7}   (arcmin, robustas)")
for lo, hi in ((0,20), (20,40), (40,60), (60,75), (75,91)):
    s = [r for r in rows if lo <= abs(r[0]) < hi]
    if not s: continue
    dl  = [r[1] for r in s]
    dlc = [r[1]*math.cos(math.radians(r[0])) for r in s]
    db  = [r[2] for r in s]
    print(f"{lo:>3}–{hi if hi<91 else 90:<3}° {len(s):>5} {rsig(dl):8.1f} {rsig(dlc):11.1f} {rsig(db):7.1f}")
