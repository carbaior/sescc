#!/usr/bin/env python3
"""
Peso efectivo de cada estrella en SESCC-pairs (longitudes) con el criterio
de vecindad angular, para actualizar la sección "Stars shared by many pairs".

Peso efectivo de un par en la mediana ponderada (Apéndice A):
    con residuos angulares (--cosb):  W_ij = dv_ij^2 / sqrt(n_i n_j)
    sin --cosb:                       W_ij = dv_ij^2 / cos(b_ij) / sqrt(n_i n_j)
donde dv_ij = |v_i - v_j| y v = mu_lambda cos(beta) es la velocidad angular
en longitud CON SIGNO (como --signed). Con --unsigned usa ||v_i| - |v_j||.
Con --nonorm quita el factor 1/sqrt(n_i n_j).

Las posiciones para formar los pares son las de Hipparcos en la eclíptica
J2000, no las del catálogo, y el movimiento propio es lineal: es una
comprobación aproximada. El número de pares debe salir cercano al de sescc.py.

Uso (desde sescc-main):
    python3 check_pairweights.py catalogs/almagest_verbunt.csv [--sep 15]
           [--maxerr 60] [--nocosb] [--unsigned] [--nonorm]
"""
import sys, math
from collections import defaultdict

argv = sys.argv[1:]
def opt(name, default):
    return float(argv[argv.index(name) + 1]) if name in argv else default
cat    = [a for i, a in enumerate(argv)
          if not a.startswith("--") and (i == 0 or argv[i-1] not in ("--sep", "--maxerr"))][0]
SEP    = opt("--sep", 15.0)
MAXERR = opt("--maxerr", 60.0)
COSB   = "--nocosb" not in argv
SIGNED = "--unsigned" not in argv
NORM   = "--nonorm" not in argv
hipfile = "hip_main.dat"

eps = math.radians(23.4392911); ce, se = math.cos(eps), math.sin(eps)
def ecl(ra, de, pa, pd):
    """lambda, beta (grados) y v_lambda = mu_lambda cos b, mu_beta (mdeg/milenio)."""
    a, d = math.radians(ra), math.radians(de)
    x  = (math.cos(d)*math.cos(a), math.cos(d)*math.sin(a), math.sin(d))
    ea = (-math.sin(a), math.cos(a), 0.0)
    ed = (-math.sin(d)*math.cos(a), -math.sin(d)*math.sin(a), math.cos(d))
    pa_r, pd_r = math.radians(pa/3.6e6), math.radians(pd/3.6e6)
    v  = [pa_r*ea[k] + pd_r*ed[k] for k in range(3)]
    rot = lambda u: (u[0], ce*u[1] + se*u[2], -se*u[1] + ce*u[2])
    X, V = rot(x), rot(v)
    lam = math.atan2(X[1], X[0]); bet = math.asin(max(-1.0, min(1.0, X[2])))
    el = (-math.sin(lam), math.cos(lam), 0.0)
    eb = (-math.sin(bet)*math.cos(lam), -math.sin(bet)*math.sin(lam), math.cos(bet))
    k = math.degrees(1.0) * 1e6                     # rad/año -> mdeg/milenio
    vl = sum(V[i]*el[i] for i in range(3)) * k
    vb = sum(V[i]*eb[i] for i in range(3)) * k
    return math.degrees(lam) % 360, math.degrees(bet), vl, vb

hip = {}
for line in open(hipfile, encoding="latin-1"):
    f = line.split("|")
    try:
        hip[int(f[1])] = (float(f[8]), float(f[9]), float(f[12]), float(f[13]))
    except (ValueError, IndexError):
        continue

stars, seen = [], set()
for line in open(cat, encoding="utf-8"):
    if line.startswith("#") or not line.strip(): continue
    f = line.rstrip("\n").split(";")
    if float(f[9]) >= MAXERR: continue
    h = int(f[0])
    if h not in hip or h in seen: continue          # sin astrometría o repetida
    seen.add(h)
    lam, bet, vl, vb = ecl(*hip[h])
    stars.append((h, f[5].strip('"'), lam, bet, vl))

pairs, npairs, npole = [], defaultdict(int), 0
N = len(stars)
for i in range(N):
    hi, _, li, bi, vi = stars[i]
    for j in range(i + 1, N):
        hj, _, lj, bj, vj = stars[j]
        dl = (li - lj + 180) % 360 - 180
        b1, b2 = math.radians(bi), math.radians(bj)
        cs = math.sin(b1)*math.sin(b2) + math.cos(b1)*math.cos(b2)*math.cos(math.radians(dl))
        if math.degrees(math.acos(max(-1.0, min(1.0, cs)))) > SEP: continue
        cb = math.cos((b1 + b2) / 2)
        if abs(dl) * cb > SEP: npole += 1; continue
        dv = abs(vi - vj) if SIGNED else abs(abs(vi) - abs(vj))
        w = dv * dv if COSB else dv * dv / max(cb, 0.02)
        pairs.append((i, j, w)); npairs[i] += 1; npairs[j] += 1

W = [(i, j, w / (math.sqrt(npairs[i]*npairs[j]) if NORM else 1.0)) for i, j, w in pairs]
tot = sum(w for _, _, w in W)
share = defaultdict(float)
for i, j, w in W:
    share[i] += w / 2; share[j] += w / 2           # mitad del par a cada estrella
neff = tot**2 / sum(w*w for _, _, w in W)

print(f"{cat}: {N} estrellas distintas, {len(pairs)} pares "
      f"(sep {SEP:g} grados; {npole} descartados cerca del polo)")
print(f"residuos {'angulares' if COSB else 'en lambda'}, velocidades "
      f"{'con signo' if SIGNED else 'sin signo'}, "
      f"{'normalizado 1/sqrt(n n)' if NORM else 'sin normalizar'}")
print(f"N_eff de pares = (sum W)^2 / sum W^2 = {neff:.1f}\n")
print("   HIP    entrada      beta   v_lambda  pares  peso %")
for i in sorted(share, key=lambda k: -share[k])[:10]:
    h, name, lam, bet, vl = stars[i]
    print(f"{h:7d}  {name:12s} {bet:6.1f} {vl:9.0f}  {npairs[i]:5d}  {100*share[i]/tot:6.1f}")
