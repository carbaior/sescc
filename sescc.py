#!/usr/bin/env python3
"""
sescc.py - SESCC: Speed-Error Signals Cross-Correlation
Dates ancient star catalogues using ecliptic latitudes (SESCC-abs).

C(T) = Σᵢ |µβᵢ| · |βᵢᶜᵃᵗ − βᵢᵐᵒᵈ(T)|

Precession does not affect ecliptic latitudes, so no reference correction
is needed. No reference star required.

For longitude dating use sescc_pairs.py instead.

Usage:
  cat catalog.csv | python3 sescc.py [OPTIONS]

Options:
  --bootstrap N   run N bootstrap resamples and report distribution (default 0)
  --subset N      random subset of N stars
  --maxmag M      exclude stars fainter than magnitude M
  --noplot        suppress plot
  -h/--help       this help

Cache: sescc_positions.pkl.xz (unified LZMA pickle)
Keys: "{hip}_{epoch_index}_lat"
  epoch_index 0  -> year 1900
  epoch_index N  -> year 1900 - N * 1

copyleft (GPLv3) 2024 Carlos Baiget Orts (asinfreedom@gmail.com)
https://github.com/carbaior/sescc
"""

import sys, os, math, random, pickle, lzma, urllib.request
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import Star, load
from skyfield.data import hipparcos

# ── Configuration ──────────────────────────────────────────────────────────────
CACHE_URL  = "https://example.com/sescc_positions.pkl.xz"
CACHE_FILE = "./sescc_positions.pkl.xz"

siglos     = 25
resolucion = 1
maxt       = siglos * 100 // resolucion
fechamax   = 1900

# ── Parse arguments ────────────────────────────────────────────────────────────
do_plot    = True
filtro     = 0
maxmag     = 100
n_boot     = 0       # bootstrap resamples (0 = disabled)

def usage():
    print(__doc__)
    sys.exit(0)

args = sys.argv[1:]
if '-h' in args or '--help' in args: usage()
if '--noplot' in args: do_plot = False; args.remove('--noplot')
if '--subset' in args:
    idx = args.index('--subset'); filtro = int(args[idx+1])
    args.pop(idx+1); args.pop(idx)
if '--maxmag' in args:
    idx = args.index('--maxmag'); maxmag = float(args[idx+1])
    args.pop(idx+1); args.pop(idx)
if '--bootstrap' in args:
    idx = args.index('--bootstrap'); n_boot = int(args[idx+1])
    args.pop(idx+1); args.pop(idx)
if args:
    print(f"Unknown arguments: {args}"); usage()

mode_label = "SESCC-abs"

# ── Cache ──────────────────────────────────────────────────────────────────────
_cache = {}
_cache_dirty = False

def load_cache():
    global _cache
    if os.path.exists(CACHE_FILE):
        print(f"Loading cache from {CACHE_FILE} ...")
        with lzma.open(CACHE_FILE, 'rb') as f:
            _cache = pickle.load(f)
        print(f"Cache loaded: {len(_cache):,} entries")
        return True
    return False

def save_cache():
    if not _cache_dirty:
        print(f"Cache unchanged — no write needed")
        return
    tmp = CACHE_FILE + '.tmp'
    print(f"Saving cache ({len(_cache):,} entries) ...")
    with lzma.open(tmp, 'wb', preset=6) as f:
        pickle.dump(_cache, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, CACHE_FILE)
    print(f"Cache saved ({os.path.getsize(CACHE_FILE)/1024/1024:.1f} MB)")

if not load_cache():
    print("Cache not found:", CACHE_FILE)
    print(f"  D = download from {CACHE_URL}")
    print(f"  N = compute locally (WARNING: several hours for full resolution)")
    try:
        with open('/dev/tty', 'r') as tty:
            while True:
                tty.write("Download cache? [Y/N]: ")
                tty.flush()
                choice = tty.readline().strip().upper()
                if choice in ('Y', 'N'):
                    break
                tty.write("Please enter Y or N\n")
                tty.flush()
    except Exception:
        # Fallback if /dev/tty not available (e.g. Windows)
        choice = 'N'
        print("Could not open terminal for input. Computing locally.")
    if choice == 'Y':
        print(f"Downloading from {CACHE_URL} ...")
        try:
            urllib.request.urlretrieve(CACHE_URL, CACHE_FILE,
                reporthook=lambda b, bs, tot:
                    print(f"  {min(b*bs,tot)/1024/1024:.1f}/{tot/1024/1024:.1f} MB", end='\r'))
            print()
            load_cache()
        except Exception as e:
            print(f"Download failed: {e}. Computing locally.")
    else:
        print("Computing positions on the fly (this may take several hours).")

# ── Astronomy setup ────────────────────────────────────────────────────────────
planets = load('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp')
earth   = planets['earth']
ts      = load.timescale()

with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

# ── Position helpers ───────────────────────────────────────────────────────────
def eclat(hip, t_idx):
    """Ecliptic latitude in milli-degrees at epoch index t_idx.
    Uses ecliptic_latlon(tt) — frame of the epoch date.
    Precession does not affect ecliptic latitudes, so this is correct."""
    key = f"{int(hip)}_{t_idx}_lat"
    if key in _cache:
        return _cache[key]
    global _cache_dirty
    S   = Star.from_dataframe(df.loc[int(hip)])
    tt  = ts.utc(fechamax - t_idx * resolucion, 1, 1)
    app = earth.at(tt).observe(S).apparent()
    lat, _, _ = app.ecliptic_latlon(tt)
    val = round(lat.degrees * 1000)
    _cache[key] = val
    _cache_dirty = True
    return val

def pmotion_lat(hip):
    """Absolute latitudinal proper motion in milli-deg/millennium."""
    S  = Star.from_dataframe(df.loc[int(hip)])
    t2 = ts.utc(1000, 1, 1)
    t1 = ts.utc(0,    1, 1)
    app2 = earth.at(t2).observe(S).apparent()
    lat2, _, _ = app2.ecliptic_latlon()
    app1 = earth.at(t1).observe(S).apparent()
    lat1, _, _ = app1.ecliptic_latlon()
    return abs(round((lat2.degrees - lat1.degrees) * 1000))

# ── Load magnitudes ────────────────────────────────────────────────────────────
hip_mag = {}
try:
    with open("./hip_main.dat") as f:
        for line in f:
            try: hip_mag[int(line[8:14])] = float(line[41:46])
            except: continue
except FileNotFoundError:
    pass

# ── Load catalog from stdin ────────────────────────────────────────────────────
print()
catalog   = []
all_lines = list(sys.stdin)
total     = sum(1 for l in all_lines
                if l.strip() and not l.startswith('#')
                and l.strip().split(';')[0].strip().isdigit()
                and int(l.strip().split(';')[0]) != 0)

for line in all_lines:
    stripped = line.strip()
    if not stripped or stripped.startswith('#'):
        continue
    parts = stripped.split(';')
    try:
        hip = int(parts[0])
    except ValueError:
        continue
    if hip == 0:
        continue
    if maxmag < 100 and hip_mag.get(hip, 0.0) > maxmag:
        continue
    try:
        lat_cat = round(float(parts[1]) * 1000)
    except (ValueError, IndexError):
        continue
    try:
        vel = pmotion_lat(hip)
    except Exception:
        continue

    entry = [hip, vel, lat_cat]
    pct = int(100 * len(catalog) / max(1, total))
    print(f"Computing positions ({pct}%)", end='\r')
    for t in range(maxt):
        entry.append(eclat(hip, t))
    catalog.append(entry)

print("Computing positions (done)    ")

if len(catalog) == 0:
    print("ERROR: no stars loaded."); sys.exit(1)

if filtro > 0:
    random.shuffle(catalog)
    catalog = catalog[:filtro]

N = len(catalog)
print(f"Catalog: {N} stars  |  Mode: {mode_label}")
print(f"Scan: {fechamax} to {fechamax - siglos*100} in steps of {resolucion} year(s)")
print()

# ── Build arrays ───────────────────────────────────────────────────────────────
cat_np      = np.array(catalog, dtype='int64')
velocidades = cat_np[:, 1]
posiciones  = cat_np[:, 2]
epoch_lats  = cat_np[:, 3:]

# ── Diagnostic ─────────────────────────────────────────────────────────────────
print(f"Diagnostic:")
print(f"  N stars          : {N}")
print(f"  Velocity max     : {velocidades.max()} milli-deg/millennium")
print(f"  Velocity median  : {np.median(velocidades):.1f}")
print(f"  Stars with vel=0 : {(velocidades==0).sum()}")
ep0_res = np.abs(posiciones - epoch_lats[:, 0])
print(f"  Mean residual at epoch 0: {ep0_res.mean()/60:.1f} arcmin")
print()

# ── Dot product ────────────────────────────────────────────────────────────────
pos_col   = posiciones[:, np.newaxis]
year_corr = []

for t in range(maxt - 1, -1, -1):
    res = np.abs(pos_col[:, 0] - epoch_lats[:, t])
    scc = np.dot(velocidades, res)
    year_corr.append([fechamax - t * resolucion, scc])

# ── Normalize ──────────────────────────────────────────────────────────────────
year_corr  = np.array(year_corr, dtype='float64')
cmin, cmax = year_corr[:, 1].min(), year_corr[:, 1].max()
if cmax > cmin:
    year_corr[:, 1] = (year_corr[:, 1] - cmin) / (cmax - cmin) * 1000
else:
    print("WARNING: flat curve.")

year_corr = year_corr[year_corr[:, 0].argsort()]

# ── Results ────────────────────────────────────────────────────────────────────
idx_min  = np.argmin(year_corr[:, 1])
year_min = int(year_corr[idx_min, 0])
print(f"Minimum at: {year_min}")
print()

sorted_yc = year_corr[year_corr[:, 1].argsort()]
print("5 lowest values:")
for row in sorted_yc[:5]:
    print(f"  year {int(row[0]):>6}: {row[1]:.1f}")
print()

below100 = year_corr[year_corr[:, 1] < 100]
if len(below100) > 0:
    y_lo = int(below100[:, 0].min())
    y_hi = int(below100[:, 0].max())
    print(f"Minimum region (C<100): {y_lo} to {y_hi} ({y_hi-y_lo} years wide)")
else:
    print("Minimum region: < 1 step wide")

# Top 5 contributing stars
t_idx = max(0, min(maxt-1, (fechamax - year_min) // resolucion))
res_min = np.abs(posiciones - epoch_lats[:, t_idx])
contribs = velocidades * res_min
top5 = np.argsort(contribs)[-5:][::-1]
print()
print(f"Top 5 contributing stars at minimum (year {year_min}):")
for idx in top5:
    hip = int(catalog[idx][0])
    print(f"  HIP{hip:>7}: vel={velocidades[idx]:>6}  "
          f"res={res_min[idx]/60:>7.1f}'  "
          f"contrib={contribs[idx]:>10}")
print()

# ── Bootstrap resampling ──────────────────────────────────────────────────────
boot_minima = []
if n_boot > 0:
    print(f"Bootstrap resampling ({n_boot} resamples) ...")
    rng = np.random.default_rng(42)
    for b in range(n_boot):
        idx_boot = rng.integers(0, N, size=N)
        v_b = velocidades[idx_boot]
        p_b = posiciones[idx_boot]
        e_b = epoch_lats[idx_boot, :]
        p_col = p_b[:, np.newaxis]
        scc_b = np.array([np.dot(v_b, np.abs(p_col[:, 0] - e_b[:, t]))
                          for t in range(maxt - 1, -1, -1)], dtype='float64')
        years_b = np.array([fechamax - t * resolucion
                            for t in range(maxt - 1, -1, -1)], dtype='float64')
        boot_minima.append(int(years_b[np.argmin(scc_b)]))
        if (b + 1) % 100 == 0 or b == n_boot - 1:
            print(f"  {b+1}/{n_boot}", end='\r')
    print()
    arr = np.array(boot_minima)
    pct16 = int(np.percentile(arr, 16))
    pct84 = int(np.percentile(arr, 84))
    pct05 = int(np.percentile(arr,  5))
    pct95 = int(np.percentile(arr, 95))
    print(f"Bootstrap results ({n_boot} resamples, seed=42):")
    print(f"  Median          : {int(np.median(arr))}")
    print(f"  68% range       : [{pct16}, {pct84}]  (width {pct84-pct16} yr)")
    print(f"  90% range       : [{pct05}, {pct95}]  (width {pct95-pct05} yr)")
    print(f"  Std deviation   : {arr.std():.0f} yr")
    print(f"  % pre-CE        : {100*(arr < 0).mean():.1f}%")
    print()

    # Save bootstrap distribution
    boot_file = "sescc_lat_abs_bootstrap.csv"
    np.savetxt(boot_file, arr, fmt='%d')
    print(f"Bootstrap distribution saved to {boot_file}")
    print()

# ── Save cache ─────────────────────────────────────────────────────────────────
save_cache()

# ── Save results ───────────────────────────────────────────────────────────────
outfile = "sescc_lat_abs.csv"
with open(outfile, 'w') as f:
    for row in year_corr:
        f.write(f"{int(row[0])};{row[1]:.2f}\n")
print(f"Results saved to {outfile}")

# ── Plot ───────────────────────────────────────────────────────────────────────
if do_plot:
    if boot_minima:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    else:
        fig, ax1 = plt.subplots(figsize=(14, 6))

    x = year_corr[:, 0]
    y = year_corr[:, 1]
    ax1.plot(x, y, color='blue', linewidth=1.5)
    ax1.axvline(year_min, color='red', linestyle='--', alpha=0.8,
               label=f'Minimum: {year_min}')
    ax1.axvline(-127, color='gray', linestyle=':', alpha=0.6,
               label='Hipparchus (−127)')
    ax1.axvline(137, color='orange', linestyle=':', alpha=0.6,
               label='Ptolemy (+137)')
    ax1.set_ylim(0, 1050)
    ax1.set_xlabel('Year')
    ax1.set_ylabel('Normalised correlation [0-1000]')
    ax1.set_title(f"SESCC dating by LATITUDES — {mode_label}  (N={N})\n"
                 f"Minimum: {year_min}")
    step = max(100, (int(x.max()) - int(x.min())) // 20)
    ax1.set_xticks(np.arange(int(x.min()), int(x.max())+1, step))
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3)
    ax1.legend()

    if boot_minima:
        arr = np.array(boot_minima)
        ax2.hist(arr, bins=40, color='steelblue', edgecolor='white', alpha=0.8)
        ax2.axvline(year_min, color='red', linestyle='--', linewidth=2,
                   label=f'Full catalogue ({year_min})')
        ax2.axvline(int(np.median(arr)), color='orange', linestyle='-',
                   linewidth=1.5, label=f'Median ({int(np.median(arr))})')
        ax2.axvline(-127, color='gray', linestyle=':', alpha=0.6,
                   label='Hipparchus (−127)')
        ax2.axvline(0, color='black', linestyle=':', alpha=0.3)
        ax2.set_xlabel('Year of minimum')
        ax2.set_ylabel('Count')
        ax2.set_title(f"Bootstrap distribution ({n_boot} resamples)\n"
                     f"Median={int(np.median(arr))}  "
                     f"68%=[{int(np.percentile(arr,16))}, {int(np.percentile(arr,84))}]  "
                     f"% pre-CE={100*(arr<0).mean():.0f}%")
        ax2.legend()
        ax2.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.show()
