#!/usr/bin/env python3
"""
sescc.py - SESCC: Speed-Error Signals Cross-Correlation
Dates ancient star catalogues using ecliptic coordinates.

Two modes:
  --lat    SESCC (latitudes):  C(T) = Σᵢ |µβᵢ| · |βᵢᶜᵃᵗ − βᵢᵐᵒᵈ(T)|
  --lon    SESCC-pairs (longitudes): Cₚ(T) = Σ₍ᵢ,ⱼ₎ |µλᵢ−µλⱼ| · |Δλᵢⱼᶜᵃᵗ − Δλᵢⱼᵐᵒᵈ(T)|

For latitudes, precession does not affect ecliptic latitudes — no reference
correction needed.
For longitudes, pairwise differences cancel any global offset including
precession, so the method is immune to systematic longitude shifts.

Usage:
  cat catalog.csv | python3 sescc.py --lat [OPTIONS]
  cat catalog.csv | python3 sescc.py --lon [OPTIONS]

Common options:
  --bootstrap N   run N bootstrap resamples (default 0)
  --noplot        suppress plot
  -h/--help       this help

Latitude options:
  --subset N      use a random subset of N stars
  --maxmag M      exclude stars fainter than magnitude M

Longitude options:
  --dlat D        max latitude difference for pairing (default 5.0 deg)
  --dlon D        max longitude difference for pairing (default 30.0 deg)
  --minvel V      min relative proper motion for a pair (default 0)
  --top N         show top N contributing pairs at minimum (default 10)

Cache: sescc_positions.pkl.xz (LZMA pickle)
  If the cache file is not found, positions are computed on the fly
  (this may take several hours for the full 2500-year scan).

copyleft (GPLv3) 2024 Carlos Baiget Orts (asinfreedom@gmail.com)
https://github.com/carbaior/sescc
"""

import sys, os, math, random, pickle, lzma
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import Star, load
from skyfield.data import hipparcos

CACHE_FILE = "./sescc_positions.pkl.xz"
siglos = 25; resolucion = 1
maxt = siglos * 100 // resolucion; fechamax = 1900
MAX_DLAT_DEFAULT = 5.0; MAX_DLON_DEFAULT = 30.0; TOP_N_DEFAULT = 10

mode = None; do_plot = True; n_boot = 0; filtro = 0; maxmag = 100
MAX_DLAT = MAX_DLAT_DEFAULT; MAX_DLON = MAX_DLON_DEFAULT
MIN_VEL = 0.0; TOP_N = TOP_N_DEFAULT

def usage():
    print(__doc__); sys.exit(0)

args = sys.argv[1:]
if '-h' in args or '--help' in args: usage()
if '--lat' in args: mode = 'lat'; args.remove('--lat')
if '--lon' in args: mode = 'lon'; args.remove('--lon')
if mode is None: print("ERROR: specify --lat or --lon"); usage()
if '--noplot' in args: do_plot = False; args.remove('--noplot')
if '--bootstrap' in args:
    idx = args.index('--bootstrap'); n_boot = int(args[idx+1])
    args.pop(idx+1); args.pop(idx)
if '--subset' in args:
    idx = args.index('--subset'); filtro = int(args[idx+1])
    args.pop(idx+1); args.pop(idx)
if '--maxmag' in args:
    idx = args.index('--maxmag'); maxmag = float(args[idx+1])
    args.pop(idx+1); args.pop(idx)
for flag, var in [('--dlat','MAX_DLAT'),('--dlon','MAX_DLON'),('--minvel','MIN_VEL')]:
    if flag in args:
        idx = args.index(flag); globals()[var] = float(args[idx+1])
        args.pop(idx+1); args.pop(idx)
if '--top' in args:
    idx = args.index('--top'); TOP_N = int(args[idx+1])
    args.pop(idx+1); args.pop(idx)
if args: print(f"Unknown arguments: {args}"); usage()

_cache = {}; _cache_dirty = False

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
        print("Cache unchanged — no write needed"); return
    tmp = CACHE_FILE + '.tmp'
    print(f"Saving cache ({len(_cache):,} entries) ...")
    with lzma.open(tmp, 'wb', preset=6) as f:
        pickle.dump(_cache, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, CACHE_FILE)
    print(f"Cache saved ({os.path.getsize(CACHE_FILE)/1024/1024:.1f} MB)")

if not load_cache():
    print(f"WARNING: cache file not found ({CACHE_FILE})")
    print("Positions will be computed on the fly.")
    print("This may take several hours for a full 2500-year scan.")
    print()

planets = load('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp')
earth = planets['earth']; ts = load.timescale()
with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

def eclat(hip, t_idx):
    key = f"{int(hip)}_{t_idx}_lat"
    if key in _cache: return _cache[key]
    global _cache_dirty
    S = Star.from_dataframe(df.loc[int(hip)])
    tt = ts.utc(fechamax - t_idx * resolucion, 1, 1)
    lat, _, _ = earth.at(tt).observe(S).apparent().ecliptic_latlon(tt)
    val = round(lat.degrees * 1000)
    _cache[key] = val; _cache_dirty = True; return val

def eclon(hip, t_idx):
    key = f"{int(hip)}_{t_idx}_lon"
    if key in _cache: return _cache[key]
    global _cache_dirty
    S = Star.from_dataframe(df.loc[int(hip)])
    tt = ts.utc(fechamax - t_idx * resolucion, 1, 1)
    _, lon, _ = earth.at(tt).observe(S).apparent().ecliptic_latlon()
    val = round(lon.degrees * 1000)
    _cache[key] = val; _cache_dirty = True; return val

def pmotion_lat(hip):
    S = Star.from_dataframe(df.loc[int(hip)])
    l2, _, _ = earth.at(ts.utc(1000,1,1)).observe(S).apparent().ecliptic_latlon()
    l1, _, _ = earth.at(ts.utc(0,1,1)).observe(S).apparent().ecliptic_latlon()
    return abs(round((l2.degrees - l1.degrees) * 1000))

def pmotion_lon(hip):
    S = Star.from_dataframe(df.loc[int(hip)])
    lat2, lon2, _ = earth.at(ts.utc(1000,1,1)).observe(S).apparent().ecliptic_latlon()
    lat1, lon1, _ = earth.at(ts.utc(0,1,1)).observe(S).apparent().ecliptic_latlon()
    beta = (lat1.degrees + lat2.degrees) / 2.0
    dlon = lon2.degrees - lon1.degrees
    if dlon > 180: dlon -= 360
    if dlon < -180: dlon += 360
    return abs(round(dlon * 1000 * math.cos(math.radians(beta))))

hip_mag = {}
try:
    with open("./hip_main.dat") as f:
        for line in f:
            try: hip_mag[int(line[8:14])] = float(line[41:46])
            except: continue
except FileNotFoundError:
    pass

pmotion_fn = pmotion_lat if mode == 'lat' else pmotion_lon
pos_fn     = eclat       if mode == 'lat' else eclon

print()
catalog = []; all_lines = list(sys.stdin)
total = sum(1 for l in all_lines
            if l.strip() and not l.startswith('#')
            and l.strip().split(';')[0].strip().isdigit()
            and int(l.strip().split(';')[0]) != 0)

for line in all_lines:
    stripped = line.strip()
    if not stripped or stripped.startswith('#'): continue
    parts = stripped.split(';')
    try: hip = int(parts[0])
    except ValueError: continue
    if hip == 0: continue
    try:
        lat_cat = float(parts[1])
        lon_cat = round(float(parts[2]) * 1000)
    except (ValueError, IndexError): continue
    if maxmag < 100 and hip in hip_mag and hip_mag[hip] > maxmag: continue
    try: vel = pmotion_fn(hip)
    except Exception: continue
    entry = [hip, vel, lat_cat, lon_cat]
    pct = int(100 * len(catalog) / max(1, total))
    print(f"Computing positions ({pct}%)", end='\r')
    for t in range(maxt):
        entry.append(pos_fn(hip, t))
    catalog.append(entry)

if mode == 'lat' and filtro > 0 and filtro < len(catalog):
    catalog = random.sample(catalog, filtro)

print(f"Computing positions (done)    ")
N = len(catalog)
if N == 0: print("ERROR: no stars loaded."); sys.exit(1)

cat_np      = np.array(catalog, dtype='float64')
hips        = [int(e[0]) for e in catalog]
velocidades = cat_np[:, 1]
lat_cat_arr = cat_np[:, 2]
lon_cat_arr = cat_np[:, 3]
epoch_pos   = cat_np[:, 4:]

def print_bootstrap(arr, n_boot):
    p16,p84 = int(np.percentile(arr,16)), int(np.percentile(arr,84))
    p05,p95 = int(np.percentile(arr,5)),  int(np.percentile(arr,95))
    print(f"Bootstrap results ({n_boot} resamples, seed=42):")
    print(f"  Median          : {int(np.median(arr))}")
    print(f"  68% range       : [{p16}, {p84}]  (width {p84-p16} yr)")
    print(f"  90% range       : [{p05}, {p95}]  (width {p95-p05} yr)")
    print(f"  Std deviation   : {arr.std():.0f} yr")
    print(f"  % pre-CE        : {100*(arr<0).mean():.1f}%")

def plot_curve_and_boot(year_corr, year_min, boot_minima, n_boot,
                        color, title_curve, title_boot_prefix):
    if boot_minima:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    else:
        fig, ax1 = plt.subplots(figsize=(14, 6))
    x, y = year_corr[:, 0], year_corr[:, 1]
    ax1.plot(x, y, color=color, linewidth=1.5)
    ax1.axvline(year_min, color='red', linestyle='--', alpha=0.8,
               label=f'Minimum: {year_min}')
    ax1.axvline(-127, color='gray', linestyle=':', alpha=0.6, label='Hipparchus (−127)')
    ax1.axvline(137,  color='orange', linestyle=':', alpha=0.6, label='Ptolemy (+137)')
    ax1.set_ylim(0, 1050); ax1.set_xlabel('Year')
    ax1.set_ylabel('Normalised correlation [0–1000]')
    ax1.set_title(title_curve)
    step = max(100, (int(x.max()) - int(x.min())) // 20)
    ax1.set_xticks(np.arange(int(x.min()), int(x.max())+1, step))
    ax1.tick_params(axis='x', rotation=45)
    ax1.grid(True, alpha=0.3); ax1.legend()
    if boot_minima:
        arr = np.array(boot_minima)
        ax2.hist(arr, bins=40, color=color, edgecolor='white', alpha=0.8)
        ax2.axvline(year_min, color='red', linestyle='--', linewidth=2,
                   label=f'Full catalogue ({year_min})')
        ax2.axvline(int(np.median(arr)), color='orange', linestyle='-',
                   linewidth=1.5, label=f'Median ({int(np.median(arr))})')
        ax2.axvline(-127, color='gray', linestyle=':', alpha=0.6, label='Hipparchus (−127)')
        ax2.axvline(0, color='black', linestyle=':', alpha=0.3)
        ax2.set_xlabel('Year of minimum'); ax2.set_ylabel('Count')
        ax2.set_title(f"{title_boot_prefix} ({n_boot} resamples)\n"
                     f"Median={int(np.median(arr))}  "
                     f"68%=[{int(np.percentile(arr,16))}, {int(np.percentile(arr,84))}]  "
                     f"% pre-CE={100*(arr<0).mean():.0f}%")
        ax2.legend(); ax2.grid(True, alpha=0.3)
    plt.tight_layout(); plt.show()

# ══════════════════════════════════════════════════════════════════════════════
if mode == 'lat':
# ══════════════════════════════════════════════════════════════════════════════
    print(f"Catalog: {N} stars  |  Mode: SESCC (latitudes)")
    print(f"Scan: {fechamax} to {fechamax-(maxt-1)*resolucion} in steps of {resolucion} yr")
    posiciones = (lat_cat_arr * 1000).astype('int64')
    ep0_res = np.abs(posiciones - epoch_pos[:, 0])
    print(f"\nDiagnostic:")
    print(f"  N stars          : {N}")
    print(f"  Velocity max     : {velocidades.max():.0f} milli-deg/millennium")
    print(f"  Velocity median  : {np.median(velocidades):.1f}")
    print(f"  Stars with vel=0 : {(velocidades==0).sum()}")
    print(f"  Mean residual at epoch 0: {ep0_res.mean()/60:.1f} arcmin\n")

    pos_col = posiciones[:, np.newaxis]
    year_corr = []
    for t in range(maxt - 1, -1, -1):
        scc = np.dot(velocidades, np.abs(pos_col[:, 0] - epoch_pos[:, t]))
        year_corr.append([fechamax - t * resolucion, scc])
    year_corr = np.array(year_corr, dtype='float64')
    cmin, cmax = year_corr[:,1].min(), year_corr[:,1].max()
    year_corr[:,1] = (year_corr[:,1]-cmin)/(cmax-cmin)*1000 if cmax>cmin else (print("WARNING: flat curve") or year_corr[:,1])
    year_corr = year_corr[year_corr[:,0].argsort()]
    year_min = int(year_corr[np.argmin(year_corr[:,1]), 0])
    print(f"Minimum at: {year_min}")
    for row in year_corr[year_corr[:,1].argsort()][:5]:
        print(f"  year {int(row[0]):>6}: {row[1]:.1f}")
    below100 = year_corr[year_corr[:,1] < 100]
    if len(below100):
        y_lo,y_hi = int(below100[:,0].min()), int(below100[:,0].max())
        print(f"Minimum region (C<100): {y_lo} to {y_hi} ({y_hi-y_lo} years wide)")
    t_idx = max(0, min(maxt-1, (fechamax-year_min)//resolucion))
    res_min = np.abs(posiciones - epoch_pos[:, t_idx])
    contribs = velocidades * res_min
    top5 = np.argsort(contribs)[-5:][::-1]
    print(f"\nTop 5 contributing stars at minimum (year {year_min}):")
    for idx in top5:
        print(f"  HIP{hips[idx]:>7}: vel={velocidades[idx]:>6.0f}  "
              f"res={res_min[idx]/60:>7.1f}'  contrib={contribs[idx]:>10.0f}")

    boot_minima = []
    if n_boot > 0:
        print(f"\nBootstrap resampling ({n_boot} resamples) ...")
        rng = np.random.default_rng(42)
        for b in range(n_boot):
            idx_b = rng.integers(0, N, size=N)
            v_b = velocidades[idx_b]; p_b = posiciones[idx_b]; e_b = epoch_pos[idx_b,:]
            p_col = p_b[:, np.newaxis]
            scc_b = np.array([np.dot(v_b, np.abs(p_col[:,0]-e_b[:,t]))
                              for t in range(maxt-1,-1,-1)], dtype='float64')
            years_b = np.array([fechamax-t*resolucion for t in range(maxt-1,-1,-1)], dtype='float64')
            boot_minima.append(int(years_b[np.argmin(scc_b)]))
            if (b+1) % 100 == 0 or b == n_boot-1: print(f"  {b+1}/{n_boot}", end='\r')
        print()
        arr = np.array(boot_minima)
        print_bootstrap(arr, n_boot)
        np.savetxt("sescc_lat_bootstrap.csv", arr, fmt='%d')
        print("Bootstrap distribution saved to sescc_lat_bootstrap.csv")

    save_cache()
    outfile = "sescc_lat.csv"
    with open(outfile, 'w') as f:
        for row in year_corr: f.write(f"{int(row[0])};{row[1]:.2f}\n")
    print(f"\nResults saved to {outfile}")

    if do_plot:
        plot_curve_and_boot(year_corr, year_min, boot_minima, n_boot,
            'steelblue',
            f"SESCC — latitudes  (N={N})\nMinimum: {year_min}",
            "Bootstrap distribution")

# ══════════════════════════════════════════════════════════════════════════════
else:  # mode == 'lon'
# ══════════════════════════════════════════════════════════════════════════════
    print(f"Catalog: {N} stars  |  Mode: SESCC-pairs (longitudes)")
    print(f"Pairing: |Δβ|<{MAX_DLAT}°, |Δλ|<{MAX_DLON}°\n")

    pairs = []
    for i in range(N):
        for j in range(i+1, N):
            if abs(lat_cat_arr[i]-lat_cat_arr[j]) > MAX_DLAT: continue
            dlon = lon_cat_arr[i]-lon_cat_arr[j]
            dlon = dlon % 360000
            if dlon > 180000: dlon -= 360000
            if abs(dlon) > MAX_DLON*1000: continue
            rel_vel = abs(velocidades[i]-velocidades[j])
            if rel_vel < MIN_VEL: continue
            pairs.append((i, j, rel_vel, dlon))

    n_pairs = len(pairs)
    print(f"Pairs found: {n_pairs}")
    if n_pairs == 0: print("ERROR: no pairs found. Try increasing --dlat or --dlon."); sys.exit(1)
    rel_vels = np.array([p[2] for p in pairs])
    print(f"Pair relative velocity: max={rel_vels.max():.0f}  median={np.median(rel_vels):.1f}  "
          f"pairs with rel_vel>0: {(rel_vels>0).sum()}\n")

    pair_i        = np.array([p[0] for p in pairs], dtype='int32')
    pair_j        = np.array([p[1] for p in pairs], dtype='int32')
    pair_vel      = np.array([p[2] for p in pairs], dtype='float64')
    pair_dlon_cat = np.array([p[3] for p in pairs], dtype='float64')

    print("Computing correlation curve ...")
    year_corr = []
    for t in range(maxt-1, -1, -1):
        dm = epoch_pos[:,t][pair_i] - epoch_pos[:,t][pair_j]
        dm = dm % 360000
        dm = np.where(dm > 180000, dm-360000, dm)
        year_corr.append([fechamax-t*resolucion, np.dot(pair_vel, np.abs(pair_dlon_cat-dm))])
        print(f"  {int(100*(maxt-1-t)/maxt)}%", end='\r')
    print("  100%   ")

    year_corr = np.array(year_corr, dtype='float64')
    cmin, cmax = year_corr[:,1].min(), year_corr[:,1].max()
    year_corr[:,1] = (year_corr[:,1]-cmin)/(cmax-cmin)*1000 if cmax>cmin else year_corr[:,1]
    year_corr = year_corr[year_corr[:,0].argsort()]
    year_min = int(year_corr[np.argmin(year_corr[:,1]), 0])
    print(f"\nMinimum at: {year_min}")
    for row in year_corr[year_corr[:,1].argsort()][:5]:
        print(f"  year {int(row[0]):>6}: {row[1]:.1f}")
    below100 = year_corr[year_corr[:,1] < 100]
    if len(below100):
        y_lo,y_hi = int(below100[:,0].min()), int(below100[:,0].max())
        print(f"Minimum region (C<100): {y_lo} to {y_hi} ({y_hi-y_lo} years wide)")

    t_min_idx = max(0, min(maxt-1, (fechamax-year_min)//resolucion))
    dm_min = epoch_pos[:,t_min_idx][pair_i] - epoch_pos[:,t_min_idx][pair_j]
    dm_min = dm_min % 360000
    dm_min = np.where(dm_min > 180000, dm_min-360000, dm_min)
    res_at_min = np.abs(pair_dlon_cat - dm_min)
    contrib    = pair_vel * res_at_min
    top_idx    = np.argsort(contrib)[-TOP_N:][::-1]
    print(f"\nTop {TOP_N} contributing pairs at minimum (year {year_min}):")
    print("  " + f"{'HIPi':>8} {'HIPj':>8} {'Dlat_deg':>9} {'rel_vel':>8} {'res_arcmin':>11} {'contrib':>14}")
    for idx in top_idx:
        i, j = int(pair_i[idx]), int(pair_j[idx])
        print(f"  {hips[i]:>8} {hips[j]:>8} {abs(lat_cat_arr[i]-lat_cat_arr[j]):>7.2f} "
              f"{pair_vel[idx]:>8.0f} {res_at_min[idx]/60:>8.1f}' {contrib[idx]:>12.0f}")

    boot_minima = []
    if n_boot > 0:
        print(f"\nBootstrap resampling ({n_boot} resamples) ...")
        rng = np.random.default_rng(42)
        for b in range(n_boot):
            idx_b = rng.integers(0, N, size=N)
            lats_b=lat_cat_arr[idx_b]; lons_b=lon_cat_arr[idx_b]
            vels_b=velocidades[idx_b]; ep_b=epoch_pos[idx_b,:]
            pi_b,pj_b,pv_b,pdlc_b=[],[],[],[]
            for i in range(N):
                for j in range(i+1,N):
                    if abs(lats_b[i]-lats_b[j])>MAX_DLAT: continue
                    dlon=lons_b[i]-lons_b[j]; dlon=dlon%360000
                    if dlon>180000: dlon-=360000
                    if abs(dlon)>MAX_DLON*1000: continue
                    pi_b.append(i); pj_b.append(j)
                    pv_b.append(abs(vels_b[i]-vels_b[j])); pdlc_b.append(dlon)
            if not pi_b: continue
            pi_b=np.array(pi_b,dtype='int32'); pj_b=np.array(pj_b,dtype='int32')
            pv_b=np.array(pv_b,dtype='float64'); pdlc_b=np.array(pdlc_b,dtype='float64')
            scc_b=[]
            for t in range(maxt-1,-1,-1):
                dm=ep_b[:,t][pi_b]-ep_b[:,t][pj_b]; dm=dm%360000
                dm=np.where(dm>180000,dm-360000,dm)
                scc_b.append(np.dot(pv_b,np.abs(pdlc_b-dm)))
            years_b=np.array([fechamax-t*resolucion for t in range(maxt-1,-1,-1)],dtype='float64')
            boot_minima.append(int(years_b[np.argmin(scc_b)]))
            if (b+1)%10==0 or b==n_boot-1: print(f"  {b+1}/{n_boot}",end='\r')
        print()
        arr = np.array(boot_minima)
        print_bootstrap(arr, n_boot)
        np.savetxt("sescc_lon_bootstrap.csv", arr, fmt='%d')
        print("Bootstrap distribution saved to sescc_lon_bootstrap.csv")

    save_cache()
    outfile = "sescc_lon.csv"
    with open(outfile, 'w') as f:
        for row in year_corr: f.write(f"{int(row[0])};{row[1]:.2f}\n")
    print(f"\nResults saved to {outfile}")

    if do_plot:
        plot_curve_and_boot(year_corr, year_min, boot_minima, n_boot,
            'seagreen',
            f"SESCC-pairs — longitudes  |Δβ|<{MAX_DLAT}°  |Δλ|<{MAX_DLON}°  N={N}  pairs={n_pairs}\nMinimum: {year_min}",
            "Bootstrap distribution")
