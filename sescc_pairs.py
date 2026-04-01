#!/usr/bin/env python3
"""
sescc_pairs.py - SESCC-pairs: dates ancient star catalogues using
pairwise longitude differences between neighbouring stars.

For each pair (i,j) where |beta_i - beta_j| < MAX_DLAT and |lon_i - lon_j| < MAX_DLON:
  signal = |mu_lon_i - mu_lon_j| * |(lon_i_cat - lon_j_cat) - (lon_i_mod(T) - lon_j_mod(T))|

Because precession adds the same offset to all longitudes, it cancels exactly
in the pairwise difference. No reference star is needed.

Usage:
  cat catalog.csv | python3 sescc_pairs.py [OPTIONS]

Options:
  --dlat D    max ecliptic latitude difference for pairing (default 5.0 deg)
  --dlon D    max longitude difference between paired stars (default 30.0 deg)
              antipodal pairs (dlon ~ 180 deg) carry no signal and are excluded
  --minvel V    min relative proper motion for a pair to contribute (default 0)
  --bootstrap N run N bootstrap resamples of the star sample (default 0)
  --noplot      suppress plot
  --top N       show top N contributing pairs at minimum (default 10)
  -h/--help     this help

Cache: sescc_positions.pkl.xz (unified LZMA pickle)
Keys: {hip}_{t}_lon (J2000 inertial frame — no precession)
"""

import sys, math, os, pickle, lzma, urllib.request
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

MAX_DLAT_DEFAULT = 5.0    # degrees
MAX_DLON_DEFAULT = 30.0   # degrees — max longitude difference between paired stars
TOP_N_DEFAULT    = 10

# ── Parse arguments ────────────────────────────────────────────────────────────
MAX_DLAT = MAX_DLAT_DEFAULT
MAX_DLON = MAX_DLON_DEFAULT
MIN_VEL  = 0.0
TOP_N    = TOP_N_DEFAULT
do_plot  = True
n_boot   = 0

def usage():
    print(__doc__)
    sys.exit(0)

args = sys.argv[1:]
if '-h' in args or '--help' in args: usage()
if '--noplot' in args: do_plot = False; args.remove('--noplot')
for flag, attr in [('--dlat', 'MAX_DLAT'), ('--dlon', 'MAX_DLON'), ('--minvel', 'MIN_VEL'), ('--top', 'TOP_N')]:
    if flag in args:
        idx = args.index(flag)
        val = float(args[idx+1]) if attr != 'TOP_N' else int(args[idx+1])
        exec(f'{attr} = {val}')
        args.pop(idx+1); args.pop(idx)
if '--bootstrap' in args:
    idx = args.index('--bootstrap'); n_boot = int(args[idx+1])
    args.pop(idx+1); args.pop(idx)

# ── Cache ──────────────────────────────────────────────────────────────────────
_cache = {}
_cache_dirty = False

def load_cache():
    global _cache
    if os.path.exists(CACHE_FILE):
        print(f"Loading cache from {CACHE_FILE} ...")
        with lzma.open(CACHE_FILE, 'rb') as f:
            _cache = pickle.load(f)
        print(f"Cache loaded: {len(_cache)} entries")
        return True
    return False

def save_cache():
    if not _cache_dirty:
        print("Cache unchanged — no write needed")
        return
    tmp = CACHE_FILE + '.tmp'
    print(f"Saving cache ({len(_cache)} entries) ...")
    with lzma.open(tmp, 'wb', preset=6) as f:
        pickle.dump(_cache, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, CACHE_FILE)
    print(f"Cache saved ({os.path.getsize(CACHE_FILE)/1024/1024:.1f} MB)")

if not load_cache():
    print("Cache not found:", CACHE_FILE)
    print(f"  Y = download from {CACHE_URL}")
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

def eclon(hip, t_idx):
    key = f"{int(hip)}_{t_idx}_lon"
    if key in _cache:
        return _cache[key]
    global _cache_dirty
    S   = Star.from_dataframe(df.loc[int(hip)])
    tt  = ts.utc(fechamax - t_idx * resolucion, 1, 1)
    app = earth.at(tt).observe(S).apparent()
    _, lon, _ = app.ecliptic_latlon()   # J2000 inertial frame — no precession
    val = round(lon.degrees * 1000)
    _cache[key] = val
    _cache_dirty = True
    return val

def eclat(hip, t_idx):
    """Ecliptic latitude in milli-degrees at epoch index t_idx.
    Uses frame of the epoch date — correct for latitudes (no precession effect)."""
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
    """Returns |µβ| in milli-deg/millennium."""
    S  = Star.from_dataframe(df.loc[int(hip)])
    t2 = ts.utc(1000, 1, 1)
    t1 = ts.utc(0,    1, 1)
    app2 = earth.at(t2).observe(S).apparent()
    lat2, _, _ = app2.ecliptic_latlon()
    app1 = earth.at(t1).observe(S).apparent()
    lat1, _, _ = app1.ecliptic_latlon()
    return abs(round((lat2.degrees - lat1.degrees) * 1000))

def pmotion_lon(hip):
    """Returns |µλ·cosβ| in milli-deg/millennium."""
    S  = Star.from_dataframe(df.loc[int(hip)])
    t2 = ts.utc(1000, 1, 1)
    t1 = ts.utc(0,    1, 1)
    app2 = earth.at(t2).observe(S).apparent()
    lat2, lon2, _ = app2.ecliptic_latlon()
    app1 = earth.at(t1).observe(S).apparent()
    lat1, lon1, _ = app1.ecliptic_latlon()
    beta_mean = (lat1.degrees + lat2.degrees) / 2.0
    cos_beta  = math.cos(math.radians(beta_mean))
    dlon = lon2.degrees - lon1.degrees
    if dlon >  180: dlon -= 360
    if dlon < -180: dlon += 360
    return abs(round(dlon * 1000 * cos_beta))

# ── Load catalog ───────────────────────────────────────────────────────────────
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
    try:
        lat_cat = float(parts[1])   # ecliptic latitude in degrees
        lon_cat = round(float(parts[2]) * 1000)   # longitude in milli-degrees
    except (ValueError, IndexError):
        continue
    try:
        vel = pmotion_lon(hip)
    except Exception:
        continue

    entry = [hip, vel, lat_cat, lon_cat]
    pct = int(100 * len(catalog) / max(1, total))
    print(f"Loading ({pct}%)", end='\r')
    for t in range(maxt):
        entry.append(eclon(hip, t))
    catalog.append(entry)

print(f"Loaded {len(catalog)} stars        ")

N = len(catalog)
if N == 0:
    print("ERROR: no stars loaded. Check catalog format and path.")
    import sys; sys.exit(1)
print(f"Catalog: {N} stars")
print(f"Method: SESCC-pairs (|Δβ|<{MAX_DLAT}°, |Δλ|<{MAX_DLON}°)")
print()

# ── Build arrays ───────────────────────────────────────────────────────────────
cat_np      = np.array(catalog, dtype='float64')
hips        = [int(e[0]) for e in catalog]
velocidades = cat_np[:, 1]   # |µλ·cosβ| for each star
lat_cat_arr = cat_np[:, 2]   # ecliptic latitude in degrees
lon_cat_arr = cat_np[:, 3]   # catalogue longitude in milli-degrees
epoch_lons  = cat_np[:, 4:]  # modern longitudes at each epoch (milli-degrees)

# ── Build pairs ────────────────────────────────────────────────────────────────
print(f"Building pairs with |Δβ|<{MAX_DLAT}°, |Δλ|<{MAX_DLON}° ...")
pairs = []   # each: (i, j, rel_vel, dcoord_cat)

for i in range(N):
    for j in range(i+1, N):
        dlat = abs(lat_cat_arr[i] - lat_cat_arr[j])
        if dlat > MAX_DLAT:
            continue
        # Longitude proximity filter (sky neighbours)
        dlon_check = lon_cat_arr[i] - lon_cat_arr[j]
        dlon_check = dlon_check % 360000
        if dlon_check > 180000: dlon_check -= 360000
        if abs(dlon_check) > MAX_DLON * 1000:
            continue
        rel_vel = abs(velocidades[i] - velocidades[j])
        if rel_vel < MIN_VEL:
            continue
        pairs.append((i, j, rel_vel, dlon_check))

n_pairs = len(pairs)
print(f"Pairs found: {n_pairs}")
if n_pairs == 0:
    print("ERROR: no pairs found. Try increasing --dlat.")
    sys.exit(1)

# Diagnostic: distribution of relative velocities
pairs_arr = np.array([(p[2], p[3]) for p in pairs], dtype='float64')
rel_vels  = pairs_arr[:, 0]
print(f"Pair relative velocity: max={rel_vels.max():.0f}  "
      f"median={np.median(rel_vels):.1f}  "
      f"pairs with rel_vel>0: {(rel_vels>0).sum()}")
print()

# ── Correlation curve ──────────────────────────────────────────────────────────
# Convert pairs to arrays for vectorised computation
pair_i      = np.array([p[0] for p in pairs], dtype='int32')
pair_j      = np.array([p[1] for p in pairs], dtype='int32')
pair_vel    = np.array([p[2] for p in pairs], dtype='float64')
pair_dlon_cat = np.array([p[3] for p in pairs], dtype='float64')  # dcoord

print("Computing correlation curve ...")
year_corr = []

for t in range(maxt - 1, -1, -1):
    # Modern longitude differences at epoch t
    lon_mod_t  = epoch_lons[:, t]
    dlon_mod   = lon_mod_t[pair_i] - lon_mod_t[pair_j]
    # Wrap to [-180000, 180000]
    dlon_mod   = dlon_mod % 360000
    dlon_mod   = np.where(dlon_mod > 180000, dlon_mod - 360000, dlon_mod)

    # Residual: how much the longitude difference changed from catalogue
    residuals  = np.abs(pair_dlon_cat - dlon_mod)

    # Correlation statistic
    scc = np.dot(pair_vel, residuals)
    year_corr.append([fechamax - t * resolucion, scc])

    pct = int(100 * (maxt - 1 - t) / maxt)
    print(f"  {pct}%", end='\r')

print("  100%   ")

# ── Normalize ──────────────────────────────────────────────────────────────────
year_corr  = np.array(year_corr, dtype='float64')
cmin, cmax = year_corr[:, 1].min(), year_corr[:, 1].max()
if cmax > cmin:
    year_corr[:, 1] = (year_corr[:, 1] - cmin) / (cmax - cmin) * 1000
else:
    print("WARNING: flat curve — no variation detected.")
    sys.exit(1)

year_corr = year_corr[year_corr[:, 0].argsort()]

# ── Results ────────────────────────────────────────────────────────────────────
idx_min  = np.argmin(year_corr[:, 1])
year_min = int(year_corr[idx_min, 0])
print(f"\nMinimum at: {year_min}")

sorted_yc = year_corr[year_corr[:, 1].argsort()]
print("5 lowest values:")
for row in sorted_yc[:5]:
    print(f"  year {int(row[0]):>6}: {row[1]:.1f}")

below100 = year_corr[year_corr[:, 1] < 100]
if len(below100) > 0:
    y_lo = int(below100[:, 0].min())
    y_hi = int(below100[:, 0].max())
    print(f"Minimum region (C<100): {y_lo} to {y_hi} ({y_hi-y_lo} years wide)")
else:
    print("Minimum region: < 1 step wide")

# ── Top contributing pairs at minimum ─────────────────────────────────────────
t_min_idx = (fechamax - year_min) // resolucion
t_min_idx = max(0, min(maxt-1, t_min_idx))
lon_mod_at_min = epoch_lons[:, t_min_idx]
dlon_mod_at_min = lon_mod_at_min[pair_i] - lon_mod_at_min[pair_j]
dlon_mod_at_min = dlon_mod_at_min % 360000
dlon_mod_at_min = np.where(dlon_mod_at_min > 180000,
                            dlon_mod_at_min - 360000, dlon_mod_at_min)
res_at_min   = np.abs(pair_dlon_cat - dlon_mod_at_min)
contrib      = pair_vel * res_at_min
top_idx      = np.argsort(contrib)[-TOP_N:][::-1]

print(f"\nTop {TOP_N} contributing pairs at minimum (year {year_min}):")
print("  " + f"{'HIPi':>8} {'HIPj':>8} {'Dlat_deg':>9} {'rel_vel':>8} {'res_arcmin':>11} {'contrib':>14}")
for idx in top_idx:
    i, j = int(pair_i[idx]), int(pair_j[idx])
    dlat = abs(lat_cat_arr[i] - lat_cat_arr[j])
    print(f"  {hips[i]:>8} {hips[j]:>8} {dlat:>7.2f} "
          f"{pair_vel[idx]:>8.0f} "
          f"{res_at_min[idx]/60:>8.1f}' "
          f"{contrib[idx]:>12.0f}")

# ── Bootstrap resampling ──────────────────────────────────────────────────────
boot_minima = []
if n_boot > 0:
    print(f"\nBootstrap resampling ({n_boot} resamples) ...")
    rng = np.random.default_rng(42)
    # Bootstrap resamples stars (not pairs — pairs are rebuilt each time)
    star_idx_all = np.arange(N)
    for b in range(n_boot):
        idx_boot = rng.integers(0, N, size=N)
        # Rebuild pairs from resampled stars
        lats_b   = lat_cat_arr[idx_boot]
        lons_b   = lon_cat_arr[idx_boot]
        vels_b   = velocidades[idx_boot]
        ep_b     = epoch_lons[idx_boot, :]

        pi_b, pj_b, pv_b, pdlc_b = [], [], [], []
        for i in range(N):
            for j in range(i+1, N):
                if abs(lats_b[i] - lats_b[j]) > MAX_DLAT: continue
                dlon = lons_b[i] - lons_b[j]
                dlon = dlon % 360000
                if dlon > 180000: dlon -= 360000
                if abs(dlon) > MAX_DLON * 1000: continue
                pi_b.append(i); pj_b.append(j)
                pv_b.append(abs(vels_b[i] - vels_b[j]))
                pdlc_b.append(dlon)

        if not pi_b:
            continue
        pi_b  = np.array(pi_b,  dtype='int32')
        pj_b  = np.array(pj_b,  dtype='int32')
        pv_b  = np.array(pv_b,  dtype='float64')
        pdlc_b= np.array(pdlc_b,dtype='float64')

        scc_b = []
        for t in range(maxt - 1, -1, -1):
            dm = ep_b[:, t][pi_b] - ep_b[:, t][pj_b]
            dm = dm % 360000
            dm = np.where(dm > 180000, dm - 360000, dm)
            scc_b.append(np.dot(pv_b, np.abs(pdlc_b - dm)))

        years_b = np.array([fechamax - t * resolucion
                            for t in range(maxt-1, -1, -1)], dtype='float64')
        boot_minima.append(int(years_b[np.argmin(scc_b)]))
        if (b+1) % 10 == 0 or b == n_boot-1:
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
    np.savetxt("sescc_lon_pairs_bootstrap.csv", arr, fmt='%d')
    print(f"Bootstrap distribution saved to sescc_lon_pairs_bootstrap.csv")

# ── Save cache ─────────────────────────────────────────────────────────────────
print()
save_cache()

# ── Save results ───────────────────────────────────────────────────────────────
outfile = "sescc_lon_pairs.csv"
with open(outfile, 'w') as f:
    for row in year_corr:
        f.write(f"{int(row[0])};{row[1]:.2f}\n")
print(f"Results saved to {outfile}")

# ── Plot ───────────────────────────────────────────────────────────────────────
if do_plot:
    x = year_corr[:, 0]
    y = year_corr[:, 1]
    if boot_minima:
        fig, (ax, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    else:
        fig, ax = plt.subplots(figsize=(14, 6))
    ax.plot(x, y, color='green', linewidth=1.5)
    ax.axvline(year_min, color='red', linestyle='--', alpha=0.8,
               label=f'Minimum: {year_min}')
    ax.axvline(-127, color='gray', linestyle=':', alpha=0.6,
               label='Hipparchus (−127)')
    ax.axvline( 137, color='orange', linestyle=':', alpha=0.6,
               label='Ptolemy (+137)')
    ax.set_ylim(0, 1050)
    ax.set_xlabel('Year')
    ax.set_ylabel('Normalised correlation [0-1000]')
    ax.set_title(f"SESCC-pairs by LONGITUDES  |Δβ|<{MAX_DLAT}°  |Δλ|<{MAX_DLON}°  N={N}  pairs={n_pairs}\n"
                 f"Minimum: {year_min}")
    step = max(100, (int(x.max()) - int(x.min())) // 20)
    ax.set_xticks(np.arange(int(x.min()), int(x.max())+1, step))
    ax.tick_params(axis='x', rotation=45)
    ax.grid(True, alpha=0.3)
    ax.legend()

    if boot_minima:
        arr = np.array(boot_minima)
        ax2.hist(arr, bins=40, color='seagreen', edgecolor='white', alpha=0.8)
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
