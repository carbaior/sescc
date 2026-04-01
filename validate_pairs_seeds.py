#!/usr/bin/env python3
# validate_pairs_seeds.py
#
# Runs SESCC-pairs on multiple synthetic catalogues generated with different
# random seeds to characterise the systematic bias and precision of the method.
#
# Generates N_SEEDS synthetic catalogues at EPOCH_YEAR with realistic noise,
# runs SESCC-pairs on each, and reports the distribution of minima.
#
# Usage:
#   python3 validate_pairs_seeds.py --epoch -127 --seeds 20 [OPTIONS]
#   python3 validate_pairs_seeds.py --epoch 137  --seeds 20 [OPTIONS]
#
# Options:
#   --epoch YEAR    epoch of synthetic catalogue (default -127)
#   --seeds N       number of random seeds to test (default 20)
#   --noise-lat L   latitude noise in arcminutes (default 23)
#   --noise-lon L   longitude noise in arcminutes (default 27)
#   --dlat D        SESCC-pairs max latitude difference (default 5.0)
#   --dlon D        SESCC-pairs max longitude difference (default 30.0)
#   --input CSV     input star list (default catalogs/almagest_verbunt.csv)
#   --exclude HIPS  comma-separated HIP numbers to exclude (default 19849,71681)
#   --noplot        suppress plot
#   -h/--help       this help

import sys, os, math, random, pickle, lzma, subprocess, tempfile
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import Star, load
from skyfield.data import hipparcos

# ── Defaults ───────────────────────────────────────────────────────────────────
EPOCH_YEAR  = -127
N_SEEDS     = 20
NOISE_LAT   = 23.0   # arcminutes
NOISE_LON   = 27.0   # arcminutes
MAX_DLAT    = 5.0
MAX_DLON    = 30.0
INPUT_CSV   = "catalogs/almagest_verbunt.csv"
EXCLUDE     = {19849, 71681}
do_plot     = True

CACHE_FILE  = "./sescc_positions.pkl.xz"

# ── Parse arguments ────────────────────────────────────────────────────────────
def usage():
    print(__doc__)
    sys.exit(0)

args = sys.argv[1:]
if '-h' in args or '--help' in args: usage()
if '--noplot'  in args: do_plot = False; args.remove('--noplot')

def getarg(flag, default, cast=str):
    if flag in args:
        idx = args.index(flag)
        val = cast(args[idx+1])
        args.pop(idx+1); args.pop(idx)
        return val
    return default

EPOCH_YEAR = getarg('--epoch',    EPOCH_YEAR, int)
N_SEEDS    = getarg('--seeds',    N_SEEDS,    int)
NOISE_LAT  = getarg('--noise-lat',NOISE_LAT,  float)
NOISE_LON  = getarg('--noise-lon',NOISE_LON,  float)
MAX_DLAT   = getarg('--dlat',     MAX_DLAT,   float)
MAX_DLON   = getarg('--dlon',     MAX_DLON,   float)
INPUT_CSV  = getarg('--input',    INPUT_CSV,  str)
excl_str   = getarg('--exclude',  '19849,71681', str)
EXCLUDE    = set(int(x) for x in excl_str.split(',') if x.strip())

print(f"Validation: epoch={EPOCH_YEAR}  seeds={N_SEEDS}  "
      f"noise_lat={NOISE_LAT}'  noise_lon={NOISE_LON}'")
print(f"SESCC-pairs: |Δβ|<{MAX_DLAT}°  |Δλ|<{MAX_DLON}°")
print(f"Excluded HIPs: {sorted(EXCLUDE)}")
print()

# ── Load cache ─────────────────────────────────────────────────────────────────
_cache = {}
if os.path.exists(CACHE_FILE):
    print(f"Loading cache from {CACHE_FILE} ...")
    with lzma.open(CACHE_FILE, 'rb') as f:
        _cache = pickle.load(f)
    print(f"Cache loaded: {len(_cache)} entries")
else:
    print("WARNING: cache not found. Will compute positions (slow).")

# ── Astronomy setup ────────────────────────────────────────────────────────────
planets = load('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp')
earth   = planets['earth']
ts      = load.timescale()

with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

siglos     = 25
resolucion = 1
maxt       = siglos * 100 // resolucion
fechamax   = 1900

def ecpos_both(hip):
    """Return (lat_mdeg, lon_mdeg) at EPOCH_YEAR for a star."""
    key_lat = f"{int(hip)}_0_lat"  # epoch index 0 = fechamax = 1900 CE
    # We need position at EPOCH_YEAR, not at 1900
    # Compute directly
    S   = Star.from_dataframe(df.loc[int(hip)])
    tt  = ts.utc(EPOCH_YEAR, 7, 1)
    app = earth.at(tt).observe(S).apparent()
    lat, lon, _ = app.ecliptic_latlon(tt)
    return lat.degrees, lon.degrees

def eclon_at_t(hip, t_idx):
    key = f"{int(hip)}_{t_idx}_lon"
    if key in _cache:
        return _cache[key] / 1000.0
    global _cache_dirty
    S   = Star.from_dataframe(df.loc[int(hip)])
    tt  = ts.utc(fechamax - t_idx * resolucion, 1, 1)
    app = earth.at(tt).observe(S).apparent()
    _, lon, _ = app.ecliptic_latlon()   # J2000 inertial frame — no precession
    val = lon.degrees
    _cache[f"{int(hip)}_{t_idx}_lon"] = round(val * 1000)
    _cache_dirty = True
    return val

def pmotion_lon(hip):
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
    return abs(dlon * 1000 * cos_beta)  # milli-deg/millennium

# ── Load star list ─────────────────────────────────────────────────────────────
print(f"Loading star list from {INPUT_CSV} ...")
hips_all = []
with open(INPUT_CSV) as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split(';')
        try:
            hip = int(parts[0])
        except ValueError:
            continue
        if hip == 0 or hip in EXCLUDE:
            continue
        hips_all.append(hip)

print(f"Stars available: {len(hips_all)}")
print()

# ── Precompute true positions and proper motions ───────────────────────────────
print("Precomputing true positions and proper motions ...")
star_data = {}   # hip -> (lat_true, lon_true, vel, epoch_lons_array)

for i, hip in enumerate(hips_all):
    print(f"  {i+1}/{len(hips_all)}", end='\r')
    try:
        lat_true, lon_true = ecpos_both(hip)
        vel = pmotion_lon(hip)
        epoch_lons = np.array([eclon_at_t(hip, t) for t in range(maxt)],
                               dtype='float64')
        star_data[hip] = (lat_true, lon_true, vel, epoch_lons)
    except Exception:
        continue

valid_hips = list(star_data.keys())
print(f"\nStars with data: {len(valid_hips)}")
print()

# ── SESCC-pairs function ───────────────────────────────────────────────────────
def run_sescc_pairs(cat_lons_deg):
    """
    Run SESCC-pairs on a catalogue with given longitudes (degrees).
    Returns the year of the minimum.
    cat_lons_deg: dict hip -> observed_longitude_deg
    """
    N = len(valid_hips)

    # Build arrays
    lats     = np.array([star_data[h][0] for h in valid_hips], dtype='float64')
    vel_arr  = np.array([star_data[h][2] for h in valid_hips], dtype='float64')
    lon_cat  = np.array([cat_lons_deg.get(h, star_data[h][1])
                         for h in valid_hips], dtype='float64') * 1000  # milli-deg
    epoch_lo = np.array([star_data[h][3] for h in valid_hips],
                         dtype='float64') * 1000  # milli-deg, shape (N, maxt)

    # Build pairs
    pair_i, pair_j, pair_vel, pair_dlon_cat = [], [], [], []
    for i in range(N):
        for j in range(i+1, N):
            if abs(lats[i] - lats[j]) > MAX_DLAT:
                continue
            rv = abs(vel_arr[i] - vel_arr[j])
            dlon = lon_cat[i] - lon_cat[j]
            dlon = dlon % 360000
            if dlon > 180000: dlon -= 360000
            if abs(dlon) > MAX_DLON * 1000:
                continue
            pair_i.append(i); pair_j.append(j)
            pair_vel.append(rv); pair_dlon_cat.append(dlon)

    if not pair_i:
        return None

    pi   = np.array(pair_i,        dtype='int32')
    pj   = np.array(pair_j,        dtype='int32')
    pv   = np.array(pair_vel,      dtype='float64')
    pdlc = np.array(pair_dlon_cat, dtype='float64')

    year_corr = []
    for t in range(maxt - 1, -1, -1):
        lon_mod_t = epoch_lo[:, t]
        dlon_mod  = lon_mod_t[pi] - lon_mod_t[pj]
        dlon_mod  = dlon_mod % 360000
        dlon_mod  = np.where(dlon_mod > 180000, dlon_mod - 360000, dlon_mod)
        residuals = np.abs(pdlc - dlon_mod)
        scc = np.dot(pv, residuals)
        year_corr.append([fechamax - t * resolucion, scc])

    yc = np.array(year_corr, dtype='float64')
    cmin, cmax = yc[:, 1].min(), yc[:, 1].max()
    if cmax <= cmin:
        return None
    yc[:, 1] = (yc[:, 1] - cmin) / (cmax - cmin) * 1000
    yc = yc[yc[:, 0].argsort()]
    return int(yc[np.argmin(yc[:, 1]), 0])

# ── Run multiple seeds ─────────────────────────────────────────────────────────
noise_lat_deg = NOISE_LAT / 60.0
noise_lon_deg = NOISE_LON / 60.0
minima = []

print(f"Running {N_SEEDS} seeds ...")
print(f"{'Seed':>6}  {'Minimum':>8}  {'Bias':>8}")
print('-' * 30)

for seed in range(N_SEEDS):
    np.random.seed(seed)
    random.seed(seed)

    # Generate noisy catalogue
    cat_lons = {}
    for hip in valid_hips:
        lat_true, lon_true, vel, _ = star_data[hip]
        lat_noise = np.random.normal(0, noise_lat_deg)
        lon_noise = np.random.normal(0, noise_lon_deg)
        lon_obs   = (lon_true + lon_noise) % 360.0
        # Skip if NaN
        import math as _math
        if _math.isnan(lon_obs) or _math.isnan(lat_true):
            continue
        # Round to nearest 1/6 degree
        lon_obs   = round(lon_obs * 6) / 6.0
        cat_lons[hip] = lon_obs

    year_min = run_sescc_pairs(cat_lons)
    bias     = (year_min - EPOCH_YEAR) if year_min is not None else None

    if year_min is not None:
        minima.append(year_min)
        print(f"{seed:>6}  {year_min:>8}  {bias:>+8}")
    else:
        print(f"{seed:>6}  {'no min':>8}")

print()
if minima:
    arr = np.array(minima)
    print(f"Results over {len(minima)} seeds:")
    print(f"  True epoch       : {EPOCH_YEAR}")
    print(f"  Mean minimum     : {arr.mean():.0f}")
    print(f"  Median minimum   : {np.median(arr):.0f}")
    print(f"  Std deviation    : {arr.std():.0f}")
    print(f"  Mean bias        : {arr.mean() - EPOCH_YEAR:+.0f} years")
    print(f"  Min / Max        : {arr.min()} / {arr.max()}")
    print(f"  % pre-CE         : {100*(arr < 0).mean():.0f}%")

# ── Plot ───────────────────────────────────────────────────────────────────────
if do_plot and minima:
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.hist(minima, bins=20, color='steelblue', edgecolor='white', alpha=0.8)
    ax.axvline(EPOCH_YEAR, color='red', linestyle='--', linewidth=2,
               label=f'True epoch ({EPOCH_YEAR})')
    ax.axvline(np.mean(minima), color='orange', linestyle='-', linewidth=1.5,
               label=f'Mean ({np.mean(minima):.0f})')
    ax.axvline(0, color='gray', linestyle=':', alpha=0.5, label='0 CE')
    ax.set_xlabel('Year of minimum')
    ax.set_ylabel('Count')
    ax.set_title(f"SESCC-pairs validation: epoch={EPOCH_YEAR}, N={N_SEEDS} seeds\n"
                 f"noise={NOISE_LAT}'/lat {NOISE_LON}'/lon  |Δβ|<{MAX_DLAT}°  |Δλ|<{MAX_DLON}°\n"
                 f"Mean={np.mean(minima):.0f}  Std={np.std(minima):.0f}  "
                 f"Bias={np.mean(minima)-EPOCH_YEAR:+.0f}")
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()
# ── Save cache ────────────────────────────────────────────────────────────────
if _cache_dirty:
    tmp = CACHE_FILE + '.tmp'
    print(f"\nSaving cache ({len(_cache):,} entries) to {CACHE_FILE} ...")
    with lzma.open(tmp, 'wb', preset=6) as f:
        pickle.dump(_cache, f, protocol=pickle.HIGHEST_PROTOCOL)
    os.replace(tmp, CACHE_FILE)   # atomic on POSIX — no corruption on interrupt
    size_mb = os.path.getsize(CACHE_FILE) / 1024 / 1024
    print(f"Cache saved ({size_mb:.1f} MB)")
else:
    print(f"\nCache unchanged ({len(_cache):,} entries) — no write needed")
