# SESCC — Speed-Error Signals Cross-Correlation

A method for dating ancient star catalogues by cross-correlating stellar
proper-motion speeds with positional residuals.

**Author:** Carlos Baiget Orts (asinfreedom@gmail.com)  
**License:** GPLv3  
**Repository:** https://github.com/carbaior/sescc

---

## Principle

A star catalogue records positions at the epoch of observation. When compared
to modern positions computed from the Hipparcos catalogue with proper motion
propagated to a trial epoch T, the residual |pos_cat − pos_mod(T)| is
minimised when T equals the true observational epoch.

At the true epoch, residuals reflect random measurement errors independent of
proper-motion speed. At any other epoch, fast-moving stars accumulate
systematically larger residuals. The SESCC statistic exploits this:

**For ecliptic latitudes (`--lat`):**

    C(T) = Σᵢ |µβᵢ| · |βᵢᶜᵃᵗ − βᵢᵐᵒᵈ(T)|

Precession does not affect ecliptic latitudes, so no correction is needed.

**For ecliptic longitudes (`--lon`, SESCC-pairs):**

    Cₚ(T) = Σ₍ᵢ,ⱼ₎ |µλᵢ·cosβᵢ − µλⱼ·cosβⱼ| · |(λᵢᶜᵃᵗ − λⱼᶜᵃᵗ) − (λᵢᵐᵒᵈ(T) − λⱼᵐᵒᵈ(T))|

Pairwise differences cancel any global longitude offset — including precession
corrections — exactly. SESCC-pairs is immune to systematic longitude shifts
by algebraic construction.

The epoch estimate is the T that minimises C(T) or Cₚ(T).

---

## Quick start

Ready-to-run scripts are provided for the three historical catalogues.
Each script shows the C(T) curve and reports the epoch estimate.

**Prerequisites:** catalogues must be prepared first (see
[Catalogue preparation](#catalogue-preparation)).

### Date the Almagest

```bash
bash date_almagest_lat.sh    # by ecliptic latitudes  → expected ~-49 BCE
bash date_almagest_lon.sh    # by ecliptic longitudes → expected ~-165 BCE
```

### Date Ulugh Beg (true epoch: 1437 CE)

```bash
bash date_ulugh_lat.sh       # by ecliptic latitudes  → expected ~1177 CE
bash date_ulugh_lon.sh       # by ecliptic longitudes → expected ~1452 CE
```

### Date Tycho Brahe (true epoch: ~1580 CE)

```bash
bash date_brahe_lat.sh       # by ecliptic latitudes  → expected ~1570 CE
bash date_brahe_lon.sh       # by ecliptic longitudes → expected ~1547 CE
```

---

## Catalogue preparation

Download the source data files and convert them to SESCC format before
running the dating scripts.

**Almagest and Ulugh Beg** — from Verbunt & van Gent (2012), A&A 544, A31:

```bash
# Download from CDS: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/544/A31
python3 verbunt2sescc.py ptolema.dat   > catalogs/almagest_verbunt.csv
python3 verbunt2sescc.py ulughbeg.dat  > catalogs/ulughbeg_verbunt.csv
```

**Tycho Brahe** — from Verbunt & van Gent (2010), A&A 516, A28:

```bash
# Download KeplerE.dat from CDS: https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/516/A28
python3 verbunt2sescc_brahe.py keplere.dat > catalogs/brahe_verbunt.csv
```

---

## Main script: `sescc.py`

All dating is done through a single script with a `--lat` or `--lon` flag.

    cat catalog.csv | python3 sescc.py --lat [OPTIONS]
    cat catalog.csv | python3 sescc.py --lon [OPTIONS]

**Common options:**

    --bootstrap N   run N bootstrap resamples and report distribution (default 0)
    --noplot        suppress the plot

**Latitude options (`--lat`):**

    --subset N      use a random subset of N stars
    --maxmag M      exclude stars fainter than magnitude M

**Longitude options (`--lon`):**

    --dlat D        max latitude difference for pairing (default 5.0 deg)
    --dlon D        max longitude difference for pairing (default 30.0 deg)
    --minvel V      min relative proper motion for a pair (default 0)
    --top N         show top N contributing pairs at minimum (default 10)

**Output files:**

    sescc_lat_abs.csv              C(T) curve (latitudes)
    sescc_lon_pairs.csv            Cₚ(T) curve (longitudes)
    sescc_lat_abs_bootstrap.csv    bootstrap distribution (if --bootstrap N)
    sescc_lon_pairs_bootstrap.csv  bootstrap distribution (if --bootstrap N)

**Important for longitude dating:** exclude stars with very large parallax
(Keid HIP 19849, α Cen HIP 71681) and filter by positional error:

```bash
# Almagest (dist column = col 10)
awk -F';' 'substr($1,1,1)=="#" || $10+0 < 60' catalogs/almagest_verbunt.csv | \
  grep -v "^19849;" | grep -v "^71681;" | \
  python3 sescc.py --lon --dlon 30

# Tycho Brahe (dist column = col 8)
awk -F';' 'substr($1,1,1)=="#" || $8+0 < 60' catalogs/brahe_verbunt.csv | \
  grep -v "^19849;" | grep -v "^71681;" | \
  python3 sescc.py --lon --dlon 30
```

---

## Catalogue converters

### `verbunt2sescc.py` — Almagest and Ulugh Beg

Converts Verbunt & van Gent (2012) fixed-width files to SESCC CSV format.

    python3 verbunt2sescc.py ptolema.dat   > catalogs/almagest_verbunt.csv
    python3 verbunt2sescc.py ulughbeg.dat  > catalogs/ulughbeg_verbunt.csv

Column 10 (`dist_arcmin`) gives the positional error for filtering.

### `verbunt2sescc_brahe.py` — Tycho Brahe

Converts KeplerE.dat (Verbunt & van Gent 2010, A&A 516, A28) to SESCC CSV format.
Coordinates are ecliptic, equinox AD 1601.0.

    python3 verbunt2sescc_brahe.py keplere.dat > catalogs/brahe_verbunt.csv

Column 8 (`dist_arcmin`) gives the positional error for filtering.

---

## Statistical validation: `validate_pairs_seeds.py`

Generates synthetic catalogues at a known epoch with realistic Gaussian noise,
runs SESCC (`--lon`) on each, and reports bias, standard deviation, and the
fraction of resamples giving a pre-Christian minimum.

    python3 validate_pairs_seeds.py --epoch -127 --seeds 20
    python3 validate_pairs_seeds.py --epoch  137 --seeds 20

Options:

    --epoch YEAR    true epoch of synthetic catalogue (default -127)
    --seeds N       number of random seeds (default 20)
    --noise-lat L   latitude noise in arcminutes (default 23)
    --noise-lon L   longitude noise in arcminutes (default 27)
    --dlat D        max latitude difference for pairing (default 5.0)
    --dlon D        max longitude difference for pairing (default 30.0)
    --exclude HIPS  comma-separated HIPs to exclude (default 19849,71681)
    --input CSV     star list (default catalogs/almagest_verbunt.csv)
    --noplot        suppress histogram plot

---

## Position cache

All computed positions are cached in `sescc_positions.pkl.xz` (LZMA-compressed
pickle, included in the repository). If the file is not found, positions are
computed on the fly — this may take several hours for a full 2500-year scan.

Key format:

    "{hip}_{epoch_index}_lat"   ecliptic latitude, milli-degrees (frame of date)
    "{hip}_{epoch_index}_lon"   ecliptic longitude, milli-degrees (J2000 inertial)

where `epoch_index` maps to year `1900 − epoch_index` (step = 1 year, range
1900 to −600 BCE). Longitude entries use the J2000 inertial frame, essential
for SESCC-pairs to cancel precession exactly. Writes are atomic.

---

## Validation results

### `--lat` (ecliptic latitudes)

| Catalogue | True epoch | Minimum | Bootstrap 68% | % pre-CE |
|---|---|---|---|---|
| Tycho Brahe | ~1580 CE | 1570 CE | — | 0% |
| Ulugh Beg | 1437 CE | 1177 CE | [1126, 1182] | 0% |
| Almagest | ? | −49 BCE | [−177, +110] | 74% |

### `--lon` (ecliptic longitudes, excl. HIP 19849 & 71681, dist < 60')

| Catalogue | True epoch | Minimum | Bootstrap 68% | Std | % pre-CE |
|---|---|---|---|---|---|
| Tycho Brahe | ~1580 CE | 1547 CE | [1535, 1582] | 27 yr | 0% |
| Ulugh Beg | 1437 CE | 1452 CE | [1186, 1500] | 187 yr | 0% |
| Synth. −127 BCE (mean, 20 seeds) | −127 BCE | −125 BCE | — | 179 yr | 85% |
| Synth. +137 CE (mean, 20 seeds) | +137 CE | +148 CE | — | 170 yr | 20% |
| Almagest | ? | −165 BCE | [−384, +176] | 235 yr | 74% |

The 74% pre-Christian bootstrap fraction for the Almagest (both coordinates
independently) is consistent with a Hipparchan origin and inconsistent with
a Ptolemaic one. Invariance to global longitude offsets verified empirically:
13 offset versions of both the Almagest and Ulugh Beg (±6°) all give
identical epoch estimates.

---

## Dependencies

    pip install skyfield numpy matplotlib scipy

Skyfield downloads the DE441 ephemeris and Hipparcos catalogue automatically
on first run (~600 MB total).

---

## References

- Verbunt, F. and van Gent, R. H. (2010). A&A 516, A28. (Tycho Brahe)
- Verbunt, F. and van Gent, R. H. (2012). A&A 544, A31. (Almagest, Ulugh Beg)
- ESA (1997). The Hipparcos and Tycho Catalogues. ESA SP-1200.
- Park et al. (2021). AJ 161, 105. (DE441 ephemeris)
- Rhodes, B. (2019). Skyfield. ASCL 1907.024.
