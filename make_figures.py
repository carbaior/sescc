#!/usr/bin/env python3
"""
make_figures.py - Regenerate the Almagest figures of the SESCC paper.

Reads the files written by sescc.py:
    sescc_lat.csv             dating curve, latitudes      (year;value)
    sescc_lon.csv             dating curve, longitudes     (year;value)
    sescc_lat_bootstrap.csv   bootstrap minima, latitudes  (one year per line)
    sescc_lon_bootstrap.csv   bootstrap minima, longitudes (one year per line)

Writes (PDF, vector; same base names used in main.tex):
    almagest_latitudes.pdf   almagest_longitudes.pdf
    fig_bootstrap_lat.pdf    fig_bootstrap_lon.pdf

Usage:
    python3 make_figures.py [--png]
"""
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXT = "png" if "--png" in sys.argv else "pdf"
HIPPARCHUS, PTOLEMY = -127, 137

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 10,
    "axes.grid": True,
    "grid.alpha": 0.3,
})


def year_label(y):
    """Signed year -> 'N BCE' / 'N CE', as in the text of the paper."""
    y = int(round(y))
    return f"{-y} BCE" if y < 0 else f"{y} CE"


def load_curve(path):
    data = np.loadtxt(path, delimiter=";")
    data = data[np.argsort(data[:, 0])]
    return data[:, 0], data[:, 1]


def plot_curve(path, title, ylabel, color, out):
    x, y = load_curve(path)
    tmin = x[np.argmin(y)]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(x, y, color=color, lw=1.4)
    ax.axvline(tmin, color="black", ls="--", lw=1,
               label=f"Min: {year_label(tmin)}")
    ax.axvline(HIPPARCHUS, color="gray", ls=":", lw=1,
               label=f"Hipparchus ({year_label(HIPPARCHUS)})")
    ax.axvline(PTOLEMY, color="tab:orange", ls=":", lw=1,
               label=f"Ptolemy ({year_label(PTOLEMY)})")
    ax.set_xlim(x.min() - 50, x.max() + 50)
    ax.set_ylim(0, 1050)
    ax.set_xlabel("Year")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"{out}: minimum at {tmin:.0f}")


def plot_bootstrap(path, title, color, out):
    m = np.loadtxt(path, dtype=int, ndmin=1)
    med = np.median(m)
    p16, p84 = np.percentile(m, [16, 84])
    pre = 100 * (m < 0).mean()
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.hist(m, bins=60, color=color, edgecolor="white", lw=0.3)
    ax.axvline(med, color="black", lw=1.2, label=f"Median: {year_label(med)}")
    ax.axvline(HIPPARCHUS, color="gray", ls="--", lw=1,
               label=f"Hipparchus ({year_label(HIPPARCHUS)})")
    ax.axvline(0, color="lightgray", lw=1, label="0 CE")
    ax.set_xlabel("Year of minimum")
    ax.set_ylabel("Count")
    ax.set_title(f"{title}\nMedian {year_label(med)}  |  "
                 f"68% = [{year_label(p16)}, {year_label(p84)}]  |  "
                 f"{pre:.0f}% pre-CE", fontsize=9)
    ax.legend(loc="upper right", fontsize=8)
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)
    print(f"{out}: N={len(m)}, median {med:.0f}, 68% [{p16:.0f}, {p84:.0f}], "
          f"{pre:.1f}% pre-CE, {(m <= -590).sum()} at scan edge")


def plot_validation(path, title, ylabel, color, true_epoch, out):
    x, y = load_curve(path)
    tmin = x[np.argmin(y)]
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(x, y, color=color, lw=1.4)
    ax.axvline(tmin, color="black", ls="--", lw=1, label=f"Min: {year_label(tmin)}")
    ax.axvline(true_epoch, color="gray", ls=":", lw=1,
               label=f"True epoch: {year_label(true_epoch)}")
    ax.set_xlim(x.min() - 50, x.max() + 50); ax.set_ylim(0, 1050)
    ax.set_xlabel("Year"); ax.set_ylabel(ylabel); ax.set_title(title)
    ax.legend(loc="upper right", fontsize=8); ax.tick_params(axis="x", rotation=45)
    fig.tight_layout(); fig.savefig(out); plt.close(fig)
    print(f"{out}: minimum at {tmin:.0f}")


# Validation curve mode:
#   python3 make_figures.py --curve CSV OUTNAME "TITLE" TRUE_EPOCH [lat|lon] [--png]
if "--curve" in sys.argv:
    k = sys.argv.index("--curve")
    csv_path, outname, title, true_ep = sys.argv[k+1:k+5]
    coord = sys.argv[k+5] if len(sys.argv) > k+5 and sys.argv[k+5] in ("lat", "lon") else "lon"
    ylabel = "Normalised $C(T)$" if coord == "lat" else "Normalised $C_p(T)$"
    color = "tab:blue" if "brahe" in outname else "tab:green"
    plot_validation(csv_path, title, ylabel, color, int(true_ep), f"{outname}.{EXT}")
    sys.exit(0)

plot_curve("sescc_lat.csv", "Almagest \u2014 latitudes",
           "Normalised $C(T)$", "tab:red", f"almagest_latitudes.{EXT}")
plot_curve("sescc_lon.csv", "Almagest \u2014 longitudes",
           "Normalised $C_p(T)$", "tab:red", f"almagest_longitudes.{EXT}")
plot_bootstrap("sescc_lat_bootstrap.csv",
               "Bootstrap distribution \u2014 Almagest latitudes (SESCC)",
               "tab:blue", f"fig_bootstrap_lat.{EXT}")
plot_bootstrap("sescc_lon_bootstrap.csv",
               "Bootstrap distribution \u2014 Almagest longitudes (SESCC-pairs)",
               "tab:green", f"fig_bootstrap_lon.{EXT}")
