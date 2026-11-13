set -e
for c in almagest ulugh brahe; do
  cp date_${c}_lon.sh date_${c}_lon_v1.sh
  cp v2d_${c}_lon.sh date_${c}_lon.sh
done
bash date_brahe_lat.sh --noplot && python3 make_figures.py --curve sescc_lat.csv brahe_latitudes "Tycho Brahe — latitudes" 1580 lat
bash date_brahe_lon.sh --noplot && python3 make_figures.py --curve sescc_lon.csv brahe_longitudes "Tycho Brahe — longitudes" 1580 lon
bash date_ulugh_lat.sh --noplot && python3 make_figures.py --curve sescc_lat.csv ulughbeg_latitudes "Ulugh Beg — latitudes" 1437 lat
bash date_ulugh_lon.sh --noplot && python3 make_figures.py --curve sescc_lon.csv ulughbeg_longitudes "Ulugh Beg — longitudes" 1437 lon
bash date_almagest_lat.sh --bootstrap 1000 --noplot
bash date_almagest_lon.sh --bootstrap 1000 --noplot
python3 make_figures.py
