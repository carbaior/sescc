#!/usr/bin/python3

import time, random, csv, sys, math

from skyfield.api import Star, load
from skyfield.data import hipparcos

def eclatlon(hip,year):
	hip=int(hip)
	S = Star.from_dataframe(df.loc[hip])
	tt = ts.utc(year,1,1)
	apparent=earth.at(tt).observe(S).apparent()
	lat, lon, distance = apparent.ecliptic_latlon(tt)
	return lat.degrees,lon.degrees

with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

#Warning: Downloads 1.5Gb
planets = load('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp') 
earth = planets['earth']

ts = load.timescale()

estrellas = list(csv.reader(sys.stdin, delimiter=';', quoting=csv.QUOTE_NONNUMERIC))

n = len(sys.argv)

precision = 0
syserror = 0
rnderror = 0

if n>1:
	year = int(sys.argv[1])
if n>2:
	precision = int(sys.argv[2]) / 60
if n>3:
	syserror = float(sys.argv[3])
if n>4:
	rnderror = float(sys.argv[4])
if n==1 or n>5:
	print ("Usage: "+sys.argv[0]+" <YEAR> [RESOLUTION] [SYSTEMATIC ERROR] [RANDOM_ERROR]")
	print ()
	print ("Year: from -1900 to 1900. 0 = same as input")
	print ("Resolution: in arcminutes. 0 = max available.")
	print ("Systematic error: in arcminutes (example: -10). 0 for no error")
	print ("Random error interval (normal): in arcminutes (example: 10). 0 for no error")
	exit(0)

for i in range(0,len(estrellas)):
	hip = int(estrellas[i][0])
	if hip==0:
		continue
	
	if year!=0:
		lat,lon=eclatlon(hip,year)
	else:
		lat=estrellas[i][1]
		lon=estrellas[i][2]
	lat+=syserror/60
	lon+=syserror/60
	
	if lat>=360:
		lat-=360
	if lon>=360:
		lon-=360

	if rnderror!=0:
		rg=rnderror/2
		lat+=random.normalvariate(0, rg/3)/60
		lon+=random.normalvariate(0, rg/3)/60
	if precision!=0:
		lats = round(round (lat / precision) * precision,3)
		lons = round(round (lon / precision) * precision,3)
		latq = round(round (lat / 0.25) * 0.25,3)
		lonq = round(round (lon / 0.25) * 0.25,3)
		if abs(lat-lats)<abs(lat-latq):
			lat=lats
		else:
			lat=latq
		if abs(lon-lons)<abs(lon-lonq):
			lon=lons
		else:
			lon=lonq	
	print (str(hip)+";"+str(round(lat,2))+";"+str(round(lon,2)))
	
	
	
