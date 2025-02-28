#!/usr/bin/env python3

#sescc is a program to date ancient catalogs compiled in ecliptical coordinates
#copyleft (GPLv3) 2024 Carlos Baiget Orts (asinfreedom@gmail.com)
#https://github.com/carbaior/sescc

import csv, sys, math, random
import numpy as np
import matplotlib.pyplot as plt
from skyfield.api import Star, load
from skyfield.data import hipparcos
from operator import itemgetter

#defaults:
maxmag = 2.6 
filtro = 0
dsource = 0
n = len(sys.argv)

def error():
	print ("Usage: "+sys.argv[0]+" [DATING_SOURCE] [RANDOM_GROUP] [MAX_MAGNITUDE] ")
	print ()
	print ("Dating Source: 0-latitudes, 1-longitudes. (default: 0)")
	print ("Magnitude: exclude stars whose magnitude < MAX_MAGNITUDE. (default: 2.6)")
	print ("Random Group: date catalog using a random group of stars of RANDOM_GROUP size  (default: all stars)")
	print ("Examples:")
	print (sys.argv[0]+" 0 70")
	print (sys.argv[0]+" 1 70 3")
	exit(0)

if n>1:
	dsource=int(sys.argv[1])
	if dsource!=0 and dsource!=1:
		error()
if n>2:
	filtro = int(sys.argv[2])
	if filtro < 0:
		error()
if n>3:
	maxmag = float(sys.argv[3])
if n>4:
	error()

siglos=30 #centuries to scan into the past
resolucion=25 #every # years (10 = every decade)
maxt=siglos*100//resolucion #number of iterations
fechamax=1900 #from year 'fechamax' to the past

#return ecliptical coordinates (lat,lon) for a star on an epoch t
def ecpos(hip,t):
	hip=int(hip)
	S = Star.from_dataframe(df.loc[hip])
	tt = ts.utc(fechamax-(t*resolucion),1,1)
	apparent=earth.at(tt).observe(S).apparent()
	lat, lon, distance = apparent.ecliptic_latlon(tt)
	return lat.degrees, lon.degrees

#compute proper motion of a star
#(not used to compare stars, so no need to adjust longitudinal speeds to latitude.)
def pmotion(hip):
	S = Star.from_dataframe(df.loc[hip])
	t2 = ts.utc(1000,1,1)
	t1 = ts.utc(0,1,1)
	apparent=earth.at(t2).observe(S).apparent()
	lat2, lon2, distance = apparent.ecliptic_latlon()
	apparent=earth.at(t1).observe(S).apparent()
	lat1, lon1, distance = apparent.ecliptic_latlon()
	vlat=round((lat2.degrees-lat1.degrees)*1000)
	vlon=round((lon2.degrees-lon1.degrees)*1000)
	return (vlat,vlon)

#load JPL ephemeris files:
#Warning: Downloads 1.5Gb! more info: https://rhodesmill.org/skyfield/planets.html
planets = load('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp') 
earth = planets['earth']
ts = load.timescale()

#load Hipparcos star database:
with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

#load magnitude for each star from Hipparcos database:
hip_mag = []
file = open("./hip_main.dat", "r")
for line in file:
	try:
		hip=line[8:14]
		hip=int(hip)
	except:
		continue
	try:
		mag=line[41:46]
		mag=float(mag)
		hip_mag.append([hip,mag])
	except:
		continue
file.close()
hip_mag = dict(hip_mag)

#load catalog to be dated from stdin:		
estrellas = list(csv.reader(sys.stdin, delimiter=';', quoting=csv.QUOTE_NONNUMERIC))
print()
almagest=[]
for i in range(1,len(estrellas)):
	print("Loading catalog ("+str(int(100*i/len(estrellas)))+"%)", end="\r")
	hip = int(estrellas[i][0])
	if hip==0:
		hip = estrellas[i][1]
		if hip!="":
			hip = int(estrellas[i][1])
			print()
			print (f"HIP{hip} excluded")
			print()
		continue
	pos = int(estrellas[i][dsource+1] * 1000)
	try:
		vel=abs(pmotion(hip)[dsource])
		#vel=pmotion(hip)[dsource]
		#faster stars were more likely to had its position updated by later astronomers
		#if you want those stars to be discarded, uncomment the following two lines (try your own values):
		#if vel>100:
		#	continue
	except ValueError:
		print(f"Not able to compute proper motion of star: HIP{hip}")
		continue
	try:
		mag = hip_mag[hip]
	except KeyError:
		print(f"Not able to get magnitude of star: HIP{hip}")
		continue
	#if mag>maxmag and hip!=35550: #HIP35550 = Delta Geminorum (mag>3), must be always included as the implicit 0,0 coordinate.
	#	continue
	almagest.append([hip,vel,pos])

print("Loading catalog (100%)",end="\r")
print()
print (f"Max. Magnitude: {maxmag}")

#option to date the catalog from a random group of stars:
if filtro!=0:
	random.shuffle(almagest)
	almagest = almagest[:filtro]
	filtro=str(filtro)
	print ("Random group of "+filtro+" stars")


n = len(almagest)
src='latitudes' if dsource==0 else 'longitudes'
print(f"Dating by {src} {n} stars of the catalog.")

#compute past position of each star on each epoch:
for i in range(n):
	print("Computing past positions ("+str(int(100*i/n))+"%)", end="\r")
	for t in range(maxt):
		almagest[i].append(round(ecpos(almagest[i][0],t)[dsource]*1000))

almagest=np.array(almagest,dtype='int64')

#stars velocities:
velocidades_catalogo=almagest[:,1]
#stars positions in the catalog:
posiciones_catalogo=almagest[:,2]

#compute the speed/error signals cross correlation (scc) for each epoch:
#SCC:
year_corr=[]
for t in range(maxt-1,-1,-1):
	year_corr.append([fechamax-t*resolucion,np.correlate(velocidades_catalogo,np.abs(np.subtract(posiciones_catalogo,almagest[:,t+3])))[0],t])
#/SCC!

#normalize correlation values, between 0 y 1000:
year_corr=np.array(year_corr,dtype='int64')
corrmin = np.min(year_corr[:,1])
corrmax = np.max(year_corr[:,1])
year_corr[:, 1] = ((year_corr[:, 1] - corrmin) / (corrmax - corrmin)) * 1000

#save results
filename="sescc_dating_"+src+".csv"

f = open(filename, "w")
for elem in year_corr:
	f.write (str(int(elem[0]))+';'+str(elem[1])+'\n')
f.close()

print ()
print (f"Output in {filename}")
print ()

#Uncomment this block to print working set: HIP star, ProperMotion, Magnitude:

print("Working set:")
print()
print(" HIP_code    ProperMotion")
almagest = sorted(almagest, key=itemgetter(1),reverse=True)

for row in almagest:
    aux=[row[0],row[1]]		
    print("{: >8} {: >8}".format(*aux))


# Plot results:
x = year_corr[:, 0]
y = year_corr[:, 1]
plt.figure(figsize=(12, 8))
color='blue' if dsource==0 else 'green'
plt.plot(x, y, color=color)
plt.ylim(0,1000)
plt.yticks([])
plt.xlabel('Year')
plt.xticks(range(int(min(x)), int(max(x))+1, 100), rotation=45)
plt.title("Almagest SESCC dating by "+src.upper())
plt.grid(True)
plt.show()
