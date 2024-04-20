#!/usr/bin/python3

import csv, sys, math
import numpy as np
from skyfield.api import Star, load
from skyfield.data import hipparcos
from operator import itemgetter

maxmag = 3 
filtro = 0
dsource = 0
n = len(sys.argv)

def error():
	print ("Usage: "+sys.argv[0]+" [DATING_SOURCE] [SPEED_FILTER] [MAX_MAGNITUDE] ")
	print ()
	print ("Dating Source: 0-latitudes, 1-longitudes. (default: 0)")
	print ("Magnitude: exclude stars whose magnitude > MAX_MAGNITUDE. (default: 3)")
	print ("Speed filter: date catalog using only the SPEED_FILTER stars. (default: no filter)")
	print ("Examples:")
	print (sys.argv[0]+" 0 20")
	print (sys.argv[0]+" 1 20 4")
	exit(0)

if n>1:
	dsource=int(sys.argv[1])
	if dsource!=0 and dsource!=1:
		error()
if n>2:
	filtro = int(sys.argv[2])
	if filtro < 2:
		error()
if n>3:
	maxmag = float(sys.argv[3])
if n>4:
	error()

siglos=30 #centuries to scan into the past
resolucion=25 #every # years (10 = every decade)
maxt=siglos*100//resolucion #*resolucion #number of iteration
fechamax=1900 #from year

def ecpos(hip,t):
	hip=int(hip)
	S = Star.from_dataframe(df.loc[hip])
	tt = ts.utc(fechamax-(t*resolucion),1,1)
	apparent=earth.at(tt).observe(S).apparent()
	lat, lon, distance = apparent.ecliptic_latlon(tt)
	return lat.degrees, lon.degrees
	
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

with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

#Warning: Downloads 1.5Gb! more info: https://rhodesmill.org/skyfield/planets.html
planets = load('https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de441_part-1.bsp') 
earth = planets['earth']
ts = load.timescale()

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
estrellas = list(csv.reader(sys.stdin, delimiter=';', quoting=csv.QUOTE_NONNUMERIC))
print()
almagest=[]
for i in range(0,len(estrellas)):
	print("Loading catalog ("+str(int(100*i/len(estrellas)))+"%)", end="\r")
	hip = int(estrellas[i][0])
	if hip==0:
		hip = int(estrellas[i][1])
		if hip==0:
			continue
		print()
		print (f"HIP{hip} excluded")
		print()
		continue
	pos = int(estrellas[i][dsource+1] * 1000)
	try:
		vel=abs(pmotion(hip)[dsource])
	except ValueError:
		print(f"Not able to compute longitudinal proper motion of star: HIP{hip}")
		continue
	try:
		mag = hip_mag[hip]
	except KeyError:
		print(f"Not able to get magnitude of star: HIP{hip}")
		continue
	if mag>maxmag and hip!=35550: #HIP35550 = Delta Geminorum (mag>3), must be always included as the implicit 0,0 coordinate.
		continue
	almagest.append([hip,vel,pos])

print("Loading catalog (100%)",end="\r")
print()
print (f"Max. Magnitude: {maxmag}")

if filtro!=0:
	#saco la referencia para ordenar la lista (la referencia ha de estar en primera posicion)
	referencia=almagest.pop(0)
	almagest = sorted(almagest, key=itemgetter(1),reverse=True)
	almagest = almagest[:filtro]
	#vuelvo a colocar la referencia en primera posicion
	almagest.insert(0,referencia)
	print ("Speed filter ON")
	filtro = "fastest & "
else:
	filtro = ""

n = len(almagest)
src='latitudes' if dsource==0 else 'longitudes'
print(f"Dating the {n} "+filtro+f"brightest stars of the catalog by {src}")

#calcular posiciones
for i in range(n):
	print("Computing past positions ("+str(int(100*i/n))+"%)", end="\r")
	for t in range(maxt):
		almagest[i].append(round(ecpos(almagest[i][0],t)[dsource]*1000))

#ajustar posiciones a la referencia (almagest[0][t+2])
for i in range(1,n):
	print("Computing past positions ("+str(int(100*i/n))+"%)", end="\r")
	for t in range(maxt+1):
		almagest[i][t+2]-=almagest[0][t+2]
		if almagest[i][t+2]<0:
			almagest[i][t+2]+=360000

almagest.pop(0)
almagest=np.array(almagest,dtype='int32')
distancias=[]
a = math.comb(n,2)
velocidades_catalogo=almagest[:,1]
posiciones_catalogo=almagest[:,2]

#scc:
year_corr=[]
def coef_epocas():
	global velocidades_catalogo
	for t in range(maxt-1,-1,-1):
		year_corr.append([fechamax-t*resolucion,np.correlate(velocidades_catalogo,np.abs(np.subtract(posiciones_catalogo,almagest[:,t+3])))[0],t])
#/scc!

coef_epocas()
year_corr=np.array(year_corr,dtype='int64')
corrmin = np.min(year_corr[:,1])
corrmax = np.max(year_corr[:,1])
pos_corrmin = np.argmin(year_corr[:,1])
e=year_corr[pos_corrmin][2]

#normalizar valores de correlacion entre 0 y 1000:
year_corr[:, 1] = ((year_corr[:, 1] - corrmin) / (corrmax - corrmin)) * 1000

filename="sescc_dating_"+src+".csv"

f = open(filename, "w")
for elem in year_corr:
	f.write (str(int(elem[0]))+';'+str(elem[1])+'\n')
f.close()

print ()
print (f"Output in {filename}")
print ()

#Uncomment to print working set: HIP star, ProperMotion, Magnitude
"""
almagest = sorted(almagest, key=itemgetter(1),reverse=True)
for elem in almagest:
	print (elem[0],elem[1],elem[2])
"""

