# sescc

## **Speed/Error Signal Cross Correlation** is a method for dating stellar catalogs in ecliptical coordinates

2024 Carlos Baiget Orts

asinfreedom@gmail.com

**WARNING: This program downloads 1.7Gb of ephemeris from NASA in its first run.** 
more info: https://rhodesmill.org/skyfield/planets.html

**sescc was made to astronomically date Ptolemy's Almagest.** 

This is the sescc dating program, you can use it to **SCIENTIFICALLY PROVE** that the Ptolemy's Almagest was initially compiled by Hipparchus and Ptolemy in the -2nd and 2nd century C.E.

sescc can reliably date the Almagest by latitudes or longitudes, as well as any other catalog compiled in ecliptic coordinates, between -1500 and 1900 C.E.

This program requires Python, Skyfield, Numpy and Panda libraries:

+ https://www.python.org/
+ https://rhodesmill.org/skyfield/
+ https://numpy.org/
+ https://pandas.pydata.org/

## some results:

![Almagest dating by latitudes.](plots/alm_lat.png)

![Almagest dating by longitudes.](plots/alm_lon.png)

## sescc.py:

**Basic usage in LINUX:**

### 0. Install libraries:

pip install pandas skyfield numpy
 
### 1.Date Almagest by latitudes:

cat catalogs/almagest.csv | ./sescc.py 0

(Load generated .csv file to spreadsheet, then graph a dispersion plot with the data.)

### 2.Date Almagest by longitudes:

cat catalogs/almagest.csv | ./sescc.py 1

### 3.Date New Chronology's 'informative kernel' by latitudes (longitudes):

cat catalogs/fkn_kernel.csv | ./sescc.py 0 (1)

### 4.Date New Chronology's 'informative kernel' without Arcturus by latitudes (longitudes):

cat catalogs/fkn_wo_arcturus.csv | ./sescc.py 0 (1)

**NOTE**: "New Chronology's informative kernel" is the **WORST** set of stars to date the Almagest. It consists of a **selection of stars** that will give a late date. 
Those stars where updated by later astronomers, Battani in 9th Cy, Regiomontanus in 15 cy.

Dissection of "New Chronology's" "informative kernel":

+ Procyon, Capella: Regiomontanus (15th century)
+ Arcturus: Al-Battani (9th century)
+ Regulus, Spica, Vega, Antares, Asellus: Hipparchus (-2nd century)

"New Chronology's" dating is dominated by Arcturus due its much higher relative proper motion:
![FKN kernel relative speeds.](plots/rel_speeds.png)

"informative kernel" w/o Arcturus is dominated by Regiomontanus stars:
![FKN kernel w/o Arcturus.](plots/fkn_lat_sin_arcturus.png)

That's how they achieved a 10th Century date.

If this error is not admitted, it's just **scientific fraud in plain sight.**

If admitted, an honest explanation be explored:

The selection criteria looking for **bright and quick** stars **had an implicit and lethal risk**: that precisely those were most carefully tracked by later astronomers who of course **noticed the error and updated the catalog.**

Results indicate that this is what happened with al-Battani, who updated Arcturus in the 9th century, and Regiomontanus, who updated Procyon and Capella in the 15th.
Those updates made their way into the canonical edition of the Almagest in the 15th century, by none other than Regiomontanus.
Later researchers thought erroneously that the whole catalog was from Ptolemy. Detected incosistencies, attributed him a fraud. 

However, the main result of SESCC is that **Star positions and established chronology of the Almagest are fully coherent.**

Clarifications:
+ **Any other consideration** about **any other information** within the Almagest **requires a separate analisys.**
+ SESCC results **do not contradict the conclusions of the book "The Crime of Claudius Ptolemy"** from R.R. Newton (1977). In fact, its a corroboration of one of the main points in it: that it contains original work from Hipparchus.
+ As SESCC shows, this part is precisely the LATITUDES part of almost all the star positions.


**Syntax for WINDOWS:**

### 0. Install python, after that, libraries like in linux (step 0)

### 1.Date Almagest by latitudes:

python sescc.py < catalogs/almagest.csv

### 2.Date Almagest by longitudes:

python sescc.py 1 < catalogs/almagest.csv

etc...

## stella.py:

Creates a catalog for any epoch, list of stars is taken from another catalog, example:

cat catalogs/fkn_kernel.csv | ./stella.py 1100 10 0 0 | ./sescc.py

This takes the list of stars detailed in 'fkn_kernel.csv' and creates a new catalog for year 1100, with 10' resolution (using the fractional system of the Almagest), with 0' of systemical error and 0' of random error.

You can feed the output to sescc.py to verify it works as expected. stella.py can also be used to prove that sescc.py doesn't date the longitudes by their magnitude. Magnitudes can be adjusted for the axial precession of any other epoch and the results will be the same.

An example of this modification is the catalog 'almagest_fake_longs.csv', where longitudes are increased by 15º

## catalogs:

Format:
star_HIPPARCOS_ code ; latitude ; longitude

'0;' before star code means a reject, which will be ignored when loading catalog

There are just three rejects:

HIP71683 - Rigil Kentaurus

HIP69673 - Arcturus

HIP97649 - Altair

You can explore the effect of uncommenting those lines in the result.

First line must be HIP35550 (Delta Geminorum), which is a surrogate of the equinox as the reference point for any epoch, rather than the one employed by the computer, which may differ from that of the ancient astronomer.
More info: https://cliolapso.blogspot.com/2024/04/fin-de-la-nueva-cronologia.html

*************************************************************************************************************
NOTA a los lectores de mi blog, el cual cerré sin previo aviso:

Abrí el blog en 2018 para generar debate sobre un tema que captó mi interés.

Nunca cumplió su objetivo: ni generó debate ni mantuvo mi interés.

Su contenido era un torro y no aportaba nada que no se pueda ver en otro lugar, haber dejado el contenido suponía hacer una revisión de todos los artículos, pues llegué a la conclusión de que debían tener errores, y dudé de lo que había asegurado.

Como digo, es un trabajo que no se puede hacer sin el interés que requiere hacerlo.

Así que lo borro todo y con ello, cierro una etapa.

Me habría despedido, pero no tenía de quién: nunca tuvo suscriptores y creo que más de la mitad de los que visitaban este blog habían acabado aquí por casualidad.

Por supuesto es posible que alguien lo eche de menos, y le parezca extraño el cierre repentino, por eso dejo esta nota. No hay nada de extraño, yo no tenía presencia en internet antes de este hobby y quiero seguir sin tenerla después de él.

Hice copia de esta nota en archive.org por si alguien llega allí buscando esta información, y lo borré. Sin embargo aparentemente no lo hice bien y guardé solo la entrada, no el blog. Ahora no la puedo encontrar.

_"El problema con el mundo es que los tontos están seguros de todo y los inteligentes están llenos de dudas"._ Bertrand Russell. 
