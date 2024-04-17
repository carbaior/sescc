# sescc
sescc (Speed/Error Signal Cross Correlation) is a method for dating stellar catalogs in ecliptical coordinates

2024 Carlos Baiget Orts
asinfreedom@gmail.com
cliolapso.blogspot.com

WARNING: This program downloads 1.5Gb of data in its first run. more info: https://rhodesmill.org/skyfield/planets.html

sescc was made to astronomically date Ptolemy's Almagest. 
more info: https://cliolapso.blogspot.com/2024/04/fin-de-la-nueva-cronologia.html

This is the sescc dating program, you can use it to PROVE the Almagest was initially compiled by Hipparchus and Ptolemy in the -2nd and 2nd century C.E.

sescc can reliably date the Almagest by latitudes or longitudes, as well as any other catalog compiled in ecliptic coordinates, between -1500 and 1900 C.E.

This program requires Python and the Skyfield and Numpy and Panda libraries:

https://www.python.org/
https://rhodesmill.org/skyfield/
https://numpy.org/
https://pandas.pydata.org/

sescc.py:

** Basic usage in LINUX: **

0. Install libraries:

 pip install pandas skyfield numpy
 
1.Date Almagest by latitudes:

cat catalogs/almagest.csv | ./sescc.py 0

(Load .csv file to spreadsheet, then graph a dispersion plot with the data.)

2.Date Almagest by longitudes:

cat catalogs/almagest.csv | ./sescc.py 1

3.Date New Chronology 'informative kernel' by latitudes (longitudes):

cat catalogs/fkn_kernel.csv | ./sescc.py 0 (1)

4.Date New Chronology 'informative kernel' without Arcturus by latitudes (longitudes):

cat catalogs/fkn_wo_arcturus.csv | ./sescc.py 0 (1)

** Syntax for WINDOWS: **

0. Install python, after that libraries like in linux (step 0)

1.Date Almagest by latitudes:

python sescc.py < catalogs/almagest.csv

2.Date Almagest by longitudes:

python sescc.py 1 < catalogs/almagest.csv

etc...

stella.py:

Creates a catalog for any epoch, list of stars is taken from another catalog, example:

cat catalogs/fkn_kernel.csv | ./stella.py 1100 10 0 0 | ./sescc.py

This takes the list of stars detailed in 'fkn_kernel.csv' and creates a new catalog for year 1100, with 10' resolution (using the fractional system of the Almagest), with 0' of systemical error and 0' of random error.

You can feed the output to sescc.py to verify it works as expected. stella.py can also be used to prove that sescc.py doesn't date the longitudes by their magnitude. Magnitudes can be adjusted for the axial precession of any other epoch and the results will be the same.

An example of this modification is the catalog 'almagest_fake_longs.csv', where longitudes are increased by 15º

catalogs:

star in HIPPARCOS code; latitude; longitude

'0;' before star code means a reject, which will be ignored when loading catalog

There are just three rejects:

HIP71683 - Rigil Kentaurus

HIP69673 - Arcturus

HIP97649 - Altair

You can explore the effect of uncommenting those lines in the result.

First line must be HIP35550, which is a surrogate of the equinox as the reference point for any epoch, rather than the one employed by the computer, which may differ from that of the ancient astronomer.
More info: https://cliolapso.blogspot.com/2024/04/fin-de-la-nueva-cronologia.html
