import os
import pandas as pd
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

print('Tesselation = ')
NSIDE = int(input())

print('Energy step = ')
estep = int(input())

print('What do you want to observe? HEALPix original: 1 - counts, 2 - exposure time, 3 - count rate')
obs = int(input())
coordsys = 0;

print('What do you want to observe coordinate system do you want (gal, ecl, rib, np)?')
coordsys = input()

dataname = os.path.join(os.path.join(os.getcwd(), "output"), "data_t{}_{}Null.txt".format(NSIDE,estep))

X = pd.read_csv(dataname, sep="\t", header=None)

print(
    "Approximate resolution at NSIDE {} is {:.2} deg".format(
        NSIDE, hp.nside2resol(NSIDE, arcmin=True) / 60
    )
)
NPIX = hp.nside2npix(NSIDE)
if obs==1:
    vals = X[1].values
if obs==2:
    vals = X[2].values
if obs==3:
    vals = X[3].values
if obs ==4:
    vals = X[1].values
print(obs)
if obs ==5:
    vals = X[1].values

m = vals

if coordsys == 'gal':
	hp.mollview(m, coord=['E','G'], cmap=plt.cm.jet, title="Mollview image RING")
if coordsys == 'ecl':
	hp.mollview(m, cmap=plt.cm.jet, title="Mollview image RING")
if coordsys == 'rib':
	hp.mollview(m, rot =(227, 35 ,0), cmap=plt.cm.jet, title="Mollview image RING")
if coordsys == 'np':
	hp.mollview(m, rot =(0, 90 ,75), cmap=plt.cm.jet, title="Mollview image RING")
hp.graticule()

plt.show()
