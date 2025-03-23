import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os


#dataname = "transformata.txt"
print('What data you want to visualise (HP2IB, HP_Ribbon_2IB, IB, IB2)?')
datname = input()
print('Tesselation = ')
tess = int(input())
print('Energy step = ')
estep = int(input())
print('What pixel resolution (angular) do you want?')
print('First parameter:')
pixres1 = int(input())
print('Second parameter:')
pixres2 = int(input())

fig = plt.figure()
ax = fig.add_subplot(111, projection='mollweide')

lon = np.linspace(-np.pi, np.pi, int(360/pixres1))
#lat = np.linspace(-np.pi/2., np.pi/2., int(180/pixres2))
lat = np.linspace(np.pi/2., -np.pi/2., int(180/pixres2))
Lon,Lat = np.meshgrid(lon,lat)
#print(np.meshgrid(lon,lat))
dataname = os.path.join(os.path.join(os.getcwd(), "output"), "normalized{}_{}_{}_{}_{}.txt".format(datname,tess,estep,pixres1,pixres2))

Y = pd.read_csv(dataname, sep="\t", header=None)
arr = Y.values

#im = ax.pcolormesh(Lon,Lat,arr, cmap=plt.cm.jet)
#im = ax.pcolormesh(Lon,Lat,arr, cmap=plt.cm.jet, vmin=0, vmax=1.69)
im = ax.pcolormesh(Lon,Lat,arr, cmap=plt.cm.jet)
fig.colorbar(im, ax=ax)

plt.show()

