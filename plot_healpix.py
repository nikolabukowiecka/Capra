import os
import pandas as pd
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

print('Tesselation = ')
NSIDE = int(input())

print('Energy step = ')
estep = int(input())

print('What do you want to observe? HEALPix original: 1 - counts, 2 - exposure time, 3 - signal, 4 - signal intensity, 5 - ENAflux')
obs = int(input())

print('What do you want to observe coordinate system do you want (gal, ecl, rib, np, nose)?')
coordsys = input()

print('Type filename, exe: data_t64_6Null.txt?')
dat = input()

dataname = os.path.join(os.path.join(os.getcwd(), "output"), f"{dat}")

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
    vals = X[4].values
if obs ==5:
    vals = X[5].values

m = vals
print(m)

if coordsys == 'gal':
    fig = plt.figure(figsize=(12, 6))

    lon_nose, lat_nose = 255.7, 5.1  # ecliptic deg

    hp.mollview(
        m,
        coord=['E','G'],        # rotate map to Galactic
        title=f"Capra, $E_{estep}$, $N_{{side}}={NSIDE}$. Galactic frame.",
        #unit="Flux [counts s$^{-1}$ sr$^{-1}$]",
        unit = "$\log_{10}(\mathrm{FluxGDF})\ [\mathrm{dimensionless}]$",
        cmap="jet",
        flip='geo',
        notext=True,
        xsize=1600,
        fig=fig#,
        #min=-0.7,max=0.7 #for GDF
    )
    hp.graticule(dpar=30, dmer=30, color='white', linestyle='--', linewidth=0.8)

    # Galactic reference points
    #hp.projscatter(0, 0, lonlat=True, coord='G', color='red', s=80, marker='*')
    #hp.projtext(0, 0, '  GC', lonlat=True, coord='G', color='red', fontsize=10)
    #hp.projscatter(0, 90, lonlat=True, coord='G', color='magenta', marker='^', s=60)
    #hp.projtext(0, 90, '  NGP', lonlat=True, coord='G', color='magenta', fontsize=9)

    # Ecliptic → Galactic poles
    rot_e2g = hp.Rotator(coord=['E','G'], deg=True)
    lon_gal, lat_gal = rot_e2g(lon_nose, lat_nose, lonlat=True)

    hp.projscatter(lon_gal, lat_gal, lonlat=True, coord='G', color='black', marker='^', s=120, edgecolor='black', linewidth=0.8)
    hp.projtext(lon_gal, lat_gal, 'Nose', lonlat=True, coord='G', color='black', fontsize=12, fontweight='bold')
    

    # Plot separately
    #hp.projscatter(lon_gal[0], lat_gal[0], lonlat=True, coord='G', color='black', marker='^', s=120, edgecolor='black', linewidth=0.8)
    #hp.projtext(lon_gal[0], lat_gal[0], 'N ecl', lonlat=True, coord='G', color='black', fontsize=12, fontweight='bold')
    #hp.projscatter(lon_gal[1], lat_gal[1], lonlat=True, coord='G', color='yellow', marker='v', s=120)
    #hp.projtext(lon_gal[1], lat_gal[1], 'S ecl', lonlat=True, coord='G', color='yellow', fontsize=10)

    plt.show()


if coordsys == 'ecl':
    fig = plt.figure(figsize=(12, 6))
    hp.mollview(
        m,
        title=f"IBEX-Hi, $E_{estep}$, $N_{{side}}={NSIDE}$. Ecliptic frame.",
        unit="Flux [counts s$^{-1}$ sr$^{-1}$]",
        #unit = "Signal rate [counts / s]",
        cmap="jet",
        notext=True,
        xsize=1600,
        #flip='geo',
        fig=fig
        # min=-1,
        # max=4
    )

    # Get the current axes
    ax = plt.gca()

    # Add white dashed graticule (healpy version)
    hp.graticule(color="white", linestyle="--", linewidth=0.8, alpha=0.7, dpar=30, dmer=30)
    lon_nose, lat_nose = 255.7, 5.1  # ecliptic coordinates

    hp.projscatter(lon_nose, lat_nose, lonlat=True, color='black', s=80, marker='o')
    hp.projtext(lon_nose, lat_nose, '  Nose', lonlat=True, color='black', fontsize=10)

    #Sanity check:
    # North and South ecliptic poles
    #hp.projscatter(0, 90, lonlat=True, color='magenta', marker='^', s=60, label='North ecl. pole')
    #hp.projscatter(0, -90, lonlat=True, color='black', marker='v', s=60, label='South ecl. pole')
    #plt.legend()

    plt.show()

if coordsys == 'nose':

    lon_nose, lat_nose = 255.7, 5.1  # ecliptic deg

    # Rotate map so nose is centered
    rot = hp.Rotator(deg=True, rot=[lon_nose - 180, 0, 0], inv=True)
    m_nose = rot.rotate_map_pixel(m)

    # Plot centered at 180° (default Mollweide center)
    hp.mollview(np.log10(np.clip(m_nose, 1e-6, None)),
                coord='E',# flip='geo',
                title=f"Capra, $E_{estep}$, $N_{{side}}={NSIDE}$. Nose-centered ecliptic frame.",
                #title=f"Theseus - HP, $E_{estep}$, $N_{{side}}={NSIDE}$. Nose-centered ecliptic frame.",
                #unit="Flux [counts s$^{-1}$ sr$^{-1}$]",
                unit = "$\log_{10}(\mathrm{Flux})\ [\mathrm{dimensionless}]$",
                #unit = "$\log_{10}(\mathrm{FluxGDF})\ [\mathrm{dimensionless}]$",
                cmap='jet', xsize=1600,
                min=1,max=1.6)
                #min=-0.3,max=0.3) #for relative flux
                #min=-0.2,max=0.8) #for GDF
                #min=-0.7,max=0.7) #for GD

    hp.graticule(dpar=30, dmer=30, color='white', linestyle='--', alpha=0.7)

    hp.projscatter(0, 0, lonlat=True, color='black', s=120)
    hp.projtext(0, 0, " Nose", lonlat=True, color='black', fontsize=11, weight='bold')

    plt.show()

    # lon_c, lat_c = rot(lon_nose, lat_nose, lonlat=True)
    # print(f"Nose in rotated frame: λ'={lon_c:.1f}°, β'={lat_c:.1f}°")

if coordsys == 'rib':
    lon_rib, lat_rib = 221.0, 39.0
    lon_nose, lat_nose = 255.7, 5.1

    def sph2vec(lon, lat):
        lon, lat = np.radians(lon), np.radians(lat)
        cl = np.cos(lat)
        return np.array([cl*np.cos(lon), cl*np.sin(lon), np.sin(lat)])

    def vec2lonlat(v):
        x,y,z = v
        lon = np.degrees(np.arctan2(y, x))
        lat = np.degrees(np.arcsin(np.clip(z, -1, 1)))
        return lon, lat

    def angsep(l1,b1,l2,b2):
        l1,b1,l2,b2 = map(np.radians, [l1,b1,l2,b2])
        return np.degrees(np.arccos(
            np.sin(b1)*np.sin(b2)+np.cos(b1)*np.cos(b2)*np.cos(l2-l1)
        ))

    # unit vectors
    er = sph2vec(lon_rib, lat_rib)                        
    k  = np.array([0.0, 0.0, 1.0])                        
    # local east/north at ribbon center
    e_east  = np.cross(k, er);  e_east  /= np.linalg.norm(e_east)
    e_north = np.cross(er, e_east); e_north /= np.linalg.norm(e_north)

    # define ribbon frame: x' = toward ribbon, y' = local north, z' = x'×y'
    xprime = er
    yprime = e_north
    zprime = np.cross(xprime, yprime)                     # right-handed
    R = np.vstack([xprime, yprime, zprime])               # rows are basis vectors

    # nose in ribbon frame
    v_nose = sph2vec(lon_nose, lat_nose)
    xp, yp, zp = R @ v_nose
    lon_p, lat_p = vec2lonlat([xp, yp, zp])

    # separations
    sep_orig = angsep(lon_rib, lat_rib, lon_nose, lat_nose)
    sep_rot  = np.degrees(np.arccos(np.clip(xp*np.cos(0)+yp*np.sin(0)+zp*0, -1, 1)))  # == arccos(xp)
    sep_rot  = np.degrees(np.arccos(np.clip(xp, -1, 1)))


    # print(f"nose (ribbon frame): lon'={lon_p:.1f}°, lat'={lat_p:.1f}°")
    # print(f"sep(ribbon,nose) orig={sep_orig:.1f}°,  from rotated={sep_rot:.1f}°")

    # --- plot by letting mollview do the same recentering visually ---
    hp.mollview(np.log10(np.clip(m,1e-6,None)),
                coord='E', rot=(lon_rib, lat_rib, 0), #flip='geo',
                title=f"Capra, $E_{estep}$, $N_{{side}}={NSIDE}$. Ribbon-centered ecliptic frame",
                unit="$\log_{10}(\mathrm{FluxGDF})\ [\mathrm{dimensionless}]$", cmap='jet', min=-0.6, max=0.4)
    hp.graticule(dpar=30, dmer=30, color='white', linestyle='--', linewidth=0.8, alpha=0.6)

    # markers (give ecliptic coords; mollview applies same rot)
    hp.projscatter(lon_rib, lat_rib, lonlat=True, coord='E', color='black', marker='*', s=140)
    hp.projtext(lon_rib, lat_rib, '  Ribbon', lonlat=True, coord='E', color='black')

    hp.projscatter(lon_nose, lat_nose, lonlat=True, coord='E', color='k', s=120)
    hp.projtext(lon_nose, lat_nose, '  Nose', lonlat=True, coord='E', color='k')

    # orientation check: North ecliptic pole
    #hp.projscatter(0, 90, lonlat=True, coord='E', color='black', marker='^', s=120)
    #hp.projtext(0, 90, '  N ecl', lonlat=True, coord='E', color='black')

    plt.show()


plt.show()
