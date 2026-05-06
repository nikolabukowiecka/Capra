# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# import os


# print('Dataname : ')
# dat = input()
# print('Energy step = ')
# estep = int(input())

# dataname = os.path.join(os.getcwd(), "output/Gaussian", f"{dat}")
# Y = pd.read_csv(dataname, sep="\t", header=None)
# arr = Y.values
# lat_size, lon_size = arr.shape


# lon = np.linspace(-180, 180, lon_size)
# lat = np.linspace(90, -90, lat_size)
# Lon, Lat = np.meshgrid(lon, lat)


# arr_log = np.log10(np.clip(arr, 1e-6, None))

# # Nose-centered rotation
# lon_nose, lat_nose = 255.7, 5.1  # ecliptic nose
# # shift longitude grid so the nose is centered
# shift_cols = int(round((lon_nose - 180) / (360 / lon_size)))
# arr_centered = np.roll(arr_log, -shift_cols, axis=1)

# # Flip to match Healpy east=left convention
# arr_centered = np.fliplr(arr_centered)

# # diag
# # print("Map centered on λ=255.7° (Ecliptic frame)")
# # print(f"Longitude bins: {lon_size}, Latitude bins: {lat_size}")
# # print("Nose-centered, Healpy-style handedness (east = left).")
# # print("North ecliptic pole is up (β increasing upward).")

# fig = plt.figure(figsize=(8, 6))
# ax = fig.add_subplot(111, projection='mollweide')

# im = ax.pcolormesh(
#     np.radians(lon),
#     np.radians(lat),
#     arr_centered,
#     cmap='jet',
#     shading='auto',
#     vmin=-0.5, vmax=0.5
# )

# ax.grid(True, color='white', linestyle='--', alpha=0.5)

# #ax.set_title(f"IBEX-Hi raw data $E_{estep}$. Nose-centered ecliptic frame.",
# ax.set_title(f"Capra gridded, $E_{estep}$. Nose-centered ecliptic frame.",
# #ax.set_title(f"Theseus, $E_{estep}$. Nose-centered ecliptic frame.",
# #ax.set_title(f"IBEX-Hi, $E_{estep}$. Nose-centered ecliptic frame.",
#              fontsize=13, color="black", y=1.03, pad=25)

# cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, fraction=0.07)
# cb.set_label("$\log_{10}(\mathrm{FluxGDF})\ [\mathrm{dimensionless}]$", fontsize=10, color="black")

# # markers
# ax.scatter(0, np.radians(lat_nose), s=100, facecolor="k", edgecolor="white", zorder=10)
# ax.text(0, np.radians(lat_nose + 3), "Nose", color="black",
#         ha="center", va="bottom", fontsize=12, weight="bold")

# # ax.scatter(0, np.radians(90), s=60, facecolor="cyan", edgecolor="black", zorder=10)
# # ax.text(0, np.radians(85), "N ecl", color="cyan",
# #         ha="center", va="top", fontsize=10, weight="bold")

# plt.show()


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import cmcrameri.cm as cmc

print('Dataname : ')
dat = input().strip()
print('Energy step = ')
estep = int(input())

dataname = os.path.join(os.getcwd(), "output/Gaussian", f"{dat}")

# Treat these as missing (add more tokens if your files use different strings)
Y = pd.read_csv(
    dataname,
    sep="\t",
    header=None,
    na_values=["Null", "NULL", "null", "", "NaN", "nan"],
    keep_default_na=True
)

Y = pd.read_csv(dataname, sep="\t", header=None, dtype=str)

Y = Y.apply(lambda col: col.str.strip())  # remove leading/trailing whitespace

Y = Y.replace({"Null": np.nan, "NULL": np.nan, "null": np.nan, "": np.nan})
Y = Y.apply(pd.to_numeric, errors="coerce")
arr = Y.to_numpy(dtype=float)



lat_size, lon_size = arr.shape

lon = np.linspace(-180, 180, lon_size)
lat = np.linspace(90, -90, lat_size)

# x0 = 5.0

# def signed_log(x, x0=5.0):
#     return np.sign(x) * np.log10(1.0 + np.abs(x) / x0)

# arr_log = np.full_like(arr, np.nan, dtype=float)
# ok = np.isfinite(arr)
# arr_log[ok] = signed_log(arr[ok], x0=x0)

# #linear
# arr_plot = np.full_like(arr, np.nan, dtype=float)
# ok = np.isfinite(arr)
# arr_plot[ok] = arr[ok]

#log10
arr_plot = np.full_like(arr, np.nan, dtype=float)
ok = np.isfinite(arr) & (arr > 0)
arr_plot[ok] = np.log10(arr[ok])

# Nose-centered rotation
lon_nose, lat_nose = 255.7, 5.1  # ecliptic nose

# shift longitude grid so the nose is centered
shift_cols = int(round((lon_nose - 180) / (360 / lon_size)))
arr_centered_plot = np.roll(arr_plot, -shift_cols, axis=1)

# Flip to match Healpy east=left convention
arr_centered_plot = np.fliplr(arr_centered_plot)

arr_centered_plot = np.ma.masked_invalid(arr_centered_plot)

#cmap = plt.cm.get_cmap('batlow').copy()
#cmap = 'turbo'
cmap = plt.colormaps['turbo'].resampled(256).copy()
#cmap.set_bad('black')

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='mollweide')

im = ax.pcolormesh(
    np.radians(lon),
    np.radians(lat),
    arr_centered_plot,
    cmap=cmap,
    shading='auto',
    vmin = 0.17837,
    #vmax = 1.509776
    # vmin=0.17837,
     vmax=1.3
)

print(np.min(arr_centered_plot))
print(np.max(arr_centered_plot))


ax.grid(True, color='white', linestyle='--', alpha=0.5)

# ax.set_title(
#     #f"IBEX gridded, $E_{estep}$. Nose-centered ecliptic frame.",
#     #f"Theseus - HP, $E_{estep}$. Nose-centered ecliptic frame.",
#     f"Capra - HP, $E_{estep}$. Nose-centered ecliptic frame.",
#     fontsize=13, color="black", y=1.03, pad=25
# )

cb = fig.colorbar(im, ax=ax, orientation="horizontal", pad=0.05, fraction=0.07)
cb.set_label(r"$\log_{10}(\mathrm{Flux})\ [\mathrm{dimensionless}]$", fontsize=10, color="black")
#cb.set_label(r"$\mathrm{Flux}$", fontsize=10, color="black")

# markers
ax.scatter(0, np.radians(lat_nose), s=100, facecolor="k", edgecolor="white", zorder=10)
ax.text(0, np.radians(lat_nose + 3), "Nose", color="black",
        ha="center", va="bottom", fontsize=12, weight="bold")

plt.show()
