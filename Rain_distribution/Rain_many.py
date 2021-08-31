# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated based on mean th and rv at cloud base cells

from sys import argv, path, maxsize
#path.insert(0,"../../local_folder/uptodate/lib/python3/dist-packages")
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
# path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")


import h5py
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from libcloudphxx import common as lcmn
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

start = time.perf_counter()
plt.rcParams.update({'font.size': 12})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i
timesteps = timesteps[1:91]
np.set_printoptions(threshold=maxsize)


paths = argv[1]
outfile = argv[2]
files = os.listdir(paths)
nr_files = len(files)


def Adia_fraction(timestep):

    rhod= [0 for i in range(nr_files)]
    p_e = [0 for i in range(nr_files)]
    dz = [0 for i in range(nr_files)]
    rl = [0 for i in range(nr_files)]
    rl_base = [0 for i in range(nr_files)]
    nc = [0 for i in range(nr_files)]
    th = [0 for i in range(nr_files)]
    rv = [0 for i in range(nr_files)]
    rr = [0 for i in range(nr_files)]
    sum_rr = [0 for i in range(nr_files)]


    for file in range(nr_files):
        # plt.clf()
        p_e[file] = h5py.File(paths+files[file] + "/const.h5", "r")["p_e"][:]
        rhod[file] = h5py.File(paths+files[file] + "/const.h5", "r")["G"][:,:]
        dz[file] = h5py.File(paths+files[file] + "/const.h5", "r").attrs["dz"]

        nx, nz = rhod[file].shape
        hght = np.arange(nz) * dz[file]
        bin = np.linspace(1,121,len(np.arange(nx)))
        sum_rr[file] = np.zeros([nx])

        filename = paths+files[file] + "/timestep" + str(int(timestep)*240).zfill(10) + ".h5"
        rl[file] = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        rl_base[file] = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        nc[file] = h5py.File(filename, "r")["cloud_rw_mom0"][:,:] * rhod[file] / 1e6; # 1 / cm^3
        rr[file] = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * rhod[file] * 4. / 3. * 3.1416 * 1e6; # kg/kg
        # cloudiness mask - as in RICO paper
        cloudy_mask = np.where(rl[file] > 1e-5, 1, 0)
        cloudy_mask_used = cloudy_mask

        # th and rv
        th[file] = h5py.File(filename, "r")["th"][:,:];
        rv[file] = h5py.File(filename, "r")["rv"][:,:];
        # T
        Vexner= np.vectorize(lcmn.exner)
        T = th[file] * Vexner(p_e[file].astype(float))
        # RH
        Vr_vs = np.vectorize(lcmn.r_vs)
        r_vs = Vr_vs(T, p_e[file].astype(float))
        RH = rv[file] / r_vs[file]
        # cloud base
        clb_idx = np.argmax(cloudy_mask_used > 0, axis=1)
        hght_2 = hght[clb_idx]
        hght_2[hght_2==0] = np.nan
        min_hght = np.nanmin(hght_2)
        max_hght = np.nanmax(hght_2)
        # print("wysokosc min",min_hght, " wysokosc max", max_hght)
        # print(i, rr[file][i])

        for j in np.arange(nx):
          if clb_idx[j] > 0:
            sum_rr[file][j] = np.sum(rr[file][j,:int(min_hght/dz[file])-1],0)
        sum_rr[file][sum_rr[file]==0] = np.nan
        # print(file)
    srednie_rr = np.nanmean(sum_rr, axis=0)
    STD_rr = np.nanstd(sum_rr, axis=0)
    plt.errorbar(bin, srednie_rr, STD_rr, fmt='o', label='average')
    plt.ylim((0,0.04))
    plt.xlim((0,120))
    plt.legend()
    plt.savefig(outfile + 'Rain_many_' + str(i*240) +'.png')
    # Od komentuj!
    plt.clf()

# bin_size = bin[2]-bin[1]
# print(bin_size)

for i in range(1, 91):
    Adia_fraction(i)
