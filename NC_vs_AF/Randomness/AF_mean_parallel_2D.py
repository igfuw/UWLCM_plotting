# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated based on mean th and rv at cloud base cells

from sys import argv, path, maxsize
#path.insert(0,"../../local_folder/uptodate/lib/python3/dist-packages")
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
# path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")
# path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")



'''
OK
HOW TO RUN:

python3 AF_mean_parallel_2D.py /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/

'''

import h5py
import numpy as np
import matplotlib
import time
import multiprocessing
import concurrent.futures
import threading
import concurrent.futures
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
import glob, os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

start = time.perf_counter()
plt.rcParams.update({'font.size': 12})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i
timesteps = timesteps[1:91]
np.set_printoptions(threshold=maxsize)

n = len(argv)
paths = argv[1]
outfile = argv[2]

const = os.listdir(paths)
# nr_files = len(files)
# const = [argv[i] for i in range(2,n)]

files = [0 for i in range(len(paths))]
series_file = [0 for i in range(len(paths))]

print(const)


def Adia_fraction(i):


    rhod = [0 for i in range(len(const))]
    p_e = [0 for i in range(len(const))]
    dz = [0 for i in range(len(const))]
    rl = [0 for i in range(len(const))]
    rl_base = [0 for i in range(len(const))]
    nc = [0 for i in range(len(const))]
    th = [0 for i in range(len(const))]
    rv = [0 for i in range(len(const))]

    for p in range(len(const)):
        rhod[p] = h5py.File(paths + const[p]+"/const.h5", "r")["G"][:,:]
        p_e[p] = h5py.File(paths + const[p]+"/const.h5", "r")["p_e"][:]
        dz[p] = h5py.File(paths + const[p]+"/const.h5", "r").attrs["dz"]
        rl[p] = (h5py.File(paths + const[p] + "/timestep" + str(240*i).zfill(10) + ".h5", "r")["actrw_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        rl_base[p] = (h5py.File(paths + const[p] + "/timestep" + str(240*i).zfill(10) + ".h5", "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        nc[p] = h5py.File(paths + const[p] + "/timestep" + str(240*i).zfill(10) + ".h5", "r")["cloud_rw_mom0"][:,:]
        th[p] = h5py.File(paths + const[p] + "/timestep" + str(240*i).zfill(10) + ".h5", "r")["th"][:,:];
        rv[p] = h5py.File(paths + const[p] + "/timestep" + str(240*i).zfill(10) + ".h5", "r")["rv"][:,:];
    rhod = np.mean(rhod,axis=0)
    p_e = np.mean(p_e,axis=0)
    dz = np.mean(dz,axis=0)
    rl = np.mean(rl,axis=0)
    rl_base = np.mean(rl_base,axis=0)
    nc = np.mean(nc,axis=0) * rhod / 1e6# 1 / cm^3
    th = np.mean(th,axis=0)
    rv = np.mean(rv,axis=0)

    nx, nz = rhod.shape
    hght = np.arange(nz) * dz
    bin = np.linspace(0,1.4,len(np.arange(nz)))
    bin_size = bin[2]-bin[1]


    # plt.clf()
    # ---- adiabatic LWC ----
    AF_min = np.zeros([nx, nz])
    adia_rl = np.zeros([nz])
    adia_rl_min = np.zeros([nz])
    clb_rv_min = np.zeros([nx])
    clb_th_min = np.zeros([nx])
    Biny = [0 for i in np.arange(nz)]
    Biny_x = [0 for i in np.arange(nz)]
    Slice = np.zeros([nx, nz])
    cloudy_mask = np.where(rl > 1e-5, 1, 0)
    cloudy_mask_used = cloudy_mask
    # print(cloudy_mask_used)
    # T
    Vexner = np.vectorize(lcmn.exner)
    T = th * Vexner(p_e.astype(float))
    # RH
    Vr_vs = np.vectorize(lcmn.r_vs)
    r_vs = Vr_vs(T, p_e.astype(float))
    RH = rv / r_vs
    # cloud base
    clb_idx = np.argmax(cloudy_mask_used > 0, axis=1)
    hght_2 = hght[clb_idx]
    hght_2[hght_2==0] = np.nan
    min_hght = np.nanmin(hght_2)
    max_hght = np.nanmax(hght_2)
    print(min_hght)
    if min_hght/dz < 10:
        min_hght = 10*dz
        d+=1

    for j in np.arange(nx):
        if clb_idx[j] > 0:
          clb_rv_min[j] = rv[j, int(min_hght/dz)]
          clb_th_min[j] = th[j, int(min_hght/dz)]

    parcel_rv_min = np.nan_to_num(np.nanmean(clb_rv_min[clb_rv_min>0]))
    parcel_th_min = np.nan_to_num(np.nanmean(clb_th_min[clb_th_min>0]))
    parcel_rl_min = 0

    for k in np.arange(nz):
      parcel_T_min = parcel_th_min * lcmn.exner(p_e.astype(float)[k])
      delta_rv_min = parcel_rv_min - lcmn.r_vs(parcel_T_min, p_e.astype(float)[k])
      if delta_rv_min <= 0:
        delta_rv_min = 0
      parcel_rv_min -= delta_rv_min
      parcel_th_min += delta_rv_min * evap_lat / lcmn.c_pd / lcmn.exner(p_e.astype(float)[k])
      parcel_rl_min += delta_rv_min
      adia_rl_min[k] = parcel_rl_min
    # print(adia_rl_min)

    for j in np.arange(nx):
        for k in np.arange(nz):
           if adia_rl_min[k] == 0:
             AF_min[j, k] = 0
           else:
             AF_min[j, k] = rl[j,k] / adia_rl_min[k]
        # print(AF_min)

  ######Biny

    AF = AF_min * cloudy_mask_used
    for k in np.arange(nz):
        for j in np.arange(nx):
            Slice[j, k] = AF[j, k]
        Slice[:,k][Slice[:,k]==0] = np.nan
        Biny[k] = np.digitize(Slice[:,k],bin)
        Biny[k] = np.bincount(Biny[k])
        Biny[k] = np.pad(Biny[k], (0, (len(bin)+1)-len(Biny[k])), mode='constant')
        Biny_x[k] = np.log10(Biny[k]/bin_size)

    AF_min[AF_min==0] = np.nan
    AF_min[cloudy_mask_used==0] = np.nan
    AF_mean_min = np.nanmean(AF_min, axis=0)
    bins = np.array(Biny_x)
    bins = bins[:,:-1]
    # print(AF_mean_min)

    fig,ax = plt.subplots(figsize=(11, 8))
    axins = inset_axes(ax,
               width="50%", # width = 10% of parent_bbox width
               height="5%", # height : 50%
               loc=2,
               bbox_to_anchor=(0.4, 0., 1, 1),
               bbox_transform=ax.transAxes,
               borderpad=0,
               )
    e = ax.contourf(bin, hght, bins, 200, cmap='gnuplot') #bin[1:]
    cbar = plt.colorbar(e, cax=axins,  orientation='horizontal', label=r"$log_{10}$($\frac{m^{3}}{\frac{unit AF}{m}}$)", format='%.2f')
    ax.set_ylabel('Height [m]')
    ax.set_xlabel('AF []')
    ax2=ax.twinx()
    ax2.plot(AF_mean_min , hght, label="AF", c='k', linewidth=2)
    ax2.set_xlim((0.01,1.4))
    ax2.set_ylim(0,10000)
    ax2.set_yticks([], [])
    ax2.set_yticks([], minor=True)
    plt.title('time = '+str(i*240)+ 's')
    plt.savefig(outfile + '/AF_average/' + str(240*i).zfill(10) +'.png')


punkty = np.intc(np.linspace(1,91,91))
with concurrent.futures.ProcessPoolExecutor() as executor:
    results = executor.map(Adia_fraction, punkty)

# for i in range(1,91):
# Adia_fraction(89)
finish = time.perf_counter()
print(f'Finished in {round(finish-start,2)} seconds(s)')
