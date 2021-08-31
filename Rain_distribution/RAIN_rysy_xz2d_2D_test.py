from sys import argv, path, maxsize
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
# path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")


'''
HOW TO RUN

python3 RAIN_rysy_xz2d_2D_test.py /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/AF/

'''

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
import glob, os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# print whole np arrays
np.set_printoptions(threshold=maxsize)
plt.rcParams.update({'font.size': 11})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i*240
timesteps = timesteps[1:91]

paths = argv[1]
outfile = argv[2]

file_names = os.listdir(paths)

pliki = [ i for i in range(len(file_names))]

for timestep in timesteps:
  plt.clf()
  rr_rows = np.zeros((121, len(file_names)))
  for p in range(len(file_names)):
    path = paths+file_names[p]
    rhod = h5py.File(path + "/const.h5", "r")["G"][:,:]
    p_e = h5py.File(path + "/const.h5", "r")["p_e"][:]
    nx, nz = rhod.shape
    rhod_cros = rhod[:,:]
    dz = h5py.File(path + "/const.h5", "r").attrs["dz"]
    dx = h5py.File(path + "/const.h5", "r").attrs["dx"]

    hght = np.arange(nz) * dz
    X = np.arange(nx)# * dx

    filename = path + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
    rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
    rr = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
    rr_c = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
    nc = h5py.File(filename, "r")["cloud_rw_mom0"][:,:] * rhod_cros / 1e6; # 1 / cm^3

    # cloudiness mask - as in RICO paper
    cloudy_mask = np.where(rl > 1e-5, 1, 0)
    cloudy_mask_used = cloudy_mask

    # th and rv
    th = h5py.File(filename, "r")["th"][:,:];
    rv = h5py.File(filename, "r")["rv"][:,:];
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

    if min_hght/dz < 10:
      min_hght = 10*int(dz)

    min_hght = np.nan_to_num(min_hght)

    # print(np.sum(rr_c[:][int(min_hght/dz)-1]>0), min_hght)
    # print(rr_c[:,int(min_hght/dz)-1].shape)
    rr_rows[:,p] = rr[:,int(min_hght/dz)-1]
    Rr_rows = np.transpose(rr_rows)

  e = plt.contourf(X, pliki, Rr_rows, vmin=7e-15, vmax=4e-4, extend='neither',levels=[7e-15, 7e-14, 7e-13, 7e-12, 7e-11, 7e-10, 7e-9, 7e-8, 7e-7, 7e-6, 7e-5, 7e-4], cmap='gnuplot')#, vmin=0, vmax=1.4, levels=[7e-15, 7e-14, 7e-13, 7e-12, 7e-11, 7e-10, 7e-9, 7e-8, 7e-7, 7e-6, 7e-5, 7e-4]) #bin[1:]
  # e = plt.pcolormesh(X, pliki, Rr_rows, vmin=7e-15, vmax=4e-4)
  plt.colorbar(e,  orientation='vertical', label=r"$r_{v}$", ticks=[7e-15, 7e-14, 7e-13, 7e-12, 7e-11, 7e-10, 7e-9, 7e-8, 7e-7, 7e-6, 7e-5, 7e-4])
  #opcjonalnie levele~!!!!!
  # plt.clim(0, 1.4)
  plt.ylabel('Simulation#')
  plt.xlabel('cell#')
  plt.title('Current time {}s '.format(int(timestep/2)))
  plt.axis(aspect='image')
  plt.savefig(outfile + '/Below_CLB_' + str(int(timestep/2)) +'.png')
