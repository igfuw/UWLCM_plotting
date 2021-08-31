# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated based on mean th and rv at cloud base cells

from sys import argv, path, maxsize
#path.insert(0,"../../local_folder/uptodate/lib/python3/dist-packages")
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")


import h5py
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

start = time.perf_counter()
# print whole np arrays
np.set_printoptions(threshold=maxsize)

plt.rcParams.update({'font.size': 12})
#plt.figure(figsize=(16,10))

evap_lat = 2.501e6 # [J/kg] latent heat of evaporation

timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i*240
timesteps = timesteps[1:91]
# timesteps = timesteps[88:91]

input_dir = argv[1]
outfile = argv[2]

rhod = h5py.File(input_dir + "/const.h5", "r")["G"][:,:]
p_e = h5py.File(input_dir + "/const.h5", "r")["p_e"][:]
nx, nz = rhod.shape
dz = h5py.File(input_dir + "/const.h5", "r").attrs["dz"]
hght = np.arange(nz) * dz
bin = np.linspace(1,121,len(np.arange(nx)))
# bin_size = bin[2]-bin[1]
# print(bin_size)


# ---- adiabatic LWC ----
AF_min = np.zeros([nx, nz])
adia_rl = np.zeros([nz])
adia_rl_min = np.zeros([nz])
clb_ql_min = np.zeros([nx])
sum_rr = np.zeros([nx])
Biny = [0 for i in np.arange(nz)]
Biny_x = [0 for i in np.arange(nz)]
Slice = np.zeros([nx, nz])


for timestep in timesteps:
  plt.clf()
  filename = input_dir + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
  rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  rl_base = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  nc = h5py.File(filename, "r")["cloud_rw_mom0"][:,:] * rhod / 1e6; # 1 / cm^3
  rr = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e6; # kg/kg

  rr = rr * rhod
  # cloudiness mask - as in RICO paper
  cloudy_mask = np.where(rl > 1e-5, 1, 0)
  cloudy_mask_used = cloudy_mask


  # print(timestep, rr)
  # print('rl>1e-5 cloudy cells: ', np.sum(cloudy_mask_used))
  # print('rl>1e-5 mean nc in cloudy cells: ', np.sum(nc * cloudy_mask_used) / np.sum(cloudy_mask_used))

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
  # print("wysokosc min",min_hght, " wysokosc max", max_hght)


  for i in np.arange(nx):
      if clb_idx[i] > 0:
      # print(timestep, np.sum(ql[i,:int(min_hght/dz)-1]),i)
        sum_rr[i] = np.sum(rr[i,:int(min_hght/dz)-1],0)


  # plt.bar(bin, sum_rr)
  plt.scatter(bin, hght_2)
  # plt.ylim((0,0.02))
  plt.title("histogram")
  plt.savefig(outfile + 'Rain_singel' + str(int(timestep)) +'.png')
  #     if clb_idx[i] > 0:
  #       clb_ql_min[i] = ql[i, int(min_hght/dz)-1]
  # for k in np.arange(nx):
  #     if clb_idx[i] > 0:
  #         sum_ql[i] = np.sum(ql[i,:int(min_hght/dz)-1],0)
