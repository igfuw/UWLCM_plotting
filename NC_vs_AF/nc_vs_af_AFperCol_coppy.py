# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated separately for each column

from sys import argv, path, maxsize
#path.insert(0,"/home/piotr/usr/local/lib/python3/dist-packages/")
path.insert(0,"../../local_folder/uptodate/lib/python3/dist-packages")
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn

# print whole np arrays
np.set_printoptions(threshold=maxsize)

plt.rcParams.update({'font.size': 20})
plt.figure(figsize=(16,10))

evap_lat = 2.5e6 # [J/kg] latent heat of evaporation

#timesteps = [12000, 36000, 72000]
#timesteps = [13200]

timesteps = np.ones(91)

for i in range(1, 91):
       timesteps[i] = i*240
timesteps = timesteps[1:91]

input_dir = argv[1]
outfile = argv[2]

rhod = h5py.File(input_dir + "/const.h5", "r")["G"][:,:,:]
p_e = h5py.File(input_dir + "/const.h5", "r")["p_e"][:]
nx, ny, nz = rhod.shape
dz = h5py.File(input_dir + "/const.h5", "r").attrs["dz"]

for timestep in timesteps:
  plt.clf()
  filename = input_dir + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
  #rl = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:] + h5py.File(filename, "r")["rain_rw_mom3"][:,:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  rl = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  nc = h5py.File(filename, "r")["cloud_rw_mom0"][:,:,:] * rhod / 1e6; # 1 / cm^3

  # cloudiness mask - as in DYCOMS paper
  cloudy_mask = np.where(nc > 20, 1, 0)

  print('nc>20 cloudy cells: ', np.sum(cloudy_mask))
  print('nc>20 mean nc in cloudy cells: ', np.sum(nc * cloudy_mask) / np.sum(cloudy_mask))

  # cloudiness mask - rl > 1e-4
  cloudy_mask = np.where(rl > 1e-4, 1, 0)

  print('rl>1e-4 cloudy cells: ', np.sum(cloudy_mask))
  print('rl>1e-4 mean nc in cloudy cells: ', np.sum(nc * cloudy_mask) / np.sum(cloudy_mask))

  # cloudiness mask - as in RICO paper
  cloudy_mask = np.where(rl > 1e-5, 1, 0)
  cloudy_mask_used = cloudy_mask

  print('rl>1e-5 cloudy cells: ', np.sum(cloudy_mask))
  print('rl>1e-5 mean nc in cloudy cells: ', np.sum(nc * cloudy_mask) / np.sum(cloudy_mask))

  hght_abv_clb = np.zeros([nx, ny, nz])
  hght = np.arange(nz) * dz

  # ---- adiabatic LWC ----
  AF = np.zeros([nx, ny, nz])
  adia_rl = np.zeros([nx, ny, nz])
  # th and rv
  th = h5py.File(filename, "r")["th"][:,:,:];
  rv = h5py.File(filename, "r")["rv"][:,:,:];

  # T
  Vexner = np.vectorize(lcmn.exner)
  T = th * Vexner(p_e.astype(float))

  # RH
  Vr_vs = np.vectorize(lcmn.r_vs)
  r_vs = Vr_vs(T, p_e.astype(float))
  RH = rv / r_vs

  # cloud base
  #clb_idx = np.argmax(RH > 1, axis=2)
  clb_idx = np.argmax(cloudy_mask_used>0, axis=2)

  # clb condition per column
  clb_rv = np.zeros([nx, ny])
  clb_th = np.zeros([nx, ny])


  # model a parcel to get an adiabatic rl
  for i in np.arange(nx):
    for j in np.arange(ny):
      parcel_rl = 0
      parcel_th = 0
      parcel_rv = 0
      if clb_idx[i,j] > 0:
        parcel_rv[i,j] = rv[i, j, clb_idx[i, j]]
        parcel_th[i,j] = th[i, j, clb_idx[i, j]]
      for k in np.arange(nz):
        parcel_T = parcel_th * lcmn.exner(p_e.astype(float)[k])
        delta_rv = parcel_rv - lcmn.r_vs(parcel_T, p_e.astype(float)[k])
        if delta_rv <= 0:
          delta_rv = 0
        parcel_rv -= delta_rv
        parcel_th += delta_rv * evap_lat / lcmn.c_pd / lcmn.exner(p_e.astype(float)[k])
        parcel_rl += delta_rv
        adia_rl[i, j, k] = parcel_rl


  for i in np.arange(nx):
    for j in np.arange(ny):
      for k in np.arange(nz):
        if adia_rl[i, j, k] == 0:
          AF[i, j, k] = 0
        else:
          AF[i, j, k] = rl[i,j,k] / adia_rl[i, j, k]
        hght_abv_clb[i, j, k] = (k - clb_idx[i,j]) * dz
  
  # set cloudy_mask=0 below cloud base and in non-cloudy columns
  # not needed? after cloud base detection based on cloudy mask and not RH?
  for i in np.arange(nx):
    for j in np.arange(ny):
      for k in np.arange(nz):
        if k < clb_idx[i,j] or clb_idx[i,j]==0:
          cloudy_mask_used[i,j,k] = 0
  # plot AF profile
  plt.clf()
  AF[AF==0] = np.nan
  AF[cloudy_mask_used==0] = np.nan
  AF_mean = np.nanmean(AF, axis=(0,1))
  AF_mean[np.isnan(AF_mean)] = 0 # set AF=0 above and below cloud
#  plt.xlim(0,1.5)
  plt.plot(AF_mean, hght)
  plt.savefig(outfile + '_AFprofile_AFperCol_' + str(int(timestep)) +'.png')
