# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated based on mean th and rv at cloud base cells

from sys import argv, path, maxsize
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
timesteps = np.ones(91)

for i in range(1, 91):
    timesteps[i] = i*240
11
timesteps = timesteps[1:91]
print(timesteps)

input_dir = argv[1]
outfile = argv[2]

rhod = h5py.File(input_dir + "/const.h5", "r")["G"][:,:,:]
p_e = h5py.File(input_dir + "/const.h5", "r")["p_e"][:]
nx, ny, nz = rhod.shape
dz = h5py.File(input_dir + "/const.h5", "r").attrs["dz"]

for timestep in timesteps:
  plt.clf()

  print(int(timestep))
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
  
  print('rl>1e-5 cloudy cells: ', np.sum(cloudy_mask))
  print('rl>1e-5 mean nc in cloudy cells: ', np.sum(nc * cloudy_mask) / np.sum(cloudy_mask))
  
  # ---- adiabatic LWC ----
  AF = np.zeros([nx, ny, nz])
  #adia_rl = np.zeros([nx, ny, nz])
  adia_rl = np.zeros([nz])
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
  clb_idx = np.argmax(RH > 1, axis=2)
  
  # clb condition per column
  clb_rv = np.zeros([nx, ny])
  clb_th = np.zeros([nx, ny])
  
  for i in np.arange(nx):
    for j in np.arange(ny):
      if clb_idx[i,j] > 0:
        clb_rv[i,j] = rv[i, j, clb_idx[i, j]]
        clb_th[i,j] = th[i, j, clb_idx[i, j]]
  
  # model a parcel to get an adiabatic rl, assume a single parcel moving starting from mean rv and th at cloud base
  parcel_rv = np.mean(clb_rv[clb_rv>0])
  parcel_th = np.mean(clb_th[clb_th>0])
  parcel_rl = 0
  
  print('parcel init: th = ', parcel_th, ' rv = ', parcel_rv)
  
  for k in np.arange(nz):
    parcel_T = parcel_th * lcmn.exner(p_e.astype(float)[k])
    delta_rv = parcel_rv - lcmn.r_vs(parcel_T, p_e.astype(float)[k])
    if delta_rv <= 0:
      delta_rv = 0
    parcel_rv -= delta_rv
    parcel_th += delta_rv * evap_lat / lcmn.c_pd / lcmn.exner(p_e.astype(float)[k])
    parcel_rl += delta_rv
    adia_rl[k] = parcel_rl
  
  print(p_e)
  print(adia_rl)
  
  #for i, j in zip(np.arange(nx), np.arange(ny)):
  #  parcel_rv = rv[i, j, clb_idx[i, j]]
  #  parcel_th = th[i, j, clb_idx[i, j]]
  #  for k in np.arange(nz):
  #    if k < clb_idx[i, j]:
  #      adia_rl[i,j,k] = 0
  #    else:
  #      parcel_T = parcel_th * lcmn.exner(p_e.astype(float)[k])
  #      delta_rv = parcel_rv - lcmn.r_vs(parcel_T, p_e.astype(float)[k])
  #      parcel_rv -= delta_rv
  #      parcel_th += delta_rv * evap_lat / lcmn.c_pd / lcmn.exner(p_e.astype(float)[k])
  #      adia_rl[i,j,k] = rl[i,j,k] / delta_rv
  
  
  #adia_rl = np.where(adia_rl > 0., adia_rl, 0)
  #print adia_rl
  
  # translate rl to AF
  #AF = np.where(adia_rl > 0., rl / adia_rl, 0)
  #AF = np.where(adia_rl > 1e-5, rl / adia_rl, 0)
  
  for i in np.arange(nx):
    for j in np.arange(ny):
      for k in np.arange(nz):
  #      if rl[i,j,k] > 0 and adia_rl[k] == 0:
  #        print 'i: ',i, ' j: ',j, ' k: ',k, 'rl: ', rl[i,j,k], 'adia_rl: ', adia_rl[k], 'nc: ', nc[i,j,k]
        if cloudy_mask[i,j,k] > 0:
          AF[i, j, k] = rl[i,j,k] / adia_rl[k]
  #        print 'i: ',i, ' j: ',j, ' k: ',k, 'rl: ', rl[i,j,k], 'adia_rl: ', adia_rl[k], 'nc: ', nc[i,j,k], 'AF: ', AF[i,j,k]
        else:
          AF[i, j, k] = 0
  
  #print cloudy_mask[cloudy_mask>0]
  #print AF[AF>0]
  #print nc[nc>0]
  
  plt.plot((AF * cloudy_mask).flatten(), (nc * cloudy_mask).flatten(), '.', markersize=1)
  plt.xlim(0,2)
  plt.ylim(0,200)
  
  plt.xlabel('AF')
  plt.ylabel('Nc [1/cc]')
  
  #plt.plot((rl*cloudy_mask).flatten(), (nc*cloudy_mask).flatten(), 'o')
  #
  #plt.xscale('log')
  #plt.yscale('log')
  #
  #plt.show()
  plt.savefig(outfile + '_NCvsAF_' + str(timestep) +'.png')
