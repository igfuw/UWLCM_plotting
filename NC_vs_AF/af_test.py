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

#timesteps = [ 8400, 10800, 13200, 17040, 18720, 20400, 21600]
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
  #rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  rl = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  # cloudiness mask - as in RICO paper
  cloudy_mask = np.where(rl > 1e-5, 1, 0)
  
  # ---- adiabatic LWC ----
  AF = np.zeros([nx, ny, nz])
  adia_rl = np.zeros([nx, ny, nz])
  #adia_rl = np.zeros([nz])
  clb_rv = np.zeros([nx, ny])
  clb_th = np.zeros([nx, ny])
  height = np.zeros(nz)
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
  '''
  for i in np.arange(nx):
    for j in np.arange(ny):
      if clb_idx[i,j] > 0:
        clb_rv[i,j] = rv[i, j, clb_idx[i, j]]
        clb_th[i,j] = th[i, j, clb_idx[i, j]]
  '''
  # model a parcel to get an adiabatic rl, assume a single parcel moving starting from mean rv and th at cloud base
  '''
  parcel_rv = np.mean(clb_rv[clb_rv>0])
  parcel_th = np.mean(clb_th[clb_th>0])
  parcel_rl = 0
  '''
  #print('parcel init: th = ', parcel_th, ' rv = ', parcel_rv)
  '''
  for k in np.arange(nz):
    parcel_T = parcel_th * lcmn.exner(p_e.astype(float)[k])
    delta_rv = parcel_rv - lcmn.r_vs(parcel_T, p_e.astype(float)[k])
    vif delta_rv <= 0:
      delta_rv = 0
    parcel_0rv -= delta_rv
    parcel_th += delta_rv * evap_lat / lcmn.c_pd / lcmn.exner(p_e.astype(float)[k])
    parcel_rl += delta_rv
    adia_rl[k] = parcel_rl
  '''
  for i in np.arange(nx):
    for j in np.arange(ny):
      parcel_rv = rv[i, j, clb_idx[i, j]]
      parcel_th = th[i, j, clb_idx[i, j]]
      parcel_rl = 0
      for k in np.arange(nz):
        if k < clb_idx[i, j] or clb_idx[i,j] == 0:
          adia_rl[i,j,k] = 0
        else:
          parcel_T = parcel_th * lcmn.exner(p_e.astype(float)[k])
          delta_rv = parcel_rv - lcmn.r_vs(parcel_T, p_e.astype(float)[k])
          if delta_rv <= 0:
            delta_rv = 0
          parcel_rv -= delta_rv
          parcel_th += delta_rv * evap_lat / lcmn.c_pd / lcmn.exner(p_e.astype(float)[k])
          parcel_rl += delta_rv
          adia_rl[i,j,k] = parcel_rl
  #print(rl[rl>0])
  #print(adia_rl[adia_rl>0])
  
  for i in np.arange(nx): 
    for j in np.arange(ny):
      for k in np.arange(nz):
        height[k] = k *50
        if cloudy_mask[i,j,k] > 0:
          if adia_rl[i, j ,k] > 0:
            AF[i, j, k] = rl[i,j,k] / (adia_rl[i, j, k] *1e3)
          else:
            AF[i, j, k] = 0
        else:
          AF[i, j, k] = 0
 
  #print(AF.shape)
  #print(AF[AF>0].shape)
  AF[AF == 0] = np.nan
  AF_yz = np.nanmean(AF, axis=0)
  AF_z = np.nanmean(AF_yz, axis=0)

  
  #AF = np.average(AF[AF>0], axis=1)
  #height = [height, height, height]
  #plt.plot((AF * cloudy_mask).flatten(), (dz * cloudy_mask).flatten(), '.', markersize=1)
 # plt.plot((AF * cloudy_mask).flatten(), height.flatten(), '.', markersize=1)
  plt.plot(AF_z , height)
  # plt.plot((AF * cloudy_mask).flatten(), height.flatten())
  plt.xlim(0,2)
  plt.ylim(0,10000)
  
  #plt.xlabel('AF')
  #plt.ylabel('Height')
  
  plt.savefig(outfile + '_NCvsAF_' + str(int(timestep)) +'.png')

  
