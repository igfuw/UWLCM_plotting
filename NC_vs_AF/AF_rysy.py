# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated based on mean th and rv at cloud base cells

from sys import argv, path, maxsize
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages") #/home/pzmij/biblioteki/local_folder/seed/lib/python3/dist-packages")

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# print whole np arrays
np.set_printoptions(threshold=maxsize)

plt.rcParams.update({'font.size': 10})
#plt.figure(figsize=(16,10))

evap_lat = 2.501e6 # [J/kg] latent heat of evaporation

#timesteps = [12960, 13200, 14400]
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i*240
timesteps = timesteps[1:90]

input_dir = argv[1]
outfile = argv[2]

rhod = h5py.File(input_dir + "/const.h5", "r")["G"][:,:,:]
p_e = h5py.File(input_dir + "/const.h5", "r")["p_e"][:]
nx, ny, nz = rhod.shape
dz = h5py.File(input_dir + "/const.h5", "r").attrs["dz"]
hght = np.arange(nz) * dz
bin = np.linspace(0,1.4,len(np.arange(nz)))
bin_size = bin[2]-bin[1]
#print(bin_size)

# ---- adiabatic LWC ----
AF_min = np.zeros([nx, ny, nz])
adia_rl = np.zeros([nz])
adia_rl_min = np.zeros([nz])
clb_rv_min = np.zeros([nx, ny])
clb_th_min = np.zeros([nx, ny])
Biny = [0 for i in np.arange(nz)]
Biny_x = [0 for i in np.arange(nz)]
Slice = np.zeros([nx,ny])

for timestep in timesteps:
  plt.clf()
  filename = input_dir + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
  rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  rl_base = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  nc = h5py.File(filename, "r")["cloud_rw_mom0"][:,:,:] * rhod / 1e6; # 1 / cm^3

  # cloudiness mask - as in RICO paper
  cloudy_mask = np.where(rl > 1e-5, 1, 0)
  cloudy_mask_used = cloudy_mask

  # print('rl>1e-5 cloudy cells: ', np.sum(cloudy_mask_used))
  # print('rl>1e-5 mean nc in cloudy cells: ', np.sum(nc * cloudy_mask_used) / np.sum(cloudy_mask_used))

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
  clb_idx = np.argmax(cloudy_mask_used > 0, axis=2)
  hght_2 = hght[clb_idx]
  hght_2[hght_2==0] = np.nan
  min_hght = np.nanmin(hght_2)
  max_hght = np.nanmax(hght_2)
  # print("wysokosc min",min_hght, " wysokosc max", max_hght)

  for i in np.arange(nx):
    for j in np.arange(ny):
      if clb_idx[i,j] > 0:
        clb_rv_min[i,j] = rv[i, j, int(min_hght/dz)]
        clb_th_min[i,j] = th[i, j, int(min_hght/dz)]

  # model a parcel to get an adiabatic rl, assume a single parcel moving starting from mean rv and th at cloud base

  parcel_rv_min = np.mean(clb_rv_min[clb_rv_min>0])
  parcel_th_min = np.mean(clb_th_min[clb_th_min>0])
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

  for i in np.arange(nx):
    for j in np.arange(ny):
      for k in np.arange(nz):
         if adia_rl_min[k] == 0:
           AF_min[i, j, k] = 0
         else:
           AF_min[i, j, k] = rl[i,j,k] / adia_rl_min[k]

######Biny

  AF = AF_min * cloudy_mask_used
  for k in np.arange(nz):
      for i in np.arange(nx):
          for j in np.arange(ny):
              Slice[i, j] = AF[i, j, k]
          Slice[i,:][Slice[i,:]==0] = np.nan
          Biny[i] = np.digitize(Slice[i,:],bin)
          Biny[i] = np.bincount(Biny[i])
          Biny[i] = np.pad(Biny[i], (0, (len(bin)+1)-len(Biny[i])), mode='constant')
      Biny = list(map(list,Biny))
      Biny_x[k] = np.log10(np.sum(Biny,0)/bin_size)

  #print(Biny_x)

  AF_min[AF_min==0] = np.nan
  AF_min[cloudy_mask_used==0] = np.nan
  AF_mean_min = np.nanmean(AF_min, axis=(0,1))
  AF_mean_min[np.isnan(AF_mean_min)] = 0
  bins = np.array(Biny_x)
  bins = bins[:,:-1]

  fig,ax = plt.subplots(figsize=(11, 8))
  axins = inset_axes(ax,
             width="60%", # width = 10% of parent_bbox width
             height="5%", # height : 50%
             loc=2,
             bbox_to_anchor=(0.3, 0., 1, 1),
             bbox_transform=ax.transAxes,
             borderpad=0,
             )
  e = ax.contourf(bin, hght, np.round(bins,2), 200, cmap='gnuplot') #bin[1:]
  cbar = plt.colorbar(e, cax=axins,  orientation='horizontal', label=r"$log_{10}$($\frac{m^{3}}{\frac{unit AF}{m}}$)")
  ax.set_ylabel('Height [m]')
  ax.set_xlabel('AF []')
  ax2=ax.twinx()
  ax2.plot(AF_mean_min , hght, label="AF", c='k', linewidth=2)
  ax2.set_xlim((0.01,1.4))
  ax2.set_ylim(0,10000)
  ax2.set_yticks([], [])
  ax2.set_yticks([], minor=True)
  plt.savefig(outfile + '_NCvsAF_' + str(int(timestep)) +'.png')
