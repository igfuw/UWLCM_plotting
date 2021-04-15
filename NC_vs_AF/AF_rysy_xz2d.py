from sys import argv, path, maxsize
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")#/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")

import h5py
import numpy as np
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# print whole np arrays
np.set_printoptions(threshold=maxsize)
plt.rcParams.update({'font.size': 11})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i*240
timesteps = timesteps[1:91]

input_dir = argv[1]
outfile = argv[2]

rhod = h5py.File(input_dir + "/const.h5", "r")["G"][:,:,:]
p_e = h5py.File(input_dir + "/const.h5", "r")["p_e"][:]
nx, ny, nz = rhod.shape
ny_cross = ny//2
rhod_cros = rhod[:,ny_cross,:]
#print(rhod_cros.shape)
dz = h5py.File(input_dir + "/const.h5", "r").attrs["dz"]
dx = h5py.File(input_dir + "/const.h5", "r").attrs["dx"]

hght = np.arange(nz) * dz
X = np.arange(nx) * dx

# ---- adiabatic LWC ----
AF_min = np.zeros([nx, nz])
adia_rl = np.zeros(nz)
adia_rl_min = np.zeros(nz)
clb_rv_min = np.zeros(nz)
clb_th_min = np.zeros(nz)

for timestep in timesteps:
  plt.clf()
  filename = input_dir + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
  rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,ny_cross,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
  nc = h5py.File(filename, "r")["cloud_rw_mom0"][:,ny_cross,:] * rhod_cros / 1e6; # 1 / cm^3

  # cloudiness mask - as in RICO paper
  cloudy_mask = np.where(rl > 1e-5, 1, 0)
  cloudy_mask_used = cloudy_mask

  # th and rv
  th = h5py.File(filename, "r")["th"][:,ny_cross,:];
  rv = h5py.File(filename, "r")["rv"][:,ny_cross,:];
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

  for i in np.arange(nz):
    if clb_idx[i] > 0:
        clb_rv_min[i] = rv[i, int(min_hght/dz)]
        clb_th_min[i] = th[i, int(min_hght/dz)]

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
      for k in np.arange(nz):
         if adia_rl_min[k] == 0:
           AF_min[i, k] = 0
         else:
           AF_min[i, k] = rl[i, k] / adia_rl_min[k]

  AF_min[AF_min==0] = np.nan
  AF_min[cloudy_mask_used==0] = np.nan
  AF_min[AF_min>1.4] = np.nan
  AF_min = np.transpose(AF_min)

  e = plt.contourf(X, hght, AF_min, 10, cmap='gnuplot', vmin=0, vmax=1.4, levels=[0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]) #bin[1:]
  plt.colorbar(e,  orientation='vertical', label=r"AF[]", ticks=[0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4])
  plt.clim(0, 1.4)
  plt.ylabel('Height [m]')
  plt.xlabel('X [m]')
  plt.title('Adiabatic Fraction [](Y=5000 m, {}s )'.format(int(timestep/2)))
  plt.savefig(outfile + '_AF_xz2d_' + str(int(timestep/2)) +'.png')
