from sys import argv, path, maxsize
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
# path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")


'''
HOW TO RUN

python3 RAIN_rysy_xz2d_2D_lvl0.py /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/AF/
UPDATE

python3 /home/pzmij/biblioteki/UWLCM_plotting/Rain_distribution/rain_histogram.py /home/pzmij/2D/PAPER/no_Piggy/SD100/classic/ /home/pzmij/2D/PAPER/Wyniki/New_hist/
Singularity> python3 /home/pzmij/biblioteki/UWLCM_plotting/Rain_distribution/rain_histogram.py /home/pzmij/2D/PAPER/Piggy/SD100/classic/ /home/pzmij/2D/PAPER/Wyniki/New_hist/Piggy/SD100/classic/

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
plt.rcParams.update({'font.size': 20})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i*240
# timesteps = timesteps[1:91]
timesteps = timesteps[1:91]
paths = argv[1]
outfile = argv[2]

file_names = os.listdir(paths)

pliki = [ i for i in range(len(file_names))]
bin = [1e-14, 1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2 ]
for timestep in timesteps:
  plt.clf()
  Hist = []
  rr_files = np.zeros((len(file_names), len(bin)-1))
  # rr_files = np.zeros((len(file_names), 121, 101))
  # rr_files = []
  # rr_rows_lvl0 = np.zeros((121, len(file_names)))
  for p in range(len(file_names)):
    path = paths+file_names[p]+'/'+file_names[p]+'_out_lgrngn'
    rhod = h5py.File(path + "/const.h5", "r")["G"][:,:]
    p_e = h5py.File(path + "/const.h5", "r")["p_e"][:]
    nx, nz = rhod.shape
    rhod_cros = rhod[:,:]
    dz = h5py.File(path + "/const.h5", "r").attrs["dz"]
    dx = h5py.File(path + "/const.h5", "r").attrs["dx"]

    X = np.arange(nx)# * dx

    filename = path + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
    rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
    rr = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
    rr_c = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
    nc = h5py.File(filename, "r")["cloud_rw_mom0"][:,:] * rhod_cros / 1e6; # 1 / cm^3

    # cloudiness mask - as in RICO paper
    cloudy_mask = np.where(rl > 1e-5, 1, 0)
    cloudy_mask_used = cloudy_mask

    RR = rr * cloudy_mask_used

    hist = np.histogram(RR, bins=bin)
    new_bin = hist[1]
    rr_files[p,:] = hist[0]
  Hist = np.mean(rr_files, axis=0)

  plt.figure(0)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  plt.step(new_bin[1:], Hist)
  plt.xscale('log')

  plt.suptitle('Current time {}s '.format(int(timestep/2)))
  plt.xlabel('q_r [kg/kg]')
  plt.ylabel('# of cells')
  plt.savefig(outfile+'New_hist'   + str(timestep) +'.png')
