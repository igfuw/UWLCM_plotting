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
plt.rcParams.update({'font.size': 20})


# print whole np arrays
np.set_printoptions(threshold=maxsize)
plt.rcParams.update({'font.size': 20})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation

# paths = argv[1]
# outfile = argv[2]

def Liczy_histogram(paths):
    timesteps = np.ones(91)
    for i in range(1, 91):
        timesteps[i] = i*240
    timesteps = timesteps[1:91]
    file_names = os.listdir(paths)
    pliki = [ i for i in range(len(file_names))]
    bin = [1e-14, 1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2, 1e-1 ]
    Hist = np.zeros((len(timesteps), len(bin)-1))
    for timestep in timesteps:
      rr_files = np.zeros((len(file_names), len(bin)-1))
      for p in range(len(file_names)):
        path = paths+file_names[p]+'/'+file_names[p]+'_out_lgrngn'
        filename = path + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
        # rl = (h5py.File(filename, "r")["actrw_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        rr = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg

        # cloudiness mask - as in RICO paper
        # cloudy_mask = np.where(rl > 1e-5, 1, 0)
        # cloudy_mask_used = cloudy_mask

        RR = rr

        hist = np.histogram(RR, bins=bin)
        new_bin = hist[1]
        rr_files[p,:] = hist[0]

      Hist[int(timestep/240)-1,:] = np.mean(rr_files, axis=0)
    Hist_data = np.sum(Hist, axis=0)
    # Hist_data = Hist_data/ np.max(Hist_data)
    return Hist_data, hist[1]


SD_order = [10, 100, 1000]
# for SD in SD_order:
# sciezka = ['/home/pzmij/2D/PAPER/Piggy/SD'+str(SD)+'/classic/','/home/pzmij/2D/PAPER/Piggy/SD'+str(SD)+'/tail/','/home/pzmij/2D/PAPER/Piggy/SD'+str(SD)+'/multi/']
outfile = '/home/pzmij/2D/PAPER/Wyniki/New_hist/Piggy/SD'#+str(SD)+'/comp/'
outfile2 = '/home/pzmij/2D/PAPER/Wyniki/New_hist/Piggy/'
  # print(SD)
  # Classic, biny = Liczy_histogram(sciezka[0])
  # print('tail')
  # tail = Liczy_histogram(sciezka[1])[0]
  # print('multi')
  # multi = Liczy_histogram(sciezka[2])[0]
  #
  # np.save(outfile+'classic_no_mask'+str(SD), Classic)
  # np.save(outfile+'tail_no_mask'+str(SD), tail)
  # np.save(outfile+'multi_no_mask'+str(SD), multi)
  # np.save(outfile+'biny_no_mask'+str(SD), biny)
#SD1
Classic_10 = np.load(outfile+'10/comp/classic_no_mask10.npy')
tail_10 = np.load(outfile+'10/comp/tail_no_mask10.npy')
multi_10 = np.load(outfile+'10/comp/multi_no_mask10.npy')
#SD10
Classic_100 = np.load(outfile+'100/comp/classic_no_mask100.npy')
tail_100 = np.load(outfile+'100/comp/tail_no_mask100.npy')
multi_100 = np.load(outfile+'100/comp/multi_no_mask100.npy')
#SD10
Classic_1000 = np.load(outfile+'1000/comp/classic_no_mask1000.npy')
tail_1000 = np.load(outfile+'1000/comp/tail_no_mask1000.npy')
multi_1000 = np.load(outfile+'1000/comp/multi_no_mask1000.npy')
biny = np.load(outfile+'10/comp/biny_no_mask10.npy')

plt.figure(0)
plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))
plt.step(biny[1:], Classic_10, c='y', label="SD10", linewidth=9)
plt.step(biny[1:], Classic_100, c='g', label='SD100', linewidth=9)
plt.step(biny[1:], Classic_1000, c='r', label='SD1000', linewidth=9)
plt.title('Classic')
plt.xscale('log')
plt.xlim((1e-13, 1e-1))
plt.grid(True)
plt.legend()
  # plt.suptitle('Current time {}s '.format(int(timestep/2)))
plt.xlabel(r'$q_r$ [kg/kg]')
plt.ylabel('# of cells')
plt.savefig(outfile2+'hist_no_mask_piggy_Classic.png')

plt.figure(0)
plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))
plt.step(biny[1:], tail_10, c='y', label="SD10", linewidth=9)
plt.step(biny[1:], tail_100, c='g', label='SD100', linewidth=9)
plt.step(biny[1:], tail_1000, c='r', label='SD1000', linewidth=9)
plt.title('Tail')
plt.xscale('log')
plt.xlim((1e-13, 1e-1))
plt.grid(True)
plt.legend()
  # plt.suptitle('Current time {}s '.format(int(timestep/2)))
plt.xlabel(r'$q_r$ [kg/kg]')
plt.ylabel('# of cells')
plt.savefig(outfile2+'hist_no_mask_piggy_Tail.png')

plt.figure(0)
plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))
plt.step(biny[1:], multi_10, c='y', label="SD10", linewidth=9)
plt.step(biny[1:], multi_100, c='g', label='SD100', linewidth=9)
plt.step(biny[1:], multi_1000, c='r', label='SD1000', linewidth=9)
plt.title('Multi')
plt.xscale('log')
plt.xlim((1e-13, 1e-1))
plt.grid(True)
plt.legend()
  # plt.suptitle('Current time {}s '.format(int(timestep/2)))
plt.xlabel(r'$q_r$ [kg/kg]')
plt.ylabel('# of cells')
plt.savefig(outfile2+'hist_no_mask_piggy_Multi.png')

  # plt.figure(1)
  # plt.rcParams.update({'font.size': 40})
  # plt.figure(figsize=(40,40))
  # plt.step(biny[1:], Classic/np.max(Classic), c='y', label="classic", linewidth=6)
  # plt.step(biny[1:], tail/np.max(tail), c='r', label='tail', linewidth=6)
  # plt.step(biny[1:], multi/np.max(multi), c='g', label='multi', linewidth=6)
  # plt.title('SD= '+str(SD)+'[#]')
  # plt.xscale('log')
  # plt.xlim((1e-13, 1e-1))
  # plt.grid(linestyle=':')
  # plt.legend()
  # # plt.suptitle('Current time {}s '.format(int(timestep/2)))
  # plt.xlabel(r'$q_r$ [kg/kg]')
  # plt.ylabel('# of cells')
  # plt.savefig(outfile+'hist_no_mask_norm_piggy_SD_'+str(SD)+'.png')
