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
import cProfile, pstats
import h5py
import functools
import numpy as np
from itertools import product
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
import glob, os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import multiprocessing
import concurrent.futures
import threading
import concurrent.futures
plt.rcParams.update({'font.size': 20})


# print whole np arrays
np.set_printoptions(threshold=maxsize)
plt.rcParams.update({'font.size': 20})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation

# paths = argv[1]
# outfile = argv[2]
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i*240
timesteps = timesteps[1:91]

def Liczy_histogram(paths, i):
    file_names = os.listdir(paths)
    
    Hist = np.zeros(( len(biny)-1))
    rr_files = np.zeros((len(file_names), len(biny)-1))
    for p in range(len(file_names)):
      path = paths+file_names[p]+'/'+file_names[p]+'_out_lgrngn'
      filename = path + "/timestep" + str(240*i).zfill(10) + ".h5"
      rhod = h5py.File(path + "/const.h5", "r")["G"][:,:]
      rr = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * rhod * 4. / 3. * 3.1416 * 1e3; # kg/kg
      hist = np.histogram(rr, bins=biny)
      rr_files[p,:] = hist[0]
    Hist = np.mean(rr_files, axis=0)
    # Hist_data = np.sum(Hist, axis=0)
    # Hist_data = Hist_data/ np.max(Hist_data)
    return Hist, hist[1]

# def Run_hist(SD, path, kroki):
    # return Liczy_histogram(main_path+str(SD)+'/'+str(path)+'/', kroki)
def Run_hist(SD, kroki):
    return Liczy_histogram(main_path+str(SD)+'/', kroki)

biny = [1e-14, 1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2, 1e-1, 1 ]
outfile = '/home/pzmij/2D/PAPER/Wyniki/New_hist/Coal/'
main_path = '/home/pzmij/2D/PAPER/Piggy/'
TITLE = 10
TITLE_save = 'testowy_SD_coal'
partial_fun_1000_VF = functools.partial(Run_hist, 'SD10_coal25')#, "VF")
partial_fun_1000_tail = functools.partial(Run_hist, 'SD10_coal50')#, "tail")
partial_fun_1000_multi = functools.partial(Run_hist, 'SD10_coal100')#, "multi")
partial_fun_1000_multi_3 = functools.partial(Run_hist, 'SD10/classic')#, "multi")
punkty = np.intc(np.linspace(1,90,91))
with concurrent.futures.ProcessPoolExecutor() as executor:
    Classic = executor.map(partial_fun_1000_VF, punkty)
    Tails = executor.map(partial_fun_1000_tail, punkty)
    Multis = executor.map(partial_fun_1000_multi, punkty)
    Multis2 = executor.map(partial_fun_1000_multi_3, punkty)
Class = np.sum([i[0] for i in Classic], axis=0)
Tail = np.sum([i[0] for i in Tails], axis=0)
Multi = np.sum([i[0] for i in Multis], axis=0)
Multi2 = np.sum([i[0] for i in Multis2], axis=0)

plt.figure(0)
plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))
plt.step(biny[1:], Multi2, c='k', label="5", linewidth=9)
plt.step(biny[1:], Class, c='y', label="25", linewidth=9)
plt.step(biny[1:], Tail, c='r', label='50', linewidth=9)
plt.step(biny[1:], Multi, c='g', label='100', linewidth=9)
plt.title('SD= '+str(TITLE)+'[#]')
plt.xscale('log')
plt.xlim((1e-13, 1e-0))
plt.grid(True)
plt.legend(title='sstp_coal')
# plt.suptitle('Current time {}s '.format(int(timestep/2)))
plt.xlabel(r'$q_r$ [kg/kg]')
plt.ylabel('# of cells')
plt.savefig(outfile+'hist_'+str(TITLE_save)+'.png')