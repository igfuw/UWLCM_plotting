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
import numpy as np
from itertools import product
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from libcloudphxx import common as lcmn
import glob, os
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams.update({'font.size': 20})
np.set_printoptions(threshold=maxsize)
plt.rcParams.update({'font.size': 20})

timesteps = [i*240 for i in range(1, 91)]
def Liczy_histogram(paths):
    file_names = os.listdir(paths)
    bin = [1e-14, 1e-13,1e-12,1e-11,1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2, 1e-1 ]
    Hist = np.zeros((len(timesteps), len(bin)-1))
    for timestep in timesteps:
      rr_files = np.zeros((len(file_names), len(bin)-1))
      for p in range(len(file_names)):
        path = paths+file_names[p]+'/'+file_names[p]+'_out_lgrngn'
        filename = path + "/timestep" + str(int(timestep)).zfill(10) + ".h5"
        rhod = h5py.File(path + "/const.h5", "r")["G"][:,:]
        rr = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * rhod * 4. / 3. * 3.1416 * 1e3; # kg/kg
        hist = np.histogram(rr, bins=bin)
        rr_files[p,:] = hist[0]
      Hist[int(timestep/240)-1,:] = np.mean(rr_files, axis=0)
    Hist_data = np.sum(Hist, axis=0)
    # Hist_data = Hist_data/ np.max(Hist_data)
    return Hist_data, hist[1]

plt.figure(0)
plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))
main_path = str(argv[1]) #'/home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD'
sciezka = argv[4::]#['VF','tail','multi']
outfile = str(argv[2]) #'/home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/'
SD_order = [argv[3]]#[1000]
for path, SD in product(sciezka, SD_order):
  Class, biny = Liczy_histogram(main_path+str(SD)+'/'+path+'/')
  plt.step(biny[1:], Class,  label=path, linewidth=9)
plt.title('SD= '+str(SD)+'[#]')
plt.xscale('log')
plt.xlim((1e-13, 1e-1))
plt.grid(True)
plt.legend()
plt.xlabel(r'$q_r$ [kg/kg]')
plt.ylabel('# of cells')
plt.savefig(outfile+'histogram_SD_'+str(SD)+'.png')
