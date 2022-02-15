from sys import argv, path, maxsize
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
# path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")




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


outfile = '/home/pzmij/2D/PAPER/Wyniki/Distribution/'
name = 'New_master_SD50_compare_1000_'
path_to_file_1 = '/home/pzmij/2D/PAPER/Distribution/Piggy_SD50_2/SD1000/tail/'
path_to_file_2 = '/home/pzmij/2D/PAPER/Distribution/Piggy_SD50_21/SD1000/tail/'
path_to_file_3 = '/home/pzmij/2D/PAPER/Distribution/Piggy_SD50_27/SD1000/tail/'
path_to_file_4 = '/home/pzmij/2D/PAPER/Distribution/Piggy_SD50_72/SD1000/tail/'



def Precipitation_rate(paths, timestep):
    files = os.listdir(paths)
    nr_files = len(files)
    rhod= [0 for i in range(nr_files)]
    precip_rate= [0 for i in range(nr_files)]
    precip_rate_STD= [0 for i in range(nr_files)]
    p = 0
    for file in range(nr_files):
        filename = paths+files[file] + '/'+ files[file]+ "_out_lgrngn"
        rhod[p] = h5py.File(filename + "/const.h5", "r")["G"][:,:]
        precip_rate[p] = (h5py.File(filename + "/timestep" + str(240*timestep).zfill(10) + ".h5", "r")["precip_rate"][:,:])
        p+=1
    rhod = np.mean(rhod,axis=0)
    precip_rate = np.mean(precip_rate,axis=0)
    precip_rate_STD = np.std(precip_rate,axis=0)

    nx, nz = rhod.shape
    hght = np.arange(nz) * 100
    return precip_rate, precip_rate_STD, hght



punkty = np.intc(np.linspace(0,90,91))
partial_fun_1 = functools.partial(Precipitation_rate, path_to_file_1)
partial_fun_2 = functools.partial(Precipitation_rate, path_to_file_2)
partial_fun_3 = functools.partial(Precipitation_rate, path_to_file_3)
partial_fun_4 = functools.partial(Precipitation_rate, path_to_file_4)
with concurrent.futures.ProcessPoolExecutor() as executor:
    Classic_1 = executor.map(partial_fun_1, punkty)
    Classic_2 = executor.map(partial_fun_2, punkty)
    Classic_3 = executor.map(partial_fun_3, punkty)
    Classic_4 = executor.map(partial_fun_4, punkty)


Class_1 = [i for i in Classic_1]
Class_2 = [i for i in Classic_2]
Class_3 = [i for i in Classic_3]
Class_4 = [i for i in Classic_4]

Class_1_array = np.array(Class_1)
Class_2_array = np.array(Class_2)
Class_3_array = np.array(Class_3)
Class_4_array = np.array(Class_4)


mean_p_rate1 = [0 for i in range(91)]
std_p_rate1 = [0 for i in range(91)]

mean_p_rate2 = [0 for i in range(91)]
std_p_rate2 = [0 for i in range(91)]

mean_p_rate3 = [0 for i in range(91)]
std_p_rate3 = [0 for i in range(91)]

mean_p_rate4 = [0 for i in range(91)]
std_p_rate4 = [0 for i in range(91)]


for i in range(91):
   
    mean_p_rate1[i] = float(0) if np.all((Class_1_array[i])[0]) == np.nan else np.mean((Class_1_array[i])[0], axis=0)
    std_p_rate1[i] = float(0) if np.all((Class_1_array[i])[1]) == np.nan else np.mean((Class_1_array[i])[1], axis=0)
    
    mean_p_rate2[i] = float(0) if np.all((Class_2_array[i])[0]) == np.nan else np.mean((Class_2_array[i])[0], axis=0)
    std_p_rate2[i] = float(0) if np.all((Class_2_array[i])[1]) == np.nan else np.mean((Class_2_array[i])[1], axis=0)
   
    mean_p_rate3[i] = float(0) if np.all((Class_3_array[i])[0]) == np.nan else np.mean((Class_3_array[i])[0], axis=0)
    std_p_rate3[i] = float(0) if np.all((Class_3_array[i])[1]) == np.nan else np.mean((Class_3_array[i])[1], axis=0)

    mean_p_rate4[i] = float(0) if np.all((Class_4_array[i])[0]) == np.nan else np.mean((Class_4_array[i])[0], axis=0)
    std_p_rate4[i] = float(0) if np.all((Class_4_array[i])[1]) == np.nan else np.mean((Class_4_array[i])[1], axis=0)

height = ((Class_1_array[i])[2])

Mean_precip_rate_1 =  np.mean(mean_p_rate1, axis=0)
STD_precip_rate_1 = np.mean(std_p_rate1, axis=0)
Mean_precip_rate_2 =  np.mean(mean_p_rate2, axis=0)
STD_precip_rate_2 = np.mean(std_p_rate2, axis=0)
Mean_precip_rate_3 =  np.mean(mean_p_rate3, axis=0)
STD_precip_rate_3 = np.mean(std_p_rate3, axis=0)
Mean_precip_rate_4 =  np.mean(mean_p_rate4, axis=0)
STD_precip_rate_4 = np.mean(std_p_rate4, axis=0)
plt.figure(figsize=(30,15))
plt.plot(Mean_precip_rate_1,height, label = 'M_2', linewidth=5)
plt.plot(Mean_precip_rate_2,height, label = 'M_21', linewidth=5)
plt.plot(Mean_precip_rate_3, height, label = 'M_27', linewidth=5)
plt.plot(Mean_precip_rate_4, height, label = 'M_72', linewidth=5)

plt.legend()
plt.title(name)
# plt.fill_betweenx(height, mean_p_rate-std_p_rate, mean_p_rate+std_p_rate, alpha=0.2)
plt.savefig(outfile+name+"_mean_"+str(91*240).zfill(10)+'.png')
plt.clf()

