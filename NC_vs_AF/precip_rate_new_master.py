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
name = 'New_master_check_'
path_to_file_1 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n1/SD100/tail/'
path_to_file_2 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n2/SD100/tail/'
path_to_file_3 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n3/SD100/tail/'
path_to_file_4 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n4/SD100/tail/'
path_to_file_5 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n5/SD100/tail/'
path_to_file_6 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n6/SD100/tail/'
path_to_file_7 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n7/SD100/tail/'
path_to_file_8 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n8/SD100/tail/'
path_to_file_9 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n9/SD100/tail/'
path_to_file_10 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n10/SD100/tail/'
path_to_file_11 = '/home/pzmij/2D/PAPER/Distribution/Piggy_n11/SD100/tail/'


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
partial_fun_5 = functools.partial(Precipitation_rate, path_to_file_5)
partial_fun_6 = functools.partial(Precipitation_rate, path_to_file_6)
partial_fun_7 = functools.partial(Precipitation_rate, path_to_file_7)
partial_fun_8 = functools.partial(Precipitation_rate, path_to_file_8)
partial_fun_9 = functools.partial(Precipitation_rate, path_to_file_9)
partial_fun_10 = functools.partial(Precipitation_rate, path_to_file_10)
partial_fun_11 = functools.partial(Precipitation_rate, path_to_file_11)
with concurrent.futures.ProcessPoolExecutor() as executor:
    Classic_1 = executor.map(partial_fun_1, punkty)
    Classic_2 = executor.map(partial_fun_2, punkty)
    Classic_3 = executor.map(partial_fun_3, punkty)
    Classic_4 = executor.map(partial_fun_4, punkty)
    Classic_5 = executor.map(partial_fun_5, punkty)
    Classic_6 = executor.map(partial_fun_6, punkty)
    Classic_7 = executor.map(partial_fun_7, punkty)
    Classic_8 = executor.map(partial_fun_8, punkty)
    Classic_9 = executor.map(partial_fun_9, punkty)
    Classic_10 = executor.map(partial_fun_10, punkty)
    Classic_11 = executor.map(partial_fun_11, punkty)


Class_1 = [i for i in Classic_1]
Class_2 = [i for i in Classic_2]
Class_3 = [i for i in Classic_3]
Class_4 = [i for i in Classic_4]
Class_5 = [i for i in Classic_5]
Class_6 = [i for i in Classic_6]
Class_7 = [i for i in Classic_7]
Class_8 = [i for i in Classic_8]
Class_9 = [i for i in Classic_9]
Class_10 = [i for i in Classic_10]
Class_11 = [i for i in Classic_11]
Class_1_array = np.array(Class_1)
Class_2_array = np.array(Class_2)
Class_3_array = np.array(Class_3)
Class_4_array = np.array(Class_4)
Class_5_array = np.array(Class_5)
Class_6_array = np.array(Class_6)
Class_7_array = np.array(Class_7)
Class_8_array = np.array(Class_8)
Class_9_array = np.array(Class_9)
Class_10_array = np.array(Class_10)
Class_11_array = np.array(Class_11)

mean_p_rate1 = [0 for i in range(91)]
std_p_rate1 = [0 for i in range(91)]

mean_p_rate2 = [0 for i in range(91)]
std_p_rate2 = [0 for i in range(91)]

mean_p_rate3 = [0 for i in range(91)]
std_p_rate3 = [0 for i in range(91)]

mean_p_rate4 = [0 for i in range(91)]
std_p_rate4 = [0 for i in range(91)]

mean_p_rate5 = [0 for i in range(91)]
std_p_rate5 = [0 for i in range(91)]

mean_p_rate6 = [0 for i in range(91)]
std_p_rate6 = [0 for i in range(91)]

mean_p_rate7 = [0 for i in range(91)]
std_p_rate7 = [0 for i in range(91)]

mean_p_rate8 = [0 for i in range(91)]
std_p_rate8 = [0 for i in range(91)]

mean_p_rate9 = [0 for i in range(91)]
std_p_rate9 = [0 for i in range(91)]

mean_p_rate10 = [0 for i in range(91)]
std_p_rate10 = [0 for i in range(91)]

mean_p_rate11 = [0 for i in range(91)]
std_p_rate11 = [0 for i in range(91)]

for i in range(91):
   
    mean_p_rate1[i] = float(0) if np.all((Class_1_array[i])[0]) == np.nan else np.mean((Class_1_array[i])[0], axis=0)
    std_p_rate1[i] = float(0) if np.all((Class_1_array[i])[1]) == np.nan else np.mean((Class_1_array[i])[1], axis=0)
    
    mean_p_rate2[i] = float(0) if np.all((Class_2_array[i])[0]) == np.nan else np.mean((Class_2_array[i])[0], axis=0)
    std_p_rate2[i] = float(0) if np.all((Class_2_array[i])[1]) == np.nan else np.mean((Class_2_array[i])[1], axis=0)
   
    mean_p_rate3[i] = float(0) if np.all((Class_3_array[i])[0]) == np.nan else np.mean((Class_3_array[i])[0], axis=0)
    std_p_rate3[i] = float(0) if np.all((Class_3_array[i])[1]) == np.nan else np.mean((Class_3_array[i])[1], axis=0)

    mean_p_rate4[i] = float(0) if np.all((Class_4_array[i])[0]) == np.nan else np.mean((Class_4_array[i])[0], axis=0)
    std_p_rate4[i] = float(0) if np.all((Class_4_array[i])[1]) == np.nan else np.mean((Class_4_array[i])[1], axis=0)

    mean_p_rate5[i] = float(0) if np.all((Class_5_array[i])[0]) == np.nan else np.mean((Class_5_array[i])[0], axis=0)
    std_p_rate5[i] = float(0) if np.all((Class_5_array[i])[1]) == np.nan else np.mean((Class_5_array[i])[1], axis=0)

    mean_p_rate6[i] = float(0) if np.all((Class_6_array[i])[0]) == np.nan else np.mean((Class_6_array[i])[0], axis=0)
    std_p_rate6[i] = float(0) if np.all((Class_6_array[i])[1]) == np.nan else np.mean((Class_6_array[i])[1], axis=0)

    mean_p_rate7[i] = float(0) if np.all((Class_7_array[i])[0]) == np.nan else np.mean((Class_7_array[i])[0], axis=0)
    std_p_rate7[i] = float(0) if np.all((Class_7_array[i])[1]) == np.nan else np.mean((Class_7_array[i])[1], axis=0)  

    mean_p_rate8[i] = float(0) if np.all((Class_8_array[i])[0]) == np.nan else np.mean((Class_8_array[i])[0], axis=0)
    std_p_rate8[i] = float(0) if np.all((Class_8_array[i])[1]) == np.nan else np.mean((Class_8_array[i])[1], axis=0)

    mean_p_rate9[i] = float(0) if np.all((Class_9_array[i])[0]) == np.nan else np.mean((Class_9_array[i])[0], axis=0)
    std_p_rate9[i] = float(0) if np.all((Class_9_array[i])[1]) == np.nan else np.mean((Class_9_array[i])[1], axis=0)

    mean_p_rate10[i] = float(0) if np.all((Class_10_array[i])[0]) == np.nan else np.mean((Class_10_array[i])[0], axis=0)
    std_p_rate10[i] = float(0) if np.all((Class_10_array[i])[1]) == np.nan else np.mean((Class_10_array[i])[1], axis=0)

    mean_p_rate11[i] = float(0) if np.all((Class_11_array[i])[0]) == np.nan else np.mean((Class_11_array[i])[0], axis=0)
    std_p_rate11[i] = float(0) if np.all((Class_11_array[i])[1]) == np.nan else np.mean((Class_11_array[i])[1], axis=0)  

height = ((Class_1_array[i])[2])

Mean_precip_rate_1 =  np.mean(mean_p_rate1, axis=0)
STD_precip_rate_1 = np.mean(std_p_rate1, axis=0)
Mean_precip_rate_2 =  np.mean(mean_p_rate2, axis=0)
STD_precip_rate_2 = np.mean(std_p_rate2, axis=0)
Mean_precip_rate_3 =  np.mean(mean_p_rate3, axis=0)
STD_precip_rate_3 = np.mean(std_p_rate3, axis=0)
Mean_precip_rate_4 =  np.mean(mean_p_rate4, axis=0)
STD_precip_rate_4 = np.mean(std_p_rate4, axis=0)
Mean_precip_rate_5 =  np.mean(mean_p_rate5, axis=0)
STD_precip_rate_5 = np.mean(std_p_rate5, axis=0)
Mean_precip_rate_6 =  np.mean(mean_p_rate6, axis=0)
STD_precip_rate_6 = np.mean(std_p_rate6, axis=0)
Mean_precip_rate_7 =  np.mean(mean_p_rate7, axis=0)
STD_precip_rate_7 = np.mean(std_p_rate7, axis=0)
Mean_precip_rate_8 =  np.mean(mean_p_rate8, axis=0)
STD_precip_rate_8 = np.mean(std_p_rate8, axis=0)
Mean_precip_rate_9 =  np.mean(mean_p_rate9, axis=0)
STD_precip_rate_9 = np.mean(std_p_rate9, axis=0)
Mean_precip_rate_10 =  np.mean(mean_p_rate10, axis=0)
STD_precip_rate_10 = np.mean(std_p_rate10, axis=0)
Mean_precip_rate_11 =  np.mean(mean_p_rate11, axis=0)
STD_precip_rate_11 = np.mean(std_p_rate11, axis=0)

plt.figure(figsize=(30,15))
plt.scatter(Mean_precip_rate_1,height, label = 'P_n1', linewidth=5)
plt.scatter(Mean_precip_rate_2,height, label = 'P_n2', linewidth=5)
plt.scatter(Mean_precip_rate_3, height, label = 'P_n3', linewidth=5)
plt.scatter(Mean_precip_rate_4, height, label = 'P_n4', linewidth=5)
plt.scatter(Mean_precip_rate_5, height, label = 'P_n5', linewidth=5)
plt.scatter(Mean_precip_rate_6, height, label = 'P_n6', linewidth=5)
plt.scatter(Mean_precip_rate_7, height, label = 'P_n7', linewidth=5)
plt.scatter(Mean_precip_rate_8, height, label = 'P_n8', linewidth=5)
plt.scatter(Mean_precip_rate_9, height, label = 'P_n9', linewidth=5)
plt.scatter(Mean_precip_rate_10, height, label = 'P_n10', linewidth=5)
plt.scatter(Mean_precip_rate_11, height, label = 'P_n11', linewidth=5)
plt.legend()
plt.title(name)
# plt.fill_betweenx(height, mean_p_rate-std_p_rate, mean_p_rate+std_p_rate, alpha=0.2)
plt.savefig(outfile+name+"_mean_"+str(91*240).zfill(10)+'.png')
plt.clf()

