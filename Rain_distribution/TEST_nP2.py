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


#########parametry do uruchomienia
outfile = '/home/pzmij/2D/PAPER/Wyniki/Distribution/no_Piggy/'
name = 'no_Piggy_TEST1_'
path_to_file_1 = '/home/pzmij/2D/PAPER/Distribution/no_Piggy/SD10000/classic/test/'
# path_to_file_2 = '/home/pzmij/2D/PAPER/Distribution/no_Piggy/SD100/tail/'
# path_to_file_3 = '/home/pzmij/2D/PAPER/Distribution/no_Piggy/SD1000/classic/'
# path_to_file_4 = '/home/pzmij/2D/PAPER/Distribution/no_Piggy/SD1000/tail/'
# path_to_file_5 = '/home/pzmij/2D/PAPER/Distribution/no_Piggy/SD10000/classic/'
# path_to_file_6 = '/home/pzmij/2D/PAPER/Distribution/no_Piggy/SD10000/tail/'


def Distro(timestep, paths, pozycja):
    files = os.listdir(paths)
    nr_files = len(files)
    rhod= [0 for i in range(nr_files)]
    dz = [0 for i in range(nr_files)]
    rl = [0 for i in range(nr_files)]
    rl_base = [0 for i in range(nr_files)]
    rr = [0 for i in range(nr_files)]
    dystrybucja_mom0 = [0 for i in range(nr_files)]
    dystrybucja_mom3 = [0 for i in range(nr_files)]
    wynik_number = [0 for i in range(nr_files)]
    wynik_mass = [0 for i in range(nr_files)]

    for file in range(nr_files):
        Name = paths+files[file] + '/'+ files[file]+ "_out_lgrngn"
        rhod[file] = h5py.File(Name + "/const.h5", "r")["G"][:,:]
        dz[file] = h5py.File(Name + "/const.h5", "r").attrs["dz"]

        filename = Name + "/timestep" + str(int(timestep)*240).zfill(10) + ".h5"
        rl[file] = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        rl_base[file] = rl[file]
        rr[file] = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) 
        dystrybucja_mom0[file] = (h5py.File(filename, "r")["rw_rng"+str(pozycja).zfill(3)+"_mom0"][:,:]) 
        dystrybucja_mom3[file] = (h5py.File(filename, "r")["rw_rng"+str(pozycja).zfill(3)+"_mom3"][:,:]) * rhod[file] * 4. / 3. * 3.1416 * 1e3; # kg/kg
        # cloudiness mask - as in RICO paper
        cloudy_mask = np.where(rl[file] > 1e-5, 1, 0)
        rainy_mask = np.where(rr[file] > 0, 1, 0)

        nx, nz = rr[file].shape
        hght = np.arange(nz) * dz[file]
        # cloud base
        clb_idx = np.argmax(cloudy_mask > 0, axis=1)
        hght_2 = hght[clb_idx]
        hght_2[hght_2==0] = np.nan
        if np.isnan(hght_2).all() == True :
            continue
        min_hght = np.nanmin(hght_2)
        if min_hght ==np.nan:
            break
        wynik_number[file] = dystrybucja_mom0[file][:,0]#* rainy_mask[:,0]
        wynik_mass [file] = dystrybucja_mom3[file][:,0]#* rainy_mask[:,0]
        wynik_number[file][wynik_number[file]==0]=np.nan
        wynik_mass[file][wynik_mass[file]==0]=np.nan
    W_n = np.array(wynik_number)
    W_m = np.array(wynik_mass)
    srednia_number = np.nanmean(W_n)
    srednia_mass = np.nanmean(W_m)
    STD = 0 #np.nanstd(wynik_number)
    return srednia_number, STD, srednia_mass

def Positions(paths, timestep):
    Results = [0 for i in range(40)] 
    for i in range(39):
        Results[i] = Distro(timestep, paths, i)
        Results[i] = np.nan_to_num(Results[i])
    return Results, timestep
    

punkty = np.intc(np.linspace(0,90,91))
partial_fun_1 = functools.partial(Positions, path_to_file_1)

with concurrent.futures.ProcessPoolExecutor() as executor:
    Classic_1 = executor.map(partial_fun_1, punkty)



Class_1 = [i for i in Classic_1]
Class_1_array = np.array(Class_1)


biny = [(1e-6*pow(10, -3+i*0.2)) for i in range(0,41)]
biny_diff = np.diff(np.log(biny[1:]))
Srednie_num1 = np.zeros((91, 39))
Srednie_mass1 = np.zeros((91, 39))
odchylenie_num1 = np.zeros((91, 39))

plt.rcParams.update({'font.size': 40})
for i in range(80,91):
    fig, (ax0, ax1) = plt.subplots(1, 2,figsize=(40,40))
    for j in range(39):
        Srednie_num1[i][j] = float(0) if np.all((((Class_1_array[i])[0])[j])[0]) != np.nan else (((Class_1_array[i])[0])[j])[0]
        odchylenie_num1[i][j] = float(0) if np.all((((Class_1_array[i])[0])[j])[1]) != np.nan else (((Class_1_array[i])[0])[j])[1]
        Srednie_mass1[i][j] = float(0) if np.all((((Class_1_array[i])[0])[j])[2]) != np.nan else (((Class_1_array[i])[0])[j])[2]
