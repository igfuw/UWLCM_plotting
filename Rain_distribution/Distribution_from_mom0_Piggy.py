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
outfile = '/home/pzmij/2D/PAPER/Wyniki/Distribution_from_mom0/Piggy_ziemia/'
name = 'histrogram_base_no_mask'

path_to_file_1 = '/home/pzmij/2D/PAPER/Piggy/SD10_distribution_2/'
path_to_file_2 = '/home/pzmij/2D/PAPER/Piggy/SD100_distribution/'
path_to_file_3 = '/home/pzmij/2D/PAPER/Piggy/SD1000_distribution/'
path_to_file_4 = '/home/pzmij/2D/PAPER/Piggy/SD10000_distribution_2/'



def Distro(timestep, paths, pozycja):
    files = os.listdir(paths)
    nr_files = len(files)
    rhod= [0 for i in range(nr_files)]
    dz = [0 for i in range(nr_files)]
    rl = [0 for i in range(nr_files)]
    rl_base = [0 for i in range(nr_files)]
    rr = [0 for i in range(nr_files)]
    dystrybucja = [0 for i in range(nr_files)]
    wynik = [0 for i in range(nr_files)]

    for file in range(nr_files):
        Name = paths+files[file] + '/'+ files[file]+ "_out_lgrngn"
        rhod[file] = h5py.File(Name + "/const.h5", "r")["G"][:,:]
        dz[file] = h5py.File(Name + "/const.h5", "r").attrs["dz"]

        filename = Name + "/timestep" + str(int(timestep)*240).zfill(10) + ".h5"
        rl[file] = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        rl_base[file] = rl[file]
        rr[file] = (h5py.File(filename, "r")["rain_rw_mom3"][:,:])
        dystrybucja[file] = (h5py.File(filename, "r")["rw_rng"+str(pozycja).zfill(3)+"_mom0"][:,:]) 
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
        wynik[file] = dystrybucja[file][:,int(min_hght/100)-1] #*rainy_mask[:,int(min_hght/100)-1]
        wynik[file][wynik[file]==0]=np.nan
    srednia = np.nanmean(wynik)
    STD = np.nanstd(wynik)
    return srednia, STD

def Positions(paths, timestep):
    Results = [0 for i in range(30)] 
    for i in range(29):
        Results[i] = Distro(timestep, paths, i)
        Results[i] = np.nan_to_num(Results[i])
    return Results, timestep
    
punkty = np.intc(np.linspace(0,90,91))
partial_fun_1 = functools.partial(Positions, path_to_file_1)
partial_fun_2 = functools.partial(Positions, path_to_file_2)
partial_fun_3 = functools.partial(Positions, path_to_file_3)
partial_fun_4 = functools.partial(Positions, path_to_file_4)

with concurrent.futures.ProcessPoolExecutor() as executor:
    Classic_1 = executor.map(partial_fun_1, punkty)
    Classic_2 = executor.map(partial_fun_2, punkty)
    Classic_3 = executor.map(partial_fun_3, punkty)
    Classic_4 = executor.map(partial_fun_4, punkty)



Class_1 = [i for i in Classic_1]
Class_2 = [i for i in Classic_2]
Class_3 = [i for i in Classic_3]
Class_4 = [i for i in Classic_4]


biny = [(1e-6*pow(10, -3+i*0.2)) for i in range(0,31)]
biny_diff = np.diff(np.log(biny[1:]))
Srednie_1 = np.zeros((91, 29))
odchylenie_1 = np.zeros((91, 29))
Srednie_2 = np.zeros((91, 29))
odchylenie_2 = np.zeros((91, 29))
Srednie_3 = np.zeros((91, 29))
odchylenie_3 = np.zeros((91, 29))
Srednie_4 = np.zeros((91, 29))
odchylenie_4 = np.zeros((91, 29))


plt.rcParams.update({'font.size': 40})
plt.rcParams.update({'font.size': 40})
plt.figure(figsize=(40,40))
for i in range(91):
    plt.clf()
    for j in range(29):
        Srednie_1[i][j] = (((Class_1[i])[0])[j])[0]
        odchylenie_1[i][j] = (((Class_1[i])[0])[j])[1]
        Srednie_2[i][j] = (((Class_2[i])[0])[j])[0]
        odchylenie_2[i][j] = (((Class_2[i])[0])[j])[1] 
        Srednie_3[i][j] = (((Class_2[i])[0])[j])[0]
        odchylenie_3[i][j] = (((Class_2[i])[0])[j])[1] 
        Srednie_4[i][j] = (((Class_4[i])[0])[j])[0]
        odchylenie_4[i][j] = (((Class_4[i])[0])[j])[1]

    plt.step(biny[1:-1], Srednie_1[i]/biny_diff[0], c='g', linewidth=9, label='SD10')
    plt.step(biny[1:-1], Srednie_2[i]/biny_diff[0], c='k', linewidth=9, label='SD100')
    plt.step(biny[1:-1], Srednie_3[i]/biny_diff[0], c='r', linewidth=9, label='SD1000')
    plt.step(biny[1:-1], Srednie_4[i]/biny_diff[0], c='b', linewidth=9, label='SD10000')

    plt.title('time= '+str(i*120)+'[s]')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(bottom=0)
    plt.grid(True)
    plt.legend()
    # plt.suptitle('Current time {}s '.format(int(timestep/2)))
    plt.xlabel(r'r [m]')
    plt.ylabel('n(ln r) [# of droplets/ d lnr]')
    plt.savefig(outfile+name+str(i*240).zfill(10)+'.png')
