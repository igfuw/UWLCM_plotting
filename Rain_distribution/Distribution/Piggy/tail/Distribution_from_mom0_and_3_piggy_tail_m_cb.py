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
outfile = '/home/pzmij/2D/PAPER/Wyniki/Distribution/Piggy/tail/mask/cb/'
name = 'Piggy_tail_mask_cb_'
path_to_file_1 = '/home/pzmij/2D/PAPER/Distribution/Piggy/SD100/tail/'
path_to_file_2 = '/home/pzmij/2D/PAPER/Distribution/Piggy/SD1000/tail/'
path_to_file_3 = '/home/pzmij/2D/PAPER/Distribution/Piggy/SD10000/tail/'



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
        wynik_number[file] = dystrybucja_mom0[file][:,int(min_hght/100)-1]* rainy_mask[:,int(min_hght/100)-1]
        wynik_mass [file] = dystrybucja_mom3[file][:,int(min_hght/100)-1]* rainy_mask[:,int(min_hght/100)-1]
        wynik_number[file][wynik_number[file]==0]=np.nan
        wynik_mass[file][wynik_mass[file]==0]=np.nan
    srednia_number = np.nanmean(wynik_number)
    srednia_mass = np.nanmean(wynik_mass)
    STD = np.nanstd(wynik_number)
    return srednia_number, STD, srednia_mass

def Positions(paths, timestep):
    Results = [0 for i in range(40)] 
    for i in range(39):
        Results[i] = Distro(timestep, paths, i)
        Results[i] = np.nan_to_num(Results[i])
    return Results, timestep
    
punkty = np.intc(np.linspace(0,90,91))
partial_fun_1 = functools.partial(Positions, path_to_file_1)
partial_fun_2 = functools.partial(Positions, path_to_file_2)
partial_fun_3 = functools.partial(Positions, path_to_file_3)
with concurrent.futures.ProcessPoolExecutor() as executor:
    Classic_1 = executor.map(partial_fun_1, punkty)
    Classic_2 = executor.map(partial_fun_2, punkty)
    Classic_3 = executor.map(partial_fun_3, punkty)


Class_1 = [i for i in Classic_1]
Class_2 = [i for i in Classic_2]
Class_3 = [i for i in Classic_3]


biny = [(1e-6*pow(10, -3+i*0.2)) for i in range(0,41)]
biny_diff = np.diff(np.log(biny[1:]))
Srednie_num1 = np.zeros((91, 39))
Srednie_mass1 = np.zeros((91, 39))
odchylenie_num1 = np.zeros((91, 39))
Srednie_num2 = np.zeros((91, 39))
Srednie_mass2 = np.zeros((91, 39))
odchylenie_num2 = np.zeros((91, 39))
Srednie_num3 = np.zeros((91, 39))
Srednie_mass3 = np.zeros((91, 39))
odchylenie_num3 = np.zeros((91, 39))


plt.rcParams.update({'font.size': 40})
for i in range(91):
    fig, (ax0, ax1) = plt.subplots(1, 2,figsize=(40,40))
    for j in range(39):
        Srednie_num1[i][j] = (((Class_1[i])[0])[j])[0]
        odchylenie_num1[i][j] = (((Class_1[i])[0])[j])[1]
        Srednie_mass1[i][j] = (((Class_1[i])[0])[j])[2]
        Srednie_num2[i][j] = (((Class_2[i])[0])[j])[0]
        odchylenie_num2[i][j] = (((Class_2[i])[0])[j])[1]
        Srednie_mass2[i][j] = (((Class_2[i])[0])[j])[2]
        Srednie_num3[i][j] = (((Class_3[i])[0])[j])[0]
        odchylenie_num3[i][j] = (((Class_3[i])[0])[j])[1]
        Srednie_mass3[i][j] = (((Class_3[i])[0])[j])[2]

    ax0.step(biny[1:-1], Srednie_num1[i]/biny_diff[0], c='g', linewidth=9, label='SD100')
    ax0.step(biny[1:-1], Srednie_num2[i]/biny_diff[0], c='k', linewidth=9, label='SD1000')
    ax0.step(biny[1:-1], Srednie_num3[i]/biny_diff[0], c='m', linewidth=9, label='SD10000')
    ax1.step(biny[2:], Srednie_mass1[i]/biny_diff[0], c='g', linewidth=9, label='SD100')
    ax1.step(biny[2:], Srednie_mass2[i]/biny_diff[0], c='k', linewidth=9, label='SD1000')
    ax1.step(biny[2:], Srednie_mass3[i]/biny_diff[0], c='m', linewidth=9, label='SD10000')

    #lub opcja ze liczba na jednym a na drugim masa
    fig.suptitle('time= '+str(i*120)+'[s]')
    ax0.set_xscale('log')
    ax1.set_xscale('log')
    ax0.set_yscale('log')
    ax1.set_yscale('log')
    # plt.yscale('log')
    # plt.xlim((1e-3, 1e-0))
    ax0.set_ylim(bottom=0)
    ax1.set_ylim(bottom=0)
    ax0.grid(True)
    ax0.legend()
    ax1.grid(True)
    ax1.legend()
    # plt.suptitle('Current time {}s '.format(int(timestep/2)))
    ax0.set_xlabel(r'r [m]')
    ax0.set_ylabel('n(ln r) [# of droplets/ d lnr]')
    ax1.set_xlabel(r'r [m]')
    ax1.set_ylabel('rr / log(bin size) [(kg/kg) / d lnr]')
    plt.savefig(outfile+name+str(i*240).zfill(10)+'.png')
    plt.clf()