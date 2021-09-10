# calculate cloud droplet conc. vs adiabatic fraction
# adiabatic rl calculated based on mean th and rv at cloud base cells

from sys import argv, path, maxsize
#path.insert(0,"../../local_folder/uptodate/lib/python3/dist-packages")
# path.insert(0,"/home/piotr/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
#path.insert(0,"/home/piotr-pc/Piotr/IGF/local_install/parcel/lib/python3/dist-packages")
path.insert(0,"/home/pzmij/biblioteki/local_folder/16_03/lib/python3/dist-packages")


'''
How to run

python3 Rain_many.py /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/dxyz100_SD100_diff_seed_NA1_Coal_2D_VF_tail_piggy_25/dxyz100_SD100_diff_seed_NA1_Coal_2D_VF_tail_piggy_25_out_lgrngn /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain

Singularity> python3 Rain_many.py /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/

python3 Rain_many_copare_time_line.py 14400 20400 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_SD1000" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_SD100"


'''

from math import exp, log, sqrt, pi, erf, cos, pow, asin, atan, acos, factorial
import h5py
from scipy.stats import moment
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
from libcloudphxx import common as lcmn
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

start = time.perf_counter()
plt.rcParams.update({'font.size': 20})
evap_lat = 2.501e6 # [J/kg] latent heat of evaporation
timesteps = np.ones(91)
for i in range(1, 91):
    timesteps[i] = i
timesteps = timesteps[1:91]
np.set_printoptions(threshold=maxsize)


def Adia_fraction(timestep, paths):
    files = os.listdir(paths)
    nr_files = len(files)
    rhod= [0 for i in range(nr_files)]
    p_e = [0 for i in range(nr_files)]
    dz = [0 for i in range(nr_files)]
    rl = [0 for i in range(nr_files)]
    rl_base = [0 for i in range(nr_files)]
    nc = [0 for i in range(nr_files)]
    th = [0 for i in range(nr_files)]
    rv = [0 for i in range(nr_files)]
    rr = [0 for i in range(nr_files)]
    sum_rr = [0 for i in range(nr_files)]


    for file in range(nr_files):
        # plt.clf()
        p_e[file] = h5py.File(paths+files[file] + '/' + files[file] + "_out_lgrngn" + "/const.h5", "r")["p_e"][:]
        rhod[file] = h5py.File(paths+files[file] + '/'+ files[file]+ "_out_lgrngn" + "/const.h5", "r")["G"][:,:]
        dz[file] = h5py.File(paths+files[file] + '/' + files[file] + "_out_lgrngn"+ "/const.h5", "r").attrs["dz"]

        nx, nz = rhod[file].shape
        hght = np.arange(nz) * dz[file]
        bin = np.linspace(1,121,len(np.arange(nx)))
        sum_rr[file] = np.zeros([nx])

        filename = paths + files[file] + '/' + files[file] + "_out_lgrngn"+ "/timestep" + str(int(timestep)*240).zfill(10) + ".h5"
        rl[file] = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        rl_base[file] = (h5py.File(filename, "r")["cloud_rw_mom3"][:,:]) * 4. / 3. * 3.1416 * 1e3; # kg/kg
        nc[file] = h5py.File(filename, "r")["cloud_rw_mom0"][:,:] * rhod[file] / 1e6; # 1 / cm^3
        rr[file] = (h5py.File(filename, "r")["rain_rw_mom3"][:,:]) * rhod[file] * 4. / 3. * 3.1416 * 1e6; # g/kg
        # cloudiness mask - as in RICO paper
        cloudy_mask = np.where(rl[file] > 1e-5, 1, 0)
        cloudy_mask_used = cloudy_mask

        # th and rv
        th[file] = h5py.File(filename, "r")["th"][:,:];
        rv[file] = h5py.File(filename, "r")["rv"][:,:];
        # T
        Vexner= np.vectorize(lcmn.exner)
        T = th[file] * Vexner(p_e[file].astype(float))
        # RH
        Vr_vs = np.vectorize(lcmn.r_vs)
        r_vs = Vr_vs(T, p_e[file].astype(float))
        RH = rv[file] / r_vs[file]
        # cloud base
        clb_idx = np.argmax(cloudy_mask_used > 0, axis=1)
        hght_2 = hght[clb_idx]
        hght_2[hght_2==0] = np.nan
        min_hght = np.nanmin(hght_2)
        max_hght = np.nanmax(hght_2)
        # print("wysokosc min",min_hght, " wysokosc max", max_hght)
        # print(i, rr[file][i])
        if min_hght/dz[file] < 10:
          min_hght = 10*int(dz[file])

        min_hght = np.nan_to_num(min_hght)

        for j in np.arange(nx):
          if clb_idx[j] > 0:
            # sum_rr[file][j] = np.sum(rr[file][j,:int(min_hght/dz[file])-1],0)
            sum_rr[file][j] = np.sum(rr[file][j,int(min_hght/dz[file])-1],0)
        sum_rr[file][sum_rr[file]==0] = np.nan
        # print(file)
    srednie_rr = np.nanmean(sum_rr, axis=0)
    STD_rr = np.nanstd(sum_rr, axis=0)
    ST_error_mean = STD_rr/sqrt(nr_files)
    error_std = np.power(1/nr_files * (moment(sum_rr,4) - (nr_files-3)/(nr_files-1)*np.power(STD_rr,4)),1/2)/(2*STD_rr)
    return(srednie_rr, STD_rr, bin, ST_error_mean, error_std)

# bin_size = bin[2]-bin[1]
# print(bin_size)
time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
outfile = argv[4]
paths = argv[5:len(argv):2]
labels = argv[6:len(argv):2]


Average = []
STD = []
Average_error = []
STD_error = []
average =np.zeros((int((time_end-time_start)/outfreq)+1,121))
average_error =np.zeros((int((time_end-time_start)/outfreq)+1,121))
std =np.zeros((int((time_end-time_start)/outfreq)+1,121))
std_error =np.zeros((int((time_end-time_start)/outfreq)+1,121))
# fig = plt.figure()
# fig.set_size_inches(18.5, 10.5)
fig, (ax0, ax1) = plt.subplots(1, 2,figsize=(30,15))
shape = ['o', 'x', '.k' ]
for path, lab in zip(paths, labels):
    b = 0
    for i in range(int(time_start/outfreq), int(time_end/outfreq)+1):
        aver, sigma , bin, st_error, sigma_error= Adia_fraction(i, path)
        average[b] = aver
        std[b] = sigma
        std_error[b] = sigma_error
        average_error[b] = st_error
        # std[b] = std
        b += 1
    Average = np.nanmean(average,axis=0)
    STD = np.nanmean(std,axis=0)
    Average_error = np.nanmean(average_error,axis=0)
    STD_error = np.nanmean(std_error,axis=0)
    ax0.plot(bin, Average, label=lab)#,  edgecolors='b'
    ax0.fill_between(bin, Average-Average_error, Average+Average_error, alpha=0.2)
    # ax1.plot(bin, STD, label=lab)#,  edgecolors='b'
    # ax1.fill_between(bin, STD-STD_error, STD+STD_error, alpha=0.2)
    ax1.errorbar(bin, STD, STD_error, fmt='o',label=lab)

    # plt.errorbar(bin, Average, Average_error, fmt='o',label=lab)#
ax0.set_ylim((0))
ax1.set_ylim((0))
ax0.set_xlim((40,90))
ax1.set_xlim((40,90))
plt.xlabel('cell#')
ax0.set_ylabel('$q_r$ (rain) [g/kg]')
ax1.set_ylabel('$q_r$ (rain) [g/kg]')
ax0.legend(title='Average')
ax1.legend(title='STD')
fig.suptitle('time range = '+str(time_start) + ' to '+str(time_end) +' s')

plt.savefig(outfile + 'Average_many_for_time range_'+str(time_start)+'_to_'+ str(time_end)+'.png')
# Od komentuj!
plt.clf()

'''
python3 Rain_many_copare_time.py 14400 20400 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ 'multi'

python3 Rain_many_copare_time.py 14400 20400 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/multi/ 'multi_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ 'multi_1c'


python3 Rain_many_copare_time.py 14400 15600 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 16800 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 18000 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 19200 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 17280 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 17760 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 17040 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'
python3 Rain_many_copare_time.py 14400 20400 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Rain/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "VF_1k" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ 'tail_1k' /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "VF_1c" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/tail/ 'tail_1c'







'''
