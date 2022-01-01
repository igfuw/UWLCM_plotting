import h5py
import numpy as np
from sys import argv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import OrderedDict
import os

'''
HOW TO run

python3 distribution_Py3_2D_mean.py 14400 20400 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "test" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/

'''

rain_data = ["rain_rw_mom3", "precip_rate", "cloud_rw_mom3"]
cloud_thresh = 1e-5 # qc at cloud base

CellVol = 100.*100*1 # hardcoded cell volume [m^3]
L_evap = 2264.76e3 # latent heat of evapporation [J/kg]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])


outfile = argv[4]
directories = argv[5:len(argv):2]
labels = argv[6:len(argv):2]

total_arr = OrderedDict()
for data in rain_data:
    total_arr[data] = OrderedDict()
    plot_labels = OrderedDict()

def Rain_distribution(directories, labels):
  RAIN, CLOUD, PRECIP = [[i for i in range(len(directories))] for i in range(3)]
  for directory, lab in zip(directories, labels):
      path_file = os.listdir(directory)
      Rain_mean = {}
      j = 0
      for path in path_file:
        path_to_file = os.listdir(directory+path)
        joined = os.path.join(directory,path)

        file = '#'+(path.split("_piggy_",-3))[-1].split("_out")[0]
        w3d = h5py.File(joined+'/'+joined.split("/",-1)[-1] +"_out_lgrngn"+"/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:]
        nx, nz = w3d.shape

        Tablica = np.zeros((nx, nz))
        for data in rain_data:
          total_arr[data][file] = np.zeros(0)
          for t in range(time_start, time_end+1, outfreq):
            filename = joined +'/'+joined.split("/",-1)[-1] +"_out_lgrngn"+ "/timestep" + str(t).zfill(10) + ".h5"
            Tablica += h5py.File(filename, "r")[data][:, : ]

          total_arr[data][file] = np.append(total_arr[data][file], Tablica/((time_end-time_start)/outfreq+1))


              # convert to typical units
          if data == "rain_rw_mom3":
            total_arr[data][file] *= 4./3. * 3.1416 * 1e3  # [kg/kg]
          if data == "cloud_rw_mom3":
            total_arr[data][file] *= 4./3. * 3.1416 * 1e3  # [kg/kg]
          if data == "precip_rate":
            total_arr[data][file] *= 4./3. * 3.1416 * 1e3 / CellVol * L_evap
      Rain = [i for i in range(len(path_file))]
      Cloud = [i for i in range(len(path_file))]
      Precip = [i for i in range(len(path_file))]
      i =0
      for path in path_file:
        file = '#'+(path.split("_piggy_",-3))[-1].split("_out")[0]
        Rain[i] =  total_arr["rain_rw_mom3"][file]
        Cloud[i] =  total_arr["cloud_rw_mom3"][file]
        Precip[i] =  total_arr["precip_rate"][file]
        i += 1
      Rain_mean = np.average(Rain, axis=0)
      Rain_Std = np.std(Rain, axis=0)
      Cloud_Std = np.std(Cloud, axis=0)
      Cloud_mean = np.average(Cloud, axis=0)
      Precip_mean = np.average(Precip, axis=0)
      RAIN[j] = Rain_mean
      CLOUD[j] = Cloud_mean
      PRECIP[j] = Precip_mean
      j += 1


  plt.figure(0)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  for k in range(len(directories)):
      # _ = plt.hist(RAIN[k], bins=np.logspace(np.log10(np.amin(RAIN[k])), np.log10(np.amax(RAIN[k])), nx),   label=labels[k], density=False, histtype='step', linewidth=6)
      _ = plt.hist(RAIN[k], bins=50,   label=labels[k], density=False, histtype='step', linewidth=6)
  plt.xscale('log')
  plt.legend(loc = 'best')
  plt.yscale('log')
  plt.xlabel('q_r [g/kg]')
  plt.ylabel('# of cells')
  plt.savefig(outfile+'A_Mean_rain_histo_'   + str(time_start) + '_' + str(time_end) +'.png')

  plt.figure(1)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  for k in range(len(directories)):
      # _ = plt.hist(PRECIP[j], bins=np.logspace(np.log10(1e-11), np.log10(np.amax(PRECIP[j])), nx), label=labels[k], density=False, histtype='step', linewidth=6)
      _ = plt.hist(PRECIP[j], bins=50, label=labels[k], density=False, histtype='step', linewidth=6)

  plt.xscale('log')
  plt.legend(loc = 'best')
  plt.yscale('log')
  plt.xlabel('precipitation flux [W / m^2]')
  plt.ylabel('# of cells')
  plt.savefig(outfile+'A_Mean_precip_rate_histo_'   + str(time_start) + '_' + str(time_end) +'.png')

  plt.figure(0)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  for k in range(len(directories)):
     # _ = plt.hist(CLOUD[k], bins=np.logspace(np.log10(1e-11), np.log10(np.amax(CLOUD[k])), nx),   label=labels[k], density=False, histtype='step', linewidth=6)
     _ = plt.hist(CLOUD[k], bins=50,   label=labels[k], density=False, histtype='step', linewidth=6)
  plt.xscale('log')
  plt.legend(loc = 'best')
  plt.yscale('log')
  plt.xlabel('q_c [kg/kg]')
  plt.ylabel('# of cells')
  plt.savefig(outfile+'A_Mean_cloud_water_histo_'   + str(time_start) + '_' + str(time_end) +'.png')


Rain_distribution(directories, labels)
