import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from collections import OrderedDict
import os

'''
HOW TO run

python3 distribution_Py3_2D.py 14400 20400 240 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/dxyz100_SD100_diff_seed_NA1_Coal_2D_VF_piggy_25_out_lgrngn "test" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/dxyz100_SD100_diff_seed_NA1_Coal_2D_VF_piggy_26_out_lgrngn "test2"

'''

#rain_data = ["u", "v", "w"]
#rain_data = ["u", "v", "w", "cloud_rw_mom3", "rv", "th", "RH", "aerosol_rw_mom3"]
#rain_data = ["cloud_rw_mom3"]
rain_data = ["rain_rw_mom3", "precip_rate"]
layer_thickness = 10
# cloud_thresh = 1e-4 # qc at cloud base
cloud_thresh = 1e-5 # qc at cloud base

CellVol = 100.*100*1 # hardcoded cell volume [m^3]
L_evap = 2264.76e3 # latent heat of evapporation [J/kg]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
#from_lvl = int(argv[4])
#to_lvl = int(argv[5])

directories = argv[4:len(argv):2]
labels = argv[5:len(argv):2]
print(directories, labels)

levels = ["ground", "cloud_base"]


def Rain_distribution(directories, labels):
  for directory, lab in zip(directories, labels):
      path_file = os.listdir(directory)
      Rain_mean = {}
      # read in nx, nz
      for lvl in levels:
        total_arr = OrderedDict()
        for data in rain_data:
          total_arr[data] = OrderedDict()

        plot_labels = OrderedDict()
        tot_cloud_base_lvl = OrderedDict()
        for path in path_file:
          path_to_file = os.listdir(directory+path)
          joined = os.path.join(directory,path)

          file = '#'+(path.split("_piggy_",-3))[-1].split("_out")[0]
          tot_cloud_base_lvl[file] = np.zeros(0)

          # print((path.split("_piggy_",-3))[-1].split("_out"))
          w3d = h5py.File(joined+"/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:]
          nx, nz = w3d.shape

          for data in rain_data:
              total_arr[data][file] = np.zeros(0)

              for t in range(time_start, time_end+1, outfreq):
                filename = joined + "/timestep" + str(t).zfill(10) + ".h5"
                # find cloud base
                # based on cloud rw
                w3d = h5py.File(filename, "r")["cloud_rw_mom3"][:,:] * 4. / 3. * 3.1416 * 1e3
                cloud_base_lvl = np.argmax(np.average(w3d, axis=(0)) > cloud_thresh)
                # based on RH
          #      w3d = h5py.File(filename, "r")["RH"][:,:,:] # * 4. / 3. * 3.1416 * 1e3
          #      cloud_base_lvl = np.argmax(np.average(w3d, axis=(0,1)) > .99)

                tot_cloud_base_lvl[file] = np.append(tot_cloud_base_lvl[file], cloud_base_lvl) # done for each data, but we dont care - wont affect average

                # print('cloud base lvl = ', cloud_base_lvl)

                if lvl == "cloud_base":
                  total_arr[data][file] = np.append(total_arr[data][file], h5py.File(filename, "r")[data][:,cloud_base_lvl-layer_thickness : cloud_base_lvl])
                if lvl == "ground":
                  total_arr[data][file] = np.append(total_arr[data][file], h5py.File(filename, "r")[data][:, 0 : layer_thickness ])

              # convert to typical units
              if data == "rain_rw_mom3":
                total_arr[data][file] *= 4./3. * 3.1416 * 1e3 * 1e3 # [g/kg]
              if data == "precip_rate":
                total_arr[data][file] *= 4./3. * 3.1416 * 1e3 / CellVol * L_evap
  A = [i for i in range(len(path_file))]
  B = [i for i in range(len(path_file))]
  i =0
  for path in path_file:
    file = '#'+(path.split("_piggy_",-3))[-1].split("_out")[0]
    A[i] =  total_arr["rain_rw_mom3"][file]
    B[i] =  total_arr["precip_rate"][file]
    i += 1
  A_mean = np.average(A, axis=0)
  B_mean = np.average(B, axis=0)


  plt.figure(0)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  _ = plt.hist(A_mean, bins=np.logspace(np.log10(1e-6), np.log10(np.amax(A_mean)), nx), label=plot_labels.values(), density=False, histtype='step', linewidth=6)

  plt.xscale('log')
  # plt.legend(loc = 'lower center')
  plt.yscale('log')
  plt.xlabel('q_r [g/kg]')
  plt.ylabel('# of cells')
  plt.savefig('rain_histo_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')

  plt.figure(1)
  plt.rcParams.update({'font.size': 40})
  plt.figure(figsize=(40,40))
  _ = plt.hist(B_mean, bins=np.logspace(np.log10(1e-3), np.log10(np.amax(B_mean)), nx), label=plot_labels.values(), density=False, histtype='step', linewidth=6)

  plt.xscale('log')
  # plt.legend(loc = 'lower center')
  plt.yscale('log')
  plt.xlabel('precipitation flux [W / m^2]')
  plt.ylabel('# of cells')
  plt.savefig('precip_rate_histo_' + lvl + '_' + str(time_start) + '_' + str(time_end) +'.png')


Rain_distribution(directories, labels)
