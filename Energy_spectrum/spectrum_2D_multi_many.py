# (C) Maciej Waruszewski

import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import os
'''
How to run

python3 spectrum_2D_multi_many.py 14400 18400 240 10 30 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Energy/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/ "SD100" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/ "SD1000" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD10000/ "SD10000"
'''
# velocities = ["u", "v", "w"]
# velocities = ["u", "v", "w", "cloud_rw_mom3", "rv", "th", "RH", "aerosol_rw_mom3"]
# velocities = ["u", "w", "cloud_rw_mom3", "rv", "th", "RH", "actrw_rw_mom3"]
# velocities = ["cloud_rw_mom3", "actrw_rw_mom3"]
# velocities = ["cloud_rw_mom3"]
velocities = ["u", "w"]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
from_lvl = int(argv[4])
to_lvl = int(argv[5])
outfile = argv[6]
directories = argv[7:len(argv):2]
labels = argv[8:len(argv):2]


def Energy(directories, labels, velocities):
    # read in nx, ny, nz
    for directory, lab in zip(directories, labels):
        # path_file = []
        path_file = os.listdir(directory)
        Exy_mean = {}
        for path in path_file:
            path_to_file = os.listdir(directory+path)
            joined = os.path.join(directory,path)

            Exy_avg_many = [0 for i in range(len(path_to_file))]
            Exy_avg = [0 for i in range(len(velocities))]
            # Exy_mean = [0 for i in range(len(velocities))]
            print(joined+ "/timestep" + str(time_start).zfill(10) + ".h5", "elo")
            for file in range(len(path_to_file)):
              w3d = h5py.File(joined+ "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:]
              print(joined+ "/timestep" + str(time_start).zfill(10) + ".h5")
              # w3d = h5py.File(joined+'/'+path_to_file[file] + "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:]
              nx, nz = w3d.shape
              for vel in range(len(velocities)):
                Exy_avg[vel] = np.zeros(int((nx+1)/2))
                # Exy_avg_many[file][vel] =np.zeros(int((nx+1)/2))

              for t in range(time_start, time_end+1, outfreq):
                filename = joined+"/timestep" + str(t).zfill(10) + ".h5"
                # filename = joined+'/'+path_to_file[file] + "/timestep" + str(t).zfill(10) + ".h5"
                # print(filename)

                for vel in range(len(velocities)):
                  if (velocities[vel] == "cloud_rw_mom3") | (velocities[vel] == "actrw_rw_mom3"):
                    w3d = h5py.File(filename, "r")[velocities[vel]][:,:] * 4. / 3. * 3.1416 * 1e3
                  else:
                    w3d = h5py.File(filename, "r")[velocities[vel]][:,:] # * 4. / 3. * 3.1416 * 1e3

                  for lvl in range(from_lvl, to_lvl+1):
                    w2d = w3d[:, lvl]
                    wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
                    Ex = (np.abs(wkx) ** 2)
                    Exy = Ex
                    Exy_avg[vel] += Exy
                Exy_avg_many[file] = Exy_avg[-len(velocities):]
            Exy_mean[path] = np.mean(Exy_avg_many[file],axis=0)

        K = np.fft.rfftfreq(nx - 1)
        lmbd = 100. / K # assume dx=100m

        for path in path_file:
            Exy_mean[path] /= (time_end - time_start) / outfreq + 1
            Exy_mean[path] /= to_lvl+1 - from_lvl
            plt.loglog(lmbd, Exy_mean[path] , linewidth=2, label=velocities+'_'+lab+"_"+path)
'''
    plt.loglog(lmbd, 2e-1* K**(-5./3.), label="-5/3" )
    plt.xlim(2*10**4,10**2)
    plt.xlabel("l[m]")
    plt.ylabel("PSD")
    plt.legend()
    plt.grid(True, which='both', linestyle='--')
        # plt.title("Mean PSD of w 322m<z<642m @3h")
    plt.savefig(outfile + 'Energy_spc.png')
    plt.show()
'''
Energy(directories, labels, velocities[0])
# Energy(directories, labels, velocities[1])
