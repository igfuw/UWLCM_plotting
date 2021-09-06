# (C) Maciej Waruszewski

import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import os
'''
How to run

python3 spectrum_2D_multi_many.py 14400 18400 240 10 30 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Energy/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/ "SD100" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/ "SD1000" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD10000/ "SD10000"

python3 spectrum_2D_multi_many_modyfikacje.py 14400 18400 240 10 30 /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Energy/ /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD100/VF/ "SD100" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/VF/ "SD100_tail" /home/piotr-pc/Piotr/WORKSHOPS/Dane_do_AF_2D/Dane/SD1000/tail/ "SD100_tail2"

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


def Energy(directories, labels):
    for directory, lab in zip(directories, labels):
        path_file = os.listdir(directory)
        Exy_mean = {}
        Exy_avg_many = {}
        Exy_avg = {}
        for path in path_file:
            path_to_file = os.listdir(directory+path)
            joined = os.path.join(directory,path)
            w3d = h5py.File(joined+ "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:]
            nx, nz = w3d.shape
            for vel in velocities:
                Exy_avg[path, vel] = np.zeros(int((nx+1)/2))
            for t in range(time_start, time_end+1, outfreq):
                filename = joined+"/timestep" + str(t).zfill(10) + ".h5"
                for vel in velocities:
                    if (vel == "cloud_rw_mom3") | (vel == "actrw_rw_mom3"):
                        w3d = h5py.File(filename, "r")[vel][:,:] * 4. / 3. * 3.1416 * 1e3
                    else:
                        w3d = h5py.File(filename, "r")[vel][:,:]
                    for lvl in range(from_lvl, to_lvl+1):
                        w2d = w3d[:, lvl]
                        wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
                        Ex = (np.abs(wkx) ** 2)
                        Exy = Ex
                        Exy_avg[path, vel] += Exy

        K = np.fft.rfftfreq(nx - 1)
        lmbd = 100. / K # assume dx=100m
        A = np.zeros((len(velocities), len(Exy_avg[path, vel]) ))
        for vel in range(len(velocities)):
            for path in path_file:
                Exy_avg[path,velocities[vel]] /= (time_end - time_start) / outfreq + 1
                Exy_avg[path,velocities[vel]] /= to_lvl+1 - from_lvl
                # Exy_mean[directory, vel] = Exy_avg[vel]
                A[vel,:] += Exy_avg[path,velocities[vel]]
                # print(Exy_avg[path,velocities[vel]])
            A = A/len(path_file)
            plt.loglog(lmbd, A[vel] , linewidth=2, label=velocities[vel]+'_'+lab)

    plt.loglog(lmbd, 2e-1* K**(-5./3.), label="-5/3" )
    plt.xlim(2*10**4,10**2)
    plt.xlabel("l[m]")
    plt.ylabel("PSD")
    plt.legend()
    plt.grid(True, which='both', linestyle='--')
        # plt.title("Mean PSD of w 322m<z<642m @3h")
    plt.savefig(outfile + 'Energy_spc.png')
    plt.show()

Energy(directories, labels)
# Energy(directories, labels, velocities[1])
