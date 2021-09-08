# (C) Maciej Waruszewski

import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt

# velocities = ["u", "v", "w"]
# velocities = ["u", "v", "w", "cloud_rw_mom3", "rv", "th", "RH", "aerosol_rw_mom3"]
# velocities = ["u", "w", "cloud_rw_mom3", "rv", "th", "RH", "actrw_rw_mom3"]
velocities = ["cloud_rw_mom3", "actrw_rw_mom3"]
# velocities = ["u", "w"]

time_start = int(argv[1])
time_end = int(argv[2])
outfreq = int(argv[3])
from_lvl = int(argv[4])
to_lvl = int(argv[5])

directories = argv[6:len(argv):2]
labels = argv[7:len(argv):2]
print(directories, labels)

# read in nx, ny, nz
for directory, lab in zip(directories, labels):
  w3d = h5py.File(directory + "/timestep" + str(time_start).zfill(10) + ".h5", "r")["u"][:,:]
  nx, nz = w3d.shape
  Exy_avg = {}
  for vel in velocities:
    Exy_avg[vel] = np.zeros(int((nx+1)/2))

  for t in range(time_start, time_end+1, outfreq):
    filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
    print(filename)

    for vel in velocities:
      if (vel == "cloud_rw_mom3") | (vel == "actrw_rw_mom3"):
        w3d = h5py.File(filename, "r")[vel][:,:]* 4. / 3. * 3.1416 * 1e3
      else:
        w3d = h5py.File(filename, "r")[vel][:,:]

      for lvl in range(from_lvl, to_lvl+1):
        w2d = w3d[:, lvl]
        wkx = 1.0 / np.sqrt(nx - 1) * np.fft.rfft(w2d, axis = 0)
        Ex = (np.abs(wkx) ** 2)
        Exy = Ex
        Exy_avg[vel] += Exy

      K = np.fft.rfftfreq(nx - 1)

      lmbd = 100. / K # assume dx=50m

    if (t == time_start and lab==labels[0]):
      plt.loglog(lmbd, 2e-1* K**(-5./3.), label="-5/3" )

  for vel in velocities:
    Exy_avg[vel] /= (time_end - time_start) / outfreq + 1
    Exy_avg[vel] /= to_lvl+1 - from_lvl
    plt.loglog(lmbd, Exy_avg[vel] , linewidth=2, label=lab+"_"+vel)

plt.xlim(2*10**4,10**2)
plt.xlabel("l[m]")
plt.ylabel("PSD")
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.show()
