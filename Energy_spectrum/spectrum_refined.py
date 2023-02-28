# (C) Maciej Waruszewski, Piotr Dziekan

import argparse
import h5py
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['figure.figsize'] = 10, 10

parser = argparse.ArgumentParser(description='Plot energy spectra from UWLCM output.')

parser.add_argument("-v", "--vars", action="extend", nargs="+", type=str, help="list of variables to be plotted", required=True)
parser.add_argument("-ts", "--time_start", type=float, required=True, help="start of the averaging period [s]")
parser.add_argument("-te", "--time_end", type=float, required=True, help="end of the averaging period [s]")
parser.add_argument("-ls", "--level_start", type=float, required=True, help="lowest level of the averaging area [m]")
parser.add_argument("-le", "--level_end", type=float, required=True, help="highest level of the averaging area [m]")
parser.add_argument("-d", "--dirs", action="extend", nargs="+", type=str, help="list of directories with the data", required=True)
parser.add_argument("-l", "--labels", action="extend", nargs="+", type=str, help="list of labels of the data (same order as --dirs)", required=True)
parser.add_argument("-of", "--outfig", help="output file name", required=True)
parser.add_argument("--outfreq", type=int, required=False, help="output frequency of the simulation [number of time steps], if not specified it will be read from const.h5 (if possible)")
args = parser.parse_args()

# directories loop
for directory, lab in zip(args.dirs, args.labels):
  Exy_avg = {}
  nx = {}
  ny = {}
  nz = {}
  dx = {}
  dz = {}
  ref = {}
  lmbd = {}
  level_start_idx = {}
  level_end_idx = {}

  # read some constant parameters
  with h5py.File(directory + "/const.h5", 'r') as consth5:
    user_params = consth5.get("user_params")
    if args.outfreq is None:
      outfreq = int(user_params.attrs["outfreq"][0])
    else:
      outfreq = args.outfreq
    advection = consth5.get("advection")
    dx_adve = advection.attrs["di"] # its the resolved dx
    dz_adve = advection.attrs["dk"] # its the resolved dx
    dt = advection.attrs["dt"]
    nx_adve = consth5["X"][:,:,:].shape[0] - 1
    nz_adve = consth5["Z"][:,:,:].shape[2] - 1
    X = dx_adve * (nx_adve-1)
    Z = dz_adve * (nz_adve-1)

  time_start_idx = int(args.time_start / dt)
  time_end_idx = int(args.time_end / dt)


  # initiliaze nx,ny,nz,dx and average energy for each variable
  for var in args.vars:
    filename = directory + "/timestep" + str(time_start_idx).zfill(10) + ".h5"
    w3d = h5py.File(filename, "r")[var][:,:,:]
    nx[var], ny[var], nz[var] = tuple(x for x in w3d.shape)
    dx[var] = X / (nx[var] - 1)
    dz[var] = Z / (nz[var] - 1) 
    Exy_avg[var] = np.zeros(int((nx[var]-1)/2 + 1))
    ref[var] = int(dx_adve / dx[var])
    assert(float(args.level_start / dz[var]).is_integer())
    assert(float(args.level_end / dz[var]).is_integer())
    level_start_idx[var] = int(args.level_start / dz[var])
    level_end_idx[var] = int(args.level_end / dz[var]) + 1

  
  # time loop
  for t in range(time_start_idx, time_end_idx+1, outfreq):
    filename = directory + "/timestep" + str(t).zfill(10) + ".h5"
    print(filename)
  
    # variables loop
    for var in args.vars:
      print(var)
  
      print(nx[var],dx[var][0])
      w3d = h5py.File(filename, "r")[var][0:nx[var]-1,0:ny[var]-1,:] # * 4. / 3. * 3.1416 * 1e3 
      
      for lvl in range(level_start_idx[var], level_end_idx[var]):
        w2d = w3d[:, :, lvl]
        #print w2d

        wkx = 1.0 / np.sqrt(nx[var] - 1) * np.fft.rfft(w2d, axis = 0)
        wky = 1.0 / np.sqrt(ny[var] - 1) * np.fft.rfft(w2d, axis = 1)
        
#        Ex = 0.5 * (np.abs(wkx) ** 2)
        Ex = (np.abs(wkx) ** 2)
        Ex = np.mean(Ex, axis = 1)
#        Ey = 0.5 * (np.abs(wky) ** 2)
        Ey = (np.abs(wky) ** 2)
        Ey = np.mean(Ey, axis = 0)
        
        Exy = 0.5 * (Ex + Ey)
#        Exy = Ex
        Exy_avg[var] += Exy
  
      K = np.fft.rfftfreq(nx[var] - 1) / dx[var] # assume dy==dx
  #    plt.loglog(K, Exy)
      #lmbd = dx[var] / K
      lmbd[var] = 1 / K
#      print(K, lmbd)
    
    if (t == time_start_idx and lab==args.labels[0]):
      L = np.array([2e2, 2e3])
      plt.loglog(L, 2e-5 * L**(5./3.), label = "-5/3" , color="black", ls='dotted')
      L = np.array([5e1, 3e2])
      plt.loglog(L, 2e-8 * L**(3.), label = "-3" , color="black", ls='dashed')
  
  for var in args.vars:
    Exy_avg[var] /= (time_end_idx - time_start_idx) / outfreq + 1
    Exy_avg[var] /= level_end_idx[var] - level_start_idx[var]

    #crudely scale
    #Exy_avg[var] /= Exy_avg[var][len(Exy_avg[var])-1]
#    Exy_avg[var] /= np.sum(Exy_avg[var])

    #plot
    plt.loglog(lmbd[var], Exy_avg[var] , linewidth=2, label=lab+"_"+var)
 
#plt.xlim(10**4,10**2)
plt.gca().invert_xaxis()
plt.xlabel("l[m]")
plt.ylabel("E")
plt.legend()
plt.grid(True, which='both', linestyle='--')
plt.title("Mean PSD z=["+str(args.level_start)+"m, "+str(args.level_end)+"m] @["+str(args.time_start)+"s, "+str(args.time_end)+"s]")
plt.savefig(args.outfig)
plt.show()
