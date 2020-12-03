from matplotlib import rc
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from matplotlib.ticker import FormatStrFormatter, NullFormatter
from matplotlib.ticker import MaxNLocator

import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")
#sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Dycoms_RF02/")
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../cases/RICO11/")

from plot_ranges import * 
from plot_profs import *

# activate latex text rendering
rc('text', usetex=True)
# font size
plt.rcParams.update({'font.size': 24})
# fig size
plt.figure(figsize=(20,10))

# profs = ["gccn_rw_cl_down"]
# nplotx=1
# nploty=1
# 
# # init the plot
# fig, axarr = plt.subplots(nplotx, nploty)
# print 'axarr type: ', type(axarr)
# print 'axarr shape: ', axarr.shape
# plot_iter=0
# 
# assert len(sys.argv) == 6
# 
# 
# n_profs=3
# file_names = []
# file_labels = []
# file_no = np.arange(1, 2 * n_profs , 2) # 2* because for each line there is data file and label
# for no in file_no:
#   print no
#   print sys.argv[no]
#   file_names.append(sys.argv[no] + "profiles_7200_18000.dat")
#   file_labels.append(sys.argv[no+1])
# 
# print file_names
# plot_iter = plot_profiles(profs, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict_series[cusc_iter], ylimdict_series, xlabel='Time [h]', file_names=file_names, file_labels=file_labels, linewidth=1.5)

linestyles = ['--', '-.', ':']
dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]

plt.xlabel('mean wet radius of droplets formed on aerosols with $r_\mathrm{dry}>2\mu\mathrm{m}$ [$\mu\mathrm{m}$]')
plt.ylabel('height / inversion height')

for no in [2,4,6]:
  file_name = sys.argv[no] + "profiles_7200_18000.dat"
  profs_file = open(file_name, "r")
  my_hgt = read_my_var(profs_file, "position")
  my_gccn_rw = read_my_var(profs_file, "gccn_rw_cl_down")
  
  profs_file.close()
  
  plt.plot(my_gccn_rw, my_hgt, label=sys.argv[no+1])

plt.legend()
plt.ylim(0,1.1)
plt.savefig(sys.argv[1])
