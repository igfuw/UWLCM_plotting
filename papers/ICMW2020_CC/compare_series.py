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
from plot_series import *

series = ["lwp", "cloud_top_height", "total_droplets_number", "acc_vol_precip"]
ylabels = {
           "lwp": "LWP [g/ m$^2$]",
           "cloud_top_height" : "cloud top [m]",
           "total_droplets_number" : "total number of cloud droplets [1]",
           "acc_vol_precip" : "accumulated precipitation volume [?]"
          }
ylims = {
         "lwp" : [0,162],
         "cloud_top_height" : [0,10000],
         "total_droplets_number" : [0,1e19],
          "acc_vol_precip" : [0,3e3]
        }

# activate latex text rendering
rc('text', usetex=True)
# font size
plt.rcParams.update({'font.size': 24})
# fig size
plt.figure(figsize=(20,10))

linestyles = ['--', '-.', ':']
dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]


file_no = (len(sys.argv) - 2.) / 2

for ser in series:

  for no in np.arange(2, 2*file_no+2, 2):
    ino = int(no)
    file_name = sys.argv[ino] + "series.dat"
    series_file = open(file_name, "r")
    my_time = read_my_var(series_file, "position")
    my_ser = read_my_var(series_file, ser)
    
    series_file.close()
    
    plt.plot(my_time / 3600., my_ser, label=sys.argv[ino+1], lw=2)
  
  plt.xlabel('time [h]')
  plt.ylabel(ylabels[ser])
  plt.legend(loc='upper left')
  plt.ylim(ylims[ser][0], ylims[ser][1])
  plt.savefig(sys.argv[1]+ser+".png")
  plt.clf()
