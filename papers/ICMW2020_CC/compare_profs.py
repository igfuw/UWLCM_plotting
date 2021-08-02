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

profs = ["cl_nc", "actrw_reff_cl", "cloud_std_dev"]
xlabels = {
           "cl_nc": "cloud droplet concentration [1/cc]",
           "actrw_reff_cl": "effective radius of cloud droplets [um]",
           "cloud_std_dev": "std. dev. of radius of cloud droplets [um]"
          }
xlims = {
           "cl_nc": [],
           "actrw_reff_cl": [],
           "cloud_std_dev": []
        }

# activate latex text rendering
rc('text', usetex=True)
# font size
plt.rcParams.update({'font.size': 24})
# fig size
plt.figure(figsize=(10,10))

linestyles = ['--', '-.', ':']
dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]

file_no = (len(sys.argv) - 2.) / 3


for prof in profs:
  print( prof)

  for no in np.arange(2, 3*file_no+2, 3):
    ino = int(no)
    time = sys.argv[ino+2]
    print (ino,time)
    file_name = sys.argv[ino] + "_" + str(time) + "_" + str(time) + ".dat"
    print (file_name)
    profs_file = open(file_name, "r")
    my_hgt = read_my_var(profs_file, "position")
    my_prof = read_my_var(profs_file, prof)
    
    profs_file.close()
    
    plt.plot(my_prof, my_hgt, label=sys.argv[ino+1]+" @ t="+str(time)+"s", lw=2)
  
  plt.ylabel('height [m]')
  plt.xlabel(xlabels[prof])
  plt.legend(loc='upper right')
#  plt.xlim(xlims[ser][0], xlims[ser][1])
  plt.ylim(0, 8000)
  plt.savefig(sys.argv[1]+prof+".png")
  plt.clf()
