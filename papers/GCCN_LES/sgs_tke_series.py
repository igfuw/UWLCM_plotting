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

# activate latex text rendering
rc('text', usetex=True)
# font size
plt.rcParams.update({'font.size': 18})
# fig size
plt.figure(figsize=(20,10))

file_name = sys.argv[1] + "series.dat"

series_file = open(file_name, "r")
my_times = read_my_var(series_file, "position")
my_sgs_tke = read_my_var(series_file, "sgs_tke")
my_sgs_tke_sd = read_my_var(series_file, "sgs_tke_sd")

# rescale time to hours
my_times = my_times / 3600.

series_file.close()

linestyles = ['--', '-.', ':']
dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]

plt.xlabel('time [h]')
plt.ylabel('SGS TKE [m$^2$ s$^{-2}$]')

plt.plot(my_times, my_sgs_tke, label='Smagorinsky')
plt.plot(my_times, my_sgs_tke_sd, label='Superdroplets')
plt.legend()
plt.savefig(sys.argv[2])
