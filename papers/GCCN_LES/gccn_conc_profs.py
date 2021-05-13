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
#plt.rcParams.update({'font.size': 24})
# fig size
#plt.figure(figsize=(20,10))

fig, axarr = plt.subplots(1,2)

linestyles = ['--', '-.', ':']
dashList = [(3,1),(1,1),(4,1,1,1),(4,2)]

# plot init conc
#hgt = np.arange(0, 1.1, 0.01)
#init_conc = hgt.copy()
#init_conc[hgt < 1] = 2.8
#init_conc[hgt >= 1] = 0
#plt.plot(init_conc, hgt)

for no in np.arange(3600, 36001, 3600):
  print(no)
  file_name = sys.argv[2] + "profiles_" + str(no) + "_" + str(no) + ".dat"
  profs_file = open(file_name, "r")
  my_hgt = read_my_var(profs_file, "position")
 # my_gccn_conc = read_my_var(profs_file, "rd_geq_0.8um_conc")
 # my_ccn_conc = read_my_var(profs_file, "rd_lt_0.8um_conc")
  my_gccn_conc = read_my_var(profs_file, "gccn_conc")
  my_ccn_conc = read_my_var(profs_file, "non_gccn_conc")
  profs_file.close()
  axarr[0].plot(my_gccn_conc, my_hgt)#, label=sys.argv[no+1])
  axarr[1].plot(my_ccn_conc, my_hgt)#, label=sys.argv[no+1])

## Dycoms rd<2um
#for no in np.arange(3600, 18001, 3600):
#  file_name = "/home/piotr/praca/GCCN_LES/wyniki/Dycoms_RF02/profs_non_gccn_conc/11_12_out_UWLCM_dycoms_sgs_GCCNx10_TurbAdve_MixLenFix_129x129x301_dt1_SstpCond10_SstpCoal10_sd100_dt1_RE1_out_lgrngn_dycoms_profiles_of_non_gccn_conc/11_12_out_UWLCM_dycoms_sgs_GCCNx10_TurbAdve_MixLenFix_129x129x301_dt1_SstpCond10_SstpCoal10_sd100_dt1_RE1_out_lgrngn_dycoms_profiles_" + str(no) + "_" + str(no) + ".dat"
#  profs_file = open(file_name, "r")
#  my_hgt = read_my_var(profs_file, "position")
#  my_non_gccn_conc = read_my_var(profs_file, "non_gccn_conc")
#  profs_file.close()
#  axarr[1].plot(my_non_gccn_conc, my_hgt)#, label=sys.argv[no+1])


#plt.legend()
#axarr[0].set_ylim(0,1.5)
#axarr[0].set_xlim(0,3)
#axarr[1].set_ylim(0,1.5)
#axarr[1].set_xlim(0,250)


#axarr[0].set_ylabel('height / inversion height')
axarr[0].set_ylabel('height [m]')
axarr[0].set_xlabel('GCCN ($r_d \geq 2 \mu m$) conc. [cm$^{-3}$]')
axarr[1].set_xlabel('CCN ($r_d < 2 \mu m$) conc. [cm$^{-3}$]')

plt.savefig(sys.argv[1])
