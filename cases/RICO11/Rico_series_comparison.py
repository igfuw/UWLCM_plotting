from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
import argparse


import os
import sys
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/../../Matplotlib_common/")

from plot_ranges import xscaledict, yscaledict, xlimdict_series, ylimdict_series
from plot_series import *
from latex_labels import labeldict
from autoscale_y import *


# activate latex text rendering
rc('text', usetex=True)

rico_vars = ["lwp", "rwp", "cloud_cover", "min_cloud_base", "inversion_height_rico", "cl_nc", "cl_nr", "surf_precip", "acc_precip", "RH_max", "wvarmax", "sd_conc", "cl_avg_cloud_meanr", "cl_avg_cloud_stddevr", "cl_avg_supersat", "cl_avg_th", "cl_avg_rv", "surf_flux_latent", "surf_flux_sensible"]#, "cl_acnv25", "cl_accr25"]# ]
#rico_vars = ["clb_bigrain_mean_rd","clb_bigrain_mean_kappa","clb_bigrain_mean_conc","clb_bigrain_mean_inclt", "cl_nr"]

# variables that need rescaling of the yrange to the limited x range of 1-6h
rescale_vars = ["lwp", "cloud_cover", "min_cloud_base", "inversion_height_rico", "cl_nc"]# rico_vars

# arguments
parser = argparse.ArgumentParser(description='Plot UWLCM series comparison for RICO simulations.')
parser.add_argument("-of", "--outfig", help="output file name", required=True)
args, extra = parser.parse_known_args()


# init the plot
nplotx = 5
nploty= 4
fig, axarr = plt.subplots(nplotx,nploty)
x_arr = np.arange(nplotx)
y_arr = np.arange(nploty)

if len(rico_vars) % nploty == 0:
  nemptyplots = 0
else:
  nemptyplots = nploty - len(rico_vars) % nploty
emptyplots = np.arange(nploty - nemptyplots, nploty)

plot_series(rico_vars, 0, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict_series, ylimdict_series, xlabel='Time [h]')

# legend font size
plt.rcParams.update({'font.size': 8})

# hide axes on empty plots
for empty in emptyplots:
  axarr[nplotx-1, empty].axis('off')

#axes = plt.gca()
#axes.tick_params(direction='in')
for x in x_arr:
  for y in y_arr:
    #tics inside
    axarr[x,y].tick_params(direction='in', which='both', top=1, right=1)
    #minor tics
    axarr[x,y].xaxis.set_minor_locator(AutoMinorLocator())
    axarr[x,y].yaxis.set_minor_locator(AutoMinorLocator())
    #labels and tics font size
    for item in ([axarr[x,y].xaxis.label, axarr[x,y].yaxis.label] + axarr[x,y].get_xticklabels() + axarr[x,y].get_yticklabels()):
      item.set_fontsize(8)
    # subplot numbering
    if y < nploty - nemptyplots or x < (nplotx - 1):
      axarr[x,y].text(0.2, 0.875, labeldict[y + x*nploty], fontsize=8, transform=axarr[x,y].transAxes)

      # rescale y range to the visible x range, note: overrides ylim!
      var = rico_vars[x*nploty + y]
      if var in rescale_vars:
        autoscale_y(axarr[x,y], margin=0.1)
      if(var == "cloud_cover" or var == "cl_nr"):
        axarr[x,y].yaxis.set_major_formatter(ticker.FormatStrFormatter('%0.2f'))

      # hide hrzntl tic labels
      if x*nploty + y < nplotx * nploty - nemptyplots - nploty:
        axarr[x,y].set_xticklabels([])


#single legend for the whole figure
handles, labels = axarr[0,0].get_legend_handles_labels()

lgd = fig.legend(handles, labels, handlelength=4, loc='lower center', bbox_to_anchor=(0.45,0))

#figure size
fig.set_size_inches(7.874, 1.5 * nplotx + (len(labels) - 2) * 0.1)# 5.214)#20.75,13.74)
#distances between subplots and from bottom of the plot
fig.subplots_adjust(bottom=0.18 + (len(labels) - 2) * 0.01, hspace=0.1, wspace=0.4)

#fig.tight_layout(pad=0, w_pad=0, h_pad=0)

#plt.show()
fig.savefig(args.outfig, bbox_inches='tight', dpi=300)#, bbox_extra_artists=(lgd,))
