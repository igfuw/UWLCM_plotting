import numpy as np
from math import floor
#from sys import argv
import argparse
import matplotlib.pyplot as plt

from latex_labels import var_labels
from read_UWLCM_arrays import *


def plot_series(var_list, plot_iter, nplotx, nploty, axarr, xscaledict, yscaledict, xlimdict, ylimdict, show_bin=False, suffix='', xlabel='', ylabeldict=var_labels, file_names=[], file_labels=[], linewidth=1):

  parser = argparse.ArgumentParser(description='Plot UWLCM series comparison')
  parser.add_argument("-ts", "--time_start", type=float, required=False, help="start of the plotted period [s] (override default)")
  parser.add_argument("-te", "--time_end", type=float, required=False, help="end of the plotted period [s] (override default)")
  parser.add_argument("-d", "--dirs", action="extend", nargs="+", type=str, help="list of directories with the data", required=True)
  parser.add_argument("-l", "--labels", action="extend", nargs="+", type=str, help="list of labels of the data (same order as --dirs)", required=True)
  args, extra = parser.parse_known_args()

  if args.time_start is not None and args.time_end is not None:
    xlimdict = {x: (args.time_start, args.time_end) for x in xlimdict}


  # if file names are not defined, read them and labels from command line
  if len(file_names)==0:
    for directory, lab in zip(args.dirs, args.labels):
      file_names.append(directory + suffix)
      file_labels.append(lab)

  for var in var_list:
    label_counter=0
    for file_name in file_names:
      print(file_name, var)
      series_file = open(file_name, "r")
      my_times = read_my_var(series_file, "position")
      my_res = read_my_var(series_file, var)
      if len(my_res) == 0: # file does not contain this type of plot
        print("skipping from: " + str(label_counter))
        label_counter+=1
        print("skipping to: " + str(label_counter))
        continue
      
      # rescale time to hours
      my_times = my_times / 3600.

      # rescale autoconv. and accr. rates to g/(m^3 * day)
      if var == "cl_accr25_rico" or var == "cl_acnv25_rico"  or var == "cl_acnv25_dycoms" or var == "cl_accr25_dycoms":
        my_res *= 3600 * 24
      
      series_file.close()
  
      linestyles = ['--', '-.', ':']
      dashList = [(3,1),(1,1),(4,1,1,1),(4,2)] 
  #    colorList = ['red', 'blue', 'green']
      colorList = plt.rcParams['axes.prop_cycle'].by_key()['color'] # default prop cycle colors

      # x label only on he lowest row
      xlabel_used = xlabel
      if plot_iter < nploty:
        xlabel_used = ''

      plot_my_array(axarr, plot_iter, my_times, my_res, nploty, xlabel=xlabel_used, ylabel=ylabeldict[var], varlabel=file_labels[label_counter], dashes = dashList[label_counter % len(dashList)], xscale=xscaledict[var], yscale=yscaledict[var], xlim=xlimdict[var], ylim=ylimdict[var], linewidth=linewidth, color = colorList[label_counter % len(colorList)])
      label_counter+=1
    plot_iter = plot_iter + 1
  return plot_iter
