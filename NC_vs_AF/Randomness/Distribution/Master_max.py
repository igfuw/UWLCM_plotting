from math import exp, log, sqrt, pi, erf, cos, pow, asin, atan, acos, factorial
import numpy as np
from scipy.stats import moment
from sys import argv
#from matplotlib.ticker import MultipleLocator
import matplotlib as plt
import matplotlib.pyplot as plt
#import matplotlib as mpl
#import matplotlib.colors as mcolors
#from matplotlib.ticker import FormatStrFormatter, LogFormatter, LogLocator, LogFormatterSciNotation, AutoMinorLocator
import glob, os
plt.rcParams.update({'font.size': 15})
 
 
 
def read_my_array(file_obj):
    arr_name = file_obj.readline()
    file_obj.readline() # discarded line with size of the array
    line = file_obj.readline()
    line = line.split(" ")
    del line[0]
    del line[len(line)-1]
    arr = list(map(float,line))
    return np.array(arr), arr_name

def read_my_var(file_obj, var_name):
    file_obj.seek(0)
    while True:
        arr, name = read_my_array(file_obj)
        if(str(name).strip() == str(var_name).strip()):
            break
    return arr[-1]

sciezki = '/home/pzmij/2D/PAPER/Master/Master_poszukiwania/times_SD100/'
Slownik = {}
dir_list = os.listdir(sciezki)
Nr_file = [(file.split("_out")[0]).split("_")[-1] for file in dir_list]

# print(Nr_file)
for i in range(len(dir_list)):
    Slownik[Nr_file[i]] = read_my_var(open(sciezki+dir_list[i], "r"), "acc_precip")

# print(Slownik)

myList = Slownik.items()
myList = sorted(myList) 
print(max(Slownik.items(), key = lambda k : k[1])) 
x, y = zip(*myList) 
plt.plot(x, y)
plt.savefig('/home/pzmij/2D/PAPER/Wyniki/Distribution/Master_SD100.png')
# for i in range(len(dir_list)):
#     Slownik[dir_list[i]] = read_my_var(open(sciezki+dir_list[i], "r"), "acc_precip")
Sorted = dict(sorted(Slownik.items(), key=lambda item: item[1]))
print([Sorted[i] for i in Sorted if Sorted[i] >= 0.1]) # prints [5]
# print(Sorted)
