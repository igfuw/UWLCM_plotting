from math import exp, log, sqrt, pi, erf, cos, pow, asin, atan, acos, factorial
import numpy as np
import pandas as pd
from scipy.stats import moment
from sys import argv
from matplotlib.ticker import MultipleLocator
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors as mcolors
from matplotlib.ticker import FormatStrFormatter, LogFormatter, LogLocator, LogFormatterSciNotation, AutoMinorLocator
import glob, os
plt.rcParams.update({'font.size': 20})

############################################################################
####  DO RECZNEGO UZUPELNIENIA#############
name = '2D_NA1'
#PATHS

path_diff_NA1_piggy_VF = '/home/piotr-pc/Piotr/WORKSHOPS/2D/NA1/time_series_diff_piggy_VF/'
path_diff_NA1_piggy_VF_tail = '/home/piotr-pc/Piotr/WORKSHOPS/2D/NA1/time_series_diff_piggy_VF_tail/'
path_diff_NA1_piggy_VF_multi = '/home/piotr-pc/Piotr/WORKSHOPS/2D/NA1/time_series_diff_piggy_VF_tail/'

outfile = ''
#OPISY
text_diff_piggy_VF = 'same velo field in all simulations'
text_diff_piggy_VF_tail = 'same velo field in all simulations \n sd conc long tail'
text_diff_piggy_VF_multi = 'same velo field in all simulations \n sd conc same multiplicity'

paths = [path_diff_NA1_piggy_VF+'SD100',path_diff_NA1_piggy_VF+'SD1000', path_diff_NA1_piggy_VF+'SD10000',path_diff_NA1_piggy_VF_tail+'SD100',path_diff_NA1_piggy_VF_tail+'SD1000', path_diff_NA1_piggy_VF+'SD10000']

label_list = ['100', '1000', '10000']
podpisy = [text_diff_piggy_VF, text_diff_piggy_VF_tail]#, text_diff_piggy_VF_multi]
width_multiplier = 0.67
##########################################################################
def Rysuj_to(sciezki, etykiety, podpisy, name):

    label = etykiety[0:int(len(sciezki)/2)]*len(podpisy)
    multi = len(podpisy)
    labels = []
    X = []
    for i in range(len(label)):
        if i < len(label)/multi :
            labels.append(podpisy[0])
            X.append(1)
        elif  i < len(label)/multi*2:
            labels.append(podpisy[1])
            X.append(2)
        else:
            labels.append(podpisy[2])
            X.append(3)

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
        return arr
    def licz_srednia(parameter_name, iter_value, sciezki):
        dl = len(series_file[iter_value])
        srednia =[0 for i in range(len(sciezki))]
        STD = [0 for i in range(len(sciezki))]
        error = [0 for i in range(len(sciezki))]
        STD_error_mean = [0 for i in range(len(sciezki))]
        Zmienna = np.zeros((len(series_file[iter_value]),len(read_my_var(series_file[iter_value][0], str(parameter_name)))))
        for j in range(len(series_file[iter_value])):
            Zmienna[j] = read_my_var(series_file[iter_value][j], str(parameter_name))
        srednia[iter_value] = Zmienna.mean(0)
        STD[iter_value] = Zmienna.std(0)
        STD_error_mean[iter_value] = Zmienna.std(0)/sqrt(dl)
        error[iter_value] = np.power(1/dl * (moment(Zmienna,4) - (dl-3)/(dl-1)*np.power(STD[iter_value],4)),1/2)/(2*STD[iter_value])
        #Error from this https://stats.stackexchange.com/questions/156518/what-is-the-standard-error-of-the-sample-standard-deviation
        return srednia[iter_value], STD[iter_value], error[iter_value], STD_error_mean[iter_value]

    def czas(iter_value):
        for j in range(len(series_file[iter_value])):
            time = read_my_var(series_file[iter_value][j], "position")
        return time
    files = [0 for i in range(len(sciezki))]
    series_file = [0 for i in range(len(sciezki))]
    for p in range(len(sciezki)):
        os.chdir(sciezki[p])
        series_file[p] = [open(file_names, "r") for file_names in glob.glob("*.dat")]
        files[p] = glob.glob("*.dat")

    colors = [ 'forestgreen', 'gold', 'forestgreen', 'gold', 'forestgreen', 'gold']
    u_init = np.linspace(-int(len(label)/len(podpisy)), int(len(label)/len(podpisy)), int(len(label)/len(podpisy)))
    u = np.tile(u_init, len(podpisy))
    width = 0.25
    colors_list = [ 'thistle', 'orchid','red', 'grey', 'green', 'orange', 'blue', 'yellow']
    colors = colors_list[0:len(etykiety)] *len(podpisy)
    multi = len(etykiety)


    fig1, ax = plt.subplots()
    fig1.set_size_inches(18.5, 10.5)
    for p in range( len(sciezki)):
        srednia = licz_srednia("acc_vol_precip", p, sciezki)[0]
        STD = licz_srednia("acc_vol_precip", p, sciezki)[1]
        error = licz_srednia("acc_vol_precip", p, sciezki)[2]
        A  = plt.bar(X[p] + u[p]*width/multi, STD[-1]/srednia[-1]*100,  color=colors[p],width =width*width_multiplier,  label=(label[p]) if (p < multi) else "")
        plt.ylabel(r"$\frac{\sigma(acc \hspace{0.5} precip)}{mean(acc \hspace{0.5} precip)}$ [%]")
        plt.xticks(X, labels, ha = 'center')
        plt.ylim(0, 2.5e2)
        plt.text(0.5, 265  , 'G', fontsize=26)
        plt.title("Standard deviation of accumulated precipitation to \n mean of accumulated precipitation in a cumulus congestus simulation",weight='bold')
        plt.legend(title="number of super-droplets per cell: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)

        for rect in A:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                '%d' % int(height) + "%", ha='center', va='bottom')
    plt.savefig(outfile+'sigma_VF_tail_multi'+name+'.png')

    fig1, ax = plt.subplots()
    fig1.set_size_inches(18.5, 10.5)
    for p in range( len(sciezki)):
        srednia = licz_srednia("acc_vol_precip", p, sciezki)[0]
        STD = licz_srednia("acc_vol_precip", p, sciezki)[1]
        error = licz_srednia("acc_vol_precip", p, sciezki)[2]
        A  = plt.bar(X[p] + u[p]*width/multi, STD[-1],  color=colors[p],width =width*width_multiplier,  yerr=error[-1], label=(label[p]) if (p < multi) else "")
        plt.ylabel(r"$\sigma$ [$m^3$]")
        plt.ylim(0, 5e-1)
        plt.text(0.5, 0.53  , 'F', fontsize=26)
        plt.xticks(X, labels, ha = 'center')
        plt.title("        Standard deviation of accumulated precipitation \n in a cumulus congestus simulation", weight='bold')
        plt.legend(title="number of super-droplets per cell: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)

        for rect in A:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                f"{height:.2e}", ha='center', va='bottom')
                # f"{height:.2e}" + "$m^3$", ha='center', va='bottom')
    plt.savefig(outfile+'STD_VF_tail_multi'+name+'.png')

    fig1, ax = plt.subplots()
    fig1.set_size_inches(18.5, 10.5)
    for p in range( len(sciezki)):
        srednia = licz_srednia("acc_vol_precip", p, sciezki)[0]
        STD = licz_srednia("acc_vol_precip", p, sciezki)[1]
        error = licz_srednia("acc_vol_precip", p, sciezki)[2]
        STD_error_mean = licz_srednia("acc_vol_precip", p, sciezki)[3]
        A  = plt.bar(X[p] + u[p]*width/multi, srednia[-1],  color=colors[p],width =width*width_multiplier,  yerr=STD_error_mean[-1], label=(label[p]) if (p < multi) else "")
        plt.ylabel(r"Mean accumulated precipitation [$m^3$]")
        plt.ylim(0, 4e-1)
        plt.text(0.5, 0.43  , 'E', fontsize=26)
        plt.xticks(X, labels, ha = 'center')
        plt.title("        Mean of accumulated precipitation in a cumulus congestus simulation", weight='bold')
        plt.legend(title="number of super-droplets per cell: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)

        for rect in A:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                f"{height:.2e}", ha='center', va='bottom')
                # f"{height:.2e}" + "$m^3$", ha='center', va='bottom')
    plt.savefig(outfile+'Mean_VF_tail_multi'+name+'.png')


Rysuj_to(paths, label_list, podpisy, name)