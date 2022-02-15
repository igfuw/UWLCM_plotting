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

############################################################################
####  DO RECZNEGO UZUPELNIENIA#############

label_list = ['100', '1000', '10000']
paths = ['/home/pzmij/2D/PAPER/Distribution/Piggy/SD100/times_classic/', '/home/pzmij/2D/PAPER/Distribution/Piggy/SD1000/times_classic/','/home/pzmij/2D/PAPER/Distribution/Piggy/SD10000/times_classic/',
        '/home/pzmij/2D/PAPER/Distribution/Piggy_2/SD100/times_classic/', '/home/pzmij/2D/PAPER/Distribution/Piggy_2/SD1000/times_classic/','/home/pzmij/2D/PAPER/Distribution/Piggy_2/SD10000/times_classic/',
        '/home/pzmij/2D/PAPER/Distribution/Piggy_3/SD100/times_classic/', '/home/pzmij/2D/PAPER/Distribution/Piggy_3/SD1000/times_classic/','/home/pzmij/2D/PAPER/Distribution/Piggy_3/SD10000/times_classic/',
        '/home/pzmij/2D/PAPER/Distribution/Piggy/SD100/times_tail/', '/home/pzmij/2D/PAPER/Distribution/Piggy/SD1000/times_tail/', '/home/pzmij/2D/PAPER/Distribution/Piggy/SD10000/times_tail/',
        '/home/pzmij/2D/PAPER/Distribution/Piggy_2/SD100/times_tail/', '/home/pzmij/2D/PAPER/Distribution/Piggy_2/SD1000/times_tail/', '/home/pzmij/2D/PAPER/Distribution/Piggy_2/SD10000/times_tail/',
        '/home/pzmij/2D/PAPER/Distribution/Piggy_3/SD100/times_tail/', '/home/pzmij/2D/PAPER/Distribution/Piggy_3/SD1000/times_tail/', '/home/pzmij/2D/PAPER/Distribution/Piggy_3/SD10000/times_tail/'] 
name = 'Piggy_3_Distribution'
text_diff_piggy = 'Classic_M1'
text_diff_piggy2 = 'Classic_M2'
text_diff_piggy3 = 'Classic_M3'
text_diff_piggy4 = 'Tail_M1'
text_diff_piggy5 = 'Tail_M2'
text_diff_piggy6 = 'Tail_M3'
podpisy = [text_diff_piggy, text_diff_piggy2, text_diff_piggy3, text_diff_piggy4, text_diff_piggy5, text_diff_piggy6]
outfile = '/home/pzmij/2D/PAPER/Wyniki/Distribution/barrs/'
width_multiplier = 0.57
##########################################################################
def Rysuj_to(sciezki, etykiety, podpisy, name):

    if len(podpisy) == 1 :
        label = etykiety
    else:
        label = etykiety[0:int(len(sciezki)/2)]*len(podpisy)
    multi = len(podpisy)
    #Y = [i+1 for i in range(int(len(label)/len(podpisy)+1))]
    Y = [i+1 for i in range(multi)]
    X = np.repeat(Y, int(len(etykiety)))
    #labels = [podpisy[i] for i in range(int(len(label)/len(podpisy)+1)) ]
    labels = np.repeat(podpisy, int(len(etykiety)))
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
        error[iter_value] = np.power(1/dl * (moment(Zmienna,4) - (dl-3)/(dl-1)*np.power(STD[iter_value],2)),1/2)/(2*np.sqrt(STD[iter_value]))
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
    width = 0.20
    colors_list = [ 'thistle', 'orchid','red', 'grey', 'green', 'orange', 'blue', 'yellow']
    colors = colors_list[0:len(etykiety)] *len(podpisy)
    multi = len(etykiety)

    print("X", X, "u", len(u), "col", len(colors))
    fig1, ax = plt.subplots()
    fig1.set_size_inches(18.5, 10.5)
    for p in range( len(sciezki)):
        srednia = licz_srednia("acc_vol_precip", p, sciezki)[0]
        STD = licz_srednia("acc_vol_precip", p, sciezki)[1]
        A  = plt.bar(X[p] + u[p]*width/multi, STD[-1]/srednia[-1]*100,  color=colors[p],width =width*width_multiplier,  label=(label[p]) if (p < multi) else "")
        plt.ylabel(r"$\frac{\sigma(acc \hspace{0.5} precip)}{mean(acc \hspace{0.5} precip)}$ [%]")
        plt.xticks(X, labels, ha = 'center')
        plt.ylim(0, 6e2) #piggy/
        # plt.ylim(0, 3.5e2) #no_piggy/
        #plt.text(0.5, 265  , 'G', fontsize=26)
        plt.title("Standard deviation of accumulated precipitation to \n mean of accumulated precipitation in a cumulus congestus simulation",weight='bold')
        # plt.legend(title="number of super-droplets per cell: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)
        plt.legend(title="SD#: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)


        for rect in A:
            height = rect.get_height()
            ax.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                '%d' % int(height) + "%", ha='center', va='bottom')
    plt.savefig(outfile+'STD_to_mean_'+name+'.png')

    fig2, ax2 = plt.subplots()
    fig2.set_size_inches(18.5, 10.5)
    for p in range( len(sciezki)):
        STD = licz_srednia("acc_vol_precip", p, sciezki)[1]
        error = licz_srednia("acc_vol_precip", p, sciezki)[2]
        A  = plt.bar(X[p] + u[p]*width/multi, STD[-1],  color=colors[p],width =width*width_multiplier,  yerr=error[-1], label=(label[p]) if (p < multi) else "")
        plt.ylabel(r"$\sigma$ [$m^3$]")
        plt.ylim(0, 1.3e-1) #piggy
        # plt.ylim(0, 8e-1) #no_piggy
        #plt.text(0.5, 0.53  , 'F', fontsize=26)
        plt.xticks(X, labels, ha = 'center')
        plt.title("        Standard deviation of accumulated precipitation \n in a cumulus congestus simulation", weight='bold')
        # plt.legend(title="number of super-droplets per cell: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)
        plt.legend(title="SD# ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)

        for rect in A:
            height = rect.get_height()
            ax2.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                f"{height:.2e}", ha='center', va='bottom')
                # f"{height:.2e}" + "$m^3$", ha='center', va='bottom')
    plt.savefig(outfile+'STD_'+name+'.png')

    fig3, ax3 = plt.subplots()
    fig3.set_size_inches(18.5, 10.5)
    for p in range( len(sciezki)):
        srednia = licz_srednia("acc_vol_precip", p, sciezki)[0]
        STD_error_mean = licz_srednia("acc_vol_precip", p, sciezki)[3]
        A  = plt.bar(X[p] + u[p]*width/multi, srednia[-1],  color=colors[p],width =width*width_multiplier,  yerr=STD_error_mean[-1], label=(label[p]) if (p < multi) else "")
        plt.ylabel(r"Mean accumulated precipitation [$m^3$]")
        plt.ylim(0, 43e-2) #piggy
        # plt.ylim(0, 65e-2) #no_piggy
        #plt.text(0.5, 0.43  , 'E', fontsize=26)
        plt.xticks(X, labels, ha = 'center')
        plt.title("        Mean of accumulated precipitation in a cumulus congestus simulation", weight='bold')
        # plt.legend(title="number of super-droplets per cell: ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)
        plt.legend(title="SD# ", prop=dict(weight='bold') ,loc='upper left', framealpha = 0.5, frameon = 0, ncol = multi)


        for rect in A:
            height = rect.get_height()
            ax3.text(rect.get_x() + rect.get_width()/2., 0.99*height,
                f"{height:.2e}", ha='center', va='bottom')
                # f"{height:.2e}" + "$m^3$", ha='center', va='bottom')
    plt.savefig(outfile+'Mean_'+name+'.png')


Rysuj_to(paths, label_list, podpisy, name)
