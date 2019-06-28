from __future__ import division
import sys
import os
import numpy as np
from matplotlib import pyplot as plt


def Plot_rate(data_directory, target_files = None):

    '''
    USAGE: 
          
          python Plot_rate <DIR> [target_files]

	  <DIR>              The full path to the directory containig the files you want to plot
	  [target_files]     Optional - full path to a file containing a list of files to plot


    Example 1: 
              python Plot_rate $(pwd)
              ## This will plot all files in the current directory. 

    Example 2:
              python Plot_rate $(pwd) /full/path/to/file/target_files.txt
              ## This will plot all files in the current directory that are listed in "target_files.txt"

    '''

    raw_data_dictionary = {}

    ## Can get data files from the directory (all files will be used)
    ## or from a file supplied with names of the data files to be used.

    if isinstance(target_files, str):
        data_files = [i.strip() for i in open(target_files, 'r').readlines()]

    elif target_files == None:
        data_files = os.listdir(data_directory)

    ## Get raw data from all the files in the given directory

    for dat_file_name in data_files:
        if "logfile_nonrec_region" in dat_file_name:
            mut_rate = dat_file_name.split(".")[2]
            recomb_rate = dat_file_name.split(".")[3]
            area = dat_file_name.split(".")[4]
            PopN = dat_file_name.split(".")[5]
            Rep = dat_file_name.split(".")[6]

            key = ".".join([mut_rate, recomb_rate, area, PopN])

            if key not in raw_data_dictionary:
                raw_data_dictionary[key] = {}

            dat_file = open("%s/%s" % (data_directory,dat_file_name), 'r').readlines()[1:]

            for line in dat_file:
                gen = int(line.split()[0])
                size = int(line.split()[1])

                if gen not in raw_data_dictionary[key]:
                    raw_data_dictionary[key][gen] = []

                raw_data_dictionary[key][gen].append(size)  
    
    ## Now average over the iterations for each parameter set
    
    Avgd_data_dict = {}

    for param_set in raw_data_dictionary:
        Avgd_data_dict[param_set] = {}

        for gen in raw_data_dictionary[param_set]:
            Avgd_data_dict[param_set][gen] = np.mean(raw_data_dictionary[param_set][gen])
    
    
    ## Make some plots
    
    plt.figure(figsize = (20,7*len(Avgd_data_dict)))
   

    subplot_index = 1
    colours = ["purple"]

    for key in sorted(Avgd_data_dict.keys()): # so that results are plotted in a meaningful order

        plt.subplot(len(Avgd_data_dict),1,subplot_index)

        gens = sorted(Avgd_data_dict[key].keys())

        x = []
        y = []

        for gen in gens:

            x.append(gen)
            y.append(Avgd_data_dict[key][gen])

        plt.plot(x,y, c = colours[0])
        plt.xlabel("Generation")
        plt.ylabel("Size (Bp) of non-recombging region around the SD")
        plt.title("Mutation rate = %s,  Recomb rate = %s,  Area of effect = %s,  Population size = %s" % 
                  (key.split(".")[0], key.split(".")[1], key.split(".")[2], key.split(".")[3]))
        subplot_index += 1


    outplot_path = "%s/%s" % (data_directory, "Non_recomb_expansion_rates.pdf")
    plt.savefig(outplot_path)
    
    print "\nDone! Your plot is here: %s\n" % outplot_path


if len(sys.argv) < 2:
    sys.exit("\nERROR: Not enough arguments\n %s" % Plot_rate.__doc__)

elif len(sys.argv) == 2:
    data_directory_cline = sys.argv[1]
    
    Plot_rate(data_directory_cline)

elif len(sys.argv) == 3:
    data_directory_cline = sys.argv[1]
    target_files_cline = sys.argv[2]

    Plot_rate(data_directory_cline, target_files_cline)




