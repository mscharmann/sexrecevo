from __future__ import division
import sys
import warnings
warnings.filterwarnings("ignore")

def Logplotter(logfile_dir, Stat_to_plot = "fst", gen_steps = 1, it_start = 1, it_stop = None):
    
    """
Logplotter - Plots outputs from Sex chromosome recombination evolution simulations

    USAGE: python Logplotter.py <logfile_dir> <Stat_to_plot> [gen_steps] [n_its]:

    logfile_dir    - The absolute path to the directory containing all logfiles from simultions to be plotted
                     (Note all simulations must have identical number of generations and windows)
    Stat_to_plot   - A string to determine the statistic to plot. The string must be in the name of the logfile and 
                     must be unique to that type of logfile. e.g. "fst", "dxy", "XX", "XY". Default = "fst"
    gen_steps      - Controls the number of generations which will be plotted - e.g. if gen_steps = 5, every 5 
                     generational measurements will be plotted. Default = 1 (i.e. All measurements will be plotted)
    it_start       - Iteration to start plotting from (Default: 1)
    it_stop        - Iteration to stop plotting to (Default: All iterations plotted)
    
    """
    
    import os
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx

    ## Get files, generation and window numbers. 

    logfiles = []
    iter_stop = 0
    
    for logfile in os.listdir(logfile_dir):
        if Stat_to_plot in logfile and "pdf" not in logfile and "png" not in logfile:
            logfiles.append("%s/%s" % (logfile_dir, logfile))
            iter_stop += 1

    if it_stop == None:
        it_stop = iter_stop

    logfiles = logfiles[it_start:it_stop]

    ## Get the number of generations and iterations
    # Mathias modified: if some runs are incomplete for the number of generations, find the largest generation that was achieved in any of the runs (files)
    generations = []
    iterations = []
    
    for f in logfiles:
#        print f
        with open(f, "r") as F:
    	    for line in F:
                if not any(["site" in line, "generation" in line]):
                    generations.append(int(line.split()[0]))    			

    generations = sorted(list(set(generations)), key = int)
	
    ## Get the windows

    for line in open(logfiles[0], 'r').readlines():
        if "generation" in line:
            windows = line.split()[2:]
    ## Note - all the files in the simulation directory must have the same number of windows!=


    ## Now get data. Could probably do this easier with pandas, but going to hack it just with base python. 

    Raw_data_by_iteration = {}
    
    for logfile in logfiles:
        if Stat_to_plot in ["dxy", "dxx"]:
            iteration = logfile.rpartition("/")[2].split("_")[1].split(".")[0]
        
        elif Stat_to_plot == "fst":
            iteration = logfile.rpartition("/")[2].split("_")[3].split(".")[0]

        Raw_data_by_iteration[iteration] = {}

        for line in open(logfile, 'r').readlines():

            if not any(["site" in line, "generation" in line]):
                generation = int(line.split()[0])
                data_fields = line.split()[2:]
                field_index = 0

                Raw_data_by_iteration[iteration][generation] = []

                for field in data_fields:

                    if field == "NA":
                        Raw_data_by_iteration[iteration][generation].append(np.nan)
                    else:
                        Raw_data_by_iteration[iteration][generation].append(float(field))
                    field_index += 1


    ### Plotting  ###

#    plt.style.use('dark_background')
    fig = plt.figure(figsize = (20,10))

    ax1 = plt.subplot()

    ## make some pretty colours

    BuPu = plt.get_cmap("BuPu")
    cNorm  = colors.Normalize(vmin=min(generations)-500, vmax=max(generations))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=BuPu)
    
    progress = 0
    for iteration in Raw_data_by_iteration:
        
        # report progress
        progress += 1
        if  progress % 10 == 0:
            print "Plotted %s iterations" % progress

        colour_index = 0
        plot_max = 0
        
        it_generations = sorted(int(i) for i in Raw_data_by_iteration[iteration].keys())
        
        for generation in it_generations: ## use generations list for for loops as its in order. Same for windows below.
            if  colour_index % gen_steps  == 0:  ## use the colour index to define the steps
                gen_data = Raw_data_by_iteration[iteration][generation]
                colorVal = scalarMap.to_rgba(generations[colour_index])
                ax1.plot(gen_data, color=colorVal, label = generation ) #, cmap = "BuPu") 
            colour_index += 1
        
            ## track highest value for manually setting the ylim on plots later. 
            if max(gen_data) > plot_max:
                plot_max = max(gen_data)

    ## Plot formatting ## 
    
    ax1.set_title(Stat_to_plot)
    
    ## remove extra axes
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.patch.set_visible(False)
    
    ## sample 50 tick labels for the x axis, Mathias modified: unless there aren't many windows anyway
    if int(np.round(len(windows)/50)) == 0:
    	xtick_sampler = 1
    else:
    	xtick_sampler = int(np.round(len(windows)/50))
    ax1.set_xticks(range(0,len(windows), xtick_sampler))
    ax1.set_xticklabels(windows[::xtick_sampler],rotation = 30) 
    
    ## turn off unwanted tickmarks on top and right axes
    ax1.tick_params(top="off", right="off") 
    
    ## set the ylim so all plots start from 0
    ax1.set_ylim(0,(plot_max+(plot_max/10)))
    
    ## add axis labels
    ax1.set_ylabel("%s for each iteration" % (Stat_to_plot))
    ax1.set_xlabel("Window")

    ## inset the axis for the colorbar legend and format it
    cb_inset = fig.add_axes([0.65,0.8,0.2,0.2])
    cb_inset.set_title("Generations")
    cb_inset.imshow((generations,generations), cmap=BuPu, extent=[min(generations),max(generations),0,1000])
    cb_inset.axes.get_yaxis().set_visible(False)    

    plt.savefig("%s/%s_summary_plot_ITS_%s-%s.pdf" % (logfile_dir, Stat_to_plot, it_start, it_stop))
    plt.savefig("%s/%s_summary_plot_ITS_%s-%s.png" % (logfile_dir, Stat_to_plot, it_start, it_stop))

    print "All done! Your plots are here: %s/%s_summary_plot_ITS_%s-%s*" % (logfile_dir, Stat_to_plot, it_start, it_stop)


if len(sys.argv) == 6:
    
    working_dir = sys.argv[1]
    stat = sys.argv[2]
    steps = int(sys.argv[3])
    its_start = int(sys.argv[4])
    its_stop = int(sys.argv[5])

    Logplotter(working_dir, stat, steps, its_start, its_stop)

elif len(sys.argv) == 3:
    working_dir = sys.argv[1]
    stat = sys.argv[2]

    Logplotter(working_dir, stat)

else:
    sys.exit("\nERROR: Incorrect number of arguments\n%s" % Logplotter.__doc__)




