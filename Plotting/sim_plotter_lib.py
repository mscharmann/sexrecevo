## For cluster only:

import matplotlib
matplotlib.use('Agg')
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)


###############################################################################################################################################################################################


def rv_to_perbase_avg(rv_file_path, n_gams = 100):

    """
    Generate average per-base recombination probability vectors for each gamete in the r output file and average across a generation.

    """

    import numpy as np

    iteration = rv_file_path.split(".")[0].rpartition("_")[2]
    print(rv_file_path)

    with open(rv_file_path, 'r') as rv_handle:

        print("Generation: ", end = " ")

        pergen_perbase_averages = {}

        for line in rv_handle:

            if line.startswith("#"): ## record by generation
                GEN = line.split()[1]
                gamcount = 0

                if int(GEN) % 100 == 0:
                    print(" %s," % GEN, end = " ")

            else:
                gamcount += 1
                blocks = line.split()
                gamete_vector = []  ## make a per-base vector for each gamete

                for block in blocks:
                    STRT = int(block.split(":")[0].split("-")[0])
                    END = int(block.split(":")[0].split("-")[1])
                    BLOCK_PROB = float(block.split(":")[1])
                    if (END - (STRT-1)) < 1:
                        PERBASE_PROB = BLOCK_PROB
                    else:
                        PERBASE_PROB = BLOCK_PROB/(END - (STRT-1))
                    
                    if PERBASE_PROB == 0:
                        PERBASE_PROB = 0.0000000001  ## this allows to plot log later while still keeping these values at essentially zero

                    gamete_vector.extend([PERBASE_PROB]*(END - (STRT-1)))

                if gamcount == 1: ## if its the first gamete of the generation
                    per_base_sum = gamete_vector # make per_base_sum the gamete vector

                else: ## if its not the first gamete, add it to the existing sum
                    per_base_sum = np.sum([per_base_sum, gamete_vector], axis = 0)

                if gamcount == n_gams:  ## divide gamete per base sum by number of gametes
                    pergen_perbase_averages[int(GEN)] = per_base_sum/n_gams

                    
        rv_avgs_outpath = "%s_AVGS.out" % rv_file_path.rpartition(".")[0]
        
        with open(rv_avgs_outpath, 'w') as rv_avgs_out:
            for gen in pergen_perbase_averages:
                rv_avgs_out.write("%s\t" % gen)
                for i in pergen_perbase_averages[gen]:
                    rv_avgs_out.write("%s " % i)
                rv_avgs_out.write("\n")
        
        print("Iteration %s done" % iteration)
        #return pergen_perbase_averages



###############################################################################################################################################################################################

def Recomb_prob_plot_per_base_per_iteration(pergen_perbase_avgs, outpath, genome_size = 100000, N_gametes = "", tick_sampler = 20):

    """
    
    Plot the per-generaton per-base averaged recombination probabilites for 
    a single simulation (iteration).
    
    <pergen_perbase_avgs>    -    The averaged (across gametes) recomb prob vectors per generation and per base for the simulation
    <outpath>                -    The name of the file to which to save plot. Suffix (e.g. ".png", ".pdf" will determine the plot format)
    <genome_size>            -    The size of the genome simulated
    <N_gametes>              -    The number of gametes in the simualtion
    <tick_sampler>           -    The number of xtick positions to plot/label
    
    
    """
    
    
    from matplotlib import pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import matplotlib.colorbar as cb
    import numpy as np

    generations = sorted([int(i) for i in pergen_perbase_avgs.keys()])
    positions = list(range(1,(genome_size+1)))

    #generations = [str(i) for i in generations]

    ### Now for the plotting  ###
    fig = plt.figure(figsize = (20,10))

    ax1 = plt.subplot()

    ## make some pretty colours

    BuPu = plt.get_cmap("BuPu")
    cNorm  = colors.Normalize(vmin=min(generations)-500, vmax=max(generations))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=BuPu)

    colour_index = 0
    plot_max = 0


    for gen in generations: ## use generations list for iteration as its in order. Same for windows below.

        gen_data = pergen_perbase_avgs[gen]

        colorVal = scalarMap.to_rgba(generations[colour_index])

        ax1.plot(gen_data, color=colorVal, label = gen ) 
        
        if gen == 0:
            ax1.plot(gen_data, color="grey", dashes = [5, 10]) ## plot starting prob as dashed line
        
        colour_index += 1

        ## track highest value for mannually setting the ylim on plots later. 
        if max(gen_data) > plot_max:
            plot_max = max(gen_data)

    ## Plot formatting ## 

    ax1.set_title("Recombination probabilities")

    ## remove extra axes
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.patch.set_visible(False)

    ## sample 50 tick labels for the x axis
    xtick_sampler = int(np.round(len(positions)/tick_sampler))
    ax1.set_xticks(range(0,len(positions), xtick_sampler))
    ax1.set_xticklabels(positions[::xtick_sampler],rotation = 30) 


    ## set the ylim so all plots start from 0
    ax1.set_ylim(0,(plot_max+(plot_max/10)))

    ## add axis labels
    ax1.set_ylabel("%s (averaged over %s gametes)" % ("Recombination probabilties", N_gametes))
    ax1.set_xlabel("Postition")

    ## inset the axis for the colorbar legend and format it
    cb_inset = fig.add_axes([0.65,0.9,0.2,0.05])
    cb_inset.set_title("Generations")
    cb_inset.imshow((generations,generations), cmap=BuPu, extent=[min(generations),max(generations),500,100])
    cb_inset.axes.get_yaxis().set_visible(False)

    ## turn off unwanted tickmarks on top and right axes
    plt.tick_params(top="off", right="off") 

    plt.savefig(outpath)
    
    #plt.show()
    plt.close()


###############################################################################################################################################################################################

def Recomb_prob_plot_per_window_per_iteration(pergen_perbase_avgs_filepath, outpath, genome_size = 100000, rc_wind = 1000, N_gametes = 50, tick_sampler = 20):

    """
    
    Plot the per-generaton per-base averaged recombination probabilites for 
    a single simulation (iteration).
    
    <pergen_perbase_avgs_filepath>    -    Filepath to the averaged (across gametes) recomb prob vectors per generation and per base for the simulation
    <outpath>                         -    The name of the file to which to save plot. Suffix (e.g. ".png", ".pdf" will determine the plot format)
    <genome_size>                     -    The size of the genome simulated
    <tick_sampler>                    -    The number of xtick positions to plot/label
    
    """
    
    
    from matplotlib import pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import matplotlib.colorbar as cb
    import numpy as np

    pergen_perbase_avgs = {}
    
    with open(pergen_perbase_avgs_filepath, 'r') as pergen_perbase_avgs_file:
        
        for line in pergen_perbase_avgs_file:
            
            gen = line.split("\t")[0]
            pergen_perbase_avgs[gen] = np.asarray([float(i) for i in line.split("\t")[1].split()])
            
    generations = sorted([int(i) for i in pergen_perbase_avgs.keys()])  ##################################
    positions = list(range(1, int(genome_size)+1))
    windows = list(range(1, int(genome_size/rc_wind)+1))
    

    ### Now for the plotting  ###
    fig = plt.figure(figsize = (20,10))

    ax1 = plt.subplot()

    ## make some pretty colours

    BuPu = plt.get_cmap("BuPu")
    cNorm  = colors.Normalize(vmin=min(generations)-500, vmax=max(generations))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=BuPu)

    colour_index = 0
    plot_max = 0

    print("Generating window averages: ")
    
    print("Generation: ", end = " ")
    
    for gen in generations: ## use generations list for iteration as its in order. Same for windows below.
        if int(gen) % 1000 == 0:
            print(" %s," % gen, end = " ")
        gen_data = pergen_perbase_avgs[gen]
        
        WIND_STRT = 0
        WIND_STOP = rc_wind-1
        
        WIND_averaged_vector = []
        
        for i in range(len(windows)):
            #if len(gen_data[WIND_STRT:WIND_STOP]) < 1:
            #    print(i,WIND_STRT, WIND_STOP)
            WIND_averaged_vector.append(np.mean(gen_data[WIND_STRT:WIND_STOP]))
            WIND_STRT += rc_wind
            WIND_STOP += rc_wind

        colorVal = scalarMap.to_rgba(generations[colour_index])

        ax1.plot(WIND_averaged_vector, color=colorVal, label = gen ) 
        
        if gen == 0:
            ax1.plot(WIND_averaged_vector, color="grey", dashes = [5, 10]) ## plot starting prob as dashed line
        
        colour_index += 1

        ## track highest value for mannually setting the ylim on plots later. 
        if max(WIND_averaged_vector) > plot_max:
            plot_max = max(WIND_averaged_vector)

    ## Plot formatting ## 

    ax1.set_title("Recombination probabilities")

    ## remove extra axes
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.patch.set_visible(False)

    ## sample 50 tick labels for the x axis
    xtick_sampler = int(np.round(len(positions)/int(len(windows)/2)))
    print(xtick_sampler)
    ax1.set_xticks(range(0,len(windows), 2)) 
    ax1.set_xticklabels(positions[::xtick_sampler],rotation = 30) 


    ## set the ylim so all plots start from 0
    ax1.set_ylim(0,(plot_max+(plot_max/10)))

    ## add axis labels
    ax1.set_ylabel("%s (averaged over %s gametes)" % ("Recombination probabilties", N_gametes))
    ax1.set_xlabel("Postition")

    ## inset the axis for the colorbar legend and format it
    cb_inset = fig.add_axes([0.65,0.9,0.2,0.05])
    cb_inset.set_title("Generations")
    cb_inset.imshow((generations,generations), cmap=BuPu, extent=[min(generations),max(generations),500,100])
    cb_inset.axes.get_yaxis().set_visible(False)

    ## turn off unwanted tickmarks on top and right axes
    plt.tick_params(top="off", right="off") 

    plt.savefig(outpath)
    
    #plt.show()
    plt.close()



###############################################################################################################################################################################################


def Recomb_prob_all_iterations(pergen_perbase_avgs_filepaths, rv_identifier = "AVGS", genome_size = 100000, rc_wind = 1000, tick_sampler = 20):

    """

    Plot the per-generaton per-base averaged recombination probabilites for
    all simualtions (iterations).

    <pergen_perbase_avgs_filepaths>     -    List of files of averaged (across gametes) recomb prob vectors per generation and per base for all simulations (iterations) to be plotted
    <rv_identifier>                     -    The string found in the rv_averages files. To be used to find all files in the directory to be plotted. 
    <genome_size>                       -    The size of the genome simulated
    <rc_wind>                           -    The value of RC_WIND used in the simulation 
    <tick_sampler>                      -    The number of xtick positions to plot/label
    """


    from matplotlib import pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx
    import matplotlib.colorbar as cb
    import numpy as np
    from collections import Counter
    import os



    positions = list(range(1, int(genome_size)+1))
    windows = list(range(1, int(genome_size/rc_wind)+1))
    
    
    N_iterations = len(pergen_perbase_avgs_filepaths)
    wd = pergen_perbase_avgs_filepaths[0].rpartition("/")[0]
    
    
    ## first pass to get the number of generations to use (cant think of a more efficient way to do this unfortunately)
    
    print("   Evaluating iterations and generations in the simulation")
    
    all_generations = []
    
    for rv_avgs_filepath in pergen_perbase_avgs_filepaths:
        
        with open(rv_avgs_filepath, 'r') as rv_avgs_file:
            
            for line in rv_avgs_file:
                
                all_generations.append(line.split()[0])
                
    generation_counter = Counter(all_generations)
    
    generations = []
    
    for generation in generation_counter:
        if generation_counter[generation] == N_iterations:
            generations.append(int(generation))

    generations = sorted(generations)
    print("   Last generation for which all (%s) iterations completed = %s" % (N_iterations, max(generations)))

    
    
    ## generate the per base averages across all simulations (iterations) - this is now the most time consuming part of the whole process. 

    per_gen_per_base_all_sim_avgs = {}  ## will contain one vector for each generation (e.g. 200 x 10k vector. Still quite big)

    
    print("\nGenerating per-base averages across %s iterations using 1 thread: " % N_iterations)
    
    perbase_pergen_all_sim_avgs_outpath = "%s/Perbase_pergen_avgs_over_%s_iterations.txt" % (wd, N_iterations)

    if os.path.exists(perbase_pergen_all_sim_avgs_outpath):
        print("\nA file containing per-base averages across %s iterations already exists - using this. Remove this file if you wish to recalculate." % N_iterations)

        with open(perbase_pergen_all_sim_avgs_outpath, 'r') as perbase_pergen_all_sim_avgs_out_handle:
            
            for line in perbase_pergen_all_sim_avgs_out_handle:
                
               per_gen_per_base_all_sim_avgs[int(line.split()[0])] = np.asarray([float(i) for i in line.split("\t")[1:]])

    else:
        print("\nGenerating per-base averages across %s iterations using 1 thread: " % N_iterations)

        with open(perbase_pergen_all_sim_avgs_outpath, 'w') as perbase_pergen_all_sim_avgs_out_handle:
    
            for iteration_file in pergen_perbase_avgs_filepaths:
            
                with open(iteration_file, 'r') as iteration:  

                    for line in iteration:

                        gen = int(line.split("\t")[0])

                        if gen in generations:  ## only use generations for which all iterations were present

                            if gen not in per_gen_per_base_all_sim_avgs:
                                per_gen_per_base_all_sim_avgs[gen] = np.asarray([float(i) for i in line.split("\t")[1].split()])
    
                            else:
                                per_gen_per_base_all_sim_avgs[gen] = np.sum([per_gen_per_base_all_sim_avgs[gen], np.asarray([float(i) for i in line.split("\t")[1].split()])], axis = 0)
                        
                for gen in per_gen_per_base_all_sim_avgs:

                    per_gen_per_base_all_sim_avgs[gen] = per_gen_per_base_all_sim_avgs[gen]/len(pergen_perbase_avgs_filepaths)  
                    perbase_pergen_all_sim_avgs_out_handle.write("%s\t%s\n" % (gen, "\t".join(str(i) for i in per_gen_per_base_all_sim_avgs[gen])))

                print("   File %s done" % iteration_file)   

 
    #######
    #######  Now reduce to the window averaged vectors for each gen across all iterations and plot ##
    #######

    ## Set up plot  ###

    fig = plt.figure(figsize = (20,10))

    ax1 = plt.subplot()

    ## make some pretty colours

    BuPu = plt.get_cmap("BuPu")
    cNorm  = colors.Normalize(vmin=min(generations)-500, vmax=max(generations))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=BuPu)

    colour_index = 0
    plot_max = 0

    ## get window averages

    print("\n\nGenerating window averages across %s iterations: " % N_iterations)
    print("Generation: ", end = " ")
    
    for gen in generations:
        if int(gen) % 1000 == 0:
            print(" %s," % gen, end = " ")

        gen_data = per_gen_per_base_all_sim_avgs[gen]
        
        WIND_STRT = 0
        WIND_STOP = rc_wind-1

        WIND_averaged_vector = []

        for i in range(len(windows)):
#            if len(gen_data[WIND_STRT:WIND_STOP]) < 1:
#                print("Window: %s - %s has no data" % (WIND_STRT, WIND_STOP))
            WIND_averaged_vector.append(np.mean(gen_data[WIND_STRT:WIND_STOP]))
            WIND_STRT += rc_wind
            WIND_STOP += rc_wind

        colorVal = scalarMap.to_rgba(generations[colour_index])


        if gen == 0:
            ax1.plot(np.log(WIND_averaged_vector), color="grey", dashes = [5, 10], label = gen) ## plot starting prob as dashed line
        else:
            ax1.plot(np.log(WIND_averaged_vector), color=colorVal, label = gen )
        colour_index += 1

        ## track highest value for mannually setting the ylim on plots later.
        if max(WIND_averaged_vector) > plot_max:
            plot_max = max(WIND_averaged_vector)


    ## Plot formatting ##

    ax1.set_title("Recombination probabilities")

    ## remove extra axes
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.patch.set_visible(False)

    ## sample tick labels for the x axis
    xtick_sampler = int(np.round(len(positions)/int(len(windows)/2)))
    print(xtick_sampler)
    ax1.set_xticks(range(0,len(windows), 2))
    ax1.set_xticklabels(positions[::xtick_sampler],rotation = 30)


    ## set the ylim so all plots start from 0
    #ax1.set_ylim(0,(plot_max+(plot_max/10)))

    ## add axis labels
    ax1.set_ylabel("%s (averaged over %s iterations)" % ("Log Recombination probabilties", N_iterations))
    ax1.set_xlabel("Postition")

    ## inset the axis for the colorbar legend and format it
    cb_inset = fig.add_axes([0.65,0.9,0.2,0.05])
    cb_inset.set_title("Generations")
    cb_inset.imshow((generations,generations), cmap=BuPu, extent=[min(generations),max(generations),500,100])
    cb_inset.axes.get_yaxis().set_visible(False)

    ## turn off unwanted tickmarks on top and right axes
    plt.tick_params(top="off", right="off")

    outpath = ""
    
    plt.savefig("%s/Recomb_probs_%s_iterations.pdf" % (wd, N_iterations))

    plt.show()
    plt.close()




###############################################################################################################################################################################################


def NRR_lengths_per_it(pergen_perbase_avgs_outpath, SD = 50000):
    
    """
    From the average per-base vectors for each generation, 
    
    """

    import numpy as np

    iteration = pergen_perbase_avgs_outpath.rpartition("/")[2].split("_")[2].split(".")[0]
    
    NRR_lengths = {}
    
    with open(pergen_perbase_avgs_outpath, 'r') as pergen_perbase_avgs:

        for line in pergen_perbase_avgs:
            
            gen = int(line.split("\t")[0])
            vector = np.asarray([float(i) for i in line.split("\t")[1].split()])
        
            left = vector[0:SD]
            right = vector[(SD+1):]

            Left_limit = (left[::-1]>0.0000000001).argmax()
            Right_limit = (right>0.0000000001).argmax()

            NRR_length = Left_limit + Right_limit

            NRR_lengths[gen] = NRR_length
            
    print("   File: %s done" % pergen_perbase_avgs_outpath)

    return (iteration, NRR_lengths)


###############################################################################################################################################################################################

def Get_last_NRR_avg(NRR_per_it):
    
    """
    Takes the dictionary of the NRR lengths per generation per iteration. Returns a dictionary the average NRR size accross all iterations in the last generation of the simulation.
    
    """
    
    import sys
    import numpy as np
    
    last_NRR_per_it = []

    last_iterations = []
    
    for i in NRR_per_it:
        last_it = max(NRR_per_it[i].keys())
        last_iterations.append(last_it)
        last_NRR_per_it.append(NRR_per_it[i][last_it])
    
    ## Check that all iterations finished - otherwise the average will not be valid.
    if len(set(last_iterations)) > 1:
        sys.exit("\nSTOPPED: Not all iterations end at the same generation!")
    elif len(set(last_iterations)) == 1:
        average_NRR_size_at_last_gen = np.mean(last_NRR_per_it)

    return average_NRR_size_at_last_gen



###############################################################################################################################################################################################

def Plot_NRR_size_per_simulation(NRR_per_it, outpath, xtick_sampler = 20):

    """
    Plots the evolution of NRR size (y) through generations (x).
    """

    import numpy as np
    from matplotlib import pyplot as plt

    ## get mean and std across iterations

    all_it_vectors = []
    for iteration in NRR_per_it:    
        it_vector = [NRR_per_it[iteration][gen] for gen in sorted(NRR_per_it[iteration].keys())]
        all_it_vectors.append(it_vector)
        
    means = np.mean(all_it_vectors, axis = 0)
    stds = np.std(all_it_vectors, axis = 0)

    y_low = means-stds
    y_high= means+stds

    fig = plt.figure(figsize = (20,10))

    ax1 = plt.subplot()

    ax1.fill_between(range(len(means)), y_low, y_high, color='blue', alpha=.05)

    
    for iteration in NRR_per_it:
        xlabs = sorted(list(NRR_per_it[iteration].keys()))
        ax1.plot(range(len(means)), [NRR_per_it[iteration][gen] for gen in sorted(NRR_per_it[iteration].keys())], color = "grey", alpha = 0.5)

    ax1.plot(range(len(means)), means, color = "black")

    ## remove extra axes
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ## sample tick labels for the x axis

    ax1.set_xticks(range(0,len(xlabs), 10))
    ax1.set_xticklabels(xlabs[::10],rotation = 30)


    plt.xlabel("Generation")
    plt.ylabel("NRR size (bp)")

    plt.title("Size of non-recombining region (NRR)")

    plt.savefig(outpath)
    #plt.show()
    plt.close()

###############################################################################################################################################################################################


def Recomb_plotter_mem(wd, num_iterations = 0, run_iterations = [], plot = True, plot_iterations = False, N_threads = 1, Gnome_size = 100000, N_gams = 100, RC_WIND = 1000, xtick_sampler = 20):

    """
    
    Plot recombination probabilities in windows along the simulated chromosome. This function can be used to plot the outputs of single simulations, or to average 
    over multiple simulations using the same parameters (refered to throughout as iterations). 
    
    The function takes recombination vector output files from the simulator. All files must be in the same directory. The script will automatically find the output files, 
    but for this the filenames must contain 1) the string "_rv_" and 2) a unique identifier for each simulation, preceded by an underscore and followed by a fullstop. E.g. m_rv_99.out"
    
    <wd>                -    The full path to the directory containing the recombination vector output files
    <num_iterations>    -    The number of iterations to use to make the plot (will be randomly chosen from all available files). Default = All used.
    <run_iterations>    -    A specific list of iterations to use (in the format ["1", "2", "3"]). Default = All used. Supersedes <num_iterations>
    <plot>              -    True (Default) or False. If False no plots will be plotted.
    <plot_iterations>   -    True or False (Default). If True, individual plots (pdf) will be produced for each processed iteration and saved to the subdirectory ./Iteration_plots/.
    <Gnome_size>        -    Size of the simulated genome
    <N_gams>            -    Number of gametes in _rv_ files (same as popualtion size used for simulation)
    <RC_WIND>           -    Window size to calculate average and plot recombination probabilities (recommended to use the same value as rc_wind in the simulation.)
    <xtick_sampler>     -    Number of positions to label on the x axis of the plots. 
    
    """

    import os
    import random
    import numpy as np
    import multiprocessing as mp


    print("Parsing rv outputs in %s" % wd)


    ## if plots for each iteration are requested, tidy them into a subdir

    if plot_iterations == True:
        if not os.path.exists("%s/Iteration_plots/" % wd):
            os.mkdir("%s/Iteration_plots/" % wd)

    rv_files = []

    for file in os.listdir(wd):
        if "_rv_" in file and "AVGS" not in file and "NRR" not in file:
            rv_files.append(file)

    ## get set of files to run, depending on args

    if all([num_iterations == 0, len(run_iterations) == 0]):
        rv_files_to_run = rv_files ## run all iterations

        print("\nUsing %s (all) iterations" % len(rv_files_to_run))

    elif all([num_iterations >= 1, len(run_iterations) == 0]):
        rv_files_to_run = random.sample(rv_files, int(num_iterations))
        print("\nRandomly chose %s iterations to use" % num_iterations)

    elif all([num_iterations == 0, len(run_iterations) >= 1]):
        rv_files_to_run = []
        for file in rv_files:
            for iteration in run_iterations:
                if "_%s." % iteration in file:
                    rv_files_to_run.append(file)
        print("\nUsing %s specified iterations" % len(rv_files_to_run))

    elif all([num_iterations >= 1, len(run_iterations) >= 1]):
        rv_files_to_run = []
        if num_iterations != len(run_iterations):
            print("\nWARNING, <num_iterations> and len(<run_iterations>) do not match, using run_iterations.")
            for file in rv_files:
                for iteration in run_iterations:
                    if "_%s." % iteration in file:
                        rv_files_to_run.append(file)
            print("\nUsing %s specified iterations" % len(rv_files_to_run))
        else:
            for file in rv_files:
                for iteration in run_iterations:
                    if "_%s." % iteration in file:
                        rv_files_to_run.append(file)

            print("\nUsing %s specified iterations" % len(rv_files_to_run))

    
    ## Check for rv_*_AVGS files already made
        
    per_base_avg_files_to_run = []
    
    for file in rv_files_to_run:
        
        if not "AVGS" in file:
            
            if not os.path.exists("%s/%s_AVGS.out" % (wd, file.rpartition(".")[0])):
                
                per_base_avg_files_to_run.append(file)
                
            else:
                print("   File: %s already has per-base recomb prob AVGS made. Not re-calculating these. Remove the existing AVGS file if you want to recalculate." % file)

    if len(per_base_avg_files_to_run) == 0:
        
        print("\n   Not calculating any per-base recomb prob averages - all necassary AVGS files already exist")

    else:
           
        ## Now run the files through
    
        ## Step 1. Convert recomb probabilities from block-wise to per base vectors
    
        if N_threads > 1:
    
            print("   Converting from blocks to per-base recomb probability vectors (using %s threads) for iterations: %s" % (N_threads, ", ".join([i.split("_")[2].split(".")[0] for i in per_base_avg_files_to_run])))
    
            #### PARALLELISED ####################
    
            pool = mp.Pool(N_threads)
    
            per_base_avg_files_to_run = sorted(per_base_avg_files_to_run) ## just make sure they are in a sensible order!
    
            rv_paths_to_run = ["%s/%s" % (wd, file) for file in per_base_avg_files_to_run]
    
            pool.map(rv_to_perbase_avg, rv_paths_to_run)
    
    
        elif N_threads == 1:
    
            print("   Converting from blocks to per-base recomb probability vectors (not parallelised) for iterations %s" % ", ".join([i.split("_")[2].split(".")[0] for i in per_base_avg_files_to_run]))
    
            #### NOT PARALLELISED ####################
    
            for file in per_base_avg_files_to_run:
    
                iteration = file.split("_")[2].split(".")[0]
    
                filepath = "%s/%s" % (wd, file)
    
                print("\n  Iteration: %s\n" % iteration)
    
                ## Step 1. Convert recomb probabilities from block-wise to per base vectors
                
                rv_to_perbase_avg(filepath)


    ## get the final list of the AVGS files now present in wd.
 
    rv_avgs_filepaths = ["%s/%s_AVGS.out" % (wd, rv_file.rpartition(".")[0]) for rv_file in rv_files_to_run]


    ## Step 2. Get NRR lenghts for each iteration
    
    ## check that NRR lengths not already calculated for each iteration

    print("\nGetting NRR lengths for %s iterations (N threads =  %s)" % (len(rv_files_to_run), N_threads))

    NRR_files_to_run = []
    NRR_per_it_unreduced = {}

    for rv_avgs_file_path in rv_avgs_filepaths:
        
        iteration_NRR_outpath = "%s_NRR.out" % (rv_avgs_file_path.rpartition("_AVGS.")[0])
        
        if os.path.exists(iteration_NRR_outpath):

            iteration = rv_avgs_file_path.rpartition("/")[2].split("_")[2].split(".")[0]
            print("   NRRs for iteration %s already exist, not calculating. Delete existing NRR file to recalculate" % iteration)

            NRR_per_it_unreduced[iteration] = {}

            with open(iteration_NRR_outpath) as iteration_out_handle:
                for line in iteration_out_handle:
                    NRR_per_it_unreduced[iteration][int(line.split()[0])] = int(line.split()[1])

        else:
            NRR_files_to_run.append(rv_avgs_file_path)

    if N_threads > 1:
         
        pool = mp.Pool(N_threads)
        
        results = pool.map(NRR_lengths_per_it, NRR_files_to_run)
        
        for i in results: ## unpack pooled results
            
            iteration = str(i[0])
            
            for rv_file in rv_files_to_run:
                if "_%s" % iteration in rv_file:

                    iteration_NRR_outpath = "%s/%s_NRR.out" % (wd, rv_file.rpartition(".")[0])

                with open(iteration_NRR_outpath, 'w') as iteration_out_handle:
            
                    NRR_per_it_unreduced[i[0]] = i[1]
                
                    for gen in i[1]:
                        iteration_out_handle.write("%s\t%s\n" % (gen, i[1][gen]))
            
        pool.close()
        pool.join()
            
        
        
    elif N_threads == 1:

        with open(iteration_NRR_outpath, 'w') as iteration_out_handle:

            for rv_avgs_filepath in NRR_files_to_run:

                rv_avgs_filepath = "%s/%s_AVGS.out" % (wd, file.rpartition(".")[0])

                iteration = rv_avgs_filepath.rpartition("/")[2].split("_")[2].split(".")[0]

                NRR_per_it_unreduced[iteration] = NRR_lengths_per_it(rv_avgs_filepath)[1]
       
                for gen in NRR_per_it_unreduced[iteration]:
                    iteration_out_handle.write("%s\t%s\n" % (gen, NRR_per_it_unreduced[iteration][gen]))

 
    print("\n   Reducing NRR sizes dataset to contain only generations completed in ALL processed iterations!")
    
    from collections import Counter
    generations = []
    generations_to_keep = []
    for iteration in NRR_per_it_unreduced:
        for gen in NRR_per_it_unreduced[iteration]:
            generations.append(gen)

    gen_count = Counter(generations)
    for gen in gen_count:
        if gen_count[gen] == len(NRR_per_it_unreduced):
            generations_to_keep.append(gen)
    max_gen_to_keep = max(generations_to_keep)

    NRR_per_it = {}

    for iteration in NRR_per_it_unreduced:
        if iteration not in NRR_per_it:
            NRR_per_it[iteration] = {}
        for gen in NRR_per_it_unreduced[iteration]:
            if gen <= max_gen_to_keep:
                NRR_per_it[iteration][gen] = NRR_per_it_unreduced[iteration][gen]

    
                

    ## Step 3. OPTIONAL. Plot each iteration

    if plot_iterations and plot == True:
        print("\nPlotting each iteration (you asked for it!)")

        for file in rv_files_to_run:

            iteration = file.split("_")[2].split(".")[0]
            rv_avgs_filepath = "%s/%s_AVGS.out" % (wd, file.rpartition(".")[0])
            outpath = "%s/Iteration_plots/%s.pdf" % (wd, file.rpartition(".")[0])
            
            Recomb_prob_plot_per_window_per_iteration(rv_avgs_filepath, outpath, genome_size = Gnome_size, rc_wind = RC_WIND, N_gametes = N_gams, tick_sampler = xtick_sampler)
            print("Iteration %s plot saved here: %s" % (iteration, outpath))


    print("\nProcessed %s rv files" % len(rv_files_to_run))

    if plot == True:

        ## Step 4. Plot
        print("\nPlotting recombination probabilties by window across all simualtions")

        Recomb_prob_all_iterations(rv_avgs_filepaths, genome_size = Gnome_size, rc_wind = RC_WIND, tick_sampler = xtick_sampler)

        ## Step 5. Plot NRR length evolution across all iteration for this simulation

        print("\nPlotting the NRR size over the course of the simulation")
        Plot_NRR_size_per_simulation(NRR_per_it, "%s/NRR_size_all_its.pdf" % wd)

    else:
        print("\nNot plotting anything")

    ## Step 6. Get the final average NRR length across all iterations for this simulation.

    Final_NRR_size = Get_last_NRR_avg(NRR_per_it)

    print("Average final size of non-recombining region (generation %s) = %s kb" % (max_gen_to_keep, np.round((Final_NRR_size/1000), 2)))

    #return Final_NRR_size




###############################################################################################################################################################################################



def Stat_plotter(logfile_dir, Stat_to_plot, verbose = False):

    """
    Stat_plotter - Plots Fst and Dxy outputs from Sex chromosome recombination evolution simulations

    USAGE: python Logplotter.py logfile_dir Stat_to_plot

    logfile_dir    - The absolute path to the directory containing all logfiles from simultions to be plotted
                     (Note all simulations must have identical number of generations and windows)
    Stat_to_plot   - A string to determine the statistic to plot. The string must be in the name of the logfile and
                     must be unique to that type of logfile. e.g. "fst", "dxy", "XX", "XY"


    """

    import os
    import numpy as np
    import matplotlib
    #matplotlib.use('Agg')
    from matplotlib import pyplot as plt
    import matplotlib.colors as colors
    import matplotlib.cm as cmx


    print("\nPlotting %s" % Stat_to_plot)

    ## Get files, generation and window numbers.

    logfiles = []

    for logfile in os.listdir(logfile_dir):
        if Stat_to_plot in logfile and "pdf" not in logfile and "png" not in logfile:
            logfiles.append("%s/%s" % (logfile_dir, logfile))

    ## Get the number of generations
    # Mathias modified: if some runs are incomplete for the number of generations, find the largest generation that was achieved in any of the runs (files)
    generations = []
    for f in logfiles:
        if verbose == True:
            print(f)
        with open(f, "r") as F:
            for line in F:
                if not any(["site" in line, "generation" in line]):
                    generations.append(int(line.split()[0]))
    g = sorted(list(set(generations)), key = int)
    generations = g

    ## Get the windows

    for line in open(logfiles[0], 'r').readlines():
        if "generation" in line:
            windows = line.split()[2:]
    ## Note - all the files in the simulation directory must have the same number of windows!=


    ## Set up the dictionary that will contain all the data

    Raw_data_by_gen = {}

    for gen in generations:
        Raw_data_by_gen[gen] = {}

        for window in windows:
            Raw_data_by_gen[gen][window] = []


    ## Now get data. Could probably do this easier with pandas, but going to hack it just with base python.

    for logfile in logfiles:
        for line in open(logfile, 'r').readlines():
            if not any(["site" in line, "generation" in line]):
                generation = int(line.split()[0])
                data_fields = line.split()[2:]
                field_index = 0


                for field in data_fields:
                    window = windows[field_index]
                    if field == "NA":
                        Raw_data_by_gen[generation][window].append(np.nan)
                    else:
                        Raw_data_by_gen[generation][window].append(float(field))
                    field_index += 1



    ## Now get the averages

    Averaged_data_by_gen = {}

    for gen in Raw_data_by_gen:
        Averaged_data_by_gen[gen] = {}
        for window in Raw_data_by_gen[gen]:
            Averaged_data_by_gen[gen][window] = np.nanmean(Raw_data_by_gen[gen][window])



    ### Now for the plotting  ###

    #plt.style.use('dark_background')
    fig = plt.figure(figsize = (20,10))

    ax1 = plt.subplot()

    ## make some pretty colours

    BuPu = plt.get_cmap("BuPu")
    cNorm  = colors.Normalize(vmin=min(generations)-500, vmax=max(generations))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=BuPu)

    colour_index = 0
    plot_max = 0
    for gen in generations: ## use generations list for iteration as its in order. Same for windows below.
        gen_data = []
        for window in windows:
            gen_data.append(Averaged_data_by_gen[gen][window])

        colorVal = scalarMap.to_rgba(generations[colour_index])

        ax1.plot(gen_data, color=colorVal, label = gen ) #, cmap = "BuPu")

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
    ax1.set_ylabel("%s (averaged over %s iterations)" % (Stat_to_plot, len(logfiles)))
    ax1.set_xlabel("Window")

    ## inset the axis for the colorbar legend and format it
    cb_inset = fig.add_axes([0.65,0.8,0.2,0.2])
    cb_inset.set_title("Generations")
    cb_inset.imshow((generations,generations), cmap=BuPu, extent=[min(generations),max(generations),0,1000])
    cb_inset.axes.get_yaxis().set_visible(False)



    plt.savefig("%s/%s_summary_plot.pdf" % (logfile_dir, Stat_to_plot))
    #plt.show()
    plt.close()

    print("\nDone! Your %s plot is here: %s/%s_summary_plot.pdf" % (Stat_to_plot, logfile_dir, Stat_to_plot))



###############################################################################################################################################################################################


def NRR_comp_plot(Multi_sim_dict, outfilepath):
    
    """
    Plots the comparative plot across multiple simualtions of the final size of the non recombining region.

    """

    from matplotlib import pyplot as plt
    import numpy as np
    
    
    Per_N_dicts = {}
    
    for sim in Multi_sim_dict:
        N = Multi_sim_dict[sim]["Params"]["N"]
        
        if N not in Per_N_dicts:
            Per_N_dicts[N] = {}
        
        Per_N_dicts[N][sim] = {}
        Per_N_dicts[N][sim]["NRR"] = Multi_sim_dict[sim]["Final_NRR"]
        Per_N_dicts[N][sim]["rp"] = Multi_sim_dict[sim]["Params"]["rp"]
        Per_N_dicts[N][sim]["u"] = Multi_sim_dict[sim]["Params"]["u"]
        
    
    scale_down = 10
    N_subplots = len(Per_N_dicts)

    fig = plt.figure(figsize = ((4*N_subplots)+3, 5))
    
    subplot_index = 1
    
    for N in sorted(Per_N_dicts.keys()):
        
        plt.subplot(1,N_subplots, subplot_index)
        
        x_labels = []
        y_labels = []

        sizes = []
        
        for sim in Per_N_dicts[N]:
            
            x_labels.append(Per_N_dicts[N][sim]["u"])
            y_labels.append(Per_N_dicts[N][sim]["rp"])
            sizes.append(Per_N_dicts[N][sim]["NRR"])
        
        size_labels = ["%s kb" % (int(i)/1000) for i in sizes]
        
        
        plt.scatter(np.log(x_labels),np.log(y_labels), s=(np.array(sizes)/scale_down), color = "black" )
        
        for i in list(range(0, len(x_labels))):
            plt.text(np.log(x_labels)[i],np.log(y_labels)[i]+(sizes[i]/10000)/2, size_labels[i] , verticalalignment='center', horizontalalignment='center')
        
        plt.xticks(np.log(x_labels), x_labels)
        plt.yticks(np.log(y_labels), y_labels)
        
        plot_lim_x_buffer = (max(np.log(x_labels))-min(np.log(x_labels)))/4
        plot_lim_y_buffer = (max(np.log(y_labels))-min(np.log(y_labels)))/4
        #plt.xlim(min(np.log(x_labels))-plot_lim_x_buffer, max(np.log(x_labels))+plot_lim_x_buffer)
        #plt.ylim(min(np.log(y_labels))-plot_lim_y_buffer, max(np.log(y_labels))+plot_lim_y_buffer)
        plt.xlabel("log u")
        
        if subplot_index % 2 != 0:
            plt.ylabel("log r")
        
        plt.title("Ne = %s" % N)
        
        subplot_index += 1
    
    plt.savefig(outfilepath)
    
    #plt.show()
    plt.close()



###############################################################################################################################################################################################


def sex_chrom_sim_plotter(Multi_sim_dir = None,
                          Single_sim_dir = None,
                          Sims_to_plot = None,
                          plots = ["fst", "dxy", "dxx", "recomb", "recomb_comp_only"],
                          parameter_file = None,
                          prefix = "sim",
                          sim_ID_pos = 0,
                          sep = "_",
                          NRR_num_iterations = 0,
                          NRR_run_iterations = [],
                          NRR_plot_iterations = False,
                          NRR_N_threads = 1,
                          NRR_Gnome_size = 100000,
                          NRR_N_gams = 100,
                          NRR_RC_WIND = 1000,
                          NRR_xtick_sampler = 20):

    """
    Wrapper for all sex chromosome simulation plots.


    Single_sim_dir     -     If plotting only outputs from a single simulation with multiple iterations all using the same parameter set. All output files must be in this directory (no subdirs).

    Multi_sim_dir      -     If plotting outputs from multiple simulations (using different parameter combinations). The script will automatically find the output files, but for this the names
                             of the directories for each simualtion must contain 1) a string found in all simulation files to be processed 2) a unique identifier for each simulation.
                             These can be set using the <prefix>, <sim_ID_pos> and <sep> arguments below. Defaults expect format ./ID_sim*/

    Sims_to_plot       -     A list specifying the IDs of a subset of simulations to process in the parent directory <Multi_sim_dir>. Default = all simualtions will be processed.

    prefix             -     A string to be found in all directory names to be included. Default = "sim"

    sim_ID_pos         -     The (zero based indexing) position in the string of the unique simulation ID (corresponding to those in the parameter file). Default = 0 (first position in the name)

    sep                -     The separator in the directory name. Default = "_"

    parameter_file     -     A tab separated file containing the parameter combination for each simulation in a multi-sim set. Format: ID Ne r u rec_prob (1st line contains headers).
                             Each in this file must correspond to those in simulation directory names.

    plots              -     A list of the plots to be produced. Choose any combination from ["Fst", "dxy", "dxx", "recomb", "recomb_comp_only"]. "recomb_comp_only" is only for multi-sim jobs.

    NRR_Gnome_size     -     Size of the simulated genome

    NRR_N_gams         -     Number of gametes in _rv_ files (same as popualtion size used for simulation)

    NRR_RC_WIND        -     Window size to calculate average and plot recombination probabilities (recommended to use the same value as rc_wind in the simulation.)

    NRR_xtick_sampler  -     Number of positions to label on the x axis of the plots.


    ## NRR Plots - these contain some potentially time consuming data processing steps, thus we include several arguments to configure how this proceeds.

    NRR_num_iterations     -    The number of iterations per simulation to use to make the plot (will be randomly chosen from all available files). Default = All used.
    NRR_run_iterations     -    A specific list of iterations to use per simulation (in the format ["1", "2", "3"]). Default = All used. Supersedes <num_iterations>
    NRR_plot_iterations    -    True or False (Default). If True, individual plots (pdf) will be produced for each processed iteration and saved to the subdirectory ./Iteration_plots/. Caution, this
                                can result in many plots and a lot of disk space usage. Only recommended as a troubleshooting guide for one or two iterations at a time.
    NRR_N_threads          -    Number of threads to use when converting recomb probabilities to per-base vectors - will significantly reduce run time. Default = 1 (i.e. not parallelised)


    """

    import os
    import sys

    #######################################
    ### Plotting a single simulation only##
    #######################################

    if Single_sim_dir:

        for plot in ["fst", "dxy", "dxx"]:
            if plot in plots:

                Stat_plotter(Single_sim_dir, plot)

        if "recomb" in plots:
            print("\nPlotting the evolution of recombination:")
            Final_NRR = Recomb_plotter_mem(Single_sim_dir,
                                       plot = True,
                                       num_iterations = NRR_num_iterations,
                                       run_iterations = NRR_run_iterations,
                                       plot_iterations = NRR_plot_iterations,
                                       Gnome_size = NRR_Gnome_size,
                                       N_threads = NRR_N_threads,
                                       RC_WIND = NRR_RC_WIND,
                                       xtick_sampler = NRR_xtick_sampler)



    ################################
    ## Plotting multiple simulations
    ################################

    elif Multi_sim_dir:

        if not parameter_file:

            sys.exit("\nERROR: Parameter file required")

        param_dict = {}


        for line in open(parameter_file, 'r').readlines()[1:]:
            ID = line.split()[0]
            param_dict[ID] = {}
            param_dict[ID]["N"] = int(line.split()[1])
            param_dict[ID]["u"] = float(line.split()[2])
            param_dict[ID]["rp"] = float(line.split()[3])


        # get simulation directories

        sim_dirs = []
        All_sim_IDs = []
        sims = {}

        for file in os.listdir(Multi_sim_dir):

            if prefix in file:
                sim_dir = "%s/%s" % (Multi_sim_dir,file)
                sim_dirs.append(sim_dir)
                sim_ID = sim_dir.rpartition("/")[2].split(sep)[sim_ID_pos]
                All_sim_IDs.append(sim_ID)


        if len(sim_dirs) == 0:
            sys.exit("\nERROR: No simulation directories found - please check format of directory names.")

        if not Sims_to_plot:
            Sims_to_plot = All_sim_IDs

        for sim_dir in sim_dirs:
            sim_ID = sim_dir.rpartition("/")[2].split(sep)[sim_ID_pos]

            if sim_ID in Sims_to_plot:
		
                if not sim_ID in param_dict:
                    sys.exit("ERROR: Sim ID not in parameter file")

                sims[sim_ID] = {}
                sims[sim_ID]["Params"] = param_dict[sim_ID]

                print("\n\nProcessing %s" % sim_dir)

                ### Plotting each simulation

                for plot in ["fst", "dxy", "dxx"]:

                    if plot in plots:

                        Stat_plotter(sim_dir, plot)

                if "recomb" in plots:
                    print("\nPlotting the evolution of recombination:")
                    sims[sim_ID]["Final_NRR"] = Recomb_plotter_mem(sim_dir,
                                                               num_iterations = NRR_num_iterations,
                                                               run_iterations = NRR_run_iterations,
                                                               plot_iterations = NRR_plot_iterations,
                                                               N_threads = NRR_N_threads,
                                                               Gnome_size = NRR_Gnome_size,
                                                               RC_WIND = NRR_RC_WIND,
                                                               xtick_sampler = NRR_xtick_sampler)


        if "recomb" in plots:
            print("\nPlotting comparative recombination supression plot:")
            NRR_comp_plot(sims, "%s/NRR_comp.pdf" % Multi_sim_dir)
    

    print("\nAll done! :) \n")
