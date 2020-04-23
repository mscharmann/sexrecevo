def NRR_plotter(wd, GENOME_SIZE = 100000, SD_POS = 50000):
    
    """
    
    USAGE:
        NRR_plotter.py </full/path/to/sim_directories/> [GENOME_SIZE = 100000] [SD_POS = 50000]
    
       </full/path/to/sim_directories/>   -   The parent dir containing all sim outputs. Child dirs have to be in the
                                              form [0-9][0-9]_*_OUTS (i.e. a two digit number starts and ends with "OUTS")
       [GENOME_SIZE]                      -   The size of the simulated genome, default is 100000 
       [SD_POS]                           -   The position in the simulated genome of the SD default is 50000
    
    
    DESCRIPTION:
    This function will plot the size of the non recombining region (NRR) over the course of simulations. 
    NRR size is calculated each generation using the distance between the two recombination events most adjacent to
    the sex determiner. This is done for each iteration and the values across all iterations are averaged for each generation.
    
    This function outputs (in the same directory as input files) a pdf and png version of a scatter plot which gives generation vs average NRR size over all iterations.
    It will also perform a linear regression for this relationship and output the slope of this as a measure of the rate at which
    the NRR grows over the course of the simulations. Note that this is linear, any non-linear relationship will not be captured by this

    
    """
    
    import numpy as np
    from sklearn.preprocessing import PolynomialFeatures
    from sklearn.linear_model import LinearRegression
    from sklearn.pipeline import make_pipeline
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import os
    
    
    sim_dirs = [i for i in os.listdir(wd) if "OUTS" in i]
    plt.figure(figsize = (20,6*len(sim_dirs)))
    plot_row_index = 0
        
    for i in sim_dirs:
        
        sim_dir = "%s/%s" % (wd, i)
        sim_number = i.split("_")[0]

        ## initate figure axes (figure will contain a subplot for each sex)
        

        for sex in ["m", "f"]: # note - lower case - must be like this in the file name

            if sex == "m":
                plot_column_index = 1
            elif sex == "f":
                plot_column_index = 2
            
            sim_files = [i for i in os.listdir(sim_dir) if sex in i and "rp" in i]

            all_recombs = {}

            for sim_file in sim_files:

                if "rp" in sim_file:

                    sim_file_lines = open("%s/%s" % (sim_dir,sim_file), 'r').readlines()

                    for line in sim_file_lines[1:]:
                        gen = line.split()[0]
                        recombs = line.split()[1:]

                        if gen not in all_recombs:
                            all_recombs[gen] = []

                        lowers = []
                        uppers = []

                        for recomb in recombs:

                            if int(recomb) < SD_POS:
                                lowers.append(SD_POS - int(recomb))
                            elif int(recomb) > SD_POS:
                                uppers.append(int(recomb) - SD_POS)

                        if len(lowers) == 0:
                            lower = 0
                        else:
                            lower = min(lowers)

                        if len(uppers) == 0:
                            upper = GENOME_SIZE
                        else:
                            upper = min(uppers)

                        all_recombs[gen].append(upper+lower)

            all_recombs_averages = {}

            for gen in all_recombs:
                all_recombs_averages[gen] = {}
                all_recombs_averages[gen]["mean"] = np.mean(all_recombs[gen])
                all_recombs_averages[gen]["std"] = np.std(all_recombs[gen])
            
            x = sorted(int(i) for i in all_recombs_averages.keys())
            y = [all_recombs_averages[str(i)]["mean"] for i in x]

            
            ## make STD margins for plot

            y_low = [all_recombs_averages[str(i)]["mean"]-all_recombs_averages[str(i)]["std"] for i in x]
            y_high=[all_recombs_averages[str(i)]["mean"]+all_recombs_averages[str(i)]["std"] for i in x]

            ## Do regression

            X = np.array(x).reshape(-1, 1)  # values converts it into a numpy array
            Y = np.array(y).reshape(-1, 1)  # -1 means that calculate the dimension of rows, but have 1 column

            linear_regressor = LinearRegression()  # create object for the class
            linear_regressor.fit(X, Y)  # perform linear regression
            Y_pred = linear_regressor.predict(X)

            ## make plot

            plt.subplot(len(sim_dirs),2,plot_row_index+plot_column_index)
            
            plt.scatter(x,y, c = "black",alpha = .7)
            plt.fill_between(x, y_low, y_high, color='black', alpha=.1)
            plt.plot(X, Y_pred, color='blue')
            plt.xlabel("Generation")
            plt.ylabel("Size of NRR (Mean of 100 iterations, bp)")

            if sex == "m":
                plt.title("Sim %s - Expansion of non-recombining region in males (slope = %s)" % (sim_number, np.round(linear_regressor.coef_[0][0],2)))
            elif sex == "f":
                plt.title("Sim %s - Expansion of non-recombining region in females (slope = %s)" % (sim_number, np.round(linear_regressor.coef_[0][0],2)))
            
        plot_row_index += 2

    plt.savefig("%s/NRR_expansion.pdf" % wd)
    plt.savefig("%s/NRR_expansion.png" % wd)
    
    
## CLINE

import sys

if len(sys.argv) <2:
    
    sys.exit("\nERROR: Incorrect number of arguments\n%s" % NRR_plotter.__doc__)

elif len(sys.argv) == 2:
    working_dir = sys.argv[1]
    NRR_plotter(working_dir)
    
elif len(sys.argv) == 3:
    working_dir = sys.argv[1]
    genome = sys.argv[2]
    NRR_plotter(working_dir, int(genome))
    
elif len(sys.argv) == 4:
    working_dir = sys.argv[1]
    genome = sys.argv[2]
    SD = sys.argv[3]
    NRR_plotter(working_dir, int(genome), int(SD))