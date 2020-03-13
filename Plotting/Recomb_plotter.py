def Recomb_plotter(data_dir, sex, iter_start = 1, iter_stop=None):
    
    """
    USAGE:   
    
          Recomb_plotter.py <full/path/to/data/directory> [start iteration=1] [stop iteration=N]
          
    OUTPUT:
    
          A histogram of all recombination events for the recorded generations
    
    """
    
    
    from matplotlib import pyplot as plt
    import os

    plt.figure(figsize = (15,10))

    if iter_stop == None:
        
        iter_stop = 0
        
        for i in os.listdir(data_dir):
            if "_rp_" in i:
                iter_stop += 1
    
    if sex == "M":
        sex = "m"
    elif sex == "F":
        sex = "f"    
    for iteration in range(iter_start,iter_stop+1):

        recombs = open("%s/%s_rp_%s.out" % (data_dir, sex, str(iteration)), 'r').readlines()
        #print "%s/f_rp_%s.out" % (wd, str(iteration))


        for line in recombs[1:]:
            gen = line.split()[0]
            recombs = [int(i) for i in line.split()[1:]]

        plt.hist(recombs, alpha = 0.2, color = "black", bins = 100)

    plt.title("Histogram of recombination locations (Sex: %s, Generations: %s, Iterations: %s)" % (sex, gen, len(range(iter_start,iter_stop))))
    plt.xlabel("Position (bp)")
    plt.savefig("%s/recombinations_%s_ITS_%s-%s.png" % (data_dir, sex, iter_start, iter_stop))
    plt.savefig("%s/recombinations_%s_ITS_%s-%s.pdf" % (data_dir, sex, iter_start, iter_stop))


## CLINE 

import sys

if len(sys.argv) == 3:
    Recomb_plotter(sys.argv[1], sys.argv[2])
elif len(sys.argv) == 5:
    Recomb_plotter(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
else:
    sys.exit("Wrong number of arguments\n %s" % Recomb_plotter.__doc__)
