import sys

if len(sys.argv) == 1:
    sys.exit("No config file given. \nSee the template.conf for help on running sim_plotter.py")

config_path = sys.argv[1]

config = open(config_path, 'r').readlines()

## set defaults

Multi_sim_dir = None
Single_sim_dir = None
Sims_to_plot = None
plots = ["fst", "dxy", "dxx", "recomb", "recomb_comp_only"]
parameter_file = None
prefix = "sim"
sim_ID_pos = 0
sep = "_"
NRR_num_iterations = 0
NRR_run_iterations = []
NRR_plot_iterations = False
NRR_N_threads = 1
NRR_Gnome_size = 100000
NRR_N_gams = 100
NRR_RC_WIND = 1000
NRR_xtick_sampler = 20


## Get and check non default params from plotting.config file

print("Non default arguments: \n")

for arg in config:
    if not arg.startswith(("#", "\n", " ")): 
        
        if len(arg.strip().split()) > 1:

            arg_name = arg.strip().split()[0]
            arg_value = "".join(arg.strip().split()[1:])
            
            if arg_name == "sim_plotter_lib_path":
                sys.path.append(arg_value)
            
            elif arg_name == "Multi_sim_dir":
                Multi_sim_dir = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "Single_sim_dir":
                Single_sim_dir = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "Sims_to_plot": 
                Sims_to_plot = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "plots":
                plots = arg_value.split(",")
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "parameter_file":
                parameter_file = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "prefix":
                prefix = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "sim_ID_pos":
                sim_ID_pos = int(arg_value)
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "sep":
                sep = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_num_iterations":
                NRR_num_iterations = int(arg_value)
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_run_iterations":
                NRR_run_iterations = arg_value.split(",")
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_plot_iterations":
                NRR_plot_iterations = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_N_threads":
                NRR_N_threads = int(arg_value)
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_Gnome_size":
                NRR_Gnome_size = int(arg_value)
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_N_gams":
                NRR_N_gams = arg_value
                print("  %s = %s" % (arg_name, arg_value))
                
            elif arg_name == "NRR_RC_WIND":
                NRR_RC_WIND = int(arg_value)
                print("  %s = %s" % arg_name, arg_value)
                
            elif arg_name == "NRR_xtick_sampler":
                NRR_xtick_sampler = int(arg_value)
                print("  %s = %s" % (arg_name, arg_value))
                
if Multi_sim_dir == None and Single_sim_dir == None:
    print("No directories supplied")

                
## Run plotting wrapper

import sim_plotter_lib as sp

sp.sex_chrom_sim_plotter(Multi_sim_dir = Multi_sim_dir,
                         Single_sim_dir = Single_sim_dir, 
                         Sims_to_plot = Sims_to_plot,
                         plots = plots,
                         parameter_file = parameter_file, 
                         prefix = prefix, 
                         sim_ID_pos = sim_ID_pos, 
                         sep = sep,
                         NRR_num_iterations = NRR_num_iterations, 
                         NRR_run_iterations = NRR_run_iterations, 
                         NRR_plot_iterations = NRR_plot_iterations, 
                         NRR_N_threads = NRR_N_threads,
                         NRR_Gnome_size = NRR_Gnome_size,
                         NRR_N_gams = NRR_N_gams,
                         NRR_RC_WIND = NRR_RC_WIND,
                         NRR_xtick_sampler = NRR_xtick_sampler)
        
        

