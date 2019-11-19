## Jörns c++ version of the simulation 

This is a rebuild of Mathias' python model in c++. Features of the model include:
	-output of male/female FST per window
	-output of dxx and dxy (differentiation of X and Y sequences in males and both Xs in females per window
	-sex specific mutation rates and recombination probabilities (for modelling achiasmy and heterochiasmy)
The software is build using the C++ 11 standard library. Unlike in Mathias script, parameters are passed through a simple config file (see example.conf).

The software can be build (under linux) using:

g++ main.cpp -std=c++11 -o schrom_sim

It is run as:

./schrom_sim example.conf

I did not test everything extensively, so please let me know if something doesn't work, if somethings behave unepectedly or if you're missing any kind of feature.
