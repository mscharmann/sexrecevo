# Neutral sex chromsome evolution model simulation

This software simulates a dioecious population (without selection) with non-overlapping generations and a fixed number of diploid individuals. The evolution of recombination surpression around the sex determining region (SDR) is simulated explicitly.

## Steps in one generation:

1. At the beginning of each generation offspring is generated from the individuals of the previous generation. Each offspring individual randomly draws a mother from the females and a father from the males of the previous generation.
Each parent supplies one haploid gamete to its offspring, which is the result of meiosis.
In each meiosis one crossing over event takes place. However, it will not necessarily take place within the simulated region, in which case the gamete is simply the copy of one of the parental haplotypes.
Within the simulated region the relative probability of a recombination event is estimated for each possible recombination position based on relative divergence between both haplotypes within a fixed window around this position and the position of the recombination event is drawn based on these probabilities.
2. The sex of each offspring individual is determined by the genotype at a single site, that has been defined as the SDR site in the config file.
3. Variants that are fixed are removed from the genotype matrix and their positions are freed for new mutations.
4. New mutations are inserted in the genotype matrix. We are assuming an infinite sites model, so mutations can only fall on sites on which there is no variant present in any individual.

## Details of the simulation

### Estimation of recombination probabilities

Recombination proababilities are based on the number of variants between both haplotypes at each meiosis. At each position at which recombination can take place recombination probabilities are calculated based on three parameters

* *rc_wind*: The area around the position in which recombination probabilities are reduced relative to the number of variants (*n_variants*) within that window.
* *min_ident*: The minimum identity $$1 - \frac{n\_variants}{rc\_wind}$$ within the window under which the recombination probability is reduced to 0
* *var_effect*: the relative degree to which a variant reduces recombination probabilities 

Both *min_ident* and *var_effect* reduce recombination probabilities, but *var_effect* is only used for a specific region if the number of variants is not sufficient toreach min_ident.
How exactly these parameters are used to calculate recombination probabilities depends on the model used (see below).

During meiosis the position of a recombination event is determined by a weighted draw, based a vector with relative recombination probabilities. This vector is made up of recombination probabilities for the simulated region and the recombination probability for the unsimulated part of the sex chromosome defined by *rec_prob*. The value for *rec_prob* is fixed and is added to the right end of the probability vector. 
As a consequence, the probability that a recombination event will fall outside the simulated region will increase with increasing recombination surpression in the simulated region.
If *rec_prob* is set to 0, recombination probabilities will be relative to other recombination probabilities within the simulated region.

#### Dynamic model

This is the standard model that is automatically run unless the *hotspot_conf* parameter is set in the config file. Here recombination probabilities are estimated for each possible recombination position in the simulated region, but probabilities of regions which have the same recombination probability and where the exact position of a recombination event doesn't matter are integrated in a dynamic way, which increases computational efficiency.

For a given site, recombination probabilities are calculated in the following way:

$$ rec\_prob = \frac{ 1 - \frac{var\_effect * n\_variants}{rc\_wind}}{genome\_size} $$ 

#### Hotspot model

This is the model that is run if the *hotspot_conf* parameter is set in the config file. In the hotspot model recombination can only take place at sites defined in the hotspot.conf file. Each mutation hotspot has a a priori recombination probability, which can be different for each hotspot and between sexes. 
Recombination probabilities are estimated in the same way as in the dynamic model: Recombination probabilities are reduced based on the number of variants within *rc_wind* around each hotspot (variants outside this region have no effect on recombination probabilities). Here the difference is that recombination probabilities are reduced relative to their a priori recombination probability.

$$ rec\_prob = a\_priori\_rec\_prob * ( 1 - \frac{var\_effect * n\_variants}{rc\_wind} )$$ 

### Mutations

This model follows the principle of the infinite sites model. At each site there can be two alleles and there can be no new mutations at a site at which two alleles are already segregating in the population. However, sites are removed when one allele gets fixed, after which a new mutation can fall on this position.
The number of mutations are defined by two parameters:

* *mut_rate*: The mutation rate, the absolute number of mutations per generation is calculated as $$ n\_muts =mut\_rate * 2 * pop\_size * genome\_size $$
* *mut_mb*: male bias of the mutation rate, this parameter can change how the absolute number of mutations is divided between males and females.
At each generation the number of mutations are drawn from a poisson distribution and are then divided between males and females.
$$ n\_muts\_f = n\_muts * \frac{1}{mut\_mb + 1} $$
$$ n\_muts\_m = n\_muts - n\_muts\_f $$

## Configuration files

Entries in configuration files are in rows and tab-separated. Empty rows and rows starting with # (comments) are ignored.
   
### General configuration file

This file determines all the parameters of the model. A parameter is defined as tab-separated parameter value pair. 

#### Parameters

* *ngens*: the number of generations the population is run

* *pop_size*: the number of individuals in the population

* *genome_size*: the size of the simulated region

* *mut_rate*: the mutation rate

* *mut_mb*: male bias of the mutation rate

* *rec_prob*: the size of the unsimulated region relative to the simulated region. *rec_prob*=2 would equal an unsimulated region two times the size of the simulated region in the dynamic model, so the probability of a recombination event falling in the simulated region (if there is no recombination surpression) would be 1/3. In the hotspot model rec_prob is relative to the sum of the a priori recombination probabilities of the hotspots, which can be different from 1.

* *var_effect*: this multiplies the effect that a single variant has on recombination probabilities.

* *sd_pos*: position of the SD position. This is a single site by which offspring sex is determined.

* *rc_wind*: window size that is used to estimate recombination probabilities around a recombination position  

* *min_ident*: the minimum identity within rc_window under which recombination probabilities are 0

* *sstat_wind*: window size for estimating summary statistics (FST, dxx, dxy)

* *sstat_gen*: estimate summary statistics every sstat_gen generations

* *fst_ofile*: filename for the FST output file

* *dxy_ofile*: filename for the dxy output file

* *dxx_ofile*: filename for the dxx output file

* *m_rp_ofile*: filename for the output file containing all male recombination positions in the generation

* *f_rp_ofile*: filename for the output file containing all female recombination positions in the generation

* *hotspot_conf*: filename of the configuration file for the hotspot model. If this parameter is set, the hotspot model will run, if not the dynamic model.

* *m_rv_ofile*: filename for the output file of male recombination probability vectors. This parameter is optional, if it is not set there is no output file. 


### Hotspot configuration file

This file defines the recombination points for the hotspot model. A hotspot is defined by its position, its relative, a priori recombination probability in male meiosis and its relative a priori recombination probability in female meiosis.
Each row corresponds to one recombination hotspot, with its position, its male and its female a priori recombination probabilities separated by tabs.

## Output files

Summary statistics are calculated every sstat_gen generations and written to output files.

### m/f Fst

Fst between males and females is calculated in sliding windows for each variable site in the window. If there are no variable sites in the window the output is NA, if there is more than one variable site it's their mean.
Fst is calculated by:

$$ fst = \frac{ \frac{Ht - ( n\_m * Hm + n\_f * Hf )}{n\_m + n\_f} }{Ht}$$

With *Ht*, *Hm* and *Hf* being the expected heterozigosity averaged over loci for the whole dataset, males and females respectively and *n_m* and *n_f* being the number of males and females.

The output file includes the generation in which FST is calculated (first column), the total number of variable sites (second column) and the mean Fst for each window (following columns).

### dxy and dxx

The mean differentiation between both haplotypes in males and females is calculated by dividing the number variable sites in males (dxy) and females (dxx) by the size of the sliding window.

The output file includes the generation in the first column and the mean dxy/dxx for each window in the following columns.

### male and female recombination points

This file records all recombination positions in a generation in males or females

Please not that:

* recombination positions are not exact in the dynamic model, because recombination probabilities are integrated over blocks. This may result in a lower number of recombination events for the first positions of the simulated region, which is not an error in the simulations, but a limitation in the way recombination events are recorded.

* recombination events that fall outside the simulated region will be recorded as genome_size+1

The output file includes the generation in the first column and a list of all recombination events in this generation in the second column.

### male recombination probability vectors

This output file contains the probability vectors that are used to draw recombination points in males. Each vector output was used to in a male meiosis create a gamete during the generation defined by *sstat_gen*. Only gametes and males which contribute to the next generation are output.  As this file can get rather large (several hundred mb!), it is optional.
Vectors are output in blocks for each generation, starting with #gen *generation* and followed by one probability vector per line.
For the dynamic model, recombination probabilities are integrated over blocks between variable sites. One block is output as: *begin*-*end*:*recombination probability*, blocks are tab-separated.
For the hotspot model recombination probabilities are output for each hotspot as: *hotspot position*:*recombination probability*. Note that recombination probabilities are for the block as a whole, to get per-base proababilities they must be divided by the block size.

## Running the code:
The software can be build (under linux) using:

g++ main.cpp -std=c++11 -o schrom_sim

It is run as:

./schrom_sim example.conf

