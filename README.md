## sexrecevo 

## 2019-06-25
- now measuring speed of expansion of the nonrec region and doing preliminary exploration of parameter space: examples in expansion_rates_examples.tar.gz
- missing: a script to estimate the slope of a linear regression, of generation against the size of the nonrec region, averaged over all replicates of a parameter combination. The final output could perhaps be a barplot, with categories (parameter combinations) along the X and the height reflecting the speed of expansion, as e.g. "basepairs per generation".

## 2019-06-05
- Dan: to code up plotting modules, use the 100 iterations in "example.100_iterations.2019-06-05.tar.gz"
Should be self-explanatory. Just one parameter combination that was tweaked to show the desired effect.
Let me know if anything needs changing!

## 2019-06-06
- Plotting function done (./Plotting/Logplotter.py). Mathias to run some more sims before we discuss with others. 


# PSEUDOCODE of the model

0. 	create population
	- diploid individuals, each with two homologous chromosomes, each chromosome is an array of 0s and 1s
	- only variable loci are tracked
	- population is dioecious, with sex determined by a single locus (SD-locus) with two alleles, located at a predefined position
 	- SD-genotype [0,1] (X,Y) is male, SD-genotype [0,0] (XX) is female
 	- population starts without any variants except at the SD-locus
 	- SD-locus itself has no effect on recombination probability

for generation in n_generations:

	1.	enforce carrying capacity of population (random culling)
	
	2.	remove fixed variants
	
	3.	produce gametes (meiosis)
	
		- male individuals produce male gametes, females produce female gametes
		
		- 4 gametes per individual (2 meioses per individual)
		
		- in each meiosis, up to 1 crossover can take place
		
		- if crossover takes place, the location is chosen randomly for each position along the chromosome
		
		- the probability to recombine per site is given by a nominal background probability
		
		- the local probability is reduced if the chromosomes in the meiosis carry different alleles, i.e. probability is zero for an "area" along the chromosome affected by the variant
		
	4.	mutate gametes
	
		- variants are thrown randomly onto the gametes
		
		- variants reduce the local recombination rate when heterozygous, but have no effect when homozygous
		
		- each variant has two properties: location along the chromosome, "area" (number of bases) affected by the reduction in recombination   
		
	5.	produce zygotes (mate)
	
		- random draws from M + F gamete pools (only M+F are allowed)
		

## Principal Parameters:

N	population size

length of chromosome

mu	mutation rate

rho	nominal background recombination rate 

"area" (number of bases) affected by the reduction in recombination (could be a single value, or a distribution)
