
"""
attempt to fix bug that nonrec-regions sometimes appeared to re-gain (some) recombination:
	- go for true infinite-sites model by excluding mutations in sites that are already variants. 
		: if re-gain of rec is caused by mutation towards of a "1" allele on the X at a position that was previously fixed heterozygous! 

read also: Ortiz-Barrientos et al. Recombination rate evolution and the origin of species.

ideas to improve algorithm:
	https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2323826/
	infinite sites+array of integers instead of sequence strings, "look-ahead" to only work on chroms that can leave offspring in future generations..

=> also here need to implement:
	- recording of Fst
	- recording of Rec-rate
	- translate allele-string into actual sequences (for exporting examples to alignment viewers or so.)

- a population of diploids, size can change but is capped at NMAX by random culling (carrying capacity) => better implementation of K: a probability of survival that is antiproportional to the number of indiv => fluctuation around NMAX should result.
- generations are non-overlapping, asexual reproduction is not considered, pre-reproductive mortality is zero
- no spatial structure / mates are drawn at random (but respecting ovule - sperm pairing)
	spatial structure implementation: X and Y coordinates as additional traits of individuals!! reflecting boundaries of grid, or "circular". Dispersal drawing two distances from dispersal kernel distribution.


principal steps of each cycle:
- cap the population size
- . . .
"""

"""
generate starting population of chromosomes by backward simulations under the coalescent for a constant-size unstructured population

mspath=/Users/scharmann/Downloads/msdir
sgpath=/Users/scharmann/Downloads/Seq-Gen-1.3.4/source

# theta = 4*Ne*mu*sites = 4*100000*7*10^-9*10000 = 28

cd /Users/scharmann/Documents/sexchrom_evol_recomb_simulations

$mspath/ms 200 1 -t 28.0 -r 28.0 10000 -T | grep '('> trees.tre

$sgpath/seq-gen -mGTR -s 0.01 -l 10000 -p $(cat trees.tre | wc -l) -f 0.15 0.35 0.15 0.35 -i 0.0 -a 5.0 -g 3 < trees.tre> start.aln

# maybe better another sequence simulator, which can evolve indels and inversions as well??
indelible - indels but no inversions

Hmm not lets not mess with these - just turn some of the SNPs into indels and inversions myself!

https://www.nature.com/articles/ng.911#s1

"Across all 80 strains, we identified 4,902,039 SNPs and 810,467 small indels" => "1- to 20-bp insertions and deletions"
=> ratio of small indels to SNPs in Ath is about 810467/4902039 = 16.5%

https://www.sciencedirect.com/science/article/pii/S0092867416306675
"After filtering, the nuclear genomes contained 10,707,430 biallelic SNPs and 1,424,879 small-scale indels (up to 40 bp)."
=> 1424879/10707430 = 13.3%



"""

import numpy, random, os, argparse


# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-m", required=True, help="mutation rate", metavar="FLOAT")
	parser.add_argument("-r", required=True, help="nominal recombination rate", metavar="FLOAT")
	parser.add_argument("-N", required=True, help="population size in number of diploids", metavar="INT")
	parser.add_argument("-g", required=True, help="number of generations", metavar="INT")
	parser.add_argument("-s", required=True, help="number of sites / length of chromosomes", metavar="INT")
	parser.add_argument("-a", required=True, help="'area' affected by reduction in recombination probability", metavar="INT")
	parser.add_argument("-logf", required=True, help="report statistics every logf generations", metavar="INT")
	
	args = parser.parse_args()

	
	return args



def make_dioecious_pop (N_per_cell, grid_xlen, grid_ylen, sexdet_site):
	
	"""
	each diploid individual has the following structure:
	[cnt, chrom1, chrom2, sexdet_site_chrom1, sexdet_site_chrom2]
	each chrom is a list of alleles: [0,0,0,1,0,0,0,...]
	
	the position of these alleles along the chromosome is given by another vector, which contains positions of all variants in the population:
	variant_idxes = [11,123,8675,98761,100001,...]
	invariant sites are NOT tracked!

	the effect sizes on recombination of each variant are stored in a third array:
	variant_effects = [[hom1,het,hom2],...]
	where each variant has three pre-defined effect sizes, for each of the three possible genotypes, e.g. no effect in either homzygote but reduction by 50% if heterozygous [0,-0.5,0]
	
	the areas affected by change in recombination of each variant are stored in a fourth array:
	variant_areas = [100,...]
	
	- Y-chromosomes have allele "1" at sexdet_site, X-chromosomes have alelle "0"
	- position of the sex-determining site is stored within each individual: 
		this is necessary if structural mutations can shift the position of this site in the alignment
	- If the sexdet site happens to be within an inversion (and only then), 
		it is possible that an individual has two chromosomes with the sexdet sites 
		at different positions in the alignment. Therefore, the sexdet site has to 
		be stored with each chromosome specifically: idv[3] and idv[4]
	"""
	variant_idxes = [sexdet_site]
	variant_effects = [[0.0,0.0,0.0]] # sex-determinator itself has NO effect on recombination rate!
	variant_areas = [0] 
	pop = []
	cnt = 0
	for x in range(grid_xlen):
		pop.append([])
		for y in range(grid_ylen):
			pop[x].append([])
			pop[x][y] = []
			for i in range(N_per_cell):
				cnt += 1
				is_male = numpy.random.choice([True,False],p=[0.5,0.5])
				if is_male:
					idv = [cnt, [1], [0], sexdet_site, sexdet_site]
				else:
					idv = [cnt, [0], [0], sexdet_site, sexdet_site]
				pop[x][y].append( idv )		
	return pop, variant_idxes, variant_effects, variant_areas		


def make_sanitycheck_pop (N_per_cell, grid_xlen, grid_ylen, sexdet_site):
	"""
	fixed nonrec Y-bloc	
	"""
	variant_idxes = [sexdet_site]
	cnt = 1
	for i in range(20):
		variant_idxes.append(sexdet_site+cnt)
		cnt += 50
	
	variant_effects = [[0.0,0.0,0.0]] + [[0.0,-1.0,0.0]]*20 # sex-determinator itself has NO effect on recombination rate!
	variant_areas = [0] + [500]*20
	pop = []
	cnt = 0
	for x in range(grid_xlen):
		pop.append([])
		for y in range(grid_ylen):
			pop[x].append([])
			pop[x][y] = []
			for i in range(N_per_cell):
				cnt += 1
				is_male = numpy.random.choice([True,False],p=[0.5,0.5])
				if is_male:
					idv = [cnt, [1]*21, [0]*21, sexdet_site, sexdet_site]
				else:
					idv = [cnt, [0]*21, [0]*21, sexdet_site, sexdet_site]
				pop[x][y].append( idv )		
	return pop, variant_idxes, variant_effects , variant_areas	


def enforce_carrying_capacity (pop, grid_xlen, grid_ylen, K_per_gridcell):	
	
	cnt_before = 0
	cnt_after = 0
	pop_culled = []
	for x in range(grid_xlen):
		pop_culled.append([])
		for y in range(grid_ylen):
			pop_culled[x].append([])
			pop_culled[x][y] = []
			gridcell_residents = pop[x][y]
			if len( gridcell_residents ) > 0:
				cnt_before += len( gridcell_residents )
				mortality_prob = (len( gridcell_residents ) - K_per_gridcell) / len( gridcell_residents )
				if mortality_prob <= 0:
					mortality_prob = 0
					pop_culled[x][y] = gridcell_residents
					cnt_after += len(gridcell_residents)					
				else:
					# now the killing:
#					survivors = []
					f = numpy.random.choice([0,1], size=len( gridcell_residents ), replace=True, p=[1.0-mortality_prob,mortality_prob])
					survivors = [c for idx,c in enumerate(gridcell_residents) if f[idx] == 0]
					pop_culled[x][y] = survivors
					cnt_after += len(survivors)	
#	print cnt_before, cnt_after
	return pop_culled

	
def measure_grid_density (thepop, grid_xlen, grid_ylen):	
	densities = []
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			densities.append( len(thepop[x][y]) )
	#print sum(densities), numpy.mean(densities), numpy.std(densities)
	return sum(densities)


def disperse_male_gametes (pop_gametes, grid_xlen, grid_ylen, scale_parameter):
	# scale_parameter = 
	# https://docs.scipy.org/doc/numpy/reference/generated/numpy.random.exponential.html
	# for x > 0 and 0 elsewhere. \beta is the scale parameter, which is the inverse of the rate parameter \lambda = 1/\beta. 
	# The rate parameter is an alternative, widely used parameterization of the exponential distribution
	#
	# grid boundaries are teleporting dispersers to the other side of the grid.. => infinite space.
	total_n_male_gametes = 0
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			total_n_male_gametes += len( pop_gametes[x][y][0] )
#	print "total male gametes produced:	", total_n_male_gametes
	x_distances = numpy.round(numpy.random.exponential(scale=scale_parameter, size=total_n_male_gametes))
	f = numpy.random.choice([-1,1], size=total_n_male_gametes, replace=True)
	x_distances = x_distances*f 
	y_distances = numpy.round(numpy.random.exponential(scale=scale_parameter, size=total_n_male_gametes))
	f = numpy.random.choice([-1,1], size=total_n_male_gametes, replace=True)
	y_distances = y_distances*f 
	
	dispersed_pop_gametes = [] # is a three-dimensional array, pop_gametes[x][y][M/F]
	for x in range(grid_xlen):
		dispersed_pop_gametes.append([ [[],pop_gametes[x][y][1]] for y in range(grid_ylen) ]) # F gametes are copied	
	cnt = 0
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			gridcell_residents = pop_gametes[x][y][0]
			if len( gridcell_residents ) > 0:
				male_gametes = []
				for i in gridcell_residents:
					x_movement = x_distances[cnt]
					y_movement = y_distances[cnt]
					if abs(x_movement) > grid_xlen: # adjust movements that traverse the grid more than 1 times!
						if x_movement < 0:
							corrmov = round( float( "-0." + str((x_movement/float(grid_xlen))).split(".")[1] )*grid_xlen )
						else:
							corrmov = round( float( "0." + str((x_movement/float(grid_xlen))).split(".")[1] )*grid_xlen )
#						print "orig_x_movement:	", x_movement, "corrected:	", corrmov
						x_movement = corrmov
					if abs(y_movement) > grid_ylen: # adjust movements that traverse the grid more than 1 times!
						if y_movement < 0:
							corrmov = round( float( "-0." + str((y_movement/float(grid_ylen))).split(".")[1] )*grid_ylen )
						else:
							corrmov = round( float( "0." + str((y_movement/float(grid_ylen))).split(".")[1] )*grid_ylen )
#						print "orig_y_movement:	", y_movement, "corrected:	", corrmov
						y_movement = corrmov
					old_coords = [x,y]
					new_coords = [old_coords[0] + x_movement, old_coords[1] + y_movement]
#					print "before:	", old_coords, new_coords
					if new_coords[0] > grid_xlen-1:
						new_coords[0] = new_coords[0] - grid_xlen
					elif new_coords[0] < 0:
						new_coords[0] = new_coords[0] + grid_xlen
					if new_coords[1] > grid_ylen-1:
						new_coords[1] = new_coords[1] - grid_ylen
					elif new_coords[1] < 0:
						new_coords[1] = new_coords[1] + grid_ylen							  
#					print "after:	", old_coords, new_coords
					dispersed_pop_gametes[int(new_coords[0])][int(new_coords[1])][0].append(i)
					cnt += 1
	# check that they are still all there!
#	cnt = 0
#	for x in range(grid_xlen):
#		for y in range(grid_ylen):
#			cnt += len(pop_gametes[x][y][0])
#			cnt += len(pop_gametes[x][y][1])
#	print "total gametes after dispersal:	", cnt
	return dispersed_pop_gametes				



def mate (dispersed_pop_gametes, grid_xlen, grid_ylen):
	# now mate
	zygotes = []
	cnt = 0
	total_f_gametes = 0
	for x in range(grid_xlen):
		zygotes.append([])
		for y in range(grid_ylen):
			zygotes[x].append([])
			zygotes[x][y] = []
			m_gametes = sorted(dispersed_pop_gametes[x][y][0], key=lambda k: random.random())
			f_gametes = sorted(dispersed_pop_gametes[x][y][1], key=lambda k: random.random())
			total_f_gametes += len(f_gametes)
#			print len(m_gametes), len(f_gametes)
			for m,f in zip(m_gametes, f_gametes):
				cnt += 1
				zygotes[x][y].append( [cnt, m[0], f[0], m[1], f[1]] )
#	print "total viable f gametes produced:	", total_f_gametes
	return zygotes
				



##
def do_meiosis(pop, variant_idxes, variant_effects, variant_areas, grid_xlen, grid_ylen, n_gametes_per_individual, chrom_length, nominal_rec_rate_array, record_rec_probs):
	"""
	for each diploid, do meiosis and return n_gametes_per_individual
	recording of XY and XX recombination probabilities is costly and thus optional; to be called only sporadically
	"""
	
	if record_rec_probs:
		X_Y_recomb_probs = numpy.zeros(chrom_length)
		X_X_recomb_probs = numpy.zeros(chrom_length)
	
	n_meioses = int(n_gametes_per_individual*0.5)
	
	idxes_of_chrom = range( chrom_length+1 )
	nominal_probs_per_site = numpy.array(nominal_rec_rate_array) 
		
	outpop = []
	pop_gametes = [] 
	# pop_gametes is a four-dimensional array, pop_gametes[x][y][[list_of_M_gametes],[list_of_F_gametes]]
	# each gamete is a list with the sequence, and second element the position of the sexdet site. list_of_M_gametes = [[seq,sexdet_site],...]
	
	for x in range(grid_xlen):
		pop_gametes.append([])
		for y in range(grid_ylen):
			pop_gametes[x].append([[],[]])
			
	male_meioses_cnt = 0
	female_meioses_cnt = 0
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for i in pop[x][y]:
				# is male or female?
				is_male = False
				is_female = False
				sexdet_alleles = i[1][variant_idxes.index(i[3])] + i[2][variant_idxes.index(i[4])] # maybe slow? could be done in abetter way!
				if sexdet_alleles == 1:
					is_male = True
					male_meioses_cnt += n_meioses
				elif sexdet_alleles == 0:
					is_female = True
					female_meioses_cnt += n_meioses
				else:
					print "is neuter, no gametes"
					continue			
				# now do meioses
				for meiosis in range(n_meioses): 
					
					if record_rec_probs:
						specific_probs_per_site = numpy.copy(nominal_probs_per_site)
					
					# rec only has consequences for the variable sites; i.e. the exact basepair on a sequence where it falls does not matter, it matters only if it fell anywhere into the interval (bloc) between pairs of consecutive variants. => use this to make code more efficient; avoid very long numpy arrays of rec. probs per each site; try to use blocs in between SNPs instead: at most a few hundred blocs per meiosis instead of 10k-millions of sites to consider.
					# there must be at least two SNPs (heteroz. sites) between the chroms in a meiosis for recombination to have any consequences at all.
					
					### get the cumulative probabilities between CONSECUTIVE pairs of SNPs in the meiosis
					# 1. find variants in this meiosis, also their effects
					snp_pos = []
					effects = [0.0] # first pos of chrom has effect 0
					the_areas = [0]
					for vidx, variant in enumerate(variant_idxes):
						genotype = i[1][vidx] + i[2][vidx] # sum of alleles can be used as index to retrieve genotye effect size, HA!
						if genotype == 1: # heterozygous in this meiosis
							snp_pos.append(variant)
						
							effects.append( nominal_probs_per_site[variant]*variant_effects[vidx][genotype] )
							the_areas.append(variant_areas[vidx])
							if record_rec_probs:
								halfarea = int(0.5*variant_areas[vidx])
								leftend = variant-halfarea
								if leftend < 0:
									leftend = 0
								rightend = variant+halfarea
								if rightend > chrom_length+1:
									rightend = chrom_length+1
								specific_probs_per_site[leftend:rightend] += (nominal_probs_per_site[variant]*variant_effects[vidx][genotype])
					
					effects.append(0.0) # last pos of chrom has effect 0 
					the_areas.append(0)
#					print "the effects", effects
					
					# record rec probs if necessary:
					if record_rec_probs:
						specific_probs_per_site[specific_probs_per_site<0] = 0.0
						if is_male:
							X_Y_recomb_probs += specific_probs_per_site
						elif is_female:
							X_X_recomb_probs += specific_probs_per_site

					
					if len(snp_pos) < 2:
						# there must be at least two SNPs (heteroz. sites) between the chroms in a
						# meiosis for recombination to have any consequences at all. If not, gametes = parental chroms.
						# print "NA"
						gamete1 = [i[1],i[3]]
						gamete2 = [i[2],i[4]]
					else:	
						# 2. get consecutive pairs of variants, also including the first and last sites of chrom (= intervals)
						intervals = [[0,snp_pos[0]]]
						for first, second in zip(snp_pos, snp_pos[1:]):
							intervals.append([first,second])
						intervals.append([snp_pos[-1],chrom_length])
						
						# 2.a construct a list of the areas_affected
						areas_affected = [[snp-int(0.5*the_areas[sidx]),snp+int(0.5*the_areas[sidx])] for sidx,snp in enumerate(snp_pos)]
						areas_affected = [[c,d] if c >= 0 else [0,d] for [c,d] in areas_affected]
						areas_affected = [[c,d] if d <= chrom_length else [c,chrom_length] for [c,d] in areas_affected]
#						print "areas", areas_affected
						
						# 3. now sum up probabilities of rec in each interval
						interval_probabilities = []
						for p in intervals:
							prob = numpy.sum(nominal_probs_per_site[p[0]:p[1]])
							
							# get effects of any variants that lie in or outside the interval, but
							# whose "area" overlaps with the interval! I worked out the following decision tree to determine the overlap of any two segments (a,p) along a line:
							cntb = 0
							for a in areas_affected:
								if a[0] < p[0]:
									if a[1] < p[1]:
										if a[1] < p[0]:
											overlap = 0
										else:
											overlap = a[1] - p[0]
									else:
										overlap = p[1] - p[0]
								else:												
									if a[1] < p[1]:
										overlap = a[1] - a[0]
									else:
										if a[0] < p[1]:
											overlap = p[1] - a[0]
										else:
											overlap = 0
								cntb += 1
#								print p, a, overlap
								effect1 = effects[cntb]
								prob += effect1*overlap
								
							if prob < 0.0:
								prob = 0.0
							interval_probabilities.append(prob)
#							print p, prob
						
						# 4. now choose recombination interval
						
						# append a dummy interval that captures all remaining prob, => probability to NOT recombine along this chromosomal segment in this meiosis						
						interval_probabilities.append(1.0-sum(interval_probabilities))
#						print interval_probabilities
						rec_interval = numpy.random.choice(range(len(interval_probabilities)), p = interval_probabilities)
#						print rec_interval
						if rec_interval == len(interval_probabilities)-1 or rec_interval == len(interval_probabilities)-2 or rec_interval == 0:
							# dummy interval chosen, or final interval chosen, or first interval chosen => no effective recombination!
							# print "NA"
							gamete1 = [i[1],i[3]]
							gamete2 = [i[2],i[4]]
						else:
							# print "recombination at ", rec_interval
							# build the recombinant gametes; only the first variant matters now for slicing the lists
							pair_of_variants_around_rec_site = intervals[rec_interval]
							idx_v1 = variant_idxes.index(pair_of_variants_around_rec_site[0]) # maybe possible fastwe with another construct?
							gamete1 = [i[1][:idx_v1+1] + i[2][idx_v1+1:],i[3]]
							gamete2 = [i[2][:idx_v1+1] + i[1][idx_v1+1:],i[4]]
		
							# now determine if the sexdet_pos switched phase:
							if idx_v1 >= variant_idxes.index(i[3]): # recombination site was downstream of sexdet site: fine as-is!
								None
							else: # recombination site was upstream of sexdet site: sexdet site has switched phase to that of other parental chrom.
								gamete1 = [gamete1[0],i[4]]
							if idx_v1 >= variant_idxes.index(i[4]): # recombination site was downstream of sexdet site: fine as-is!
								None
							else: # recombination site was upstream of sexdet site: sexdet site has switched phase to that of other parental chrom.
								gamete2 = [gamete2[0],i[4]]					

					# finally append gametes to male, female gamete pools of that grid cell	
					if is_male:
						pop_gametes[x][y][0].append(gamete1)
						pop_gametes[x][y][0].append(gamete2)
					elif is_female:
						pop_gametes[x][y][1].append(gamete1)
						pop_gametes[x][y][1].append(gamete2)
#					print len(gamete1[0]), len(gamete2[0])
	
	if record_rec_probs:
		avg_X_Y_rec_probs = X_Y_recomb_probs / male_meioses_cnt 
		avg_X_X_rec_probs = X_X_recomb_probs / female_meioses_cnt
	else:
		avg_X_Y_rec_probs = "" 
		avg_X_X_rec_probs = ""
		
	return pop_gametes, avg_X_Y_rec_probs, avg_X_X_rec_probs


def disperse_pop (inpop, grid_xlen, grid_ylen, scale_parameter): 
	
	### NEEDS also seed dispersal because empty gridcells will otherwise forever remain empty!
	total_n = 0
	dispersed_pop = [] # is a three-dimensional array, pop_gametes[x][y][the individuals]
	for x in range(grid_xlen):
		dispersed_pop.append([])
		for y in range(grid_ylen):
			dispersed_pop[x].append([])
			dispersed_pop[x][y] = []
			total_n += len( inpop[x][y] )
#	print "total offspring produced:	", total_n
	x_distances = numpy.round(numpy.random.exponential(scale=scale_parameter, size=total_n))
	f = numpy.random.choice([-1,1], size=total_n, replace=True)
	x_distances = x_distances*f 
	y_distances = numpy.round(numpy.random.exponential(scale=scale_parameter, size=total_n))
	f = numpy.random.choice([-1,1], size=total_n, replace=True)
	y_distances = y_distances*f 
	cnt = 0
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			gridcell_residents = inpop[x][y]
			if len( gridcell_residents ) > 0:
				for i in gridcell_residents:
					x_movement = x_distances[cnt]
					y_movement = y_distances[cnt]
					if abs(x_movement) > grid_xlen: # adjust movements that traverse the grid more than 1 times!
						if x_movement < 0:
							corrmov = round( float( "-0." + str((x_movement/float(grid_xlen))).split(".")[1] )*grid_xlen )
						else:
							corrmov = round( float( "0." + str((x_movement/float(grid_xlen))).split(".")[1] )*grid_xlen )
#						print "orig_x_movement:	", x_movement, "corrected:	", corrmov
						x_movement = corrmov
					if abs(y_movement) > grid_ylen: # adjust movements that traverse the grid more than 1 times!
						if y_movement < 0:
							corrmov = round( float( "-0." + str((y_movement/float(grid_ylen))).split(".")[1] )*grid_ylen )
						else:
							corrmov = round( float( "0." + str((y_movement/float(grid_ylen))).split(".")[1] )*grid_ylen )
#						print "orig_y_movement:	", y_movement, "corrected:	", corrmov
						y_movement = corrmov
					old_coords = [x,y]
					new_coords = [old_coords[0] + x_movement, old_coords[1] + y_movement]
					if new_coords[0] > grid_xlen-1:
						new_coords[0] = new_coords[0] - grid_xlen
					elif new_coords[0] < 0:
						new_coords[0] = new_coords[0] + grid_xlen
					if new_coords[1] > grid_ylen-1:
						new_coords[1] = new_coords[1] - grid_ylen
					elif new_coords[1] < 0:
						new_coords[1] = new_coords[1] + grid_ylen
#					print x,y, new_coords							  
					dispersed_pop[int(new_coords[0])][int(new_coords[1])].append(i)
					cnt += 1
	return dispersed_pop				



def random_choice_with_exceptions (sample_from, sample_n, exclusion_set):
	
	sampled_items = []
#	cnt = 0
	while len(sampled_items) < sample_n:
#		cnt += 1
		rc = numpy.random.choice(sample_from, sample_n)
		sampled_items = [i for i in rc if not i in exclusion_set]
	
#	print "number of calls:	", cnt
	return numpy.array(sampled_items)

def mutate_gametes(ingametes, grid_xlen, grid_ylen, mu, chrom_length, variant_idxes, variant_effects, variant_areas, new_area_size):
		
	idxes_of_chrom = range( chrom_length )
	
	n_gametes = 0
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for s in [0,1]:
				n_gametes += len(ingametes[x][y][s])	
	
	# draw all mutations; if total expected < 1: mutate only 1 but with probability = total_expected_mutations.
	total_expected_mutations = n_gametes*chrom_length*mu
	if total_expected_mutations > 1.0:
		total_expected_mutations = int(total_expected_mutations) # so this is not entirely correct!
		proceed = True
	else:
		proceed = numpy.random.choice([True, False], p = [total_expected_mutations, 1.0-total_expected_mutations])
		total_expected_mutations = 1
	if proceed:
		mutating_sites = random_choice_with_exceptions(idxes_of_chrom, total_expected_mutations, set(variant_idxes))
	#	heteroz_effect_sizes = numpy.random.normal(loc = 0.0, scale = 0.1, size = total_expected_mutations)
	#	heteroz_effect_sizes = numpy.array(heteroz_effect_sizes, dtype = "float64")
	#	heteroz_effect_sizes = -1*numpy.random.beta(0.3, 1.0, size = total_expected_mutations)
	#	heteroz_effect_sizes = -1*numpy.random.beta(1, 0.01, size = total_expected_mutations)
		heteroz_effect_sizes = numpy.array([-1.0]*total_expected_mutations)
	#	print "mutation effects", heteroz_effect_sizes

	#	heteroz_effect_areas = [int(fu) for fu in numpy.ceil(numpy.random.beta(0.5, 10.0, size = total_expected_mutations)*chrom_length)]
	#	heteroz_effect_areas = [int(fu) for fu in numpy.ceil(numpy.random.beta(0.1, 30.0, size = total_expected_mutations)*chrom_length)]
		## R-histogram of this distribution: hist(rbeta(10000,0.5,10)*100000, breaks = 100)
	#	heteroz_effect_areas = heteroz_effect_areas.astype(int)
		heteroz_effect_areas = [new_area_size]*total_expected_mutations
	#	print heteroz_effect_areas
	
		# draw gridcell and sex where new alleles arise
		random_x_coordinates = numpy.random.choice(range(grid_xlen), total_expected_mutations)
		random_y_coordinates = numpy.random.choice(range(grid_ylen), total_expected_mutations)
		random_sex = numpy.random.choice([0,1], total_expected_mutations)
	
		old_gamete_length = len(variant_idxes)
	
		addsites = [0]*total_expected_mutations
		
		# create new outgametes
		outgametes = []
		for x in range(grid_xlen):
			outgametes.append([])
			for y in range(grid_ylen):
				outgametes[x].append([])
				for s in [0,1]:
					outgametes[x][y].append([])
				
					for ingamete in ingametes[x][y][s]:
						seq = ingamete[0][:] ## have to COPY the list, otherwise will change the list in-place and cause havoc!!
						seq += addsites
						outgametes[x][y][s].append([seq,ingamete[1]])

		# now throw the mutations onto the gametes!
		cnt = 0
		for mut in sorted(mutating_sites):
	
			idx_of_gamete_in_sexpool_of_gridcell = random.choice(range(len( outgametes[random_x_coordinates[cnt]][random_y_coordinates[cnt]][random_sex[cnt]] )))
			outgametes[random_x_coordinates[cnt]][random_y_coordinates[cnt]][random_sex[cnt]][idx_of_gamete_in_sexpool_of_gridcell][0][old_gamete_length+cnt] = 1
			cnt += 1

		# update variant_effects, variant_idxes, variant_areas
		variant_idxes += sorted(mutating_sites)
		variant_effects += [[0.0,x,0.0] for x in heteroz_effect_sizes]
		variant_areas += heteroz_effect_areas
	#	print variant_idxes
	
		# sort everything
		# from https://stackoverflow.com/questions/7851077/how-to-return-index-of-a-sorted-list
		indexes_after_sorting = sorted(range(len(variant_idxes)), key=variant_idxes.__getitem__)

		variant_idxes = [variant_idxes[i] for i in indexes_after_sorting[:]]
		variant_effects = [variant_effects[i] for i in indexes_after_sorting[:]]
		variant_areas = [variant_areas[i] for i in indexes_after_sorting[:]]
	#	print variant_idxes
	
		outgametes_newly_sorted = []
		for x in range(grid_xlen):
			outgametes_newly_sorted.append([])
			for y in range(grid_ylen):
				outgametes_newly_sorted[x].append([])
				for s in [0,1]:
					outgametes_newly_sorted[x][y].append([])
				
					for ingamete in outgametes[x][y][s]:
						sortedseq = [ingamete[0][i] for i in indexes_after_sorting]
						outgametes_newly_sorted[x][y][s].append([sortedseq,ingamete[1]])
	
	#	print outgametes[0][0][0][0]	
	#	print outgametes_newly_sorted[0][0][0][0]	
		return outgametes_newly_sorted, variant_idxes, variant_effects, variant_areas
	else:
		return ingametes, variant_idxes, variant_effects, variant_areas		

def clean_from_invariant_sites (pop, variant_idxes, variant_effects, variant_areas, grid_xlen, grid_ylen):
	
	## BUT I should not remove recombination-effectors if they go to fixation!!? How to deal with this shit?? 
	## => Modify the nominal-rec rate??
	## => have to track an array with the nominal rec-rate, because it can change over time and change specifically for each site.
	## BUT if only heterozygous state has an effect, IT IS APPROPRIATE to remove the effects of alleles that become fixed! :)
	
	# find invariant sites by checking all chromosomes in pop:
	all_chroms = []
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for idv in pop[x][y]:
				all_chroms.append( idv[1] )
				all_chroms.append( idv[2] )
	
	total_chrom_cnt = len(all_chroms)
	
	invariant_idxes = []
	for idx in range(len(all_chroms[0])):
		allsum = sum([chr[idx] for chr in all_chroms])
		if allsum == 0:
			invariant_idxes .append(idx)
		elif allsum == 2*total_chrom_cnt:
			invariant_idxes .append(idx)
		else:
			None
	
	variant_idxes = del_list_numpy(variant_idxes, invariant_idxes)
	variant_effects = del_list_inplace(variant_effects, invariant_idxes)
	variant_areas = del_list_inplace(variant_areas, invariant_idxes)
	
	cleaned_pop = []
	for x in range(grid_xlen):
		cleaned_pop.append([])
		for y in range(grid_ylen):
			cleaned_pop[x].append([])
			for idv in pop[x][y]:
				chr1 = del_list_numpy(idv[1], invariant_idxes)
				chr2 = del_list_numpy(idv[2], invariant_idxes)
				cleaned_pop[x][y].append([idv[0],chr1,chr2,idv[3],idv[4]])
#				print variant_idxes
#				print chr1
#				print chr2
	 	
	return cleaned_pop, variant_idxes, variant_effects, variant_areas


def del_list_numpy(l, id_to_del):
	
	# according to https://stackoverflow.com/questions/11303225/how-to-remove-multiple-indexes-from-a-list-at-the-same-time
	# this is the fastest way to delete many items from a list!
	arr = numpy.array(l, dtype='int32')
	return list(numpy.delete(arr, id_to_del))


def del_list_inplace(l, id_to_del):
	# according to https://stackoverflow.com/questions/11303225/how-to-remove-multiple-indexes-from-a-list-at-the-same-time
	# this is the 2nd fastest way to delete many items from a list; good for non-integer items
    for i in sorted(id_to_del, reverse=True):
        del(l[i])
    return l

def report_window_rec_probs (generation, window_size, step_size, avg_rec_rates, XY_or_XX):
	
	header = ["generation"]
	window_avgs = [generation]
	left_idx = 0
	right_idx = left_idx+window_size
	while right_idx <= len(avg_rec_rates):
		w_avg = numpy.mean(avg_rec_rates[left_idx:right_idx])
		window_avgs.append(w_avg)	
		header.append(left_idx+int(0.5*window_size))
		left_idx += step_size
		right_idx += step_size

	if os.path.isfile("logfile_" + XY_or_XX + "_rec_probs.txt"): 
		with open("logfile_" + XY_or_XX + "_rec_probs.txt", "a") as F:
			F.write("\t".join([str(x) for x in window_avgs])+ "\n")
	else:
		with open("logfile_" + XY_or_XX + "_rec_probs.txt", "w") as A:
			A.write("\t".join(["","site"]) + "\n")
			A.write("\t".join([str(x) for x in header]) + "\n")
			A.write("\t".join([str(x) for x in window_avgs])+ "\n")
		

def report_window_male_female_Fst_and_dxy (inpop, variant_idxes, grid_xlen, grid_ylen, sample_size, generation, window_size, step_size, chrom_length):

	# random sample from population
	all_indivs = []
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			all_indivs += inpop[x][y]
	
	males = [i for i in all_indivs if i[1][variant_idxes.index(i[3])]+i[2][variant_idxes.index(i[4])] == 1]
	females = [i for i in all_indivs if i[1][variant_idxes.index(i[3])]+i[2][variant_idxes.index(i[4])] == 0]
	
	msamples = numpy.random.choice(range(len(males)), size=int(sample_size*0.5))
	msamples = [[males[x][1],males[x][2]] for x in msamples]
	
	fsamples = numpy.random.choice(range(len(females)), size=int(sample_size*0.5))
	fsamples = [[females[x][1],females[x][2]] for x in fsamples]
	
	msamples_flat = [y for x in msamples for y in x]
	fsamples_flat = [y for x in fsamples for y in x]
	
	samples_flat = msamples_flat + fsamples_flat
		
	n_chroms_sampled = float(sample_size*2)
	
	# find bi-allelic SNPs and calculate their freqs
	SNP_pos = []
	SNP_freqs = []
	SNP_freqs_males = []
	SNP_freqs_females = []
	Fst_list = []
	dxy_list = []
	for idx, site in enumerate(variant_idxes):
		totsum = sum([x[idx] for x in samples_flat])
		if not totsum == 0 and not totsum == n_chroms_sampled: # is not fixed in the sample
			SNP_pos.append( site )
			freq_global = totsum / n_chroms_sampled
			freq_m = sum([x[idx] for x in msamples_flat]) / (n_chroms_sampled*0.5)
			freq_f = sum([x[idx] for x in fsamples_flat]) / (n_chroms_sampled*0.5)
			
			if freq_global > 0.5:
				maf_g = 1.0-freq_global
				freq_m = 1.0-freq_m
				freq_f = 1.0-freq_f
			else:
				maf_g = freq_global
			SNP_freqs.append(maf_g)
			SNP_freqs_males.append(freq_m)
			SNP_freqs_females.append(freq_f)
			
			pi_overall = maf_g * (1-maf_g)
			pi_m = freq_m * (1-freq_m)
			pi_f = freq_f * (1-freq_f)
			avg_pi_within = (pi_m + pi_f) / 2.0
			Fst = ( pi_overall - avg_pi_within ) / pi_overall 
			Fst_list.append(Fst)
			dxy = freq_m*(1-freq_f) + freq_f*(1-freq_m)
			dxy_list.append( dxy )
			
	#	

	# Fst in sliding windows
	header = ["generation","N_variants"]
	left_idx = 0
	right_idx = left_idx+window_size
	window_Fsts = [generation, len(Fst_list)]
	while right_idx < chrom_length :
		SNPs_in_window = [x for x in SNP_pos if x >= left_idx and x < right_idx]
#		print SNPs_in_window
		header.append(left_idx+int(0.5*window_size))
		left_idx += step_size
		right_idx += step_size
		mean_Fst = numpy.mean([ Fst_list[SNP_pos.index(x)] for x in SNPs_in_window ])
		window_Fsts.append( mean_Fst )

	if os.path.isfile("logfile_window_Fsts.txt"): 
		with open("logfile_window_Fsts.txt", "a") as F:
			F.write("\t".join([str(x) for x in window_Fsts])+ "\n")
	else:
		with open("logfile_window_Fsts.txt", "w") as A:
			A.write("\t".join(["","site"]) + "\n")
			A.write("\t".join([str(x) for x in header]) + "\n")
			A.write("\t".join([str(x) for x in window_Fsts])+ "\n")


	# dxy in sliding windows
	header = ["generation","N_variants"]
	left_idx = 0
	right_idx = left_idx+window_size
	window_dxys = [generation, len(dxy_list)]
	while right_idx < chrom_length :
		SNPs_in_window = [x for x in SNP_pos if x >= left_idx and x < right_idx]
#		print SNPs_in_window
		header.append(left_idx+int(0.5*window_size))
		left_idx += step_size
		right_idx += step_size
		mean_dxy = numpy.sum([ dxy_list[SNP_pos.index(x)] for x in SNPs_in_window ]) / float(window_size)
		window_dxys.append( mean_dxy )

	if os.path.isfile("logfile_window_dxy.txt"): 
		with open("logfile_window_dxy.txt", "a") as F:
			F.write("\t".join([str(x) for x in window_dxys])+ "\n")
	else:
		with open("logfile_window_dxy.txt", "w") as A:
			A.write("\t".join(["","site"]) + "\n")
			A.write("\t".join([str(x) for x in header]) + "\n")
			A.write("\t".join([str(x) for x in window_dxys])+ "\n")


def find_consecutive_zeros (sequence):
	
	zero_idxes = []
	for idx,x in enumerate(sequence):
		if x == 0.0:
			zero_idxes.append(idx)

	if not len(zero_idxes) == 0:
		consec_zero_blocs = []
		cnt = 0
		blk = []
		while cnt < len(zero_idxes)-1:
			if zero_idxes[cnt]+1 == zero_idxes[cnt+1]:
	#			print "consec"
				blk.append(zero_idxes[cnt])
			else:
	#			print "broke"
				blk.append(zero_idxes[cnt])
				consec_zero_blocs.append(blk)
				blk = []
			cnt += 1
		blk.append(zero_idxes[cnt])
		consec_zero_blocs.append(blk)
	else:
		consec_zero_blocs = [[]]
	
	return consec_zero_blocs


def report_size_of_nonrec_region (generation, avg_rec_rates, sexdet_site ):
		
	consec_zero_blocs = find_consecutive_zeros (avg_rec_rates)
	nonrec_length = 0
	if len(consec_zero_blocs[0]) == 0:
		SD_bloc = False
	else:
		# is SD-site among them?
		SD_bloc = False
		for bloc in consec_zero_blocs:
			start = bloc[0]
			end = bloc[-1]
			if start <= sexdet_site and end >= sexdet_site:
				nonrec_length = len(bloc)

	if os.path.isfile("logfile_nonrec_region_size.txt"): 
		with open("logfile_nonrec_region_size.txt", "a") as F:
			F.write("\t".join([str(generation), str(nonrec_length)])+ "\n")
	else:
		with open("logfile_nonrec_region_size.txt", "w") as A:
			A.write("\t".join(["generation","size_of_nonrec_region"]) + "\n")
			A.write("\t".join([str(generation), str(nonrec_length)])+ "\n")

			

###################### MAIN

args = get_commandline_arguments ()

gridlength_x = 1
gridlength_y = 1
N_per_cell = float(args.N)
exponential_decay_shape_param = 2.0


sexdet_site = int(0.5*float(args.s)) ## always placed at mid of the chromosome

# 2 events on a 10 Mb chromosome = 2e-7 per site
nominal_overall_rec_rate_per_site = float(args.r)

mu = float(args.m)


# length of chromosome
chrom_length = int(args.s)

logf = int(args.logf)

pop, variant_idxes, variant_effects, variant_areas = make_dioecious_pop (int(N_per_cell), gridlength_x , gridlength_y, sexdet_site)
#pop, variant_idxes, variant_effects, variant_areas = make_sanitycheck_pop (int(N_per_cell), gridlength_x , gridlength_y, sexdet_site)
nominal_rec_rate_array = numpy.array([nominal_overall_rec_rate_per_site]*chrom_length)

for generation in range(int(args.g)+1):
#	measure_chrom_lengths (pop, gridlength_x, gridlength_y)
	pop = enforce_carrying_capacity (pop, gridlength_x, gridlength_y, N_per_cell) # carrying capacity must be a float!
#	print "cleaning"
	pop, variant_idxes, variant_effects, variant_areas = clean_from_invariant_sites (pop, variant_idxes, variant_effects, variant_areas, gridlength_x, gridlength_y)
	n_total = measure_grid_density (pop, gridlength_x, gridlength_y)
	print "gen = ", generation, ", N = ", n_total, ", variants:	", len(variant_idxes), ", largest area of any variant: ", max(variant_areas)
	if n_total == 0.0:
		print "extinction at generation	", generation
		break
#	print "doing meiosis"
	if (float(generation)/logf).is_integer(): # every x generations report something
		pop_gametes, avg_X_Y_rec_probs, avg_X_X_rec_probs = do_meiosis(pop, variant_idxes, variant_effects, variant_areas, gridlength_x, gridlength_y, 4, chrom_length, nominal_rec_rate_array, True)
		report_window_rec_probs (generation, 500, 100, avg_X_Y_rec_probs, "XY")
		report_window_rec_probs (generation, 500, 100, avg_X_X_rec_probs, "XX")
		report_window_male_female_Fst_and_dxy (pop, variant_idxes, gridlength_x, gridlength_y, 20.0, generation, 500, 100, chrom_length)
		report_size_of_nonrec_region (generation, avg_X_Y_rec_probs, sexdet_site)
	else:
		pop_gametes, avg_X_Y_rec_probs, avg_X_X_rec_probs = do_meiosis(pop, variant_idxes, variant_effects, variant_areas, gridlength_x, gridlength_y, 4, chrom_length, nominal_rec_rate_array, False)
#	print "dispersing"
	dispersed_pop_gametes = disperse_male_gametes (pop_gametes, gridlength_x, gridlength_y, exponential_decay_shape_param)
#	print "mutating"
	dispersed_pop_gametes, variant_idxes, variant_effects, variant_areas = mutate_gametes(dispersed_pop_gametes, gridlength_x, gridlength_y, mu , chrom_length, variant_idxes, variant_effects, variant_areas, int(args.a))
#	print "mating"
	pop = mate (dispersed_pop_gametes, gridlength_x, gridlength_y)
#	print "dispersing again"
	pop = disperse_pop (pop, gridlength_x, gridlength_y, exponential_decay_shape_param)
#	print variant_idxes
#	print [x[1] for x in variant_effects]

print "Done!"



