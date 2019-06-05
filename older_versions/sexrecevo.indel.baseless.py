
"""
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
- there is no selfing -- outcrossing is already mandatory and hence not subject to evolution (leaving only sexual conflict and random drift as evolutionary factors)
- each individual is granted the same glucose budget upon birth
- budget can be spent on ovules, sperm, or both
- ovules and sperm have fixed costs (e.g. one ovule == 2 sperm)
- zygote development until birth (gestation) incurs a cost to mothers (maternal care / provisioning): fixed XY glucose per zygote => if a mother has more ovules fertilised than can be reared given its glucose budget reamining after gamete production, the excess zygotes will be randomly culled (aborted) . This should mean that mothers (irrespective of the % of ovules/sperm they made) are penalised for making too many gametes and leaving too few resources for gestation.
- allocation to ovule number and sperm number are principally independent, quantitative traits, each controlled by X bi-allelic loci with pure additivity 
- individual allocation to ovules and sperm may evolve.

- we observe / record: histogram of "% resources invested into sperm production"
	(bimodal at 0,1 : dioecy ; bimodal at c. 0,intermediate : gynodioecy ; bimodal at c. intermediate,1 : androdioecy, unimodal: only hermaphrodites. Also interesting: width of the distribution/variance is sex allocation in pop!!)
	total resources = start budget = sperm + ovules + postzygotic care
	



principal steps of each cycle:
- cap the population size
- calculate individual traits from genotypes: N_ovules, N_sperm, N_gestation
- 
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

import numpy, random


def write_phylip (concat_dict, outfile):
	
	with open(outfile, "w") as OUTFILE:
		ntaxa = len(concat_dict.keys())
		len_align = len(concat_dict[concat_dict.keys()[0]]) # gets value of first element in dictionary -> this OK since all seqs have same length
		header = str(ntaxa) + " " + str(len_align) + "\n"
		OUTFILE.write(header)
		for sample in sorted(concat_dict.keys()):
			out_line = sample + "    " + concat_dict[sample] + "\n"
			OUTFILE.write(out_line)
	OUTFILE.close()
	

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
	return pop, variant_idxes, variant_effects 		


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
def do_meiosis(pop, variant_idxes, variant_effects, grid_xlen, grid_ylen, n_gametes_per_individual, chrom_length, nominal_rec_rate_array, sliding_window_size):
	"""
	for each diploid, do meiosis and return n_gametes_per_individual
	- local_rec_rate = (events_per_meiosis/sites) * sequence_identity_in_window
	"""
	wlenfloat = float(sliding_window_size)
	n_meioses = int(n_gametes_per_individual*0.5)
	
	half_of_windowsize = int(0.5*sliding_window_size)
	
	idxes_of_chrom = range( chrom_length+1 )
	t = list(nominal_rec_rate_array)
	t.append(1.0)
	nominal_probs_per_site = numpy.array(t) # last site is a dummy to absorb all the remaining probability; no recombination if sampling this "pseudo-site"
		
	outpop = []
	pop_gametes = [] 
	# pop_gametes is a four-dimensional array, pop_gametes[x][y][[list_of_M_gametes],[list_of_F_gametes]]
	# each gamete is a list with the sequence, and second element the position of the sexdet site. list_of_M_gametes = [[seq,sexdet_site],...]
	
	for x in range(grid_xlen):
		pop_gametes.append([])
		for y in range(grid_ylen):
			pop_gametes[x].append([[],[]])
			
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for i in pop[x][y]:
				# is male or female?
				is_male = False
				is_female = False
				sexdet_alleles = i[1][variant_idxes.index(i[3])] + i[2][variant_idxes.index(i[4])] # maybe slow? could be done in abetter way!
				if sexdet_alleles == 1:
					is_male = True
				elif sexdet_alleles == 0:
					is_female = True
				else:
					print "is neuter, no gametes"
					continue			
				# now do meioses
				for meiosis in range(n_meioses): 
					# find site where recombination occurs, if any				
					# thus take the nominal_probs_per_site array and 
					# modify it by multiplication with the cumulative effect sizes of the variants in windows along it!
					specific_probs_per_site = nominal_probs_per_site
					cnt = -1
					for variant, effects in zip(variant_idxes, variant_effects):
						cnt += 1
#						print variant, effects
						genotype = i[1][cnt] + i[2][cnt] # sum of alleles can be used as index to retrieve genotye effect size, HA!
						effect = nominal_probs_per_site[variant]*effects[genotype]
						specific_probs_per_site[variant-half_of_windowsize:variant+half_of_windowsize] -= effect
					
					# if rec rate became negative, set to zero:
					specific_probs_per_site[specific_probs_per_site<0] = 0.0
					# let dummy-site have all remaining probability, such that it sums to 1:
					specific_probs_per_site[-1] = 1.0-numpy.sum(specific_probs_per_site[:len(specific_probs_per_site)-1])

					rec_site = numpy.random.choice(idxes_of_chrom, p = specific_probs_per_site)
					if rec_site == chrom_length: # dummy site, no recombination event
						gamete1 = [i[1],i[3]]
						gamete2 = [i[2],i[4]]
					else:
						print "recombination at ", rec_site
						for idx, variant in enumerate(variant_idxes):
							if variant > rec_site:
								gamete1 = [i[1][:idx] + i[2][idx:],i[3]]
								gamete2 = [i[2][:idx] + i[1][idx:],i[4]]
								break
							else: # recombination without consequence, because haplotypes identical beyond this point
								gamete1 = [i[1],i[3]]
								gamete2 = [i[2],i[4]]				
						# now determine if the sexdet_pos switched phase:
						if rec_site >= i[3]: # recombination site was downstream of sexdet site: fine as-is!
							None
						else: # recombination site was upstream of sexdet site: sexdet site has switched phase to that of other parental chrom.
							gamete1 = [gamete1[0],i[4]]
						if rec_site >= i[4]: # recombination site was downstream of sexdet site: fine as-is!
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
		
	return pop_gametes


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



def measure_chrom_lengths (inpop, grid_xlen, grid_ylen):
	
	chrom_lengths = []
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for i in inpop[x][y]:
				chrom_lengths.append(len(i[1]))
				chrom_lengths.append(len(i[2]))
	print max(chrom_lengths), min(chrom_lengths)
	

def pop_to_phylip (inpop, grid_xlen, grid_ylen, sample_size_M, sample_size_F, outfilename):

	# random sample
	all_indivs = []
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			all_indivs += inpop[x][y]
	
	males = [i for i in all_indivs if i[1][i[3]] == "A" or i[2][i[4]] == "A"]
	females = [i for i in all_indivs if i[1][i[3]]+i[2][i[4]] == "TT"]

	msamples = numpy.random.choice(range(len(males)), size=sample_size_M)
	msamples = [males[x] for x in msamples]
	fsamples = numpy.random.choice(range(len(females)), size=sample_size_F)
	fsamples = [females[x] for x in fsamples]
	
	samples = msamples+fsamples
	
	out_seq_dict = {}
	cnt = 100
	for s in samples:
		cnt += 1
		out_seq_dict[str(cnt)+"_1"] = s[1]
		out_seq_dict[str(cnt)+"_2"] = s[2]
	
	write_phylip (out_seq_dict, outfilename)




def mutate_gametes(ingametes, grid_xlen, grid_ylen, mu, chrom_length, variant_idxes, variant_effects):
		
	idxes_of_chrom = range( chrom_length )
	
	n_gametes = 0
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for s in [0,1]:
				n_gametes += len(ingametes[x][y][s])	
	
	# draw all mutations
	total_expected_mutations = int(n_gametes*chrom_length*mu)
	mutating_sites = numpy.random.choice(idxes_of_chrom, total_expected_mutations)
	heteroz_effect_sizes = numpy.random.normal(loc = 0.0, scale = 1, size = total_expected_mutations)
#	heteroz_effect_sizes = numpy.array(heteroz_effect_sizes, dtype = "float64")
	
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

	# update variant_effects and variant_idxes
	variant_idxes += sorted(mutating_sites)
	variant_effects += [[0.0,x,0.0] for x in heteroz_effect_sizes]
	
	return outgametes, variant_idxes, variant_effects
				

def clean_from_invariant_sites (pop, variant_idxes, variant_effects, grid_xlen, grid_ylen):
	
	## BUT I should not remove recombination-effectors if they go to fixation!! How to deal with this shit?? 
	## => Modify the nominal-rec rate??
	## => have to track an array with the nominal rec-rate, because it can change over time and change specifically for each site.
	
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
	
	cleaned_pop = []
	for x in range(grid_xlen):
		cleaned_pop.append([])
		for y in range(grid_ylen):
			cleaned_pop[x].append([])
			for idv in pop[x][y]:
				chr1 = del_list_numpy(idv[1], invariant_idxes)
				chr2 = del_list_numpy(idv[2], invariant_idxes)
				cleaned_pop[x][y].append([idv[0],chr1,chr2,idv[3],idv[4]])
	
	
	print chr1 
	return cleaned_pop, variant_idxes, variant_effects


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


######################

########

gridlength_x = 1
gridlength_y = 1
N_per_cell = 40
exponential_decay_shape_param = 2.0


sexdet_site = 1000 ## males = AT, females = TT ; other genotypes are neuters = no gametes! 

# 1 events on a 10 Mb chromosome = 2e-7 per site
nominal_overall_rec_rate_per_site = 0.9999e-7 # can be at most 1 event; => max(nominal_overall_rec_rate_per_site) = 1/sequence_length
sliding_window_size = 100

# mutation rate
mu = 2.5e-6


# length of chromosome
chrom_length = 10000

pop, variant_idxes, variant_effects = make_dioecious_pop (N_per_cell, gridlength_x , gridlength_y, sexdet_site)
nominal_rec_rate_array = numpy.array([nominal_overall_rec_rate_per_site]*chrom_length)
print nominal_rec_rate_array

for generation in range(500):
#	measure_chrom_lengths (pop, gridlength_x, gridlength_y)
	pop = enforce_carrying_capacity (pop, gridlength_x, gridlength_y, 40.0) # carrying capacity must be a float!
	pop, variant_idxes, variant_effects = clean_from_invariant_sites (pop, variant_idxes, variant_effects, gridlength_x, gridlength_y)
	n_total = measure_grid_density (pop, gridlength_x, gridlength_y)
	print generation, n_total
	if n_total == 0.0:
		print "extinction at generation	", generation
		break
	pop_gametes = do_meiosis(pop, variant_idxes, variant_effects, gridlength_x, gridlength_y, 4, chrom_length, nominal_rec_rate_array, sliding_window_size)
	dispersed_pop_gametes = disperse_male_gametes (pop_gametes, gridlength_x, gridlength_y, exponential_decay_shape_param)
	dispersed_pop_gametes, variant_idxes, variant_effects = mutate_gametes(dispersed_pop_gametes, gridlength_x, gridlength_y, mu , chrom_length, variant_idxes, variant_effects)
	pop = mate (dispersed_pop_gametes, gridlength_x, gridlength_y)
	pop = disperse_pop (pop, gridlength_x, gridlength_y, exponential_decay_shape_param)
	print variant_idxes
	print [round(x[1],3) for x in variant_effects]

print pop
#pop_to_phylip (pop, gridlength_x, gridlength_y, 10, 10, "start_XY_ident.Ne40.generation_"+ str(generation) +".aln")




	
#	# find bi-allelic SNPs and calculate their freq
#	SNP_pos = []
#	SNP_freqs = []
#	for idx in range(len(inseqs.values()[0])):
#		col = [x[idx] for x in inseqs.values()]
#		if len(set(col)) > 1:
#			is_indel = numpy.random.choice([True,False],p=[indel_rate,1.0-indel_rate])
#			if is_indel:
#				indel_pos.append(idx)
#				SNP_alleles_to_become_indel.append( numpy.random.choice(list(set(col)))) 






