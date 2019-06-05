
"""
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


def read_phylip(INFILE):
	
	indict = {}
	with open(INFILE, "r") as infile:
		infile.readline() # remove header
		for line in infile:
			fields = line.strip("\n").split()
			indict[fields[0]] = fields[1]
	return indict



def replace_SNPs_with_indels (inseqs, indel_rate, indel_lengths_scale):
	
	# find SNPs that will become indels
	indel_pos = []
	SNP_alleles_to_become_indel = []
	for idx in range(len(inseqs.values()[0])):
		col = [x[idx] for x in inseqs.values()]
		if len(set(col)) > 1:
			is_indel = numpy.random.choice([True,False],p=[indel_rate,1.0-indel_rate])
			if is_indel:
				indel_pos.append(idx)
				SNP_alleles_to_become_indel.append( numpy.random.choice(list(set(col)))) 
	
	# choose indel lengths
	n_indels = len(indel_pos)
	indel_lengths = numpy.round(numpy.random.exponential(scale=indel_lengths_scale, size=n_indels)+1)
	
	# delete at indel positions
	out_dict = {}
	for id,seq in inseqs.items():
		print id
		outseq = seq
		for i in range(len(indel_pos)):
			idx = indel_pos[i]
			allele = SNP_alleles_to_become_indel[i]
			indel_length = int(indel_lengths[i])
			if seq[idx] == allele:
				outseq = outseq[:idx] + "-"*indel_length + outseq[idx+indel_length:]
		out_dict[id] = outseq
	return out_dict
	


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
	

def make_dioecious_pop_imported_startseqs (N_per_cell, grid_xlen, grid_ylen, start_seqs, sexdet_site):
	
	"""
	Y-chromosomes have allele "A" at sexdet_site, X-chromosomes have alelle "T"
	- position of the sex-determining site is stored within each individual: 
		this is necessary if structural mutations can shift the position of this site in the alignment
	- If the sexdet site happens to be within an inversion (and only then), 
		it is possible that an individual has two chromosomes with the sexdet sites 
		at different positions in the alignment. Therefore, the sexdet site has to 
		be stored with each chromosome specifically: idv[3] and idv[4]
	"""
	
	pop = []
	cnt = 0
	seq_idx = 0
	for x in range(grid_xlen):
		pop.append([])
		for y in range(grid_ylen):
			pop[x].append([])
			pop[x][y] = []
			for i in range(N_per_cell):
				cnt += 1
				is_male = numpy.random.choice([True,False],p=[0.5,0.5])
				if is_male:
					t = list(start_seqs[seq_idx])
					t[sexdet_site] = "A"
					Y = "".join(t)
					t = list(start_seqs[seq_idx+1])
					t[sexdet_site] = "T"
					X = "".join(t)
					idv = [cnt, Y, X, sexdet_site, sexdet_site]
				else:
					t = list(start_seqs[seq_idx])
					t[sexdet_site] = "T"
					X1 = "".join(t)
					t = list(start_seqs[seq_idx+1])
					t[sexdet_site] = "T"
					X2 = "".join(t)
					idv = [cnt, X1, X2, sexdet_site, sexdet_site]
				seq_idx += 2
				pop[x][y].append( idv )		
	return pop		


def make_dioecious_pop_single_random_startseq (N_per_cell, grid_xlen, grid_ylen, nsites, sexdet_site):
	
	"""
	Y-chromosomes have allele "A" at sexdet_site, X-chromosomes have alelle "T"
	- position of the sex-determining site is stored within each individual: 
		this is necessary if structural mutations can shift the position of this site in the alignment
	- If the sexdet site happens to be within an inversion (and only then), 
		it is possible that an individual has two chromosomes with the sexdet sites 
		at different positions in the alignment. Therefore, the sexdet site has to 
		be stored with each chromosome specifically: idv[3] and idv[4]
	"""
	
	dna = ["A","G","C","T"]
	seq = ''.join(random.choice(dna) for i in range(nsites))
	
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
					t = list(seq)
					t[sexdet_site] = "A"
					Y = "".join(t)
					t = list(seq)
					t[sexdet_site] = "T"
					X = "".join(t)
					idv = [cnt, Y, X, sexdet_site, sexdet_site]
				else:
					t = list(seq)
					t[sexdet_site] = "T"
					X1 = "".join(t)
					t = list(seq)
					t[sexdet_site] = "T"
					X2 = "".join(t)
					idv = [cnt, X1, X2, sexdet_site, sexdet_site]
				pop[x][y].append( idv )		
	return pop		


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

def do_recombination(pop, grid_xlen, grid_ylen, nominal_overall_rec_rate_per_site, sliding_window_size):
	
	"""
	for each diploid, do meiosis and return diploid
	- local_rec_rate = (events_per_meiosis/sites) * sequence_identity_in_window
	"""
	wlenfloat = float(sliding_window_size)
	
	outpop = []
	for x in range(grid_xlen):
		outpop.append([])
		for y in range(grid_ylen):
			outpop[x].append([])
			for i in pop[x][y]:	
				does_recombine = numpy.random.choice([True,False],p=[nominal_overall_rec_rate_per_site*len(i[1]),1-nominal_overall_rec_rate_per_site*len(i[1])])
				if not does_recombine:
					outpop[x][y].append(i)
					continue
				else:				
					# per-site measure sequence identity
					pairwise_identity_per_site = []
					idx = 0
					while idx < len( i[1] ):
						if i[1][idx] == i[2][idx]:
							if i[1][idx] == "-":
								pairwise_identity_per_site.append(0)
							else:
								pairwise_identity_per_site.append(1)
						else:
							pairwise_identity_per_site.append(0)
						idx += 1
					# now window-slide; at ends extrapolate from value of ends' windows
					left_idx = 0
					center_idx = int(left_idx+sliding_window_size*0.5)
					right_idx = left_idx+sliding_window_size
				
					pairwise_identity_per_site_window_avg = []
					while right_idx < len( i[1] ):
						window_idenity_avg = sum(pairwise_identity_per_site[left_idx:right_idx])/wlenfloat
						pairwise_identity_per_site_window_avg.append(window_idenity_avg)
						left_idx += 1
						center_idx += 1
						right_idx += 1
					# add the ends:
					pairwise_identity_per_site_window_avg = [pairwise_identity_per_site_window_avg[0]]*int(0.5*sliding_window_size) + pairwise_identity_per_site_window_avg + [pairwise_identity_per_site_window_avg[-1]]*int(0.5*sliding_window_size)
	#				print len(pairwise_identity_per_site_window_avg), pairwise_identity_per_site_window_avg			
					recombination_probabilities_per_site = numpy.array(pairwise_identity_per_site_window_avg)/sum(pairwise_identity_per_site_window_avg)
					# now choose recombination sites:
					rec_site = numpy.random.choice(range( len( i[1] ) ), p = recombination_probabilities_per_site)
					print "recombination at ", rec_site
					recombinant_1 = i[1][:rec_site] + i[2][rec_site:] 
					recombinant_2 = i[2][:rec_site] + i[1][rec_site:]
					outpop[x][y].append([i[0],recombinant_1,recombinant_2])
	return outpop




def make_gametes(pop, grid_xlen, grid_ylen, sexdet_site):
	
	pop_gametes = [] # is a three-dimensional array, pop_gametes[x][y][[list_of_M_gametes],[list_of_F_gametes]]
	for x in range(grid_xlen):
		pop_gametes.append([ [[],[]] for y in range(grid_ylen) ])
	
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for i in pop[x][y]:
				if set([chrom[sexdet_site] for chrom in i[1:]]) == set(["A","T"]):
#					print "is male"
					pop_gametes[x][y][0].append(i[1])
					pop_gametes[x][y][0].append(i[2])
				elif set([chrom[sexdet_site] for chrom in i[1:]]) == set(["T","T"]):
#					print "is female"
					pop_gametes[x][y][1].append(i[1])
					pop_gametes[x][y][1].append(i[2])
				else:
					print "is neuter, no gametes"
					continue
	return pop_gametes
	

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
def do_meiosis(pop, grid_xlen, grid_ylen, n_gametes_per_individual, nominal_overall_rec_rate_per_site, sliding_window_size):
	"""
	for each diploid, do meiosis and return n_gametes_per_individual
	- local_rec_rate = (events_per_meiosis/sites) * sequence_identity_in_window
	"""
	wlenfloat = float(sliding_window_size)
	n_meioses = int(n_gametes_per_individual*0.5)
	
	outpop = []
	pop_gametes = [] 
	# pop_gametes is a four-dimensional array, pop_gametes[x][y][[list_of_M_gametes],[list_of_F_gametes]]
	# each gamete is a list with the sequence, and second element the position of the sexdet site. list_of_M_gametes = [[seq,sexdet_site],...]
		
	for x in range(grid_xlen):
		pop_gametes.append([ [[],[]] for y in range(grid_ylen) ])
	
	mset = set(["A","T"])
	fset = set("T")
	
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			for i in pop[x][y]:
				# is male or female?
				is_male = False
				is_female = False
				sexdet_alleles = set([i[1][i[3]],i[2][i[4]]])
				if sexdet_alleles == mset:
					is_male = True
				elif sexdet_alleles == fset:
					is_female = True
				else:
					print "is neuter, no gametes"
					continue			
				# now do meioses
				for meiosis in range(n_meioses): 
					does_recombine = numpy.random.choice([True,False],p=[nominal_overall_rec_rate_per_site*len(i[1]),1-nominal_overall_rec_rate_per_site*len(i[1])])
					if not does_recombine: # gametes have same haplotypes, sexdet_pos as parent
						if is_male:
							pop_gametes[x][y][0].append([i[1],i[3]])
							pop_gametes[x][y][0].append([i[2],i[4]])
						elif is_female:
							pop_gametes[x][y][1].append([i[1],i[3]])
							pop_gametes[x][y][1].append([i[2],i[4]])
					else: # find site where recombination occurs				
						# per-site measure sequence identity
						pairwise_identity_per_site = []
						idx = 0
						while idx < len( i[1] ):
							if i[1][idx] == i[2][idx]:
								if i[1][idx] == "-":
									pairwise_identity_per_site.append(0)
								else:
									pairwise_identity_per_site.append(1)
							else:
								pairwise_identity_per_site.append(0)
							idx += 1
						# now window-slide; at the left and right ends (too short for window) 
						# 	copy value of first resp. last window
						left_idx = 0
						right_idx = left_idx+sliding_window_size
				
						pairwise_identity_per_site_window_avg = []
						while right_idx < len( i[1] ):
							window_idenity_avg = sum(pairwise_identity_per_site[left_idx:right_idx])/wlenfloat
							pairwise_identity_per_site_window_avg.append(window_idenity_avg)
							left_idx += 1
							right_idx += 1
						# add the ends:
						pairwise_identity_per_site_window_avg = [pairwise_identity_per_site_window_avg[0]]*int(0.5*sliding_window_size) + pairwise_identity_per_site_window_avg + [pairwise_identity_per_site_window_avg[-1]]*int(0.5*sliding_window_size)
		#				print len(pairwise_identity_per_site_window_avg), pairwise_identity_per_site_window_avg			
						recombination_probabilities_per_site = numpy.array(pairwise_identity_per_site_window_avg)/sum(pairwise_identity_per_site_window_avg)
						# now choose recombination sites:
						rec_site = numpy.random.choice(range( len( i[1] ) ), p = recombination_probabilities_per_site)
#						print "recombination at ", rec_site
						recombinant_1 = i[1][:rec_site] + i[2][rec_site:] 
						recombinant_2 = i[2][:rec_site] + i[1][rec_site:]
						
						# now determine if the sexdet_pos switched phase:
						if rec_site >= i[3]: # recombination site was downstream of sexdet site: fine as-is!
							gamete1 = [recombinant_1,i[3]]
						else: # recombination site was upstream of sexdet site: sexdet site has switched phase to that of other parental chrom.
							gamete1 = [recombinant_1,i[4]]
						if rec_site >= i[4]: # recombination site was downstream of sexdet site: fine as-is!
							gamete2 = [recombinant_2,i[4]]
						else: # recombination site was upstream of sexdet site: sexdet site has switched phase to that of other parental chrom.
							gamete2 = [recombinant_2,i[3]]						
						
						if is_male:
							pop_gametes[x][y][0].append(gamete1)
							pop_gametes[x][y][0].append(gamete2)
						elif is_female:
							pop_gametes[x][y][1].append(gamete1)
							pop_gametes[x][y][1].append(gamete2)
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


def window_LD (inpop, grid_xlen, grid_ylen, sample_size, window_length, step_size):

	window_LD_per_site = []
	
	# random sample
	all_indivs = []
	for x in range(grid_xlen):
		for y in range(grid_ylen):
			all_indivs += inpop[x][y]
	
	males = [i for i in all_indivs if i[1][i[3]] == "A" or i[2][i[4]] == "A"]
	females = [i for i in all_indivs if i[1][i[3]]+i[2][i[4]] == "TT"]
	msamples = numpy.random.choice(range(len(males)), size=sample_size*0.5)

	msamples = [[males[x][1],males[x][2]] for x in msamples]
	fsamples = numpy.random.choice(range(len(females)), size=sample_size*0.5)
	fsamples = [[females[x][1],females[x][2]] for x in fsamples]
	
	samples_flat = [y for x in msamples+fsamples for y in x]

	n_chroms_sampled = float(sample_size*2)
	
	# find bi-allelic SNPs and calculate their freqs
	SNP_pos = []
	SNP_freqs = []
	for idx in range(len(samples_flat[0])):
		allconcat = "".join([x[idx] for x in samples_flat])
		alleles = sorted(list(set(allconcat)))
		if not "-" in alleles:
			if len(alleles) == 2:
				SNP_pos.append( idx )
				f = allconcat.count(alleles[0]) / n_chroms_sampled 
				if f > 0.5:
					maf = 1.0-f
				else:
					maf = f
				SNP_freqs.append(maf)
	
	# pairs of SNPs in sliding windows
	left_idx = 0
	right_idx = left_idx+window_length
	while right_idx < len( samples_flat[0]) :
		SNPs_in_window = [x for x in SNP_pos if x >= left_idx and x < right_idx]
#		print SNPs_in_window
		left_idx += step_size
		right_idx += step_size
		# get LD from max 50 randomly selected pairs of SNPs per window
		possible_pairs = set()
		for a in SNPs_in_window:
			for b in SNPs_in_window:
				if a != b:	
					pair = ".".join(sorted([str(a),str(b)]))
					possible_pairs.add(pair)
		try:
			selected_pairs = numpy.random.choice(list(possible_pairs), 50, replace = False)
		except ValueError:
			selected_pairs = list(possible_pairs)
#		print selected_pairs	
			
		LD_list = []
		for p in selected_pairs:
			[a,b] = [int(x) for x in p.split(".")]
			allconcat_a = "".join([x[a] for x in samples_flat])
			alleles_a = sorted(list(set(allconcat_a)))

			allconcat_b = "".join([x[b] for x in samples_flat])
			alleles_b = sorted(list(set(allconcat_b)))
					
			a_allele = alleles_a[0]
			b_allele = alleles_b[0]
			pa = allconcat_a.count(a_allele)/float(len(allconcat_a))
			pb = allconcat_b.count(b_allele)/float(len(allconcat_b))
			if not pa <= 0.5: # ensure that these are minor allele frequencies
				pa = 1.0 - pa
				a_allele = alleles_a[1]
			if not pb <= 0.5:
				pb = 1.0 - pb
				b_allele = alleles_b[1]
#			print pa, pb			
			pab = len([x for x in samples_flat if x[a] == a_allele and x[b] == b_allele])/float(len(allconcat_a))
			r2 = ((pab - (pa*pb))**2)/float( pa*(1-pa)*pb*(1-pb) )
			# r2_max, conditional on where in allele frequency space we are:
			if pa >= pb: # section S6
				r2_max = ((1-pa)*pb)/(pa*(1-pb))
			else: # section S7
				r2_max = (pa*(1-pb))/((1-pa)*pb)
						
	#		print pa, pb, pab, r2
			LD_list.append( r2/r2_max )
	#		print r2/r2_max, pa, pb, pab, r2, r2_max
		window_LD_per_site.append(numpy.mean( LD_list))	
	print window_LD_per_site[:10], window_LD_per_site[20:25]
			
			
#						pairwise_identity_per_site_window_avg = []
#						while right_idx < len( i[1] ):
#							window_idenity_avg = sum(pairwise_identity_per_site[left_idx:right_idx])/wlenfloat
#							pairwise_identity_per_site_window_avg.append(window_idenity_avg)
#							left_idx += 1
#							center_idx += 1
#							right_idx += 1

	
	
	
				
#				pp = genotypes.count(alleles[0]+alleles[0])
#				pq = genotypes.count(alleles[0]+alleles[1])
#				qq =
				
	ssites = len(SNP_pos)
	pi = numpy.mean(numpy.array(SNP_freqs)*(1.0-numpy.array(SNP_freqs)))/len(samples_flat[0])
	print ssites, pi			


def mutate_gametes_SNPs (ingametes, grid_xlen, grid_ylen, mu):
		
	# now can do multiple mutations per each chromosome!
	outgametes = []
	for x in range(grid_xlen):
		outgametes.append([])
		for y in range(grid_ylen):
			outgametes[x].append([[],[]])
			for s in [0,1]:
				for gamete in ingametes[x][y][s]:
					chrom = gamete[0]
					gaps = chrom.count("-")
					total_mutations = mu*(len(chrom)-gaps)
					if total_mutations <= 1.0:	# at most one mutation expected along entire chromosome				
						mutate = numpy.random.choice([True,False], p = [total_mutations,1.0-total_mutations])
						if mutate:
							char = "-" 
							while char == "-":
								mutated_site = numpy.random.choice(range(len(chrom)))
								char = chrom[mutated_site]
							new_allele = numpy.random.choice([b for b in ["A","T","G","C"] if not b == char])
							lchrom = list(chrom)
							lchrom[mutated_site] = new_allele
							outgametes[x][y][s].append(["".join(lchrom),gamete[1]]) # preserve also sexdet position
						else:
							outgametes[x][y][s].append(gamete)
					else: # more than 1 mutation expected along entire chromosome: iterate!
						print "more than one mutation" 
						focal_chrom = chrom
						remaining_mutations = total_mutations
						while remaining_mutations > 1.0:
							print remaining_mutations
							remaining_mutations -= 1.0						
							char = "-" 
							while char == "-":
								mutated_site = numpy.random.choice(range(len(focal_chrom)))
								char = focal_chrom[mutated_site]
							new_allele = numpy.random.choice([b for b in ["A","T","G","C"] if not b == char])
							lchrom = list(focal_chrom)
							lchrom[mutated_site] = new_allele
							focal_chrom = "".join(lchrom)
						# finally, draw if mutate with remaining probability:
						print remaining_mutations
						mutate = numpy.random.choice([True,False], p = [remaining_mutations,1.0-remaining_mutations])
						if mutate:
							char = "-" 
							while char == "-":
								mutated_site = numpy.random.choice(range(len(focal_chrom)))
								char = focal_chrom[mutated_site]
							new_allele = numpy.random.choice([b for b in ["A","T","G","C"] if not b == char])
							lchrom = list(focal_chrom)
							lchrom[mutated_site] = new_allele
							outgametes[x][y][s].append(["".join(lchrom),gamete[1]]) # preserve also sexdet position
						else:
							outgametes[x][y][s].append([focal_chrom,gamete[1]])
	return outgametes


def mutate_gametes_indel (ingametes, grid_xlen, grid_ylen, mu, indel_lengths_scale):
	
	## because insertion mutations change the length of EVERY sequence in the aligned form, need a two-pass algorithm:
	# 1. find mutating sites + indel lengths for all sequences
	# 2. apply indels, thereby resolving total & same length for all sequences
	
	## DIFFICULT!
		
	# now can do multiple mutations per each chromosome!
	outgametes = []
	for x in range(grid_xlen):
		outgametes.append([])
		for y in range(grid_ylen):
			outgametes[x].append([[],[]])
			for s in [0,1]:
				for gamete in ingametes[x][y][s]:
					chrom = gamete[0]
					gaps = chrom.count("-")
					total_mutations = mu*(len(chrom)-gaps)
					if total_mutations <= 1.0:	# at most one mutation expected along entire chromosome				
						mutate = numpy.random.choice([True,False], p = [total_mutations,1.0-total_mutations])
						if mutate:
							char = "-" 
							while char == "-":
								mutated_site = numpy.random.choice(range(len(chrom)))
								char = chrom[mutated_site]
							new_allele = numpy.random.choice([b for b in ["A","T","G","C"] if not b == char])
							lchrom = list(chrom)
							lchrom[mutated_site] = new_allele
							outgametes[x][y][s].append(["".join(lchrom),gamete[1]]) # preserve also sexdet position
						else:
							outgametes[x][y][s].append(gamete)
					else: # more than 1 mutation expected along entire chromosome: iterate!
						print "more than one mutation" 
						focal_chrom = chrom
						remaining_mutations = total_mutations
						while remaining_mutations > 1.0:
							print remaining_mutations
							remaining_mutations -= 1.0						
							char = "-" 
							while char == "-":
								mutated_site = numpy.random.choice(range(len(focal_chrom)))
								char = focal_chrom[mutated_site]
							new_allele = numpy.random.choice([b for b in ["A","T","G","C"] if not b == char])
							lchrom = list(focal_chrom)
							lchrom[mutated_site] = new_allele
							focal_chrom = "".join(lchrom)
						# finally, draw if mutate with remaining probability:
						print remaining_mutations
						mutate = numpy.random.choice([True,False], p = [remaining_mutations,1.0-remaining_mutations])
						if mutate:
							char = "-" 
							while char == "-":
								mutated_site = numpy.random.choice(range(len(focal_chrom)))
								char = focal_chrom[mutated_site]
							new_allele = numpy.random.choice([b for b in ["A","T","G","C"] if not b == char])
							lchrom = list(focal_chrom)
							lchrom[mutated_site] = new_allele
							outgametes[x][y][s].append(["".join(lchrom),gamete[1]]) # preserve also sexdet position
						else:
							outgametes[x][y][s].append([focal_chrom,gamete[1]])
	return outgametes

####### DELETE THIS
def replace_SNPs_with_indels (inseqs, indel_rate, indel_lengths_scale):
	
	# find SNPs that will become indels
	indel_pos = []
	SNP_alleles_to_become_indel = []
	for idx in range(len(inseqs.values()[0])):
		col = [x[idx] for x in inseqs.values()]
		if len(set(col)) > 1:
			is_indel = numpy.random.choice([True,False],p=[indel_rate,1.0-indel_rate])
			if is_indel:
				indel_pos.append(idx)
				SNP_alleles_to_become_indel.append( numpy.random.choice(list(set(col)))) 
	
	# choose indel lengths
	n_indels = len(indel_pos)
	indel_lengths = numpy.round(numpy.random.exponential(scale=indel_lengths_scale, size=n_indels)+1)
	
	# delete at indel positions
	out_dict = {}
	for id,seq in inseqs.items():
		print id
		outseq = seq
		for i in range(len(indel_pos)):
			idx = indel_pos[i]
			allele = SNP_alleles_to_become_indel[i]
			indel_length = int(indel_lengths[i])
			if seq[idx] == allele:
				outseq = outseq[:idx] + "-"*indel_length + outseq[idx+indel_length:]
		out_dict[id] = outseq
	return out_dict
###### DEL HERE
	
######################

########

inseqs = read_phylip("start.aln")

with_indels = replace_SNPs_with_indels (inseqs, 0.15, 1)

#write_phylip (with_indels, "start_with_indels.aln")

# OK works! next: inversions and/or larger insertions (TEs!)
# see "Inversions are bigger on the X chromosome" Cheng & Kirkpatrick


gridlength_x = 1
gridlength_y = 1
N_per_cell = 40
exponential_decay_shape_param = 2.0


sexdet_site = 1000 ## males = AT, females = TT ; other genotypes are neuters = no gametes! 

# 2 events on a 10 Mb chromosome = 2e-7 per site
nominal_overall_rec_rate_per_site = 1e-8 # can be at most 1 event; => max(nominal_overall_rec_rate_per_site) = 1/sequence_length
sliding_window_size = 100

# SNP mutation rate
mu_snps = 2.5e-6

# indel mutation rate
mu_indel = 0.15*mu_snps

# TE indel rates: see here:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5447328/
TE_insertion_rate = "2.11*10^-9"
TE_deletion_rate = "1.37*10^-10"
TE_length = 2000
dna = ["A","G","C","T"]
TE_seq = ''.join(random.choice(dna) for i in range(TE_length))
#print TE_seq



#pop = make_dioecious_pop_imported_startseqs (N_per_cell, gridlength_x , gridlength_y, with_indels.values(), sexdet_site) # OK

pop = make_dioecious_pop_single_random_startseq (N_per_cell, gridlength_x , gridlength_y, 10000, sexdet_site)


for generation in range(20):
#	measure_chrom_lengths (pop, gridlength_x, gridlength_y)
	pop = enforce_carrying_capacity (pop, gridlength_x, gridlength_y, 70.0) # carrying capacity must be a float!
	if (float(generation)/100).is_integer():
		window_LD (pop, gridlength_x, gridlength_y, 20, 1000, 200)
		pop_to_phylip (pop, gridlength_x, gridlength_y, 10, 10, "start_XY_ident.Ne40.generation_"+ str(generation) +".aln")
	n_total = measure_grid_density (pop, gridlength_x, gridlength_y)
	print generation, n_total
	if n_total == 0.0:
		print "extinction at generation	", generation
		break
	pop_gametes = do_meiosis(pop, gridlength_x, gridlength_y, 4, nominal_overall_rec_rate_per_site, sliding_window_size)
	dispersed_pop_gametes = disperse_male_gametes (pop_gametes, gridlength_x, gridlength_y, exponential_decay_shape_param)
	dispersed_pop_gametes = mutate_gametes_SNPs(dispersed_pop_gametes, gridlength_x, gridlength_y,mu_snps )
#	print dispersed_pop_gametes
	pop = mate (dispersed_pop_gametes, gridlength_x, gridlength_y)
	pop = disperse_pop (pop, gridlength_x, gridlength_y, exponential_decay_shape_param)


pop_to_phylip (pop, gridlength_x, gridlength_y, 10, 10, "start_XY_ident.Ne40.generation_"+ str(generation) +".aln")




	
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






