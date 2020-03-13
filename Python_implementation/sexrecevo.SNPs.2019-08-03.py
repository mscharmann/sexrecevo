
"""

to avoid slow-down of the simulations by large numbers of mutations which have no further effect:
once a region becomes fully non-recombining, no further mutations are allowed at its core

"""

import numpy, random, os, argparse, datetime


# parses command line arguments
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("-m", required=True, help="mutation rate", metavar="FLOAT")
	parser.add_argument("-r", required=True, help="nominal recombination rate", metavar="FLOAT")
	parser.add_argument("-rws", required=True, help="recombination window size, an even integer smaller than 0.5 * s to get at least two windows", metavar="INT")
	parser.add_argument("-min_ident", required=True, help="minimum sequence identity in a recombination window, proportion <= 1.0 below which recombination probability is zero", metavar="FLOAT")
	parser.add_argument("-N", required=True, help="population size in number of diploids", metavar="INT")
	parser.add_argument("-g", required=True, help="number of generations", metavar="INT")
	parser.add_argument("-s", required=True, help="number of sites / length of chromosomes", metavar="INT")
	parser.add_argument("-logf", required=True, help="report statistics every logf generations", metavar="INT")
	parser.add_argument("-logw", required=True, help="report statistics in windows of size w bases, sliding with a step-size of 0.5*w", metavar="INT")
	
	args = parser.parse_args()

	
	return args



def make_dioecious_pop (N, sexdet_site):
	
	"""
	pop is a list of two lists: [[males],[females]]
	
	each diploid individual has the following structure:
	[chrom1, chrom2]
	each chrom is a list of alleles: [0,0,0,1,0,0,0,...]
	
	the position of these alleles along the chromosome is given by another vector, which contains positions of all variants in the population:
	variant_idxes = [11,123,8675,98761,100001,...]
	invariant sites are NOT tracked!
		
	- Y-chromosomes have allele "1" at sexdet_site, X-chromosomes have alelle "0"
	"""
	variant_idxes = [sexdet_site]
	pop = [[],[]]
	cnt = 0
	for i in range(int(0.5*N)):
		cnt += 1
		idv = [[1], [0]]
		pop[0].append( idv )
	for i in range(int(0.5*N)):
		cnt += 1
		idv = [[0], [0]]
		pop[1].append( idv )	
	return pop, variant_idxes


def make_sanitycheck_pop (N, sexdet_site):
	"""
	fixed nonrec Y-bloc	with 20% sequence divergence (4 mutations in 20 bases), a total of at least 200 bases long
	=> suitable to stop recombination in at least 1 window of up to 100 bases length (parameter rws), and if identity to stop recombination between 100 and 80%
	=> 
	"""
	variant_idxes = [sexdet_site]
	cnt = 1
	for i in range(26):
		variant_idxes.append(sexdet_site-cnt)
		cnt += 5
	cnt = 1
	for i in range(26):
		variant_idxes.append(sexdet_site+cnt)
		cnt += 5

	
	pop = []
	for s in [0,1]:
		pop.append([])
	for i in range(int(N/2.0)):
		idv = [[1]*53, [0]*53]
		pop[0].append( idv )
	for i in range(int(N/2.0)):
		idv = [[0]*53, [0]*53]
		pop[1].append( idv )
	return pop, variant_idxes

	

def random_choice_with_exceptions (sample_from, sample_n, exclusion_set):
	
	sampled_items = []
#	cnt = 0
	while len(sampled_items) < sample_n:
#		cnt += 1
		rc = numpy.random.choice(sample_from, sample_n)
		sampled_items = [i for i in rc if not i in exclusion_set]
	
#	print "number of calls:	", cnt
	return numpy.array(sampled_items)

def mutate_gametes(ingametes, mu, chrom_length, variant_idxes, banned_start, banned_end):
		
	idxes_of_chrom = range( chrom_length )
	
	n_gametes = len(ingametes[0])
	n_gametes += len(ingametes[1])
	
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
		
		# purge mutations that fell into the banned region, don't allow them to occur.
		if banned_start >= 0:
			if banned_start < banned_end:
#				print "before", mutating_sites
				mutating_sites = [site for site in mutating_sites if not (site > banned_start and site < banned_end)  ]
#				print "after", mutating_sites
				if len(mutating_sites) == 0:
					return ingametes, variant_idxes
				
		# draw sex where new alleles arise
		random_sex = numpy.random.choice([0,1], total_expected_mutations)
	
		old_gamete_length = len(variant_idxes)
	
		addsites = [0]*total_expected_mutations
		
		# create new outgametes
		outgametes = []
		for s in [0,1]:
			outgametes.append([])
			
			for ingamete in ingametes[s]:
				seq = ingamete[:] ## have to COPY the list, otherwise will change the list in-place and cause havoc!!
				seq += addsites
				outgametes[s].append(seq)

		# now throw the mutations onto the gametes!
		cnt = 0
		for mut in sorted(mutating_sites):
	
			idx_of_gamete_in_sexpool = random.choice(range(len( outgametes[random_sex[cnt]] )))
			outgametes[random_sex[cnt]][idx_of_gamete_in_sexpool][old_gamete_length+cnt] = 1
			cnt += 1

		# update variant_idxes, 
		variant_idxes += sorted(mutating_sites)
	#	print variant_idxes
	
		# sort everything
		# from https://stackoverflow.com/questions/7851077/how-to-return-index-of-a-sorted-list
		indexes_after_sorting = sorted(range(len(variant_idxes)), key=variant_idxes.__getitem__)

		variant_idxes = [variant_idxes[i] for i in indexes_after_sorting[:]]
	
		outgametes_newly_sorted = []
		for s in [0,1]:
			outgametes_newly_sorted.append([])
				
			for ingamete in outgametes[s]:
				sortedseq = [ingamete[i] for i in indexes_after_sorting]
				outgametes_newly_sorted[s].append(sortedseq)
	
	#	print outgametes[0][0][0][0]	
	#	print outgametes_newly_sorted[0][0][0][0]	
		return outgametes_newly_sorted, variant_idxes
	else:
		return ingametes, variant_idxes	


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
		

def report_window_male_female_Fst_and_dxy (inpop, variant_idxes, sample_size, generation, window_size, step_size, chrom_length):

	# random sample from population
	all_indivs = inpop[0] + inpop[1]
	
	males = inpop[0]
	females = inpop[1]
	
	msamples = numpy.random.choice(range(len(males)), size=int(sample_size*0.5))
	msamples = [[males[x][0],males[x][1]] for x in msamples]
	
	fsamples = numpy.random.choice(range(len(females)), size=int(sample_size*0.5))
	fsamples = [[females[x][0],females[x][1]] for x in fsamples]
	
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
	start = -1
	end = -1
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
				break
			else:
				start = -1
				end = -1
							
	if os.path.isfile("logfile_nonrec_region_size.txt"): 
		with open("logfile_nonrec_region_size.txt", "a") as F:
			F.write("\t".join([str(generation), str(nonrec_length)])+ "\n")
	else:
		with open("logfile_nonrec_region_size.txt", "w") as A:
			A.write("\t".join(["generation","size_of_nonrec_region"]) + "\n")
			A.write("\t".join([str(generation), str(nonrec_length)])+ "\n")
	return nonrec_length, start, end
			

def sample_parents (N, pop):
	
	fathers_idxes = list(numpy.random.choice(range(len(pop[0])),N))
	mothers_idxes = list(numpy.random.choice(range(len(pop[1])),N))
	fathers_unique = list(set(fathers_idxes))
	mothers_unique = list(set(mothers_idxes))
	fathers_offspring_counts = [fathers_idxes.count(x) for x in fathers_unique]
	mothers_offspring_counts = [mothers_idxes.count(x) for x in mothers_unique]	
	return [fathers_idxes, mothers_idxes, fathers_unique, mothers_unique, fathers_offspring_counts, mothers_offspring_counts]

def make_gametes_meioses (pop, N, sampled_parents, variant_idxes, chrom_length, nominal_overall_rec_rate_per_site, sexdet_site, rec_window_size, min_ident, record_rec_probs):

	if record_rec_probs:
		X_Y_recomb_probs = numpy.zeros(chrom_length)
		X_X_recomb_probs = numpy.zeros(chrom_length)
	else:
		X_Y_recomb_probs = ""
		X_X_recomb_probs = ""
		
		
	male_gametes = []
	n_meioses_M_total = 0
	for idx,father in enumerate(sampled_parents[2]):
		father_genotype = pop[0][father]
		n_meioses = sampled_parents[4][idx]/2.0
		if not n_meioses.is_integer():
			n_meioses += 0.5
		n_meioses = int(n_meioses)
		n_meioses_M_total += n_meioses
		for meiosis in range(n_meioses):
			g1, g2, rec_probs = simple_meiosis (father_genotype, variant_idxes, chrom_length,  nominal_overall_rec_rate_per_site, sexdet_site, rec_window_size, min_ident, record_rec_probs)
			male_gametes += [g1,g2]
			X_Y_recomb_probs += rec_probs
			
	female_gametes = []
	n_meioses_F_total = 0
	for idx,mother in enumerate(sampled_parents[3]):
		mother_genotype = pop[1][mother]
		n_meioses = sampled_parents[5][idx]/2.0
		if not n_meioses.is_integer():
			n_meioses += 0.5
		n_meioses = int(n_meioses)
		n_meioses_F_total += n_meioses
		for meiosis in range(n_meioses):
			g1, g2, rec_probs = simple_meiosis (mother_genotype, variant_idxes, chrom_length,  nominal_overall_rec_rate_per_site, sexdet_site, rec_window_size, min_ident, record_rec_probs)
			female_gametes += [g1,g2]
			X_X_recomb_probs += rec_probs	
			
	if record_rec_probs: 
		avg_X_Y_rec_probs = X_Y_recomb_probs / n_meioses_M_total 
		avg_X_X_rec_probs = X_X_recomb_probs / n_meioses_F_total
	else:
		avg_X_Y_rec_probs = "" 
		avg_X_X_rec_probs = ""
		
	# random shuffle the gamete pools
	# If gamete pools not randomised, divergence along chromosome can become skewed, because predominantly brother-sister mating!
	random.shuffle(male_gametes)
	
	return male_gametes, female_gametes, avg_X_Y_rec_probs, avg_X_X_rec_probs
	
def simple_meiosis (parent_genotype, variant_idxes, chrom_length, nominal_overall_rec_rate_per_site, sexdet_site, rec_window_size, min_ident, record_rec_probs):
	
	effect = True # if mutations have no effect on rec probs; sanity-checking; place random 
		
	# measure sequence identity 
	idents = [0 if parent_genotype[0][idx] != parent_genotype[1][idx] else 1 for idx in range(len(variant_idxes)) ]
	if not effect:	
		idents = [1]*len(variant_idxes) # always fully identical... if mutations have no effect on rec probs; sanity-checking; place random 
	# set identity at sexdet_site to 1, even if heterozyogous (XY male meiosis)
	idents[variant_idxes.index(sexdet_site)] = 1
#	print idents[variant_idxes.index(chrom_length/2)]
	if effect:
		if sum(idents) >= len(idents)-1:
			# identical chromosomes or only one SNP, don't try to place breakpoints
			if record_rec_probs:
				return parent_genotype[0], parent_genotype[1], numpy.array([nominal_overall_rec_rate_per_site]*chrom_length)
			else:
				return parent_genotype[0], parent_genotype[1], ""
				
	# dirty: now translate to non-overlapping windows, each meiosis with a random starting index!
	window_size = rec_window_size # an even integer smaller than 0.5 * chrom_length to get at least two windows
	wsize_float = float(window_size)
	
	remain = int((chrom_length/wsize_float % 1)*wsize_float)
	n_windows = (chrom_length/window_size)
	if remain == 0:
		remain = window_size
		n_windows = (chrom_length/window_size)-1
#	print "remaining", remain
#	window_start_idx = 1
	window_start_idx = int(numpy.random.choice(range(1,remain),1))
	windows_end_idx = window_start_idx + window_size*n_windows
#	print window_start_idx, windows_end_idx
#	print n_windows
	
	window_idents = [window_size]*n_windows 
	# print n_windows
	idx = -1
	for variant in variant_idxes:
		idx += 1
		if idents[idx] == 0:
			# only consider variants inside the windows, exclude leading/trailing edge variants
			if variant > window_start_idx:
				if variant < windows_end_idx:
					# substract window_start_idx from variant
					# divide by window size
					# results in index for windows_right_borders_nozero AND window_idents
					variant_window_idx = (variant-window_start_idx)/window_size
					
#					print "start, var and right border", window_start_idx, variant, windows_right_borders_nozero[variant_window_idx]
#					print windows_right_borders_nozero
#					print "before", window_idents
					window_idents[variant_window_idx] -= 1
#					print "after", window_idents	
		
	# set hard limit "min_ident" of sequence divergence beyond which rec probability drops to zero:
	window_idents = [0 if g/wsize_float < min_ident else g for g in window_idents]		
#	print len(window_idents)
			
	# now seq identity to probabilities; first and last windows (= outer edges) get additional sizes, which were not covered by seq-identity measurement, extrapolated from nearest inner windows
		
	rec_probs = [(window_start_idx+window_size)*(window_idents[0]/wsize_float)*nominal_overall_rec_rate_per_site] + [g*nominal_overall_rec_rate_per_site for g in window_idents[1:-1]]+ [(window_size + (chrom_length-windows_end_idx))*(window_idents[-1]/wsize_float)*nominal_overall_rec_rate_per_site]
#	print len(rec_probs), sum(rec_probs)
	if record_rec_probs:
#		print "recording"
		# expand windows to full per-base sequence of probabilities!
		per_base_rec_probs = [(window_idents[0]/wsize_float)*nominal_overall_rec_rate_per_site]*(window_start_idx+window_size)
		for w in rec_probs[1:-1]:
#			print w
			per_base_rec_probs += [w/wsize_float]*window_size
		per_base_rec_probs += [rec_probs[-1]/(window_size + (chrom_length-windows_end_idx))]*(window_size + (chrom_length-windows_end_idx))
#		print len(per_base_rec_probs)
		per_base_rec_probs_arr = numpy.array(per_base_rec_probs)
	rec_probs.append( 1.0-sum(rec_probs)) # a dummy window is appended which represents the not-simulated remainder of the chromosome!
		
	# now choose rec. window
#	rec_window = numpy.random.choice(range(len(rec_probs)), p = rec_probs)
	rec_window = numpy.random.choice(range(len(rec_probs)), p = rec_probs)
#	rec_window = multidimensional_shifting(1, 1, range(len(rec_probs)), rec_probs)			
	
#	print "v_inW"
#	print window_start_idx, variant_idxes_window_space			
#	print variant_idxes
	if rec_window == len(rec_probs)-1:
		# dummy window chosen => recombination falls outside of the simulated section of chromosome
		# print "NA"
		if record_rec_probs:
			return parent_genotype[0], parent_genotype[1], per_base_rec_probs_arr
		else:
			return parent_genotype[0], parent_genotype[1], ""
	else:
#		print rec_window, variant_idxes_window_space
#		rec_idx_among_variants = max([t for t in variant_idxes_window_space if t < rec_window])
#		print "LA"
#		print window_idents, rec_window
#		print variant_idxes
#		print rec_window, rec_idx_among_variants, windows_right_borders[rec_window]
		# build the recombinant gametes			
#		print parent_genotype[0][:rec_idx_among_variants+1], parent_genotype[1][rec_idx_among_variants+1:]
#		print parent_genotype[1][:rec_idx_among_variants+1], parent_genotype[0][rec_idx_among_variants+1:]
		# back-translate which variants the windows encompass
		windows_right_borders = [0] + [window_size*n+window_start_idx for n in range(1,n_windows+1)]
#		windows_right_borders.append(chrom_length)

		right_border_of_rec_window = windows_right_borders[rec_window]
		try:
			rec_idx_among_variants = max([cidx for cidx,c in enumerate(variant_idxes) if c < right_border_of_rec_window])
		except ValueError:
			rec_idx_among_variants = 0
		
		# randomly choose if breakpoint placed before or after the chosen window:
#		before_or_after = numpy.random.randint(2)
#		if before_or_after == 1:
		# breakpoint AFTER the chosen window
		gamete_1 = parent_genotype[0][:rec_idx_among_variants+1] + parent_genotype[1][rec_idx_among_variants+1:]
		gamete_2 = parent_genotype[1][:rec_idx_among_variants+1] + parent_genotype[0][rec_idx_among_variants+1:]
#		else:
#			# breakpoint BEFORE the chosen window
#			gamete_1 = parent_genotype[0][:rec_idx_among_variants] + parent_genotype[1][rec_idx_among_variants:]
#			gamete_2 = parent_genotype[1][:rec_idx_among_variants] + parent_genotype[0][rec_idx_among_variants:]
		
		
					
		if record_rec_probs:
			return gamete_1, gamete_2, per_base_rec_probs_arr
		else:
			return gamete_1, gamete_2, ""

def make_offspring_pop (N, male_gametes, female_gametes, variant_idxes, sexdet_site):
	
	idx_of_sex =  variant_idxes.index(sexdet_site)
	
	new_pop = [[],[]]
	for i in range(N):
		if male_gametes[i][idx_of_sex] == 1:
			new_pop[0].append( [male_gametes[i], female_gametes[i]] )
		else:
			new_pop[1].append( [male_gametes[i], female_gametes[i]] )

	return new_pop

def clean_from_invariant_sites (pop, variant_idxes):
	
	## BUT I should not remove recombination-effectors if they go to fixation!!? How to deal with this shit?? 
	## => Modify the nominal-rec rate??
	## => have to track an array with the nominal rec-rate, because it can change over time and change specifically for each site.
	## BUT if only heterozygous state has an effect, IT IS APPROPRIATE to remove the effects of alleles that become fixed! :)
	
	# find invariant sites by checking all chromosomes in pop:
	all_chroms = []
	for s in [0,1]:
		for idv in pop[s]:
			all_chroms.append( idv[0] )
			all_chroms.append( idv[1] )
	
	total_chrom_cnt = len(all_chroms)
	
	invariant_idxes = []
	for idx in range(len(all_chroms[0])):
		allsum = sum([chr[idx] for chr in all_chroms])
		if allsum == 0:
			invariant_idxes.append(idx)
		elif allsum == 2*total_chrom_cnt:
			invariant_idxes.append(idx)
		else:
			None
	
	variant_idxes = del_list_numpy(variant_idxes, invariant_idxes)
	
	cleaned_pop = []
	for s in [0,1]:
		cleaned_pop.append([])
		for idv in pop[s]:
			chr1 = del_list_numpy(idv[0], invariant_idxes)
			chr2 = del_list_numpy(idv[1], invariant_idxes)
			cleaned_pop[s].append([chr1,chr2])
	 	
	return cleaned_pop, variant_idxes


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



def report_duration(start_time):
	
	end_time = datetime.datetime.now()
	duration = end_time - start_time                         # For build-in functions
	duration_in_s = duration.total_seconds()

	days    = divmod(duration_in_s, 86400)        # Get days (without [0]!)
	hours   = divmod(days[1], 3600)               # Use remainder of days to calc hours
	minutes = divmod(hours[1], 60)                # Use remainder of hours to calc minutes
	seconds = divmod(minutes[1], 1)               # Use remainder of minutes to calc seconds

	with open("logfile_duration.txt", "w") as F:
		F.write("start: " + str(start_time) + "\n")
		F.write("end: " + str(end_time)+ "\n")
		F.write( "total duration: %d days, %d hours, %d minutes and %d seconds" % (days[0], hours[0], minutes[0], seconds[0]) + "\n")
	

###################### MAIN

start_time = datetime.datetime.now()

args = get_commandline_arguments ()


sexdet_site = int(0.5*float(args.s)) ## always placed at mid of the chromosome

# 2 events on a 10 Mb chromosome = 2e-7 per site
nominal_overall_rec_rate_per_site = float(args.r)

mu = float(args.m)


# length of chromosome
chrom_length = int(args.s)

logf = int(args.logf)
logw = int(args.logw)

rec_window_size = int(args.rws)

N = int(args.N)

min_ident = float(args.min_ident)

pop, variant_idxes = make_dioecious_pop (N, sexdet_site)



#pop, variant_idxes = make_sanitycheck_pop (N, sexdet_site)
nominal_rec_rate_array = numpy.array([nominal_overall_rec_rate_per_site]*chrom_length)

for generation in range(int(args.g)+1):
	if (float(generation)/10).is_integer(): # every x generations clean from monomorphic sites
		pop, variant_idxes = clean_from_invariant_sites (pop, variant_idxes)
	sampled_parents = sample_parents (N, pop)
	if (float(generation)/logf).is_integer(): # every x generations report something
		male_gametes, female_gametes, avg_X_Y_rec_probs, avg_X_X_rec_probs = make_gametes_meioses (pop, N, sampled_parents, variant_idxes, chrom_length, nominal_overall_rec_rate_per_site, sexdet_site, rec_window_size, min_ident, True)
		report_window_male_female_Fst_and_dxy (pop, variant_idxes, 20.0, generation, logw, int(0.5*logw), chrom_length)
		report_window_rec_probs (generation, logw, int(0.5*logw), avg_X_Y_rec_probs, "XY")
		report_window_rec_probs (generation, logw, int(0.5*logw), avg_X_X_rec_probs, "XX")
		print generation, len(variant_idxes)
		nonrec_length, nonrec_start, nonrecend = report_size_of_nonrec_region (generation, avg_X_Y_rec_probs, sexdet_site )
#		print nonrec_length, nonrec_start, nonrecend
		# we stop the simulation if the nonrec_size reaches the full length of the chromosome, because there is no way back anyway:
		if nonrec_length == chrom_length:
			print "entire chrom fragment reached zero prob of recomb in male meiosis, exiting"
			with open("logfile_nonrec_region_size.txt", "a") as F:
				for gi in range(generation, int(args.g)+1, logf)[1:]:
					F.write("\t".join([str(gi), str(chrom_length)])+ "\n")
			report_duration(start_time)
			exit()

	else:
		male_gametes, female_gametes, avg_X_Y_rec_probs, avg_X_X_rec_probs = make_gametes_meioses (pop, N, sampled_parents, variant_idxes, chrom_length, nominal_overall_rec_rate_per_site, sexdet_site, rec_window_size, min_ident, False)
	mutated_gametes, variant_idxes = mutate_gametes([male_gametes, female_gametes], mu, chrom_length, variant_idxes, nonrec_start+rec_window_size, nonrecend-rec_window_size)
	pop = make_offspring_pop (N, mutated_gametes[0], mutated_gametes[1], variant_idxes, sexdet_site)
#	pop, variant_idxes, sexdet_site = invert_sequence (pop, variant_idxes, sexdet_site, chrom_length)

report_duration(start_time)
print "Done!"

