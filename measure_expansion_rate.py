#!/usr/local/bin/python
# Python 2.7.6
# 
# June 2019
# Mathias Scharmann



import argparse
import os

########################## HEAD

# parses arguments from command line
def get_commandline_arguments ():
	
	parser = argparse.ArgumentParser()
	
	parser.add_argument("--d", required=True)

 	args = parser.parse_args()
		
	# finish
	return args

################################## CORE

def read (INDIR, SD_site, windowsize, stepsize):

	infiles = [x for x in os.listdir(INDIR) if "logfile_XY_rec_probs.txt" in x]
	
	outlines = []
	with open(INDIR+"/"+infiles[0], "r") as F:
		F.readline()
		window_sites = [int(x) for x in F.readline().strip("\n").split("\t")[1:]]
		for line in F:
			fields = line.strip("\n").split("\t")
			generation = fields[0]
			outlines.append([generation])
		
	window_sites_with_SD_idxes = set([window_sites.index(x) for x in window_sites if x >= SD_site-windowsize and x <= SD_site+windowsize])
	
	header = ["generation"]
	cnt = 0
	for f in infiles:
		cnt += 1
		header.append("replicate_"+ str(cnt))
		with open(INDIR+"/"+f, "r") as F:
			F.readline()
			F.readline()
			linecnt = 0
			for line in F:
				fields = line.strip("\n").split("\t")
				generation = fields[0]
				consec_zero_blocs = find_consecutive_zeros (fields[1:])
				nonrec_length = 0
				if len(consec_zero_blocs[0]) == 0:
					SD_bloc = False
				else:
					# is SD-site among them?
					SD_bloc = False
					for bloc in consec_zero_blocs:
						for idx in bloc:
							if idx in window_sites_with_SD_idxes:
								SD_bloc = True
								nonrec_length = windowsize+(len(bloc)-1)*stepsize
				if SD_bloc == True:
					print f, generation, nonrec_length, SD_bloc
					outlines[linecnt].append(str(nonrec_length))
				else:
					print f, generation, nonrec_length
					outlines[linecnt].append(str(nonrec_length))
				linecnt += 1
	
	with open("nonrec_sizes_through_time.txt", "w") as F:
		F.write("\t".join( header ) + "\n")
		for l in outlines:
			F.write("\t".join( l ) + "\n")
	

def find_consecutive_zeros (sequence):
	
	zero_idxes = []
	for idx,x in enumerate(sequence):
		if x == "0.0":
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
				
################################## MAIN


	
args = 	get_commandline_arguments ()
	
#print args

read(args.d, 50000, 500, 100)
	
print "Done!"
	


