#!/usr/bin/python

# matrix_scorer.py
# Jeremy Horst, 04/06/2010

#####################################
# this program takes as input       #
# [1]- a FASTA file                 #
# [2]- a directory of FASTA files   #
# [3]- a scoring matrix             #
# [4,5]- gap penalties              #
# calculates total similarity score #
#####################################

# break up query protein into all possible fragments 
# that match the size of sequences in the database
# score position by TSS of all fragments covering the position
# normalize by amount of fragments considered at position,
# perhaps multiplied by the size of the fragment itself
# test by recapture


# details:
#	- all fastas must have equivalent name.fasta & first line >name, with at least white space after >name


import os
from sys import argv
from os import mkdir
from os import listdir
from random import shuffle
from random import randrange
from subprocess import call
from subprocess import Popen
from subprocess import PIPE



##################
# matrix details #
value_limit = 20
ggsearch_dir = './fasta-35.4.11/'
if '-gg' in argv:
	ggsearch_dir = argv[argv.index('-gg')+1]
#matrix_dir   = ggsearch_dir + 'data/'
matrix_dir   = './'
AA_order = "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X".split()
matrix_size = len(AA_order)
min_value = -1.4
max_value = 5.5
##################


#make process dir
try: mkdir('tmp')
except: nada=True
try: mkdir('tmp_mat')
except: nada=True


#############################
def run_ggsearch(gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library):
	command = "%s/bin/ggsearch35 -q -p -m 4 -f %s -g %s -s %s -O %s %s %s"\
	% (ggsearch_dir, gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library)
#	print command
	crap = Popen(command.split(),stderr=PIPE,stdout=PIPE).communicate()[0]
#############################
def TSS(fastasA, setA, fastasB, setB, gap_open, gap_extension, substitution_matrix):
	# TSS_A-B([A]_NA - [B]_NB) = 
	# = (1/ (NA*delta_AB)) sum{1toNA}[ sum{1toNB}[ PSS_ij(1 - delta_ij * delta_AB)]]
	
	######################
	# calculate first_term
	# Kronecker delta: check if sets are the same
	delta_AB = 0
	if setA == setB:  delta_AB = 1
	first_term = 1/float(len(setA)*(len(setB)-delta_AB))
	
	######################
	# for second term iterate through sets
	second_term = 0
	for protein_a in fastasA:
		seq_A = open(protein_a).readlines()[1].strip()
		seq_length = float(len(seq_A))
#@#		seq_length = 1
		################
		# run ggsearch35
		a_name = protein_a.split('/')[-1].split('.')[0]
		ggsearch_output = 'tmp/'+ a_name +'.ggsearch_output'
		run_ggsearch(gap_open, gap_extension, substitution_matrix, ggsearch_output, protein_a, setB)
		################
		# parse ggsearch output for each protein in library set
		protein_scores = {}
		start = False
		for line in open(ggsearch_output).readlines():
			# collect protein scores
			if start:
				if not line.strip():  start = False
				#sorry, this next line is a parsing hack, for ggsearch output
				elif ')' in line and not ':' in line:
					prot_name = line.split()[0]
					score     = int(line.split(')')[1].split()[0])
					if not protein_scores.has_key(prot_name):
						protein_scores[prot_name] = score
			# initiate collection
			if 'The best scores are:' in line:  start=True
		######################
		# score protein a vs setB
		for protein_b in fastasB:
			# Kronecker delta: check if proteins are the same
			seq_B = open(protein_b).readlines()[1].strip()
			delta_ab = 0
			if seq_A == seq_B:  delta_ab = 1
			elif THIS_protein in protein_b:  delta_ab = 1 #; print "caught same same!"
			
			# strip '.fasta' in name for indexing
			name_b = protein_b.split('/')[-1].split('.')[0]
			# if not present, assume PSS_ab = 0
			PairSimScore_ab = 0
			if protein_scores.has_key(name_b):
				PairSimScore_ab = protein_scores[name_b]
			######################
			
			###########################
			# normalize to seq length #
			# sum to second term #
			second_term += PairSimScore_ab * (1 - delta_ab) / seq_length
			######################
#	print 'first_term, second_term:',first_term, second_term
	return first_term * second_term
############################


###########################################################
def directory_to_library(fasta_set,library):
	writer = open(library,'w')
	for f in fasta_set:
		for line in open(f).readlines():
			writer.write(line.strip()+'\n')
	writer.close()
#############################
def dir_to_fastalist_N_lib(db_dir,db_lib):
	db_fastas = []
	for f in listdir(db_dir):
		if f.endswith('.fasta'): db_fastas += [db_dir+f]
	# make searchable library from database
	directory_to_library(db_fastas,db_lib)
	return db_fastas
#############################
#############################
#############################
def db_seq_lengths(db_dir):
	lengths = []
	for f in listdir(db_dir):
		if f.endswith('.fasta'):
			length = len(''.join(open(db_dir+f).readlines()[1:]).strip().replace('\n','').replace('-',''))
			if not length in lengths: lengths+=[length]
	return lengths
#############################
def db_seq_lengths_all(db_dir):
	lengths = []
	for f in listdir(db_dir):
		if f.endswith('.fasta'):
			length = len(''.join(open(db_dir+f).readlines()[1:]).strip().replace('\n','').replace('-',''))
			lengths+=[length]
	return lengths
#############################
def fragment_protein(sequence, lengths):
	seq_length = len(sequence)
	fragments = []
	for db_frag_size in lengths:
		# sequence needs to be at least as long as the fragment
		if seq_length == db_frag_size:
			# collect & index
			fragments += [[sequence,0]]
		elif seq_length >  db_frag_size:
			# collect from NH3 terminus to end
			for i in range( seq_length-db_frag_size+1):
				# collect & index
				fragments += [[sequence[i:i+db_frag_size],i]]
	return fragments
#############################
def window_protein(sequence, window_length):
	seq_length = len(sequence)
	fragments = []
	for i in range(len(sequence)-window_length):
		fragments += [[sequence[i:i+window_length],i]]
	return fragments
#############################
#############################
def write_fasta(sequence,fasta_filename):
	writer=open(fasta_filename,'w')
	writer.write('>'+fasta_filename.split('.')[0]+'\n'+sequence+'\n')
	writer.close()
#############################
def write_fasta_lib(sequences,fasta_filename):
	writer=open(fasta_filename,'w')
	index=1
	for sequence in sequences:
		writer.write('>fragment_'+str(index)+'\n'+sequence+'\n')
	writer.close()
##############################################
##############################################


		#########
		# START #
		#########


if __name__=='__main__':
	#####################################
	# this program takes as input       #
	# [1]- a FASTA file                 #
	# [2]- a directory of FASTA files + #
	# [3]- a directory of FASTA files - #
	# [4]- a scoring matrix             #
	# [5,6]- gap penalties (gop, gep)   #
	# calculates total similarity score #
	#####################################
	
	# break up query protein into all possible fragments 
	# that match the size of sequences in the database
	# score position by TSS of all fragments covering the position
	# normalize by amount of fragments considered at position,
	# perhaps multiplied by the size of the fragment itself
	# test by recapture
	
	try:
		######################
		# prepare input sets #
		######################
		query_file       = argv[1]
		
#		db_strong_dir    = argv[2]+'/'
#		db_weak_dir      = argv[3]+'/'
		db_strong_dir    = 'binders/'
		db_weak_dir      = 'nonbinders/'
		matrix_file_name = argv[2]
		gap_open         = argv[3]
		gap_extension    = argv[4]
		
		
		#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
		# there are 2 conceptually distinct ways of doing this #
		# with the full matching compliment of fragments       #
		# 1. TSS of all fragments possible at the position     #
		# 2. highest TSS of the applicable fragments           #
		#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
		scoring_method = 1
		if '-mean' in argv:  scoring_method = 1
		elif '-max' in argv:  scoring_method = 2
		#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
		# there are 2 conceptually distinct ways of doing this #
		# with a sliding window of fragments                   #
		# 3. TSS of all fragments possible at the position     #
		# 4. highest TSS of the applicable fragments           #
		#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
		frags_built_to_match_distribution = True
		if '-window' in argv:  frags_built_to_match_distribution = False
		
		# report scores of each residue (default) or whole protein
		residues = True
		if '-w' in argv:
			residues = False
		
		# grab query sequence
		query_seq=''
		for line in open(query_file).readlines():
			if not line.startswith('>') and not line.startswith('#') and line.strip():
				query_seq+=line.strip()
		THIS_protein = query_file.split('/')[-1].split('.')[0].split('_')[0]
		#@#	print 'the query sequence is:',query_seq
		
		#################
		### strong db ###
		#################
		# input database directory
		db_strong_lib = "lib_db_strong.fasta"
		db_strong_fastas = dir_to_fastalist_N_lib(db_strong_dir,db_strong_lib)
		
		#################
		###  weak db  ###
		#################
		# input database directory
		db_weak_lib = "lib_db_weak.fasta"
		db_weak_fastas = dir_to_fastalist_N_lib(db_weak_dir,db_weak_lib)
		
		
		###############################################
		
		if frags_built_to_match_distribution:
			# measure the variability in db sequence length
			db_strong_lengths = db_seq_lengths(db_strong_dir)
			db_strong_lengths.sort()
			# fragment protein by sizes in db
			fragments_matching_distribution = fragment_protein(query_seq, db_strong_lengths)
		else:
			# measure the variability in db sequence length
			lengths = db_seq_lengths_all(db_strong_dir)
			window_length = int(round(sum(lengths) / float(len(lengths))))
			# fragment protein by windows, sized by avg in db
			fragments_matching_distribution = window_protein(query_seq, window_length)
		# fragments[i] = ['SEQUENCE',start position]
		###############################################
		
		query_dir    = 'tmp_mat/'
		query_lib    = query_dir+'lib_query.fasta'
		try: os.mkdir(query_dir)
		except: do_nothing = True
		
		######################################
		# run TSS on all fragments
		name = query_dir+'fragment.fasta'
		for fragment in range(len(fragments_matching_distribution)):
			sequence = fragments_matching_distribution[fragment][0]
			write_fasta(sequence,name)
			########################################        ###
			# calculate TSS for each fragment in the strong set
			TSS_s = TSS([name], name, db_strong_fastas, db_strong_lib, gap_open, gap_extension, matrix_file_name)
			########################################      ###
			# calculate TSS for each fragment in the weak set
			TSS_w = TSS([name], name, db_weak_fastas, db_weak_lib, gap_open, gap_extension, matrix_file_name)
			fragments_matching_distribution[fragment] += [TSS_s, TSS_w]
		
		if not residues:
			all_scores = []
		for position in range(len(query_seq)):
			TSS_strongs = []
			TSS_weaks   = []
			#########################################
			# find fragments relevant to the position
			for fragment in fragments_matching_distribution:
				# fragment = ['SEQUENCE',start_position, TSS_s, TSS_w]
				if position >= fragment[1] and position < fragment[1]+len(fragment[0]):
					TSS_strongs += [fragment[2]]
					TSS_weaks   += [fragment[3]]
			##########################
			# print the position score
			if scoring_method == 1:  #mean
				try:
					score = -1000 * ( (sum(TSS_strongs) / len(TSS_strongs)) - (sum(TSS_weaks) / len(TSS_weaks)) )
					norm_score = 100 * (max_value - score) / (max_value - min_value)
					if residues:
						print score
					else:
						all_scores += [score]
				except: print 0
			if scoring_method == 2:  #max
				try: 
					score = -1000 * (max(TSS_strongs) - max(TSS_weaks))
					norm_score = 100 * (max_value - score) / (max_value - min_value)
					if residues:
						print score
					else:
						all_scores += [score]
				except: print 0
		if not residues:
			print sum(all_scores) / len(all_scores)
	#@#			try: print -1000 * ((max(TSS_strongs)/len(TSS_strongs)) - (max(TSS_weaks) / len(TSS_weaks)))
	#@#			try: print -1000 * (min(TSS_strongs) - min(TSS_weaks))
				########################################################
	
	except:
		# [1]- a FASTA file                 #
		# [2]- a directory of FASTA files + #
		# [3]- a directory of FASTA files - #
		# [4]- a scoring matrix             #
		# [5,6]- gap penalties (gop, gep)   #
#		print "Usage: MATRIX_SCORER.py <fasta file> <dir of fxn> <dir of nonfxn> <scoring matrix> <gop> <gep>"
		print "Usage: MATRIX_SCORER.py <fasta file> <scoring matrix> <gop> <gep>"
		print "\toptions:\t-w\twhole protein score"
		print "\t\t\t-max\tscore by maximum similarity value for each residue rather than the mean\n"


