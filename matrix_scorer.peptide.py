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
min_value = -0.01371074746
max_value = 0.01116456299
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
	
	try:
		# break up query protein into all possible fragments 
		# that match the size of sequences in the database
		# score position by TSS of all fragments covering the position
		# normalize by amount of fragments considered at position,
		# perhaps multiplied by the size of the fragment itself
		# test by recapture
		
		######################
		# prepare input sets #
		######################
		query_file       = argv[1]
	#	db_strong_dir    = argv[2]+'/'
	#	db_weak_dir      = argv[3]+'/'
		db_strong_dir    = 'binders/'
		db_weak_dir      = 'nonbinders/'
		matrix_file_name = argv[2]
		try:
			gap_open         = argv[3]
			gap_extension    = argv[4]
		except:
			see_below = True
		
		# normalization
		if matrix_file_name == 'b3DAli.trained.mat':
			max_value = 0.011164563
			min_value = -0.013710747
			gap_open = -5; gap_extension = -2
		elif matrix_file_name == 'b3DAli.mat':
			max_value = 0.002823529
			min_value = -0.005347594
			gap_open = -5; gap_extension = -2
		elif matrix_file_name == 'pam250.mat':
			max_value = -0.019898452
			min_value = 0.025114926
			gap_open = -5; gap_extension = -1
		elif matrix_file_name == 'pam250.trained.mat':
			max_value = -0.018220846
			min_value = 0.066990337
			gap_open = -4; gap_extension = -4
		else:
			max_value = False; min_value = False
		
		
		# grab query sequence
		query_seq=''
		for line in open(query_file).readlines():
			if not line.startswith('>') and not line.startswith('#') and line.strip():
				query_seq+=line.strip()
		#@#	print 'the query sequence is:',query_seq
		
		#################
		### strong db ###
		#################
		# input database directory
		db_strong_lib = "lib_db_strong.fasta"
		db_strong_fastas = dir_to_fastalist_N_lib(db_strong_dir,db_strong_lib)
		
		#################
		### weak db ###
		#################
		# input database directory
		db_weak_lib = "lib_db_weak.fasta"
		db_weak_fastas = dir_to_fastalist_N_lib(db_weak_dir,db_weak_lib)
		
		########################################        ###
		# calculate TSS for each fragment in the strong set
		TSS_s = TSS([query_file], query_file, db_strong_fastas, db_strong_lib, gap_open, gap_extension, matrix_file_name)
		
		########################################      ###
		# calculate TSS for each fragment in the weak set
		TSS_w   = TSS([query_file], query_file, db_weak_fastas, db_weak_lib, gap_open, gap_extension, matrix_file_name)
		#################################################
		
		#################################################################
	#	print "Predictions"
		score = TSS_s - TSS_w
		if max_value:
			norm_score = 100*(max_value - score) / (max_value - min_value)
			print norm_score
		else:
			print score
		

##############################################
	except:
		# [1]- a FASTA file                 #
		# [2]- a directory of FASTA files + #
		# [3]- a directory of FASTA files - #
		# [4]- a scoring matrix             #
		# [5,6]- gap penalties (gop, gep)   #
		print "Usage: matrix_scorer.peptide.py <fasta file> <scoring matrix> <gop> <gep>"

##############################################
