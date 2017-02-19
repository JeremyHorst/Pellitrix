#!/usr/bin/python

# matrix_mutator.py
# Jeremy A Horst, 04/04/2010

################################################################
# this program takes as input two directories of FASTA files,  #
# trains both the gap opening and extension penalty variables, #
# and trains a scoring matrix to maximize the score difference #
# requirements:
# - python2.5
# - FASTA software package
################################################################



# details:
#	- all fastas must have equivalent name.fasta & first line >name, with white space after ">name".
# 	- if matrices are to be run as integers rather than float point values, switch "#@@" comment from each line to the corresponding (next) line.

##################
# matrix details #

#@@value_limit = 20
value_limit = 2.00

ggsearch_dir = './fasta-35.4.11/'

matrix_dir   = './'
AA_order = "A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X".split()
matrix_size = len(AA_order)
##################


from sys import argv
from os import mkdir
from os import listdir
from random import shuffle
from random import randrange
from numpy import argmax
from subprocess import call
from subprocess import Popen
from subprocess import PIPE


#make process dir
try: mkdir('tmp')
except: nada=True


#############################
#############################
def directory_to_library(fasta_set,library):
	writer = open(library,'w')
	for f in fasta_set:
		for line in open(f).readlines():
			writer.write(line.strip()+'\n')
	writer.close()
#############################
def train_gap_penalties(matrix_file_name):
	###############################################
	# gop = gap open penalties:       -16, 0, +1  #
	# gep = gap extension penalties:   -8, 8, +1  #
	###############################################
	
	# Emre's paper says in the text (section 3.1):
	# we also acquiesse on the positive thing.
	# but I don't think we want to minimize BOTH
	# rather we should maximize the difference!!!
	
	report_count=0
	first = True
	gap_open_n_extend = [0,0]
	for gap_open in range(-16,0):
		for gap_extension in range(-8,0):
			gap_open      = str(gap_open)
			gap_extension = str(gap_extension)
			# run TSSs
			if   TSS_min_diff:   this_combination = TSS_minimum_difference(gap_open, gap_extension, matrix_file_name)
			elif TSS_distr_mult: this_combination = TSS_distribution_mult(gap_open, gap_extension, matrix_file_name)
			elif TSS_distr_add:  this_combination = TSS_distribution_add(gap_open, gap_extension, matrix_file_name)
			
			if first:
				best_combination = this_combination
				gap_open_n_extend = [gap_open, gap_extension]
				first = False
			if this_combination > best_combination:
				best_combination = this_combination
				gap_open_n_extend = [gap_open, gap_extension]
			#
			report_count+=1
			if gap_open_n_extend[0]:
				print 'gap training... TSS full calculation iterations complete:',report_count,'\t',gap_open_n_extend[0],gap_open_n_extend[1],'\t',gap_open, gap_extension,' '*(4-len(str(gap_extension))),'\tbest_combination:   ',best_combination,'\tthis_combination',this_combination
			else:
				print 'gap training... TSS full calculation iterations complete:',report_count,'\t',gap_open_n_extend[0],gap_open_n_extend[1],'\t',gap_open, gap_extension
	return gap_open_n_extend
#############################
#############################
#############################
def TSS_minimum_difference(gap_open, gap_extension, matrix_file_name):
	# the goal is to separate the strong from weak binders
	# yet maximizing between worst scoring strong & best scoring weak
	# will skew to this peptides AND
	# AND allow for the inconsistencies of ggsearch
	# so... maximize between the THIRD worst strong & THIRD best weak
	############################
	TSS_strong = []
	for high_fasta in high_fastas:
		TSS_ss = TSS([high_fasta], high_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_sw = TSS([high_fasta], high_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_strong += [  TSS_ss - TSS_sw  ]
	
#@@	TSS_strong.sort()
#@@	third_smallest_strong_TSS = TSS_strong[2]
	TSS_strong.sort(reverse=True)
	third_BIGGEST_strong_TSS = TSS_strong[2]
	############################
	TSS_weak = []
	for low_fasta in low_fastas:
		TSS_ws = TSS([low_fasta], low_fasta, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
		TSS_ww = TSS([low_fasta], low_fasta, low_fastas,  low_lib,  gap_open, gap_extension, matrix_file_name)
		TSS_weak += [  TSS_ws - TSS_ww  ]
	
#@@	TSS_weak.sort(reverse=True)
#@@	third_biggest_weak_TSS = TSS_weak[2]
	TSS_weak.sort()
	third_SMALLEST_weak_TSS = TSS_weak[2]
	############################
	
#@@	return (third_smallest_strong_TSS - third_biggest_weak_TSS)   * 1000
	return (third_SMALLEST_weak_TSS   - third_BIGGEST_strong_TSS) * 1000

#############################
def TSS_distribution_mult(gap_open, gap_extension, matrix_file_name):
	TSS_ss = TSS(high_fastas, high_lib, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_sw = TSS(high_fastas, high_lib,  low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ws = TSS(low_fastas,  low_lib,  high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ww = TSS(low_fastas,  low_lib,   low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	
#@@	TSScore = (TSS_ss * TSS_ww) / (TSS_sw * TSS_ws)
	TSScore = (TSS_sw * TSS_ws) / (TSS_ss * TSS_ww)
	
	return TSScore
#############################
def TSS_distribution_add(gap_open, gap_extension, matrix_file_name):
	TSS_ss = TSS(high_fastas, high_lib, high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_sw = TSS(high_fastas, high_lib,  low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ws = TSS(low_fastas,  low_lib,  high_fastas, high_lib, gap_open, gap_extension, matrix_file_name)
	TSS_ww = TSS(low_fastas,  low_lib,   low_fastas,  low_lib, gap_open, gap_extension, matrix_file_name)
	
#@@	TSScore = TSS_ss + TSS_ww - TSS_sw - TSS_ws
	TSScore = TSS_sw + TSS_ws - TSS_ss - TSS_ww
	
	return TSScore
#############################
#############################
#############################
def TSS(fastasA, setA, fastasB, setB, gap_open, gap_extension, substitution_matrix):
	# TSS_A-B([A]_NA - [B]_NB) = 
	# = (1/ (NA*delta_AB)) sum{1toNA}[ sum{1toNB}[ PSS_ij(1 - delta_ij * delta_AB)]]
	# we add normalization of the second term by the sequence length (not considered before)
	
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
		seq_length = float(len(open(protein_a).readlines()[1].strip()))
		################
		# run ggsearch35
		ggsearch_output = 'tmp/'+protein_a.split('/')[-1].split('.')[0]+'.ggsearch_output'
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
			delta_ab = 0
			if protein_a == protein_b:  delta_ab = 1 #; print 'caught self'
			
			# strip '.fasta' in name for indexing
			protein_b = protein_b.split('/')[-1].split('.')[0]
			# if not present, assume PSS_ab = 0
			PairSimScore_ab = 0
			if protein_scores.has_key(protein_b):
				PairSimScore_ab = protein_scores[protein_b]
			###########################
			###########################
			# normalize to seq length #
			#   sum to second term    #
			second_term += PairSimScore_ab * (1 - delta_ab * delta_AB) / seq_length
			###########################
			###########################
	return first_term * second_term
############################
def run_ggsearch(gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library):
	command = "%s/bin/ggsearch35 -q -p -m 4 -f %s -g %s -s %s -O %s %s %s"\
	% (ggsearch_dir, gap_open, gap_extension, substitution_matrix, ggsearch_output, fasta, library)
	crap = Popen(command.split(),stderr=PIPE,stdout=PIPE).communicate()[0]
#############################
def matrix_evolution(base_substitution_matrix, gap_open, gap_extension):
	
	###################################
	#              CONCEPT            #
	###################################
	# columns & rows mutated en mass to catch large tendancies of amino acid types
	# columns mutated before rows to catch features of the compared target
	#         i.e. the S & W of ss & sw
	# order of individual columns/rows/cells to follow:
	# first follow greedy trajectory by checking each at each step
	# then random local perturbations
	# then monte carlo perturbations
	###################################
	
	# load matrix file into dictionary
	matrix_dict = matrix_file_to_dict(base_substitution_matrix)
	# write matrix dictionary into file
	matrix_file_name = base_substitution_matrix.split('/')[-1].split('.mat')[0]+'_revised.mat'
	write_matrix_from_dict(matrix_dict, matrix_file_name)
	write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
	
	###########################
	# calculate base accuracy #
	###########################
	if   TSS_min_diff:   TSS_base = TSS_minimum_difference(gap_open, gap_extension, matrix_file_name)
	elif TSS_distr_mult: TSS_base = TSS_distribution_mult(gap_open, gap_extension, matrix_file_name)
	elif TSS_distr_add:  TSS_base = TSS_distribution_add(gap_open, gap_extension, matrix_file_name)
	TSS_laststep = TSS_base
	###########################
	
	
	
	#########################
	#@# mutate dictionary #@#
	#########################
	
	#########################################
	# start with a deterministic projection #  i.e. GREEDY
	#########################################
	print 'greedy path...'
	if columns_then_cells_OR_cells_only == 0:
		# mutate columns in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 0)
		# mutate rows in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 1)
		# mutate all individual cells in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
		# (this includes local minimization, but if running columns & rows, can do again)
	elif columns_then_cells_OR_cells_only == 1:
		# mutate all individual cells in order of benefit
		matrix_dict, TSS_laststep = deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
		# (this includes local minimization, but if running columns & rows, can do again)
	
	##########################################
	# proceed with random local minimization #
	##########################################
	if columns_then_cells_OR_cells_only == 0:
		print 'local minimization...'
		# mutate columns in random order
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 0)
		# mutate rows in random order
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 1)
	
	########################################
	# finish with Monte Carlo minimization #
	########################################
	print 'monte carlo minimization...'
	# repeat with Monte Carlo search
	if columns_then_cells_OR_cells_only == 0:
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 0)
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 1)
		# mutate all individual cells in random order #
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
	elif columns_then_cells_OR_cells_only == 1:
		matrix_dict, TSS_laststep = mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, 2)
	########################################
	return matrix_file_name, TSS_laststep
#########################################
def matrix_file_to_dict(base_substitution_matrix):
	# instantiate dictionary to hold matrix
	matrix_dict = {}
	for AA in AA_order:
		matrix_dict[AA]={}
		for aa in AA_order:
			matrix_dict[AA][aa]=0
	# load in matrix
	for line in open(base_substitution_matrix).readlines():
		if not line.startswith('#') and line.strip():
			if line[0].strip():
				wt_AA = line.split()[0]
				if wt_AA in AA_order:
					for mut_aa in range(matrix_size):
#@@						matrix_dict[wt_AA][AA_order[mut_aa]] = int(line.split()[mut_aa+1])
						matrix_dict[wt_AA][AA_order[mut_aa]] = float(line.split()[mut_aa+1])
	return matrix_dict
#############################
def write_matrix_from_dict(matrix_dict, matrix_file_name):
	spacer = '   '
	writer = open(matrix_file_name,'w')
	writer.write('    ' + spacer.join(AA_order)+'\n')
	for AA in AA_order:
		printline = AA+' '
		for aa in AA_order:
#@@			value = str(  matrix_dict[AA][aa]  )[:4]
			value = str(  round(matrix_dict[AA][aa],2)  )[:4]
			if value[0]=='-':  str(  round(matrix_dict[AA][aa],2)  )[:5]
#@@			printline += ' '*(3-len( value )) + value + ' '
			printline += ' '*(5-len( value )) + value + ' '
		writer.write(printline.strip()+'\n')
	writer.close()
#############################
def deterministic_mutation(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, row_or_column_or_cell):
	
	# try mutating all columns/rows/cells
	# return to original
	# follow greedy path
	type_crc = ['column','row','cell']
	
	count_iterations = 0
	if row_or_column_or_cell < 2:
		aa=False
		
		AAs_done = []
		for i in range(matrix_size):
			# find best AA to improve
			# follow trajectory until end of path
			# find next best AA to improve, repeat
			# go through all AAs, only ONCE (AAs_done)
			
			increment = 0.04  #this is similar to increments of 1 with min=-19 and max=19
			
			greedy_TSS        = []
			greedy_trajectory = []
			
			# find best AA to improve
			for AA in AA_order:
				if not AA in AAs_done:
					#
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					TSS_up,matrix_dict = check_row_column_cell_change(+1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					TSS_down,matrix_dict = check_row_column_cell_change(-1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					#
					trajectory = 0;  TSS_nextstep = 0
					if   TSS_up   > TSS_down:  trajectory = +1;  TSSdif = TSS_up
					elif TSS_down >   TSS_up:  trajectory = -1;  TSSdif = TSS_down
					elif TSS_down == TSS_up:   trajectory = -1;  TSSdif = TSS_down
					greedy_TSS        += [TSSdif]
					greedy_trajectory += [trajectory]
					
#					print AA,TSSdif,'vs',TSS_laststep
				else:
					greedy_TSS        += [0]
					greedy_trajectory += [0]
			
			
			if len(AAs_done) < matrix_size:
				# follow trajectory until end of path
				AA           = AA_order[argmax(greedy_TSS)]
				trajectory   = greedy_trajectory[argmax(greedy_TSS)]
				TSS_nextstep = greedy_TSS[argmax(greedy_TSS)]
				print AA,'is the best next',type_crc[row_or_column_or_cell],' at step',len(AAs_done)+1,':',TSS_nextstep,'vs',TSS_laststep
				
				# still we only follow the trajectory if there is improvement
				if TSS_nextstep > TSS_laststep:
					
					crap,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, False)
					
					while TSS_nextstep > TSS_laststep:
						# save
						write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
						TSS_laststep = TSS_nextstep
						
						# move along trajectory
						TSS_nextstep,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						
						count_iterations += 1
						print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',TSS_laststep
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					AAs_done += [AA]
				else:
					# skip the rest; there will be no improvement
					AAs_done += [AA]
					for aaaa in AA_order:
						if not aaaa in AAs_done:
							AAs_done += [aaaa]
	
	# instance of one amino acid pair at a time
	elif row_or_column_or_cell == 2:
		print 'cell time'
		AA_combos = []
		for AA in AA_order:
			for aa in AA_order:
				if not [aa,AA] in AA_combos:
					AA_combos += [[AA,aa]]
		
		AAs_done = []
		how_many_we_will_do = len(AA_combos)
		
		for i in range(how_many_we_will_do):
			if len(AAs_done) < how_many_we_will_do:
				# find best AA combination to improve
				# follow trajectory until end of path
				# find next best AA to improve, repeat
				# go through all AAs, only ONCE (AAs_done)
				
#@@				increment = 1
				increment = 0.04  #this is most similar to -19 to 19
				
				greedy_TSS        = []
				greedy_trajectory = []
				
				# find best AA to improve
				for AA,aa in AA_combos:
					if not [AA,aa] in AAs_done and not [aa,AA] in AAs_done:
						
						matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
						TSS_up,matrix_dict = check_row_column_cell_change(+1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
						TSS_down,matrix_dict = check_row_column_cell_change(-1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
						
						trajectory = 0;  TSS_nextstep = 0
						if   TSS_up   > TSS_down:  trajectory = +1;  TSSdif = TSS_up
						elif TSS_down >   TSS_up:  trajectory = -1;  TSSdif = TSS_down
						elif TSS_down == TSS_up:   trajectory = -1;  TSSdif = TSS_down
						if TSSdif > TSS_laststep:
							print AA,aa,TSS_laststep,TSSdif
						greedy_TSS        += [TSSdif]
						greedy_trajectory += [trajectory]
					else:
						greedy_TSS        += [0]
						greedy_trajectory += [0]
				
				if not len(AA_combos) == len(greedy_TSS) or not len(AA_combos) == len(greedy_trajectory):
					print "len(AA_combos) == len(greedy_TSS), len(AA_combos) == len(greedy_trajectory)"
					print len(AA_combos) == len(greedy_TSS), len(AA_combos) == len(greedy_trajectory)
					print "len(AA_combos), len(greedy_TSS), len(AA_combos), len(greedy_trajectory)"
					print len(AA_combos), len(greedy_TSS), len(AA_combos), len(greedy_trajectory)
				
				# follow trajectory until end of path
				AA,aa        = AA_combos[        argmax(greedy_TSS)]
				trajectory   = greedy_trajectory[argmax(greedy_TSS)]
				TSS_nextstep = greedy_TSS[       argmax(greedy_TSS)]
				print AA,aa,'is the best at step:',len(AAs_done)+1,TSS_laststep,TSS_nextstep
				
				# still we only follow the trajectory if there is improvement
				if TSS_nextstep > TSS_laststep:
					crap,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, False)
					while TSS_nextstep > TSS_laststep:
						# save
						write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
						TSS_laststep = TSS_nextstep
						
						TSS_nextstep,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
						count_iterations += 1
						print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',TSS_laststep
					matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
					AAs_done += [[AA,aa]]
					if not AA==aa:  AAs_done += [[aa,AA]]
				else:
					# skip the rest; there will be no improvement
					AAs_done += [[AA,aa]]
					if not AA==aa:  AAs_done += [[aa,AA]]
					for AAAaaa in AA_combos:
						if not AAAaaa in AAs_done:
							AAs_done += [AAAaaa]
	return matrix_dict, TSS_laststep
############################
def check_row_column_cell_change(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, run_TSS):
	# columns
	if row_or_column_or_cell == 0:
		for aa in AA_order:
			matrix_dict[aa][AA] = check_to_maxmin(matrix_dict[aa][AA], increment)
	# row
	elif row_or_column_or_cell == 1:
		for aa in AA_order:
			matrix_dict[AA][aa] = check_to_maxmin(matrix_dict[AA][aa], increment)
	# cells
	elif row_or_column_or_cell == 2:
		matrix_dict[AA][aa] = check_to_maxmin(matrix_dict[AA][aa], increment)
		if not AA==aa:
			matrix_dict[aa][AA] = check_to_maxmin(matrix_dict[aa][AA], increment)
	# calc
	write_matrix_from_dict(matrix_dict, matrix_file_name)
	if run_TSS:
		if   TSS_min_diff:   TSS = TSS_minimum_difference(gap_open, gap_extension, matrix_file_name)
		elif TSS_distr_mult: TSS = TSS_distribution_mult(gap_open, gap_extension, matrix_file_name)
		elif TSS_distr_add:  TSS = TSS_distribution_add(gap_open, gap_extension, matrix_file_name)
	else: TSS = 0
	return TSS, matrix_dict
#############################
def check_to_maxmin(valueA, increment):
	out = 0
	if valueA+increment > value_limit:
		out = value_limit
	elif valueA+increment < -1*value_limit:
		out = -1*value_limit
	else:
		out = valueA+increment
	return out
############################
def mutate_row_or_column_or_cell(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, row_or_column_or_cell):
	
	increment = 0.04  #this is similar to increments of 1 with min=-19 and max=19
	
	count_iterations = 0
	# mutate rows in random order
	if row_or_column_or_cell < 2:
		row_order  = []
		row_order += AA_order
		shuffle(row_order)
		for AA in row_order:
			matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, False)
	elif row_or_column_or_cell == 2:
		AA_matrix = []
		for AA in AA_order:
			for aa in AA_order:
				if not [aa,AA] in AA_matrix:
					AA_matrix += [[AA,aa]]
		shuffle(AA_matrix)
		for (AA,aa) in AA_matrix:
			matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, aa)
	return matrix_dict, TSS_laststep
############################
def mutate_row_or_column_or_cell_MonteCarlo(TSS_laststep, matrix_dict, matrix_file_name, gap_open, gap_extension, row_or_column_or_cell):
	count_iterations = 0
	
	###############################
	# mutate rows in random order #
	###############################
	if row_or_column_or_cell < 2:
		row_order  = []
		row_order += AA_order
		shuffle(row_order)
		for AA in row_order:
		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
		#$$$$$  Bounce out to random distance  $$$$$#
		#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$#
			for i in range(5):
				# increment should not be zero, but pos & neg are vital
				
			#@@	increment = 1
			#@@	increment = 0.04
			#@@	increment = randrange(10)
				increment = randrange(43)/10.0	#this is most similar to -19 to 19 / 10
				
				posneg = randrange(2)
				if posneg: increment = -1*increment
				# run...
				matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, False)
	elif row_or_column_or_cell == 2:
		AA_matrix = []
		for AA in AA_order:
			for aa in AA_order:
				if not [aa,AA] in AA_matrix:
					AA_matrix += [[AA,aa]]
		shuffle(AA_matrix)
		for (AA,aa) in AA_matrix:
			for i in range(5):
				# increment should not be zero, but pos & neg are vital
				
			#@@	increment = randrange(10)
				increment = randrange(43)/10.0	#this is most similar to -19 to 19 / 10
				
				posneg = randrange(2)
				if posneg: increment = -1*increment
				# run...
				matrix_dict, TSS_laststep, count_iterations = increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, aa)
	return matrix_dict, TSS_laststep
############################
def increment_row_or_column_or_cell(increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, count_iterations, TSS_laststep, aa):
	#@# first, get trajectory
	
	type_crc = ['column','row','cell']
	###################
	# save dictionary #
	###################
	matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	###################
	
	##################
	# try adding one #
	##################
	TSS_up,matrix_dict = check_row_column_cell_change(+1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
	matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	##################
	
	#######################
	# try subtracting one #
	#######################
	TSS_down,matrix_dict = check_row_column_cell_change(-1*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
	matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
	
	###########
	# compare #
	###########
	trajectory = 0;  TSS_nextstep = 0
	if   TSS_up   > TSS_down:  trajectory = +1;  TSS_nextstep = TSS_up
	elif TSS_down >   TSS_up:  trajectory = -1;  TSS_nextstep = TSS_down
	elif TSS_down ==  TSS_up:  trajectory = -1;  TSS_nextstep = TSS_down
	
	count_iterations += 1
	print AA,aa,count_iterations,'iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',TSS_laststep,TSS_nextstep
	###########################################################
	#@# now, move along trajectory unil no more improvement #@#
	###########################################################
	if TSS_nextstep > TSS_laststep:
		
		# push back to improved trajectory; the next line will proceed... n+1+1-1=n+1
		matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
		crap,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, False)
		
		while TSS_nextstep > TSS_laststep:
			# save
			write_matrix_from_dict(matrix_dict, 'save.'+matrix_file_name)
			TSS_laststep = TSS_nextstep
			
			TSS_nextstep,matrix_dict = check_row_column_cell_change(trajectory*increment, AA, row_or_column_or_cell, matrix_dict, matrix_file_name, gap_open, gap_extension, aa, True)
			count_iterations += 1
			print AA,aa,count_iterations,'IN iterations of',type_crc[row_or_column_or_cell],'mutations','\tTSS_laststep:',TSS_laststep
		
		# revert to most recent BTTR form
		matrix_dict = matrix_file_to_dict('save.'+matrix_file_name)
		
	return matrix_dict, TSS_laststep, count_iterations
############################
#########################################


		#########
		# START #
		#########


if __name__=='__main__':
#	if 1==1:
	try:
		######################
		# prepare input sets #
		######################
		# input directories
		high_dir = argv[1]+'/'
		low_dir = argv[2]+'/'
		# input fasta lists
		high_fastas=[]
		for f in listdir(high_dir):
			if f.endswith('.fasta'): high_fastas+=[high_dir+f]
		low_fastas=[]
		for f in listdir(low_dir):
			if f.endswith('.fasta'): low_fastas+=[low_dir+f]
		# make searchable libraries
		high_lib = "lib_high.fasta"
		low_lib  = "lib_low.fasta"
		directory_to_library(high_fastas,high_lib)
		directory_to_library(low_fastas,low_lib)
		
		###########
		# options #
		###########
		columns_then_cells_OR_cells_only = 0
		if '-cells' in argv:  columns_then_cells_OR_cells_only = 1
		matrix = 'pam250.mat'
		if '-matrix' in argv:  matrix = argv[argv.index('-matrix')+1]
		
		TSS_distr_add  = True; TSS_min_diff = False; TSS_distr_mult = False
		if   '-multiply' in argv: TSS_distr_mult = True; TSS_distr_add = False
		elif '-maximize_diff' in argv: TSS_distr_mult = True; TSS_distr_add  = False
		
		gaps_only = False
		if '-gaps_only' in argv: gaps_only = True  #train gap open values ONLY
		gap_open_pam=False; gap_extension_pam=False
		if '-gop_gep' in argv: # set gap open & extend values
			gap_open_pam = int(argv[argv.index('-gop_gep')+1])
			gap_extension_pam = int(argv[argv.index('-gop_gep')+2])
			
		###############################################################
		# the above names & lists will be used throughout the process #
		###############################################################
		
		
		
		####################################################################
		substitution_matrix = matrix_dir + matrix
		if not gap_open_pam or not gap_extension_pam:
			gap_open_pam, gap_extension_pam = train_gap_penalties(substitution_matrix)
		print matrix,'  gap open penalty:',gap_open_pam,  'gap extend penalty:',gap_extension_pam
		####################################################################
		
		
		
		#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
		# mutate each matrix to maximize (TSS_strong-strong - TSS_strong-weak) #
		#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@##@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
		
		if not gaps_only:
			pam_matrix_file_name,    TSS_laststep_pam    = matrix_evolution(matrix_dir+matrix,   gap_open_pam,    gap_extension_pam)
			print 'PAM acheives TSS difference of:',TSS_laststep_pam
			print 'PAM matrix file:',pam_matrix_file_name
			print matrix,'  gap open penalty:',gap_open_pam,  'gap extend penalty:',gap_extension_pam
#	if 1==2:
	except:
		print "\nUsage: matrix_mutator.py <desired_peptides/> <unwanted_peptides/>"
		print "Options:"
		print "\t-cells\t\t\tmutate cell by cell; trajectory & then monte carlo"
		print "\t\t\t\tdefault is trajectory by column,row,cell; then monte carlo by column,row"
		print "\t-matrix <matrix>\tbase substitution matrix; pam250 is the default matrix\n\n"
		print "\t-multiply\t\t\tmatrix mutation driving force is ss*ww/(sw*ws)"
		print "\t-maximize_diff\t\t\tmaximize separation of 3rd highest scoring weak & 3rd lowest scoring strong peptide"
		print "\t\t\t\tdefault matrix mutation driving force is ss+ww-sw-ws"

