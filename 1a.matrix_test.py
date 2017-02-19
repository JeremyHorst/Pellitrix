#!/usr/bin/python

# matrix_test.py
# run matrix similarity on both binders & nonbinders
# normalize output score range to 100-0

maxvalue = 100
minvalue = 0

from sys import argv
from subprocess import Popen
from subprocess import PIPE
import os

matrix_method_dir =	'./'
bindir = 		'./binders/'
nonbindir = 		'./nonbinders/'

matrix = argv[1]
gop = argv[2]
gep = argv[3]

binders = os.listdir(bindir)
nonbinders = os.listdir(nonbindir)
binder_scores = []
for binder in binders:
	command = "%s/1b.matrix_scorer.py %s/%s %s %s %s %s %s" % (matrix_method_dir, bindir,binder, bindir, nonbindir, matrix, gop, gep)
	binder_scores += [float(Popen(command.split(),stdout=PIPE).communicate()[0].strip())]
nonbinder_scores = []
for nonbinder in nonbinders:
	command = "%s/1b.matrix_scorer.py %s/%s %s %s %s %s %s" % (matrix_method_dir, nonbindir,nonbinder, bindir, nonbindir, matrix, gop, gep)
	nonbinder_scores += [float(Popen(command.split(),stdout=PIPE).communicate()[0].strip())]

most = max( max(binder_scores) , max(nonbinder_scores) )
least = min( min(binder_scores) , min(nonbinder_scores) )
longest = max( len(binder_scores) , len(nonbinder_scores) )

def norm(instance,highest,lowest):
	return ((instance-lowest)/(highest-lowest))*(maxvalue-minvalue) + minvalue

for score in range(longest):
	printline = ''
	try: printline += str(    norm(binder_scores[score],most,least)    )+'\t'
	except: printline += '\t'
	try: printline += str(    norm(nonbinder_scores[score],most,least)    )+'\t'
	except: printline += '\t'
	print printline

