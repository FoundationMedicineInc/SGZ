#!/usr/bin/env python

import os, sys, glob
from numpy import *
from scipy.stats import binom_test
from scipy.stats import binom
import csv
import argparse
# import pdb

__author__  = 'Yuting He'
__contact__ = 'yhe@foundationmedicine.com'
__version__ = '1.0.0'
__doc__     = '''Basic method to predict CNA origin and zygosity.'''

def call_SGZ(short_variants, pathology_purity, alpha, high_purity_tresh):
	sv_sgz = []
	for sv in short_variants:
		if sv.chr_ == 'chrX':continue # Ignore sv from chrX
		sv.zygosity, sv.SG_prediction = core_SGZ(sv.frequency, sv.depth, pathology_purity, alpha, high_purity_tresh)
		sv_sgz.append(sv)

	return sv_sgz

def core_SGZ(frequency, depth, pathology_purity, alpha = 0.05, high_purity_tresh = 0.9):
	if frequency > 0.95:
		zygosity = 'homozygous'
		if pathology_purity == 'NA':
			SG_prediction = 'ambiguous'
		elif pathology_purity < high_purity_tresh:
			SG_prediction = 'germline'
		else:
			SG_prediction = 'ambiguous'
	else:
		pvalue = binom_test(frequency*depth, depth, p=0.5)
		if pvalue< alpha:
			SG_prediction = 'somatic'
			zygosity  = 'ambiguous'
		else:
			SG_prediction = 'germline'
			zygosity  = 'het'

	return zygosity, SG_prediction

def read_mut_aggr_full(fname):
	short_variants = []
	with open(fname, 'r') as fin:
		firstline = fin.readline()
		firstline_notes = True  if 'Notes' in firstline else False

	with open(fname,'r') as fin:
		if firstline_notes:
			next(fin)
		reader = csv.DictReader(fin, dialect = 'excel-tab')
		for line in reader:
			sv = SV(line['mutation'], line['frequency'], line['depth'], line['pos'])
			short_variants.append(sv)
	return short_variants

def read_pathology_purity_file(pathology_purity_file):
	if not os.path.isfile(pathology_purity_file):
		pathology_purity = 'NA'
		return pathology_purity

	with open(pathology_purity_file) as fin:
		for line in fin:
			pathology_purity = line.strip('\r\n').strip('\n')
			break
	try:
		pathology_purity = float(pathology_purity)/100
	except ValueError:
		pathology_purity = 'NA' 

	return pathology_purity

def _arg_parser():
	"""
	Handle input arguments.

	Returns:
	  Argparse parser object.
	"""
	script_version = globals().get('__version__')
	script_description = globals().get('__doc__')
	script_epilog = None
	script_usage = '''%(prog)s [options] aggregated_mutations_file 
							%(prog)s [-h|--help]
							%(prog)s [--version]'''

	parser = argparse.ArgumentParser(usage=script_usage,
												description=script_description,
												epilog=script_epilog,
												add_help=False)

	g = parser.add_argument_group('Program Help')
	g.add_argument('-h', '--help', action='help',
						help='show this help message and exit')
	g.add_argument('--version', action='version', version=script_version)

	g = parser.add_argument_group('Required arguments')
	g.add_argument('aggregated_mutations_file', action='store', type=str,
						help='Full aggregated mutations text file')

	g = parser.add_argument_group('Optional arguments')
	g.add_argument('-f', dest='pathology_purity_file', action='store', type =str,
						help ='file contains percent of tumor nuclei estimated by pathologists', default = 'NA')

	g.add_argument('-a', dest='alpha', action='store', type=float, default=0.05,
						help='Significance threshold in the bionomial hypothesis testing in core_SGZ method') 
	g.add_argument('-p', dest='pthresh', action='store', type=float, default = 0.9,
						help='If sample purity is higher than pthresh, the somatic/germline prediction \
						of high frequency short-variants cannot be obtained.')
	g.add_argument('-o', dest ='output_header', action='store', type=str, default='',
						help='Output file folder and header')

	return parser


def main(args):

	fname_muts = args.aggregated_mutations_file  # vars/sample11.mut_aggr.full.txt
	pathology_purity_file = args.pathology_purity_file
	alpha = args.alpha
	pthresh = args.pthresh

	out_header = args.output_header
	if len(out_header)==0:
		out_header = fname_muts[0:-4]
	sgz_out = out_header + '.basic.sgz.txt'

	# 1) Read mutations file
	short_variants = read_mut_aggr_full(fname_muts)
	pathology_purity = read_pathology_purity_file(pathology_purity_file)

	# 2) Compute SGZ
	sv_sgz = call_SGZ(short_variants,pathology_purity, alpha, pthresh)
	# 3) Write output
	fout = open(sgz_out, 'w')   
	header = ['mutation', 'pos', 'depth', 'frequency', 'germline/somatic']
	fout.write('\t'.join(header) + '\n')

	for sv in sv_sgz:
		out_list = [sv.mutation,
						'chr%s:%d' %(sv.chr_, sv.position),
						'%d' %sv.depth,
						'%0.2f' %sv.frequency,
						sv.SG_prediction]
		fout.write('\t'.join(out_list)+'\n')

	fout.close()

class SV(object):
	"""lightweight structure for short variants"""
	def __init__(self, mutation, frequency, depth, position):
		super(SV, self).__init__()
		self.mutation = mutation
		self.frequency = float(frequency)
		if self.frequency>1: self.frequency = 1.0

		self.depth = float(depth)

		self.chr_ = position.split(':')[0].strip('chr')
		self.position = int(position.split(':')[1])

if __name__ == "__main__":
	args = _arg_parser().parse_args(sys.argv[1:])
	sys.exit(main(args))

