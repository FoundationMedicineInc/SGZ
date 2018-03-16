#!/usr/bin/env python
import os, csv, glob
import sys, argparse

# import pdb
import filecmp

def main(args):
	basicSGZ = args.basicSGZ
	fmiSGZ   = args.fmiSGZ

	test_samples_folder = args.test_samples_folder
	outputfolder        = args.outputfolder
	expectedfolder      = args.expected_output_folder

	os.system(' '.join(['mkdir', '-p', outputfolder]))

	test_sample_list_file = args.test_sample_list_file

	cmp_results = []
	with open(test_sample_list_file) as f:
		reader = csv.DictReader(f, dialect = 'excel-tab')
		for line in reader:
			sample = line['sampleID']
			aggregated_mutations_file = glob.glob(os.path.join(test_samples_folder, '%s*.mut_aggr.full.txt' %sample))[0]
			pathology_purity_file = glob.glob(os.path.join(test_samples_folder, '%s*.pathology_purity.txt' %sample))[0]
			cna_model_file = glob.glob(os.path.join(test_samples_folder, '%s*.cna_calls.txt' %sample))[0]

			# run basicSGZ algorithm
			cmd = ['python', basicSGZ, aggregated_mutations_file, '-f', pathology_purity_file, '-o', os.path.join(outputfolder, sample)]
			os.system(' '.join(cmd))

			# run FMI_SGZ algorithm
			cmd = ['python', fmiSGZ, aggregated_mutations_file, cna_model_file, '-o', os.path.join(outputfolder, sample)]
			os.system(' '.join(cmd))

			## compare the results
			expected_files = [os.path.join(expectedfolder, '%s.basic.sgz.txt' %sample), os.path.join(expectedfolder, '%s.fmi.sgz.txt' %sample), os.path.join(expectedfolder, '%s.fmi.sgz.full.txt' %sample)]
			result_files = [os.path.join(outputfolder, '%s.basic.sgz.txt' %sample), os.path.join(outputfolder, '%s.fmi.sgz.txt' %sample), os.path.join(outputfolder, '%s.fmi.sgz.full.txt' %sample)]

			for e,r in zip(expected_files, result_files):
				cmp_result = filecmp.cmp(e, r)
				cmp_results.append(cmp_result)
				if not cmp_result:
					print 'Actual output file %s is different from expected output file %s.\n' %(r, e)

	if all(cmp_results):
		print '\n-------Test succeeded.-------\n'
	else:
		print '\n-------Test failed.-------\n'

def _arg_parser():
	"""
	Handle input arguments.

	Returns:
	  Argparse parser object.
	"""
	script_version = globals().get('__version__')
	script_description = globals().get('__doc__')
	script_epilog = None
	script_usage = '''%(prog)s [options] 
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

	g = parser.add_argument_group('Optional arguments')
	g.add_argument('-s', '--test_sample_list_file', action='store', type=str, default = os.path.join(os.getcwd(),'test', 'test_sample_list.txt'),
						help='file contains the ID of the example set.(default=%(default)s)')

	g.add_argument('-f','--fmiSGZ', action='store', type =str, default = os.path.join(os.getcwd(), 'fmiSGZ.py'),
						help ='path to the FMI SGZ script. (default=%(default)s)')
	g.add_argument('-b','--basicSGZ', action='store', type=str, default=os.path.join(os.getcwd(), 'basicSGZ.py'),
						help='path to the basic SGZ script. (default=%(default)s)')

	g.add_argument('-t', '--test_samples_folder', action='store', type=str, default = os.path.join(os.getcwd(), 'test', 'test_samples'),
						help='path to the test sample set provided. (default=%(default)s)')
	g.add_argument('-o','--outputfolder', action='store', type=str, default=os.path.join(os.getcwd(), 'test_results'),
						help='Output results folder. (default=%(default)s)')

	g.add_argument('-e','--expected_output_folder', action='store', type=str, default=os.path.join(os.getcwd(), 'test', 'expected_output'),
						help='Expected results folder. For test purpose. (default=%(default)s)')

	return parser


if __name__ == "__main__":
	args = _arg_parser().parse_args(sys.argv[1:])
	sys.exit(main(args))