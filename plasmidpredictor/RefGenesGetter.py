# From ARIBA https://github.com/sanger-pathogens/ariba/blob/master/ariba/ref_genes_getter.py
# bfd6cc9a409828ecc940b554164fabc8e6cc6b9a
# Author: Martin Hunt
# Copyright (c) 2015 - 2018 by Genome Research Ltd.
# GNU GPL version 3

class Error (Exception): pass

import os
import re
import shutil
import pyfastaq
import time
import sys
import subprocess

class RefGenesGetter:
	def __init__(self, verbose=False):
		self.ref_db='plasmidfinder'
		self.verbose = verbose
		self.genetic_code = 11
		self.max_download_attempts = 3
		self.sleep_time = 2

	def _get_from_plasmidfinder(self, outprefix):
		outprefix = os.path.abspath(outprefix)
		final_fasta = outprefix + '.fa'
		final_tsv = outprefix + '.tsv'
		tmpdir = outprefix + '.tmp.download'
		current_dir = os.getcwd()

		try:
			os.mkdir(tmpdir)
			os.chdir(tmpdir)
		except:
			raise Error('Error mkdir/chdir ' + tmpdir)

		cmd = 'curl -o enterobacteriaceae.fsa https://bitbucket.org/genomicepidemiology/plasmidfinder_db/raw/master/enterobacteriaceae.fsa'
		print('Downloading data with:', cmd, sep='\n')
		subprocess.check_call(cmd, shell=True)
		
		cmd = 'curl -o gram_positive.fsa https://bitbucket.org/genomicepidemiology/plasmidfinder_db/raw/master/gram_positive.fsa'
		print('Downloading data with:', cmd, sep='\n')
		subprocess.check_call(cmd, shell=True)

		print('Combining downloaded fasta files...')
		fout_fa = pyfastaq.utils.open_file_write(final_fasta)
		fout_tsv = pyfastaq.utils.open_file_write(final_tsv)
		name_count = {}

		for filename in os.listdir(tmpdir):
			if filename.endswith('.fsa'):
				print('   ', filename)
				file_reader = pyfastaq.sequences.file_reader(os.path.join(tmpdir, filename))
				for seq in file_reader:
					original_id = seq.id
					seq.id = seq.id.replace('_', '.', 1)
					seq.seq = seq.seq.upper()
					if seq.id in name_count:
						name_count[seq.id] += 1
						seq.id = seq.id + '.' + str(name_count[seq.id])
					else:
						name_count[seq.id] = 1

					print(seq, file=fout_fa)
					print(seq.id, '0', '0', '.', '.', 'Original name was ' + original_id, sep='\t', file=fout_tsv)

		pyfastaq.utils.close(fout_fa)
		pyfastaq.utils.close(fout_tsv)
		print('\nFinished combining files\n')
		os.chdir(current_dir)
		if not self.verbose:
			shutil.rmtree(tmpdir)
		print('Finished. Final files are:', final_fasta, final_tsv, sep='\n\t', end='\n\n')
		print('If you use this downloaded data, please cite:')
		print('"PlasmidFinder and pMLST: in silico detection and typing of plasmids", Carattoli et al 2014, PMID: 24777092\n')

	def run(self, outprefix):
		exec('self._get_from_' + self.ref_db + '(outprefix)')
		

